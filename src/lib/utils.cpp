/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <filesystem>
#include <sasi/utils.hpp>

namespace sasi::utils {

void trim_whitespace(std::string& str) {
    str.erase(std::remove_if(str.begin(), str.end(),
                             [](char x) { return std::isspace(x); }),
              str.end());
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("trim_whitespace") {
    std::string spaces{"  multiple  - spaces   "};
    std::string feed{"\fform feed\f"};
    std::string new_line{"\nnew line \n"};
    std::string c_return{"carriage \rreturn\r"};
    std::string tab{"\ttab\t"};
    std::string vertical{"  vertical\v space\n"};

    sasi::utils::trim_whitespace(spaces);
    sasi::utils::trim_whitespace(feed);
    sasi::utils::trim_whitespace(new_line);
    sasi::utils::trim_whitespace(c_return);
    sasi::utils::trim_whitespace(tab);
    sasi::utils::trim_whitespace(vertical);

    CHECK(spaces == "multiple-spaces");
    CHECK(feed == "formfeed");
    CHECK(new_line == "newline");
    CHECK(c_return == "carriagereturn");
    CHECK(tab == "tab");
    CHECK(vertical == "verticalspace");
}
// GCOVR_EXCL_STOP

file_type_t extract_file_type(std::string path) {
    constexpr auto npos = std::string::npos;

    // trim whitespace
    trim_whitespace(path);

    // Format ext:path
    auto colon = path.find_first_of(':');
    if(colon != npos && colon > 1) {
        auto filepath = path.substr(colon + 1);
        auto ext = "." + path.substr(0, colon);
        return {std::move(filepath), std::move(ext)};
    }
    std::filesystem::path fpath{path};
    return {std::move(path), fpath.extension()};
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("extract_file_type") {
    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_extract_ft = [](const std::string& input, const std::string& path,
                              // NOLINTNEXTLINE(misc-unused-parameters)
                              const std::string& ftype) {
        sasi::file_type_t result = sasi::utils::extract_file_type(input);
        CHECK(result.path == path);
        CHECK(result.type_ext == ftype);
    };

    test_extract_ft("test.fasta", "test.fasta", ".fasta");
    test_extract_ft(" test2.fasta ", "test2.fasta", ".fasta");
    test_extract_ft("test3.phy", "test3.phy", ".phy");
    test_extract_ft("phy : test4 ", "test4", ".phy");
    test_extract_ft("fas:test5.mid ", "test5.mid", ".fas");
    test_extract_ft("", "", "");
}
// GCOVR_EXCL_STOP

sasi::args_t set_cli_options(CLI::App& app) {
    sasi::args_t args;

    // Commands - 1 required: gap & sequence
    args.gap = app.add_subcommand("gap", "Gap information");
    args.seq = app.add_subcommand("sequence", "Sequence information");
    app.require_subcommand(1);

    // Add input positional argument
    args.seq->add_option("input", args.input, "Input file(s) (FASTA format)")
        ->take_all()
        ->check(CLI::ExistingFile);
    args.gap->add_option("input", args.input, "Input file(s) (FASTA format)")
        ->take_all()
        ->check(CLI::ExistingFile);

    // Gap subcommands - 1 required: frameshift, count, position, phase
    args.gap->add_subcommand("frameshift",
                             "Count gaps with length not multiple of 3");
    args.gap->add_subcommand("count", "Count number of gaps");
    args.gap->add_subcommand("position", "Position of gaps");
    args.gap->add_subcommand("phase", "Distribution of gap phases");
    args.gap->require_subcommand(1);

    // Seq subcommands - 1 required: stop, frameshift, ambiguous
    auto* stop = args.seq->add_subcommand("stop", "Count early stop codons");
    args.seq->add_subcommand("frameshift",
                             "Count sequences with length not multiple of 3");
    args.seq->add_subcommand("ambiguous", "Count ambiguous nucleotides");
    args.seq->require_subcommand(1);

    // Command & subcommand specific options & flags
    stop->add_option("-i,--information", args.stop_inf,
                     "Stop codons: total = 0, file = 1, sequence = 2");
    args.seq->add_flag("-g,--discard-gaps", args.discard_gaps,
                       "Remove gaps before analysis");

    return args;
}
}  // namespace sasi::utils
