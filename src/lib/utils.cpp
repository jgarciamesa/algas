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
    auto test_extract_ft = [](const std::string& input, const std::string& path,
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

int write_histogram(std::vector<size_t>& counts, bool zeros,
                    const std::string& out_file) {
    // if no file specified or file is "-" write to stdout
    std::ostream* pout(nullptr);
    std::ofstream outfile;
    if(out_file.empty() || out_file == "-") {
        pout = &std::cout;
    } else {
        outfile.open(out_file);
        if(!outfile) {
            throw std::invalid_argument("Opening output file " + out_file +
                                        "failed.");
        }
        pout = &outfile;
    }
    std::ostream& out = *pout;

    // if no gaps, print 1 gap of length zero
    auto any_gaps = std::find_if(begin(counts), end(counts),
                                 [](size_t s) { return s > 0; });
    if(any_gaps == std::end(counts)) {
        out << 0 << "," << 0 << std::endl;
        return EXIT_SUCCESS;
    }

    // print all gaps counts that are not zero
    for(size_t i = 1; i < counts.size(); i++) {
        if(zeros || (counts[i] > 0)) {
            out << i << "," << counts[i] << std::endl;
        }
    }
    return EXIT_SUCCESS;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write_histogram") {
    auto test_histogram = [](std::vector<size_t>& counts, std::string& expected,
                             bool zeros = false) {
        CHECK_EQ(
            sasi::utils::write_histogram(counts, zeros, "test_histogram.csv"),
            EXIT_SUCCESS);
        // read output
        std::ifstream infile("test_histogram.csv");
        std::string result;
        infile.seekg(0, std::ios::end);
        result.resize(infile.tellg());
        infile.seekg(0, std::ios::beg);
        infile.read(&result[0], result.size());

        CHECK(result == expected);
        CHECK(std::filesystem::remove("test_histogram.csv"));
    };

    SUBCASE("All zeros") {
        std::vector<size_t> counts = {0, 0, 0, 0, 0};
        std::string expected{"0,0\n"};
        test_histogram(counts, expected);
    }

    SUBCASE("Empty") {
        std::vector<size_t> counts = {};
        std::string expected{"0,0\n"};
        test_histogram(counts, expected);
    }

    SUBCASE("No zeros") {
        std::vector<size_t> counts = {0, 2, 3, 4, 0, 0, 5, 6, 1000};
        std::string expected{"1,2\n2,3\n3,4\n6,5\n7,6\n8,1000\n"};
        test_histogram(counts, expected);
    }

    SUBCASE("Include zeros") {
        std::vector<size_t> counts = {0, 2, 3, 4, 0, 0, 5, 6, 1000};
        std::string expected{"1,2\n2,3\n3,4\n4,0\n5,0\n6,5\n7,6\n8,1000\n"};
        test_histogram(counts, expected, true);
    }
}
// GCOVR_EXCL_STOP
}  // namespace sasi::utils
