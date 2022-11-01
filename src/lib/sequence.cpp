/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/sequence.hpp>

namespace sasi::seq {

/**
 * @brief Number of sequences with length not multiple of 3.
 *
 * @return size_t count.
 */
std::pair<size_t, size_t> frameshift(const sasi::args_t& args) {
    size_t frm{0};
    size_t total{0};
    // for each fasta file in input
    for(const auto& file : args.input) {
        sasi::data_t data = sasi::fasta::read_fasta(file);

        for(const std::string& seq : data.seqs) {
            total++;
            if(seq.length() % 3 != 0) {
                frm++;
            }
        }
    }
    std::pair<size_t, size_t> count{frm, total};
    return count;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("sequence_frameshift") {
    sasi::args_t args;
    auto test = [](const sasi::args_t& args,
                   const std::vector<std::string>& seqs,
                   const std::pair<size_t, size_t>& expected = {}) {  // NOLINT
        // create input files
        std::ofstream out;
        REQUIRE(args.input.size() == seqs.size());
        for(size_t i = 0; i < seqs.size(); ++i) {
            out.open(args.input[i]);
            REQUIRE(out);
            out << seqs[i] << std::endl;
            out.close();
        }

        CHECK_EQ(sasi::seq::frameshift(args), expected);
        for(const auto& file : args.input) {  // NOLINT
            REQUIRE(std::filesystem::remove(file));
        }
    };
    SUBCASE("one file - no frameshift") {
        args.input = {"test-frm-1.fa"};
        std::vector<std::string> seqs = {">1\nAAA AAA AAA"};
        std::pair<size_t, size_t> expected{0, 1};
        test(args, seqs, expected);
    }
    SUBCASE("one file - frameshift phase 1") {
        args.input = {"test-frm-1.fa"};
        std::vector<std::string> seqs = {">1\nAAA AAA A"};
        std::pair<size_t, size_t> expected{1, 1};
        test(args, seqs, expected);
    }
    SUBCASE("one file - frameshift phase 2") {
        args.input = {"test-frm-1.fa"};
        std::vector<std::string> seqs = {">1\nAAA AAA AA"};
        std::pair<size_t, size_t> expected{1, 1};
        test(args, seqs, expected);
    }
    SUBCASE("multiple file - 2 frm / 3 total") {
        args.input = {"test-frm-1.fa", "test-frm-2.fa", "test-frm-3.fa"};
        std::vector<std::string> seqs = {">1\nAA- -A- --A --- --- --- -AA",
                                         ">2\nAA- -AC CCA --- T-- --G -A",
                                         ">3\nAA- -A- --A --- TTT --- -"};
        std::pair<size_t, size_t> expected{2, 3};
        test(args, seqs, expected);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Count **early** stop codons.
 *
 * @return std::vector<std::string> number of early stop codons, either by file,
 * by sequence, or total count.
 */
std::vector<std::string> stop_codons(const sasi::args_t& args) {
    std::vector stop_codons{"TAA", "TAG", "TGA"};

    std::vector<size_t> counts{0, 0, 0};  // total, file, seq
    std::vector<std::string> file_counts, seq_counts;
    file_counts.reserve(args.input.size());
    seq_counts.reserve(args.input.size());
    file_counts.emplace_back("filename,stop_codons");
    seq_counts.emplace_back("filename,seqname,stop_codons");

    // for each fasta file in input
    for(const auto& file : args.input) {
        sasi::data_t data = sasi::fasta::read_fasta(file);

        // find early stop codons on each sequence
        for(size_t i = 0; i < data.seqs.size(); ++i) {
            std::string seq{data.seqs[i]};

            if(args.discard_gaps) {
                seq.erase(std::remove(seq.begin(), seq.end(), '-'), seq.end());
            }

            // remove last codon or nucleotides (1 or 2) if sequence length
            // not multiple of 3
            if(seq.length() % 3 == 0 && !args.stop_keep_last) {
                seq = seq.substr(0, seq.length() - 3);
            } else {
                seq = seq.substr(0, seq.length() - (seq.length() % 3));
            }

            for(size_t pos = 0; pos < seq.length(); pos += 3) {
                std::string codon = seq.substr(pos, 3);
                if(std::find(stop_codons.cbegin(), stop_codons.cend(), codon) !=
                   stop_codons.cend()) {
                    counts[0]++;
                    counts[1]++;
                    counts[2]++;
                }  // if found stop codon
            }      // for position in sequence

            // save sequence counts
            if(args.stop_inf == info_detail::SEQ && counts[2] > 0) {
                seq_counts.emplace_back(file + "," + data.names[i] + "," +
                                        std::to_string(counts[2]));
                counts[2] = 0;
            }
        }  // for each sequence
        if(args.stop_inf == info_detail::FILE && counts[1] > 0) {
            file_counts.emplace_back(file + "," + std::to_string(counts[1]));
            counts[1] = 0;
        }
    }  // for each file
    if(args.stop_inf == info_detail::FILE) {
        if(file_counts.size() == 1) {
            file_counts.emplace_back("files,0");
        }
        return file_counts;
    }
    if(args.stop_inf == info_detail::SEQ) {
        if(seq_counts.size() == 1) {
            seq_counts.emplace_back("files,sequences,0");
        }
        return seq_counts;
    }
    return {"stop_codons\n" + std::to_string(counts[0])};
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("sequence_stop_codons") {
    sasi::args_t args;
    auto test = [](const sasi::args_t& args,
                   const std::vector<std::string>& seqs,
                   const std::vector<std::string>& expected) {  // NOLINT
        // create input files
        std::ofstream out;
        REQUIRE(args.input.size() == seqs.size());
        for(size_t i = 0; i < seqs.size(); ++i) {
            out.open(args.input[i]);
            REQUIRE(out);
            out << seqs[i] << std::endl;
            out.close();
        }

        CHECK_EQ(sasi::seq::stop_codons(args), expected);
        for(const auto& file : args.input) {  // NOLINT
            REQUIRE(std::filesystem::remove(file));
        }
    };
    SUBCASE("no stops - total") {
        args.input = {"test-stop-1.fa"};
        std::vector<std::string> seqs = {">1\nAAA AAA AAA"};
        std::vector<std::string> expected{"stop_codons\n0"};
        test(args, seqs, expected);
    }
    SUBCASE("no stops - file") {
        args.input = {"test-stop-1.fa"};
        args.stop_inf = sasi::info_detail::FILE;
        std::vector<std::string> seqs = {">1\nAAA AAA AAA"};
        std::vector<std::string> expected{"filename,stop_codons", "files,0"};
        test(args, seqs, expected);
    }
    SUBCASE("no stops - sequence") {
        args.input = {"test-stop-1.fa"};
        args.stop_inf = sasi::info_detail::SEQ;
        std::vector<std::string> seqs = {">1\nAAA AAA AAA"};
        std::vector<std::string> expected{"filename,seqname,stop_codons",
                                          "files,sequences,0"};
        test(args, seqs, expected);
    }
    SUBCASE("end stop but not early") {
        args.input = {"test-stop-1.fa"};
        std::vector<std::string> seqs = {">1\nAAA AAA TAA"};
        std::vector<std::string> expected{"stop_codons\n0"};
        test(args, seqs, expected);
    }
    SUBCASE("keep last true") {
        args.input = {"test-stop-1.fa"};
        args.stop_keep_last = true;
        std::vector<std::string> seqs = {">1\nAAA AAA TAG"};
        std::vector<std::string> expected{"stop_codons\n1"};
        test(args, seqs, expected);
    }
    SUBCASE("by file") {
        args.input = {"test-stop-1.fa", "test-stop-2.fa", "test-stop-3.fa"};
        args.stop_inf = sasi::info_detail::FILE;
        std::vector<std::string> seqs = {
            ">1\nAAA TAA AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TTA TAA"};
        std::vector<std::string> expected{
            "filename,stop_codons", "test-stop-1.fa,2", "test-stop-2.fa,1"};
        test(args, seqs, expected);
    }
    SUBCASE("by sequence") {
        args.input = {"test-stop-1.fa", "test-stop-2.fa", "test-stop-3.fa"};
        args.stop_inf = sasi::info_detail::SEQ;
        std::vector<std::string> seqs = {
            ">1\nAAA TAA AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TTA TAA"};
        std::vector<std::string> expected{
            "filename,seqname,stop_codons", "test-stop-1.fa,1,1",
            "test-stop-1.fa,3,1", "test-stop-2.fa,3,1"};
        test(args, seqs, expected);
    }
    SUBCASE("by file - keep last true") {
        args.input = {"test-stop-1.fa", "test-stop-2.fa", "test-stop-3.fa"};
        args.stop_inf = sasi::info_detail::FILE;
        args.stop_keep_last = true;
        std::vector<std::string> seqs = {
            ">1\nAAA TAA AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TTA TAA"};
        std::vector<std::string> expected{
            "filename,stop_codons", "test-stop-1.fa,5", "test-stop-2.fa,4",
            "test-stop-3.fa,3"};
        test(args, seqs, expected);
    }
    SUBCASE("by sequence - keep last true") {
        args.input = {"test-stop-1.fa", "test-stop-2.fa", "test-stop-3.fa"};
        args.stop_inf = sasi::info_detail::SEQ;
        args.stop_keep_last = true;
        std::vector<std::string> seqs = {
            ">1\nAAA TAA AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TGA TAA",
            ">1\nAAA TAC AAA TAG\n>2\nAAA AAA AAA TAA\n>3\nAAA AAA TTA TAA"};
        std::vector<std::string> expected{
            "filename,seqname,stop_codons", "test-stop-1.fa,1,2",
            "test-stop-1.fa,2,1",           "test-stop-1.fa,3,2",
            "test-stop-2.fa,1,1",           "test-stop-2.fa,2,1",
            "test-stop-2.fa,3,2",           "test-stop-3.fa,1,1",
            "test-stop-3.fa,2,1",           "test-stop-3.fa,3,1"};
        test(args, seqs, expected);
    }
    SUBCASE("with gaps") {
        args.input = {"test-stop-1.fa"};
        std::vector<std::string> seqs = {">1\n--T AAT --- TAG"};
        std::vector<std::string> expected{"stop_codons\n0"};
        test(args, seqs, expected);
    }
    SUBCASE("discard gaps") {
        args.input = {"test-stop-1.fa"};
        args.discard_gaps = true;
        std::vector<std::string> seqs = {">1\n--T AAT --- TAG"};
        std::vector<std::string> expected{"stop_codons\n1"};
        test(args, seqs, expected);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Count number of ambiguous nucleotides.
 *
 * @return std::size_t count.
 */
std::size_t ambiguous(const sasi::args_t& args) {
    size_t n_amb{0};
    const std::string amb{"ryswkmbdhvnRYSWKMBDHVN"};
    // for each fasta file in input
    for(const auto& file : args.input) {
        sasi::data_t data = sasi::fasta::read_fasta(file);
        // for sequence in file
        for(const std::string& seq : data.seqs) {
            n_amb += std::count_if(seq.begin(), seq.end(), [amb](auto s) {
                return std::any_of(amb.begin(), amb.end(),
                                   [s](auto c) { return s == c; });
            });
        }
    }
    return n_amb;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("sequence_ambiguous") {
    auto test =
        [](const std::vector<std::string>& files,
           const std::vector<std::string>& fnames,
           // NOLINTNEXTLINE(misc-unused-parameters)
           std::size_t expected) {  // NOLINT(clang-diagnostic-unusedvariable)
            std::ofstream out;
            REQUIRE(files.size() == fnames.size());
            size_t size = files.size();

            // write out files
            for(size_t i = 0; i < size; ++i) {
                out.open(fnames[i]);
                REQUIRE(out);
                out << files[i];
                out.close();
            }

            sasi::args_t args;
            args.input = fnames;

            // NOLINTNEXTLINE(clang-diagnostic-unused-variable)
            std::size_t result = sasi::seq::ambiguous(args);
            // NOLINTNEXTLINE(clang-dianostic-unused-variable
            for(const auto& name : fnames) {
                REQUIRE(std::filesystem::remove(name));
            }

            CHECK(result == expected);
        };

    SUBCASE("no ambiguous nucleotides") {
        std::string file{">seq1\nAAAAAAAAAA\n>seq2\nCCCCCCCCC"};
        test({file}, {"test.fasta"}, 0);
    }
    SUBCASE("one file - uppercase") {
        std::string file{">seq1\nAAANNAAAAA\n>seq2\nCCWCCCYCC"};
        test({file}, {"test.fasta"}, 4);
    }
    SUBCASE("multiple files - lowercase") {
        std::string file1{">n\nnannwwyccgtwrk\n"};
        std::string file2{">1\nbaaandnnhwk\n>2\ncccccacacagt"};
        test({file1, file2}, {"test1.fasta", "test2.fasta"}, 17);
    }
}
// GCOVR_EXCL_STOP

}  // namespace sasi::seq
