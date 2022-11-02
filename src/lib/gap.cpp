/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/gap.hpp>

namespace sasi::gap {

/**
 * @brief Gap frequency.
 *
 * @param[in] args sasi::args_t contains name of sequence files.
 *
 * @return std::vector<size_t> gap frequency vector - position is length and
 * value is count;
 */

std::vector<std::pair<size_t, size_t>> frequency(const sasi::args_t& args) {
    // gap counts vector -  each position is the number of gaps with its length
    // (e.g. value  at position 1 is number of gaps of size 1)
    std::vector<size_t> counts;

    // for each fasta file in input
    for(const auto& file : args.input) {
        sasi::data_t data = sasi::fasta::read_fasta(file);

        // if capacity of vector is smaller than size of aln, resize
        if(counts.capacity() < data.len()) {
            counts.resize(data.len() + 1);
        }

        // for each sequence find next gap starting at 0
        // if we reach the end, next sequence
        // otherwhise count the length of the gap
        for(const std::string& seq : data.seqs) {
            size_t pos{seq.find(GAP, 0)};
            while(pos != std::string::npos) {
                size_t lcount{0};
                // add current gap count until a nucleotide is found
                while(seq.at(pos) == GAP) {
                    ++pos;
                    ++lcount;
                    if(pos >= seq.size()) {
                        break;
                    }
                }
                counts[lcount]++;
                // look for next gap
                pos = seq.find(GAP, pos + 1);
            }
        }
    }

    // remove zero counts and create vector of pairs <length, count>
    std::vector<std::pair<size_t, size_t>> freqs;
    for(size_t i = 0; i < counts.size(); ++i) {
        if(counts[i] > 0) {
            freqs.emplace_back(i, counts[i]);
        }
    }

    return freqs;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("gap_frequency") {
    sasi::args_t args;
    auto test =
        [](const sasi::args_t& args, const std::vector<std::string>& seqs,
           // NOLINTNEXTLINE
           const std::vector<std::pair<size_t, size_t>>& expected = {}) {
            // create input files
            std::ofstream out;
            REQUIRE(args.input.size() == seqs.size());
            for(size_t i = 0; i < seqs.size(); ++i) {
                out.open(args.input[i]);
                REQUIRE(out);
                out << seqs[i] << std::endl;
                out.close();
            }

            CHECK_EQ(sasi::gap::frequency(args), expected);
            for(const auto& file : args.input) {  // NOLINT
                REQUIRE(std::filesystem::remove(file));
            }
        };
    SUBCASE("one file") {
        args.input = {"test-freq-1.fa"};
        std::vector<std::string> seqs = {">1\nAA- -A- --A --- --- --- -AA"};
        std::vector<std::pair<size_t, size_t>> expected = {
            {2, 1}, {3, 1}, {10, 1}};
        test(args, seqs, expected);
    }
    SUBCASE("multiple file") {
        args.input = {"test-freq-1.fa", "test-freq-2.fa", "test-freq-3.fa"};
        std::vector<std::string> seqs = {">1\nAA- -A- --A --- --- --- -AA",
                                         ">2\nAA- -AC CCA --- TTT --G -AA",
                                         ">3\nAA- -A- --A --- TTT --- -AA"};
        std::vector<std::pair<size_t, size_t>> expected = {
            {1, 1}, {2, 4}, {3, 4}, {4, 1}, {10, 1}};
        test(args, seqs, expected);
    }
    SUBCASE("no gaps") {
        args.input = {"test-freq-1.fa"};
        std::vector<std::string> seqs = {">1\nAAA AAA"};
        std::vector<std::pair<size_t, size_t>> expected = {};
        test(args, seqs, expected);
    }
    SUBCASE("gap at beginning") {
        args.input = {"test-freq-1.fa"};
        std::vector<std::string> seqs = {">1\n--- --A --- --- --- -AA"};
        std::vector<std::pair<size_t, size_t>> expected = {{5, 1}, {10, 1}};
        test(args, seqs, expected);
    }
    SUBCASE("gap at end") {
        args.input = {"test-freq-1.fa", "test-freq.2.fa"};
        std::vector<std::string> seqs = {">1\nAA- -A- --A --- --- --- ---",
                                         ">2\nAA- -A- --A AAA AAA AAA AA-"};
        std::vector<std::pair<size_t, size_t>> expected = {
            {1, 1}, {2, 2}, {3, 2}, {12, 1}};
        test(args, seqs, expected);
    }
    SUBCASE("entire seq is gap") {
        args.input = {"test-freq-1.fa"};
        std::vector<std::string> seqs = {">1\n--- --- --- --- ---"};
        std::vector<std::pair<size_t, size_t>> expected = {{15, 1}};
        test(args, seqs, expected);
    }
    SUBCASE("trigger resize") {  // seq ending with gap short then long seq
        args.input = {"test-freq-1.fa", "test-freq.2.fa"};
        std::vector<std::string> seqs = {">1\nAA- -A-",
                                         ">2\nAA- -A- --A AAA A-A AAA AA-"};
        std::vector<std::pair<size_t, size_t>> expected = {
            {1, 3}, {2, 2}, {3, 1}};
        test(args, seqs, expected);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Position of gaps in sequence.
 *
 * @details Relative position (percentage) of gaps in sequence.
 * Each value [0,100] is the number of gaps at that position
 * (e.g. value at position 0 is number of gaps at beginning of sequence).
 *
 * @returns std::vector<size_t> with number of gaps per position [0,100].
 */

std::vector<size_t> position(const sasi::args_t& args) {
    // gap position vector
    std::vector<size_t> gaps(101, 0);

    for(const auto& file : args.input) {
        // read fasta file
        sasi::data_t data = sasi::fasta::read_fasta(file);

        // find gaps on each sequence
        for(const std::string& seq : data.seqs) {
            size_t pos{0};
            while(pos < seq.length()) {
                // look for next gap
                pos = seq.find(GAP, pos);
                if(pos > seq.length()) {
                    break;
                }
                // store current gap position
                auto percentage = static_cast<float>(pos) /
                                  static_cast<float>(seq.length() - 1) * 100;
                ++gaps[static_cast<size_t>(percentage)];
                // skip length of current gap (only beginning is reported)
                while(pos < seq.length() && seq.at(pos) == GAP) {
                    ++pos;
                }
            }
        }
    }

    return gaps;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("gap_position") {
    SUBCASE("end with gap") {
        std::ofstream out;
        out.open("test-positions.fa");
        out << ">seqA\nAAA---\n";
        out.close();

        sasi::args_t args;
        args.input = {"test-positions.fa"};
        std::vector<size_t> expected(101, 0);
        expected[60] = 1;
        auto pos = position(args);
        for(size_t i = 0; i < expected.size(); ++i) {
            CHECK_EQ(pos[i], expected[i]);
        }
        REQUIRE(std::filesystem::remove("test-positions.fa"));
    }
    SUBCASE("gap length 1 at end") {
        std::ofstream out;
        out.open("test-positions.fa");
        out << ">seqA\nAAA-\n";
        out.close();

        sasi::args_t args;
        args.input = {"test-positions.fa"};
        std::vector<size_t> expected(101, 0);
        expected[100] = 1;
        auto pos = position(args);
        for(size_t i = 0; i < expected.size(); ++i) {
            CHECK_EQ(pos[i], expected[i]);
        }
        REQUIRE(std::filesystem::remove("test-positions.fa"));
    }
    SUBCASE("start with gap") {
        std::ofstream out;
        out.open("test-positions.fa");
        out << ">seqA\n--AA\n";
        out.close();

        sasi::args_t args;
        args.input = {"test-positions.fa"};
        std::vector<size_t> expected(101, 0);
        expected[0] = 1;
        auto pos = position(args);
        for(size_t i = 0; i < expected.size(); ++i) {
            CHECK_EQ(pos[i], expected[i]);
        }
        REQUIRE(std::filesystem::remove("test-positions.fa"));
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Count gaps that occur within codons (position not multiple of 3).
 *
 * @details Run `sasi::gap::frequency` and count how many gaps have length
 * not multiple of 3.
 *
 * @return std::pair<size_t, size_t> Number of frameshifts and total number
 * of gaps for comparison.
 */
std::pair<size_t, size_t> frameshift(
    const std::vector<std::pair<size_t, size_t>>& counts) {
    std::pair<size_t, size_t> total{0, 0};  // total_frameshifts, total_gaps
    if(any_of(counts.begin(), counts.end(), [](std::pair<size_t, size_t> pair) {
           return pair.second <= 0;
       })) {
        throw std::runtime_error("Counts vector contains frequencies of zero");
    }
    for(const auto& pair : counts) {
        total.second += pair.second;
        if(pair.first % 3 != 0) {
            total.first += pair.second;
        }
    }

    return total;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("gap_frameshift") {
    SUBCASE("no gaps") {
        const std::vector<std::pair<size_t, size_t>> counts{};
        CHECK_EQ(frameshift(counts), std::pair<size_t, size_t>(0, 0));
    }
    SUBCASE("one gap") {
        const std::vector<std::pair<size_t, size_t>> counts{{1, 1}};
        CHECK_EQ(frameshift(counts), std::pair<size_t, size_t>(1, 1));
    }
    SUBCASE("no frameshifts") {
        const std::vector<std::pair<size_t, size_t>> counts{
            {3, 7}, {6, 3}, {9, 2}};
        CHECK_EQ(frameshift(counts), std::pair<size_t, size_t>(0, 12));
    }
    SUBCASE("only frameshifts") {
        const std::vector<std::pair<size_t, size_t>> counts{
            {1, 15}, {5, 9}, {8, 4}, {11, 1}};
        CHECK_EQ(frameshift(counts), std::pair<size_t, size_t>(29, 29));
    }
    SUBCASE("mix case") {
        const std::vector<std::pair<size_t, size_t>> counts{
            {1, 30}, {2, 15}, {3, 18}, {5, 8}, {6, 10},
            {7, 3},  {9, 4},  {10, 1}, {12, 1}};
        CHECK_EQ(frameshift(counts), std::pair<size_t, size_t>(57, 90));
    }
    SUBCASE("zeros - fail") {
        const std::vector<std::pair<size_t, size_t>> counts{
            {1, 1}, {2, 0}, {3, 1}};
        CHECK_THROWS_AS(frameshift(counts), std::runtime_error);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Number of gaps of each phase {0, 1, 2}.
 *
 * @return std::vector<float> Vector of size 3 with number of gaps
 * of phase 0, 1, and 2.
 */
std::vector<size_t> phase(const sasi::args_t& args) {
    // phase counts vector
    std::vector<size_t> phase = {0, 0, 0};

    // for each fasta file in input
    for(const auto& file : args.input) {
        // read fasta file
        sasi::data_t data = sasi::fasta::read_fasta(file);

        // find gaps on each sequence
        for(const std::string& seq : data.seqs) {
            size_t pos{seq.find(GAP, 0)};
            while(pos != std::string::npos) {
                // add current gap
                phase[pos % 3]++;
                do {
                    pos++;
                } while((pos + 1) < seq.length() && seq.at(pos + 1) == GAP);
                // look for next gap
                pos = seq.find(GAP, pos + 1);
            }
        }
    }

    return phase;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("gap_phase") {
    sasi::args_t args;
    auto test = [](const sasi::args_t& args,
                   const std::vector<std::string>& seqs,
                   const std::vector<size_t>& expected = {}) {  // NOLINT
        // create input files
        std::ofstream out;
        REQUIRE(args.input.size() == seqs.size());
        for(size_t i = 0; i < seqs.size(); ++i) {
            out.open(args.input[i]);
            REQUIRE(out);
            out << seqs[i] << std::endl;
            out.close();
        }

        CHECK_EQ(sasi::gap::phase(args), expected);
        for(const auto& file : args.input) {  // NOLINT
            REQUIRE(std::filesystem::remove(file));
        }
    };
    SUBCASE("one file") {
        args.input = {"test-phase-1.fa"};
        std::vector<std::string> seqs = {">1\nAA- -A- --A --- --- --- -AA"};
        std::vector<size_t> expected = {1, 0, 2};
        test(args, seqs, expected);
    }
    SUBCASE("multiple file") {
        args.input = {"test-phase-1.fa", "test-freq-2.fa", "test-freq-3.fa"};
        std::vector<std::string> seqs = {">1\nAA- -A- --A --- --- --- -AA",
                                         ">2\nAA- -AC CCA --- T-- --G -AA",
                                         ">3\nAA- -A- --A --- TTT --- -AA"};
        std::vector<size_t> expected = {5, 1, 5};
        test(args, seqs, expected);
    }
    SUBCASE("no gaps") {
        args.input = {"test-phase-1.fa"};
        std::vector<std::string> seqs = {">1\nAAA AAA"};
        std::vector<size_t> expected = {0, 0, 0};
        test(args, seqs, expected);
    }
    SUBCASE("gap at beginning") {
        args.input = {"test-phase-1.fa"};
        std::vector<std::string> seqs = {">1\n--- --A --- --- --- -AA"};
        std::vector<size_t> expected = {2, 0, 0};
        test(args, seqs, expected);
    }
    SUBCASE("gap at end") {
        args.input = {"test-phase-1.fa", "test-phase.2.fa"};
        std::vector<std::string> seqs = {">1\nAA- -AA A-- --- AAA --- --A",
                                         ">2\nAAA AA- -AA AAA --- AAA AA-"};
        std::vector<size_t> expected = {2, 1, 3};
        test(args, seqs, expected);
    }
    SUBCASE("entire seq is gap") {
        args.input = {"test-phase-1.fa"};
        std::vector<std::string> seqs = {">1\n--- --- --- --- ---"};
        std::vector<size_t> expected = {1, 0, 0};
        test(args, seqs, expected);
    }
}
// GCOVR_EXCL_STOP
}  // namespace sasi::gap
