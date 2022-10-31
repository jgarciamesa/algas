/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/gap.hpp>

namespace sasi::gap {

std::vector<size_t> count(const sasi::args_t& args) {
    // gap counts vector
    //  each position is the number of gaps with its length (e.g. value
    //  at position 1 is number of gaps of size 1)
    std::vector<size_t> counts;

    // for each fasta file in input
    for(const auto& file : args.input) {
        // read fasta file
        sasi::data_t data = sasi::fasta::read_fasta(file);

        // if capacity of vector is smaller than size of aln, resize
        if(counts.capacity() < data.len()) {
            counts.resize(data.len());
        }

        // find gaps on each sequence
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

    return counts;
}

/// @private
// GCOVR_EXCL_START
// TEST_CASE("gap count") {
//     auto test = [](int num, char* files[]) {

//     }
// seq ending with gap
// short then long seq to trigger `resize`
// }
// GCOVR_EXCL_STOP

std::vector<size_t> position(const sasi::args_t& args) {
    // gap position vector
    //  each value [0,100] is the number of gaps at that position (e.g. value
    //  at position 0 is number of gaps at beginning of sequence)
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
TEST_CASE("position") {
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

std::pair<size_t, size_t> frameshift(const std::vector<size_t>& counts) {
    std::pair<size_t, size_t> total{0, 0};  // total_frameshifts, total_gaps
    for(size_t pos = 0; pos < counts.size(); ++pos) {
        if(counts[pos] > 0) {
            total.second += counts[pos];
            if(pos % 3 != 0) {
                total.first += counts[pos];
            }
        }
    }

    return total;
}

std::vector<float> phase(const sasi::args_t& args) {
    // phase counts vector
    std::vector<float> phase = {0, 0, 0};

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

    auto total =
        static_cast<float>(std::accumulate(phase.begin(), phase.end(), 0));
    phase[0] /= total;
    phase[1] /= total;
    phase[2] /= total;

    return phase;
}
}  // namespace sasi::gap
