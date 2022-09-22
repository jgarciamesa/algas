/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/gap.hpp>

namespace sasi {

int gap(int argc, char* argv[]) {
    if((argc < 2) || (std::strcmp(argv[1], "help") == 0)) {
        std::cout << "Usage sasi gap command [options]" << std::endl
                  << std::endl;
        std::cout << "Commands available:   help - display this message"
                  << std::endl;
        std::cout << "                      frameshift - number of gaps with "
                     "length not multiple of 3"
                  << std::endl;
        std::cout << "                      count - count number of gaps"
                  << std::endl;
        std::cout << "                      phase - distribution of gaps phases"
                  << std::endl;
        std::cout << "                      position - position of gaps"
                  << std::endl;
        return EXIT_SUCCESS;
    }

    if(std::strcmp(argv[1], "count") == 0) {
        if((argc < 3) || (std::strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    sasi gap count fasta(s)" << std::endl;
        }

        std::vector<size_t> gap_counts = count(argc - 2, argv + 2);

        // write counts to stdout
        return sasi::utils::write_histogram(gap_counts);

    } else if(std::strcmp(argv[1], "frameshift") == 0) {
        if((argc < 3) || (std::strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    sasi gap frameshift fasta(s)" << std::endl;
        }

        std::pair<size_t, size_t> frameshifts =
            frameshift(count(argc - 2, argv + 2));

        std::cout << "number of gaps with length not multiple of 3: "
                  << frameshifts.first << " (" << frameshifts.first << "/"
                  << frameshifts.second << ")" << std::endl;

        return EXIT_SUCCESS;

    } else if(std::strcmp(argv[1], "phase") == 0) {
        if((argc < 3) || (std::strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    sasi gap phase fasta(s)" << std::endl;
            return EXIT_SUCCESS;
        }

        std::vector<float> phase_counts = phase(argc - 2, argv + 2);

        // write counts to stdout
        std::cout << "phase 0:" << phase_counts[0] << '\n';
        std::cout << "phase 1:" << phase_counts[1] << '\n';
        std::cout << "phase 2:" << phase_counts[2] << '\n';
        return EXIT_SUCCESS;
    } else if(std::strcmp(argv[1], "position") == 0) {
        if((argc < 3) || (std::strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    sasi gap position fasta(s)" << std::endl;
            return EXIT_SUCCESS;
        }

        std::vector<size_t> gap_positions = position(argc - 2, argv + 2);

        // write counts to stdout
        return sasi::utils::write_histogram(gap_positions, true);
    }

    std::cout << "Command " << argv[1] << " not supported." << std::endl;
    return EXIT_FAILURE;
}

/// @private
// GCOVR_EXCL_START
// TEST_CASE("gap") {
//     CHECK(sasi::gap("help"))
//     CHECK(sasi::gap("histogram"))
//     CHECK(sasi::gap("frameshift"))
//     CHECK(sasi::gap("phase"))
//     // fail
//     CHECK(!sasi::gap("phhase"))
//     CHECK(!sasi::gap("other"))
// }
// GCOVR_EXCL_STOP

std::vector<size_t> count(int num_files, char* files[]) {
    // gap counts vector
    //  each position is the number of gaps with its length (e.g. value
    //  at position 1 is number of gaps of size 1)
    std::vector<size_t> counts;

    // for each fasta file in input
    for(int file = 0; file < num_files; ++file) {
        // read fasta file
        sasi::data_t data = sasi::fasta::read_fasta(files[file]);

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

std::vector<size_t> position(int num_files, char* files[]) {
    // gap position vector
    //  each value [0,100] is the number of gaps at that position (e.g. value
    //  at position 0 is number of gaps at beginning of sequence)
    std::vector<size_t> gaps(101, 0);

    for(int file = 0; file < num_files; ++file) {
        // read fasta file
        sasi::data_t data = sasi::fasta::read_fasta(files[file]);

        // find gaps on each sequence
        for(const std::string& seq : data.seqs) {
            size_t pos{seq.find(GAP, 0)};
            while(pos != std::string::npos) {
                // store current gap position
                ++gaps[pos / seq.length()];
                // skip length of current gap (only beginning is reported)
                while(seq.at(pos) == GAP) {
                    ++pos;
                }
                // look for next gap
                pos = seq.find(GAP, pos + 1);
            }
        }
    }

    return gaps;
}

std::pair<size_t, size_t> frameshift(std::vector<size_t> counts) {
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

std::vector<float> phase(int num_files, char* files[]) {
    // phase counts vector
    std::vector<float> phase = {0, 0, 0};

    // for each fasta file in input
    for(int file = 0; file < num_files; ++file) {
        // read fasta file
        sasi::data_t data = sasi::fasta::read_fasta(files[file]);

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

    float total = std::accumulate(phase.begin(), phase.end(), 0);
    phase[0] /= total;
    phase[1] /= total;
    phase[2] /= total;

    return phase;
}
}  // namespace sasi
