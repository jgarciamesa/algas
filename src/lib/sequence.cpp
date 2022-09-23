/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/sequence.hpp>

namespace sasi {

int sequence(int argc, char* argv[]) {
    if((argc < 2 || (std::strcmp(argv[1], "help") == 0))) {
        std::cout << "Usage sasi sequence command [options]" << std::endl
                  << std::endl;
        std::cout << "Commands available:   help - display this message"
                  << std::endl;
        std::cout << "                      stop - find early stop codons"
                  << std::endl;
        std::cout << "                      frameshift - number of sequences"
                     " with length not multiple of 3"
                  << std::endl;
        std::cout << "                      ambiguous - number of ambiguos "
                     "nucleotides"
                  << std::endl;
        return EXIT_SUCCESS;
    }

    // stop
    if(std::strcmp(argv[1], "stop") == 0) {
        if((argc < 3) || (std::strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    sasi sequence stop fasta(s) [options]"
                      << std::endl
                      << std::endl;
            std::cout << " -i\t\t Level of information "
                         "detail. 0: total count, 1: count per file, 2: count "
                         "per file per sequence"
                      << std::endl;
            std::cout << " -g\t\t Discard gaps" << std::endl;
            return EXIT_SUCCESS;
        }

        enum struct info_detail { TOTAL = 0, FILE = 1, SEQ = 2 };
        info_detail info = info_detail::TOTAL;

        std::vector<std::string> files(argv + 2, argv + argc);
        bool gaps{true};

        for(int i = 0; i < argc; i++) {
            if(std::strcmp(argv[i], "-i") == 0) {
                if(i == argc - 1) {
                    std::cout << "Missing value for -i option, default value "
                                 "(0) used."
                              << std::endl;
                    // remove '-i' argument from argv
                    files.erase(files.begin() + i - 2);

                } else {
                    ++i;
                    switch(std::stoi(argv[i])) {
                    case 0:
                        break;  // default value
                    case 1:
                        info = info_detail::FILE;
                        break;
                    case 2:
                        info = info_detail::SEQ;
                        break;
                    default:
                        throw std::invalid_argument("-i value must be {0,1,2}");
                    }

                    // remove '-i' argument and its value
                    files.erase(files.begin() + i - 2, files.begin() + i);
                }
            } else if(std::strcmp(argv[i], "-g") == 0) {
                gaps = false;
                // remove '-d' argument
                files.erase(files.begin() + i - 2);
            }
        }

        std::vector stop_codons{"TAA", "TAG", "TGA"};

        std::vector<size_t> counts{0, 0, 0};  // total, file, seq
        std::vector<std::pair<std::string, size_t>> file_counts;
        file_counts.reserve(files.size());
        std::vector<std::string> seq_counts;
        seq_counts.reserve(files.size());

        // for each fasta file in input
        for(const auto& file : files) {
            sasi::data_t data = sasi::fasta::read_fasta(file);

            // find early stop codons on each sequence
            for(size_t i = 0; i < data.seqs.size(); i++) {
                std::string seq{data.seqs[i]};

                if(!gaps) {
                    seq.erase(std::remove(seq.begin(), seq.end(), '-'),
                              seq.end());
                }

                // remove last codon or nucleotides (1 or 2) if sequence length
                // not multiple of 3
                if(seq.length() % 3 == 0) {
                    seq = seq.substr(0, seq.length() - 3);
                } else {
                    seq = seq.substr(0, seq.length() - (seq.length() % 3));
                }

                // TODO(juanjo): think about using template for each info_detail
                // type and create function `add_count()`
                for(size_t pos = 0; pos < seq.length(); pos += 3) {
                    std::string codon = seq.substr(pos, 3);
                    if(std::find(stop_codons.begin(), stop_codons.end(),
                                 codon) != std::end(stop_codons)) {
                        counts[0]++;
                        counts[1]++;
                        counts[2]++;
                    }  // if found stop codon
                }      // for position in sequence

                // save sequence counts
                if(info == info_detail::SEQ) {
                    seq_counts.emplace_back(file + "," + data.names[i] + "," +
                                            std::to_string(counts[2]));
                    counts[2] = 0;
                }
            }  // for each sequence
            if(info == info_detail::FILE) {
                file_counts.emplace_back(std::pair(file, counts[1]));
                counts[1] = 0;
            }
        }  // for each file
        switch(info) {
        case info_detail::TOTAL:
            std::cout << "total count" << ',' << counts[0] << std::endl;
            break;
        case info_detail::FILE:
            for(const auto& pair : file_counts) {
                std::cout << pair.first << ',' << pair.second << '\n';
            }
            break;
        case info_detail::SEQ:
            for(const auto& seq_count : seq_counts) {
                std::cout << seq_count << std::endl;
            }
            break;
        }
        return EXIT_SUCCESS;
    }
    if(std::strcmp(argv[1], "frameshift") == 0) {
        if((argc < 3) || (std::strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    sasi sequence frameshift fasta(s)"
                      << std::endl;
            return EXIT_SUCCESS;
        }

        size_t total_count{0};
        // for each fasta file in input
        for(int file = 2; file < argc; file++) {
            sasi::data_t data = sasi::fasta::read_fasta(argv[file]);

            for(const std::string& seq : data.seqs) {
                if(seq.length() % 3 != 0) {
                    total_count++;
                }
            }
        }
        std::cout << "total_count" << ',' << total_count << std::endl;
        return EXIT_SUCCESS;
    }
    if(std::strcmp(argv[1], "ambiguous") == 0) {
        if((argc < 3) || (std::strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    sasi sequence ambiguous fasta(s)"
                      << std::endl;
            return EXIT_SUCCESS;
        }

        std::size_t num_ambiguous = ambiguous(argc - 2, argv + 2);
        std::cout << "total number of ambiguous nucleotides: " << num_ambiguous
                  << std::endl;
        return EXIT_SUCCESS;
    }
    std::cout << "Command " << argv[1] << " not supported." << std::endl;
    return EXIT_FAILURE;
}

std::size_t ambiguous(int num_files, char* files[]) {
    size_t n_amb{0};
    const std::string amb{"ryswkmbdhvnRYSWKMBDHVN"};
    // for each fasta file in input
    for(int file = 0; file < num_files; ++file) {
        sasi::data_t data = sasi::fasta::read_fasta(files[file]);
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
    auto test = [](const std::vector<std::string>& files,
                   const std::vector<std::string>& fnames,
                   std::size_t expected) {
        std::ofstream out;
        REQUIRE(files.size() == fnames.size());
        size_t size = files.size();
        char* cstrings[size];

        for(size_t i = 0; i < size; ++i) {
            out.open(fnames[i]);
            REQUIRE(out);
            out << files[i];
            out.close();
            cstrings[i] = const_cast<char*>(fnames[i].c_str());
        }

        std::size_t result = sasi::ambiguous(size, cstrings);
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
        std::string file2{">1\nbaaandnnhwk\n>2cccccacacagt"};
        test({file1, file2}, {"test1.fasta", "test2.fasta"}, 17);
    }
}
// GCOVR_EXCL_STOP

}  // namespace sasi
