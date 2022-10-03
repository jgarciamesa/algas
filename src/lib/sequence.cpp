/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/sequence.hpp>

namespace sasi {

int sequence(const sasi::args_t& args, const CLI::App& app) {
    // stop
    if(args.seq.stop) {
        std::vector stop_codons{"TAA", "TAG", "TGA"};

        std::vector<size_t> counts{0, 0, 0};  // total, file, seq
        std::vector<std::pair<std::string, size_t>> file_counts;
        file_counts.reserve(args.input.size());
        std::vector<std::string> seq_counts;
        seq_counts.reserve(args.input.size());

        // for each fasta file in input
        for(const auto& file : args.input) {
            sasi::data_t data = sasi::fasta::read_fasta(file);

            // find early stop codons on each sequence
            for(size_t i = 0; i < data.seqs.size(); i++) {
                std::string seq{data.seqs[i]};

                if(args.seq.discard_gaps) {
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
                if(args.seq.stop_inf == info_detail::SEQ) {
                    seq_counts.emplace_back(file + "," + data.names[i] + "," +
                                            std::to_string(counts[2]));
                    counts[2] = 0;
                }
            }  // for each sequence
            if(args.seq.stop_inf == info_detail::FILE) {
                file_counts.emplace_back(std::pair(file, counts[1]));
                counts[1] = 0;
            }
        }  // for each file
        switch(args.seq.stop_inf) {
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
    if(args.seq.frameshift) {
        size_t total_count{0};
        // for each fasta file in input
        for(auto file : args.input) {
            sasi::data_t data = sasi::fasta::read_fasta(file);

            for(const std::string& seq : data.seqs) {
                if(seq.length() % 3 != 0) {
                    total_count++;
                }
            }
        }
        std::cout << "total_count" << ',' << total_count << std::endl;
        return EXIT_SUCCESS;
    }
    if(args.seq.ambiguous) {
        std::size_t num_ambiguous = ambiguous(args);
        std::cout << "total number of ambiguous nucleotides: " << num_ambiguous
                  << std::endl;
        return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}

std::size_t ambiguous(const sasi::args_t& args) {
    size_t n_amb{0};
    const std::string amb{"ryswkmbdhvnRYSWKMBDHVN"};
    // for each fasta file in input
    for(auto file : args.input) {
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
            std::size_t result = sasi::ambiguous(args);
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
