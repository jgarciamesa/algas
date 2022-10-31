/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/sequence.hpp>

namespace sasi::seq {

size_t frameshift(const sasi::args_t& args) {
    size_t total_count{0};
    // for each fasta file in input
    for(const auto& file : args.input) {
        sasi::data_t data = sasi::fasta::read_fasta(file);

        for(const std::string& seq : data.seqs) {
            if(seq.length() % 3 != 0) {
                total_count++;
            }
        }
    }
    return total_count;
}

// Count **early** stop codons
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
            if(seq.length() % 3 == 0) {
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
            if(args.stop_inf == info_detail::SEQ) {
                seq_counts.emplace_back(file + "," + data.names[i] + "," +
                                        std::to_string(counts[2]));
                counts[2] = 0;
            }
        }  // for each sequence
        if(args.stop_inf == info_detail::FILE) {
            file_counts.emplace_back(file + "," + std::to_string(counts[1]));
            counts[1] = 0;
        }
    }  // for each file
    if(args.stop_inf == info_detail::FILE) {
        return file_counts;
    }
    if(args.stop_inf == info_detail::SEQ) {
        return seq_counts;
    }
    return {"stop_codons\n" + std::to_string(counts[0])};
}

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
        std::string file2{">1\nbaaandnnhwk\n>2cccccacacagt"};
        test({file1, file2}, {"test1.fasta", "test2.fasta"}, 17);
    }
}
// GCOVR_EXCL_STOP

}  // namespace sasi::seq
