/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <salsa/sequence.hpp>

namespace salsa {

int sequence(int argc, char* argv[]) {
    if((argc < 2 || (strcmp(argv[1], "help") == 0))) {
        std::cout << "Usage salsa sequence command [options]" << std::endl
                  << std::endl;
        std::cout << "Commands available:   help - display this message"
                  << std::endl;
        std::cout << "                      stop - find early stop codons"
                  << std::endl;
        return EXIT_SUCCESS;
    }

    // stop
    if(strcmp(argv[1], "stop") == 0) {
        if((argc < 3) || (strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    salsa sequence stop fasta(s)" << std::endl;
            return EXIT_SUCCESS;
        }

        // TOOD: add feature specifying level of detail/information:
        //       0: total count, 1: count per file, 2: count per seq per file

        std::vector stop_codons{"TAA", "TAG", "TGA"};

        // // vector storing number of early stop codons for each input fasta
        // file std::vector<std::pair<std::string, size_t>> stop_counts(argc -
        // 2);

        size_t total_count{0};
        // for each fasta file in input
        for(int file = 2; file < argc; file++) {
            salsa::data_t data = salsa::fasta::read_fasta(argv[file]);

            // find early stop codons on each sequence
            for(std::string& seq : data.seqs) {
                // remove last codon or nucleotides (1 or 2) if sequence length
                // not multiple of 3
                if(seq.length() % 3 == 0) {
                    seq = seq.substr(0, seq.length() - 3);
                } else {
                    seq = seq.substr(0, seq.length() - (seq.length() % 3));
                }

                for(size_t pos = 0; pos < seq.length(); pos += 3) {
                    std::string codon = seq.substr(pos, 3);
                    if(std::find(stop_codons.begin(), stop_codons.end(),
                                 codon) != std::end(stop_codons)) {
                        total_count++;
                    }
                }
            }
        }
        std::cout << "total_count" << ',' << total_count << std::endl;
        return EXIT_SUCCESS;
    }

    if(strcmp(argv[1], "frameshift") == 0) {
        if((argc < 3) || (strcmp(arv[2], "help") == 0)) {
            std::cout << "Usage:    salsa sequence frameshift fasta(s)"
                      << std::endl;
            return EXIT_SUCCESS;
        }

        size_t total_count{0};
        // for each fasta file in input
        for(int file = 2; file < argc; file++) {
            salsa::data_t data = salsa::fasta::read_fasta(argv[file]);

            for(const std::string& seq : data.seqs) {
                if(seq.length() % 3 != 0) {
                    total_count++;
                }
            }
        }
        std::cout << "total_count" << ',' << total_count << std::endl;
        return EXIT_SUCCESS;
    }

    std::cout << "Command " << argv[1] << " not supported." << std::endl;
    return EXIT_FAILURE;
}
}  // namespace salsa
