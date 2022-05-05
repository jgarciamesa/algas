/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <salsa/gap.hpp>

namespace salsa {

int gap(int argc, char* argv[]) {
    if((argc < 2) || (strcmp(argv[1], "help") == 0)) {
        std::cout << "Usage salsa gap command [options]" << std::endl
                  << std::endl;
        std::cout << "Commands available:   help - display this message"
                  << std::endl;
        std::cout << "                      histogram" << std::endl;
        std::cout << "                      phase" << std::endl;
        return EXIT_SUCCESS;
    }

    // histogram
    if(strcmp(argv[1], "histogram") == 0) {
        if((argc < 3) || (strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    salsa gap histogram fasta(s)" << std::endl;
            return EXIT_SUCCESS;
        }

        // gap counts vector
        //  each position is the number of gaps with its length (e.g. value at
        //  position 1 is number of gaps of size 1).
        //  assuming no sequence is entirely formed of gaps, we consider gap
        //  lengths from 0 to (longest seq size-1)
        std::vector<size_t> counts;

        // for each fasta file in input
        for(int file = 2; file < argc; ++file) {
            // read fasta file
            salsa::data_t data = salsa::fasta::read_fasta(argv[file]);

            // if capacity of vector is smaller than size of aln, resize
            if(counts.capacity() < data.len()) {
                counts.resize(data.len());
            }

            // find gaps on each sequence
            for(const std::string& seq : data.seqs) {
                size_t pos{seq.find(GAP, 0)};
                while(pos != std::string::npos) {
                    size_t count{1};
                    // add current gap count until a nucleotide is found
                    do {
                        pos++;
                        count++;
                        if(pos + 1 >= seq.size()) {
                            break;
                        }
                    } while(seq.at(pos + 1) == GAP);
                    counts[count]++;
                    // look for next gap
                    pos = seq.find(GAP, pos + 1);
                }
            }
        }

        // write counts to stdout
        return salsa::utils::write_histogram(counts);
    }
    // phase
    if(strcmp(argv[1], "phase") == 0) {
        if((argc < 3) || (strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    salsa gap phase fasta(s)" << std::endl;
            return EXIT_SUCCESS;
        }

        // phase counts vector
        std::vector<float> phase = {0, 0, 0};

        // for each fasta file in input
        for(int file = 2; file < argc; ++file) {
            // read fasta file
            salsa::data_t data = salsa::fasta::read_fasta(argv[file]);

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

        // write counts to stdout
        float total = std::accumulate(phase.begin(), phase.end(), 0);
        std::cout << "phase 0:" << phase[0] / total << '\n';
        std::cout << "phase 1:" << phase[1] / total << '\n';
        std::cout << "phase 2:" << phase[2] / total << '\n';

        return EXIT_SUCCESS;
    }

    std::cout << "Command " << argv[1] << " not supported." << std::endl;
    return EXIT_FAILURE;
}
}  // namespace salsa
