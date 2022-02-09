/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <algas/fasta.hpp>
#include <algas/utils.hpp>

#define GAP '-'

int main(int argc, char* const argv[]) {
    if((argc < 2) || (strcmp(argv[1], "help") == 0)) {
        std::cout << "Usage:    algas command [options]" << std::endl
                  << std::endl;
        std::cout << "Commands available:   help" << std::endl;
        std::cout << "                      histogram" << std::endl;
        return EXIT_SUCCESS;
    }

    // histogram
    if(strcmp(argv[1], "histogram") == 0) {
        if((argc < 3) || (strcmp(argv[2], "help") == 0)) {
            std::cout << "Usage:    algas histogram fasta" << std::endl;
        }
        // read fasta file
        algas::data_t data = algas::fasta::read_fasta(argv[2]);

        // gap counts vector
        //  each position is the number of gaps with its length (e.g. value at
        //  position 1 is number of gaps of size 1).
        //  assuming no sequence is entirely formed of gaps, we consider gap
        //  lengths from 0 to (longest seq size-1)
        std::vector<size_t> counts(data.len());

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

        // write counts to stdout
        return algas::utils::write_histogram(counts);
    }

    std::cout << "Command not supported." << std::endl;
    return EXIT_FAILURE;
}
