/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <salsa/gap.hpp>
#include <salsa/sequence.hpp>

int main(int argc, char* argv[]) {
    if((argc < 2) || (strcmp(argv[1], "help") == 0)) {
        std::cout << "Usage:    salsa command [options]" << std::endl
                  << std::endl;
        std::cout << "Commands available:   help" << std::endl;
        std::cout << "                      gap" << std::endl;
        std::cout << "                      sequence" << std::endl;
        return EXIT_SUCCESS;
    }

    if(strcmp(argv[1], "gap") == 0) {
        return salsa::gap(argc - 1, argv + 1);
    }
    if(strcmp(argv[1], "sequence") == 0) {
        return salsa::sequence(argc - 1, argv + 1);
    }

    std::cout << "Command " << argv[1] << " not supported." << std::endl;
    return EXIT_FAILURE;
}
