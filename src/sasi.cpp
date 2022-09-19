/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <sasi/gap.hpp>
#include <sasi/sequence.hpp>

int main(int argc, char* argv[]) {
    if((argc < 2) || (strcmp(argv[1], "help") == 0)) {
        std::cout << "Usage:    sasi command [options]" << std::endl
                  << std::endl;
        std::cout << "Commands available:   help" << std::endl;
        std::cout << "                      gap" << std::endl;
        std::cout << "                      sequence" << std::endl;
        return EXIT_SUCCESS;
    }

    // TODO(JJ): use CLI11 to parse cli arguments

    if(strcmp(argv[1], "gap") == 0) {
        return sasi::gap(argc - 1, argv + 1);
    }
    if(strcmp(argv[1], "sequence") == 0) {
        return sasi::sequence(argc - 1, argv + 1);
    }

    std::cout << "Command " << argv[1] << " not supported." << std::endl;
    return EXIT_FAILURE;
}
