/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <sasi/gap.hpp>
#include <sasi/sequence.hpp>
#include <sasi/utils.hpp>

int main(int argc, char* argv[]) {
    CLI::App app{"SASi - simple sequence alignment statistics - v0.1.9000"};

    try {
        sasi::args_t args = sasi::utils::set_cli_options(argc, argv, app);
        CLI11_PARSE(app, argc, argv);
        if(app.got_subcommand("gap")) {
            return sasi::gap(args, app);
        }
        if(app.got_subcommand("sequence")) {
            return sasi::sequence(args, app);
        }

        std::cout << app.help() << std::endl;
    } catch(std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
