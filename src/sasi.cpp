/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <sasi/gap.hpp>
#include <sasi/output.hpp>
#include <sasi/sequence.hpp>
#include <sasi/utils.hpp>

int main(int argc, char* argv[]) {
    CLI::App app{"SASi - simple sequence alignment statistics - v0.1.9000"};

    try {
        sasi::args_t args = sasi::utils::set_cli_options(app);
        CLI11_PARSE(app, argc, argv);
        if(app.got_subcommand("gap")) {
            if(args.gap->got_subcommand("count")) {
                sasi::gap::output::count(sasi::gap::count(args));

            } else if(args.gap->got_subcommand("frameshift")) {
                sasi::gap::output::frameshift(
                    sasi::gap::frameshift(sasi::gap::count(args)));

            } else if(args.gap->got_subcommand("phase")) {
                sasi::gap::output::phase(sasi::gap::phase(args));

            } else if(args.gap->got_subcommand("position")) {
                sasi::gap::output::position(sasi::gap::position(args));
            }
            return EXIT_SUCCESS;
        }

        if(app.got_subcommand("sequence")) {
            if(args.seq->got_subcommand("ambiguous")) {
                sasi::seq::output::ambiguous(sasi::seq::ambiguous(args));

            } else if(args.seq->got_subcommand("frameshift")) {
                sasi::seq::output::frameshift(sasi::seq::frameshift(args));

            } else if(args.seq->got_subcommand("stop")) {
                sasi::seq::output::stop_codons(sasi::seq::stop_codons(args));
            }
            return EXIT_SUCCESS;
        }
    } catch(std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }

    return EXIT_FAILURE;
}
