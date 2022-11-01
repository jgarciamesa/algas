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

        // set output stream pointer
        std::ostream* pout(nullptr);
        std::ofstream outfile;
        if(args.output.empty()) {
            pout = &std::cout;
        } else {
            outfile.open(args.output);
            pout = &outfile;
        }
        std::ostream& out = *pout;

        // gap command
        if(app.got_subcommand("gap")) {
            if(args.gap->got_subcommand("frequency")) {
                sasi::gap::output::frequency(sasi::gap::frequency(args), out);

            } else if(args.gap->got_subcommand("frameshift")) {
                sasi::gap::output::frameshift(
                    sasi::gap::frameshift(sasi::gap::frequency(args)), out);

            } else if(args.gap->got_subcommand("phase")) {
                sasi::gap::output::phase(sasi::gap::phase(args), out);

            } else if(args.gap->got_subcommand("position")) {
                sasi::gap::output::position(sasi::gap::position(args), out);
            }
            return EXIT_SUCCESS;
        }

        // sequence command
        if(app.got_subcommand("sequence")) {
            if(args.seq->got_subcommand("ambiguous")) {
                sasi::seq::output::ambiguous(sasi::seq::ambiguous(args), out);

            } else if(args.seq->got_subcommand("frameshift")) {
                sasi::seq::output::frameshift(sasi::seq::frameshift(args), out);

            } else if(args.seq->got_subcommand("stop")) {
                sasi::seq::output::stop_codons(sasi::seq::stop_codons(args),
                                               out);
            }
            return EXIT_SUCCESS;
        }
    } catch(std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }

    return EXIT_FAILURE;
}
