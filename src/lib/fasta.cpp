/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <salsa/fasta.hpp>

namespace salsa::fasta {

salsa::data_t read_fasta(const std::string& f_path) {
    salsa::data_t fasta(f_path);

    // set input pointer and file type
    std::istream* pin(nullptr);
    std::ifstream infile;  // input file
    salsa::file_type_t in_type = salsa::utils::extract_file_type(f_path);
    if(in_type.path.empty() || in_type.path == "-") {
        pin = &std::cin;  // set to stdin
        in_type.path = "-";
    } else {
        infile.open(f_path);
        if(!infile || !infile.good()) {
            throw std::invalid_argument("Opening input file " + f_path +
                                        " failed.");
        }
        pin = &infile;  // set to file
        in_type = salsa::utils::extract_file_type(f_path);
    }
    std::istream& in = *pin;

    std::string line, name, content;
    while(in.good()) {
        getline(in, line);
        if(line.empty()) {
            continue;  // omit empty lines
        }
        if(line[0] == ';') {  // omit comment lines
            continue;
        }
        if(line[0] == '>') {  // Identifier marker
            if(!name.empty()) {
                fasta.seqs.push_back(content);
                name.clear();
            }
            // Add name of sequence
            name = line.substr(1);
            fasta.names.push_back(name);
            content.clear();
        } else if(!name.empty()) {
            // Remove spaces
            line.erase(remove_if(line.begin(), line.end(), ::isspace),
                       line.end());
            content += line;
        }
    }
    if(!name.empty()) {  // Add last sequence FSA if needed
        fasta.seqs.push_back(content);
    }

    return fasta;
}

}  // namespace salsa::fasta
