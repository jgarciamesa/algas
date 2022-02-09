/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef FASTA_HPP
#define FASTA_HPP

#include <filesystem>
#include <string>
#include <vector>

#include "structs.hpp"
#include "utils.hpp"

namespace algas::fasta {

algas::data_t read_fasta(const std::string& f_path);
bool write_fasta(algas::data_t& fasta);

}  // namespace algas::fasta

#endif
