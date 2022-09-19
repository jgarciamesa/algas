/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef GAP_HPP
#define GAP_HPP

#include <cstring>

#include "fasta.hpp"
#include "utils.hpp"

namespace sasi {
int gap(int argc, char* argv[]);
std::vector<size_t> count(int num_files, char* files[]);
std::pair<size_t, size_t> frameshift(std::vector<size_t> counts);
std::vector<float> phase(int num_files, char* files[]);
std::vector<size_t> position(int num_files, char* files[]);
}  // namespace sasi
#endif
