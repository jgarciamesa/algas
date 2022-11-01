/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <algorithm>
#include <cstring>

#include "fasta.hpp"

namespace sasi::seq {
enum struct verb { STOP = 0, FRMST = 1, AMB = 2 };
std::size_t ambiguous(const sasi::args_t& args);
std::pair<size_t, size_t> frameshift(const sasi::args_t& args);
std::vector<std::string> stop_codons(const sasi::args_t& args);
}  // namespace sasi::seq
#endif
