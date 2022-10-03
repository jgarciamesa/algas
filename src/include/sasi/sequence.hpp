/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <algorithm>
#include <cstring>

#include "fasta.hpp"

namespace sasi {
int sequence(const sasi::args_t& args, const CLI::App& app);
std::size_t ambiguous(const sasi::args_t& args);
}  // namespace sasi
#endif
