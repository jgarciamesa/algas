/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef GAP_HPP
#define GAP_HPP

#include <CLI11.hpp>
#include <cstring>

#include "fasta.hpp"
#include "utils.hpp"

namespace sasi {
int gap(const sasi::args_t& args, const CLI::App& app);
std::vector<size_t> count(const sasi::args_t& args);
std::pair<size_t, size_t> frameshift(const std::vector<size_t>& counts);
std::vector<float> phase(const sasi::args_t& args);
std::vector<size_t> position(const sasi::args_t& args);
}  // namespace sasi
#endif
