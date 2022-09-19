/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef UTILS_HPP
#define UTILS_HPP

#define GAP '-'

#include <boost/algorithm/string.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

#include "structs.hpp"

namespace sasi::utils {

void trim_whitespace(std::string& str);
file_type_t extract_file_type(std::string path);
int write_histogram(std::vector<size_t>& counts, bool zeros = false,
                    const std::string& out_file = "");

}  // namespace sasi::utils
#endif
