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

file_type_t extract_file_type(std::string path);

int write_histogram(std::vector<size_t>& count);

}  // namespace sasi::utils
#endif
