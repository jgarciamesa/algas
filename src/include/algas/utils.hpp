/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef UTILS_HPP
#define UTILS_HPP

#include <boost/algorithm/string.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "structs.hpp"

namespace algas::utils {

file_type_t extract_file_type(std::string path);

int write_histogram(std::vector<size_t>& count);

}  // namespace algas::utils
#endif
