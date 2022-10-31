/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "structs.hpp"

namespace sasi::gap::output {
void count(const std::vector<size_t>& counts);
void frameshift(const std::pair<size_t, size_t>& gaps);
void phase(const std::vector<float>& phases);
void position(const std::vector<size_t>& positions);
}  // namespace sasi::gap::output

namespace sasi::seq::output {
void ambiguous(const size_t count);
void frameshift(const size_t count);
void stop_codons(const std::vector<std::string>& count);
}  // namespace sasi::seq::output
#endif
