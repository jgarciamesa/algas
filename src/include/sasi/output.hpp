/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "structs.hpp"

namespace sasi::gap::output {
void frequency(const std::vector<std::pair<size_t, size_t>>& counts,
               std::ostream& out);
void frameshift(const std::pair<size_t, size_t>& gaps, std::ostream& out);
void phase(const std::vector<size_t>& phases, std::ostream& out);
void position(const std::vector<size_t>& positions, std::ostream& out);
}  // namespace sasi::gap::output

namespace sasi::seq::output {
void ambiguous(const size_t count, std::ostream& out);
void frameshift(const std::pair<size_t, size_t> count, std::ostream& out);
void stop_codons(const std::vector<std::string>& count, std::ostream& out);
}  // namespace sasi::seq::output
#endif
