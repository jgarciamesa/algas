/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/output.hpp>

namespace sasi::gap::output {

void count(const std::vector<size_t>& counts) {
    // if no gaps, print 1 gap of length zero
    auto any_gaps = std::find_if(counts.rbegin(), counts.rend(),
                                 [](size_t s) { return s > 0; });
    if(any_gaps == counts.rend()) {
        std::cout << "Gap_length,count" << std::endl;
        std::cout << 0 << "," << 0 << std::endl;
    } else {
        // print all gaps
        auto last = std::distance(counts.rbegin(), any_gaps);
        // for(auto it = counts.begin(); it != any_gaps; ++it) {
        for(std::int64_t i = 1; i < last; ++i) {
            std::cout << i << "," << counts[i] << std::endl;
        }
    }
}

void frameshift(const std::pair<size_t, size_t>& gaps) {
    std::cout << "number of gaps with length not multiple of 3: " << gaps.first
              << " (" << gaps.first << "/" << gaps.second << ")" << std::endl;
}

void phase(const std::vector<float>& phases) {
    // write counts to stdout
    std::cout << "phase, count" << std::endl;
    std::cout << "phase 0," << phases[0] << std::endl;
    std::cout << "phase 1," << phases[1] << std::endl;
    std::cout << "phase 2," << phases[2] << std::endl;

}  // namespace sasi::gap::output

void position(const std::vector<size_t>& positions) {
    // if no gaps, print 1 gap of length zero
    auto any_gaps = std::find_if(begin(positions), end(positions),
                                 [](size_t s) { return s > 0; });
    if(any_gaps == std::end(positions)) {
        std::cout << "Gap_position,count" << std::endl;
        std::cout << 0 << "," << 0 << std::endl;
    } else {
        // print all gaps
        for(size_t i = 1; i < positions.size(); i++) {
            std::cout << i << "," << positions[i] << std::endl;
        }
    }
}

}  // namespace sasi::gap::output

namespace sasi::seq::output {

void ambiguous(const size_t count) {
    std::cout << "ambiguous_nucleotides" << std::endl;
    std::cout << count << std::endl;
}

void frameshift(const size_t count) {
    std::cout << "frameshifts" << std::endl;
    std::cout << count << std::endl;
}

void stop_codons(const std::vector<std::string>& count) {
    for(const auto& line : count) {
        std::cout << line << std::endl;
    }
}

// void stop;
}  // namespace sasi::seq::output
