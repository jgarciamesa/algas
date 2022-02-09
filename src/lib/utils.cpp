/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <algas/utils.hpp>
#include <filesystem>

namespace algas::utils {

file_type_t extract_file_type(std::string path) {
    constexpr auto npos = std::string::npos;

    // trim whitespace
    boost::algorithm::trim(path);

    // Format ext:path
    auto colon = path.find_first_of(':');
    if(colon != npos && colon > 1) {
        auto filepath = path.substr(colon + 1);
        auto ext = "." + path.substr(0, colon);
        return {std::move(filepath), std::move(ext)};
    }
    std::filesystem::path fpath{path};
    return {std::move(path), fpath.extension()};
}

int write_histogram(std::vector<size_t>& counts) {
    // if no gaps, print 1 gap of length zero
    auto any_gaps = std::find_if(begin(counts), end(counts),
                                 [](size_t s) { return s > 0; });
    if(any_gaps == std::end(counts)) {
        std::cout << 0 << "," << 1 << std::endl;
        return EXIT_SUCCESS;
    }

    // print all gaps counts that are not zero
    for(size_t i = 0; i < counts.size(); i++) {
        if(counts[i] > 0) {
            std::cout << i << "," << counts[i] << std::endl;
        }
    }
    return EXIT_SUCCESS;
}
}  // namespace algas::utils
