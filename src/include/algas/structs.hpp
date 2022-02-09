/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#ifndef STRUCTS_HPP
#define STRUCTS_HPP

#include <algorithm>
#include <filesystem>

namespace algas {

// extracts extension and filename from both file.foo and ext:file.foo
struct file_type_t {
   public:
    std::string path;
    std::string type_ext;
};

struct data_t {
   public:
    std::filesystem::path path;     /*!< path to input file */
    std::vector<std::string> names; /*!< names of fasta sequences */
    std::vector<std::string> seqs;  /*!< fasta sequences */

    data_t() = default;
    explicit data_t(std::filesystem::path p, std::vector<std::string> n = {},
                    std::vector<std::string> s = {})
        : path{std::move(p)}, names{std::move(n)}, seqs{std::move(s)} {}

    /** \brief Return number of names/sequences */
    std::size_t size() {
        if(names.size() != seqs.size()) {
            throw std::invalid_argument(
                "Different number of sequences and names.");
        }
        return names.size();
    }

    /** \brief Return number of names/sequences */
    [[nodiscard]] size_t size() const {
        if(names.size() != seqs.size()) {
            throw std::invalid_argument(
                "Different number of sequences and names.");
        }
        return names.size();
    }

    /** \brief Return length of sequence on position index */
    size_t len(size_t index) { return seqs[index].length(); }
    size_t len() {
        return std::max_element(seqs.begin(), seqs.end(),
                                [](const auto& s1, const auto& s2) {
                                    return s1.size() < s2.size();
                                })
            ->size();
    }
};

}  // namespace algas
#endif
