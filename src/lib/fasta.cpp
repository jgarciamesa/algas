/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <filesystem>
#include <sasi/fasta.hpp>

namespace sasi::fasta {

sasi::data_t read_fasta(const std::string& f_path, bool ignore) {
    sasi::data_t fasta(f_path);

    // set input pointer and file type
    std::istream* pin(nullptr);
    std::ifstream infile;  // input file
    sasi::file_type_t in_type = sasi::utils::extract_file_type(f_path);
    if(in_type.path.empty() || in_type.path == "-") {
        pin = &std::cin;  // set to stdin
        in_type.path = "-";
    } else {
        infile.open(f_path);
        if(!infile || !infile.good()) {
            throw std::invalid_argument("Opening input file " + f_path +
                                        " failed.");
        }
        pin = &infile;  // set to file
        in_type = sasi::utils::extract_file_type(f_path);
    }
    std::istream& in = *pin;

    std::string line, name, content;
    while(in.good()) {
        getline(in, line);
        if(line.empty()) {
            continue;  // omit empty lines
        }
        if(line[0] == ';') {  // omit comment lines
            continue;
        }
        if(line[0] == '>') {  // Identifier marker
            if(!name.empty()) {
                if(content.empty()) {  // If we have a name with no seq, remove
                    fasta.names.pop_back();
                } else {
                    fasta.seqs.push_back(content);
                    name.clear();
                }
            }
            // Add name of sequence
            name = line.substr(1);
            fasta.names.push_back(name);
            content.clear();
            continue;
        }
        // Remove spaces
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
        content += line;
    }
    if(!content.empty()) {
        fasta.seqs.push_back(content);  // Add last sequence if needed
    }

    if(fasta.seqs.size() == 0 && !ignore) {
        throw std::invalid_argument("Input file " + f_path + " is empty");
    }

    if(fasta.seqs.size() != fasta.names.size()) {
        throw std::invalid_argument(
            "Different number of sequences and names in " + f_path + ".");
    }

    return fasta;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("read_fasta") {
    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test = [](const std::string_view file, const sasi::data_t& expected) {
        std::ofstream out;
        out.open("test.fasta");
        REQUIRE(out);
        out << file;
        out.close();

        sasi::data_t result = sasi::fasta::read_fasta("test.fasta");
        REQUIRE(std::filesystem::remove("test.fasta"));

        CHECK(result.names == expected.names);
        CHECK(result.seqs == expected.seqs);
    };

    SUBCASE("Read fasta") {
        std::string file{"; comment line\n>1\n\nCTCTGGATAGTC\n>2\nCTATAGTC\n"};
        sasi::data_t expected("test.fasta", {"1", "2"},
                              {"CTCTGGATAGTC", "CTATAGTC"});
        test(file, expected);
    }
    SUBCASE("File not found") {
        REQUIRE_THROWS_AS(read_fasta("test-not-found.fasta"),
                          std::invalid_argument);
    }
    SUBCASE("One sequence, multiple lines") {
        std::string file{"; comment line\n>1\nNTNTGGATAGTC\nACGTACGTACGT\n"};
        sasi::data_t expected("test.fasta", {"1"},
                              {"NTNTGGATAGTCACGTACGTACGT"});
        test(file, expected);
    }
    SUBCASE("Two sequences, spaces between nucleotides") {
        std::string file{
            "; comment line\n>1\n\nCTC TGG ATA GTC\n>2\nCTA TAG TC\n"};
        sasi::data_t expected("test.fasta", {"1", "2"},
                              {"CTCTGGATAGTC", "CTATAGTC"});
        test(file, expected);
    }
    SUBCASE("Empty lines at end of file") {
        std::string file{
            "; comment line\n>nombre\nNTNTGGATAGTC\n>name2\n\n\nAACG"};
        sasi::data_t expected("test.fasta", {"nombre", "name2"},
                              {"NTNTGGATAGTC", "AACG"});
        test(file, expected);
    }
    SUBCASE("Empty sequence") {
        std::string file{">nombre\nNTNTGGATAGTC\n>name2\n>name3\nAACG\n"};
        sasi::data_t expected("test.fasta", {"nombre", "name3"},
                              {"NTNTGGATAGTC", "AACG"});
        test(file, expected);
    }
    SUBCASE("Only name") {
        std::ofstream out;
        out.open("test-seq.fasta");
        REQUIRE(out);
        out << ">name1\n" << std::endl;
        out.close();
        CHECK_THROWS_AS(read_fasta("test-empty-seq.fasta"),
                        std::invalid_argument);
        REQUIRE(std::filesystem::remove("test-seq.fasta"));
    }
    SUBCASE("Diff number of names and sequences") {
        std::ofstream out;
        out.open("test-seq.fasta");
        REQUIRE(out);
        out << "AAA\nAAA" << std::endl;
        out.close();
        CHECK_THROWS_AS(read_fasta("test-seq.fasta"), std::invalid_argument);
        REQUIRE(std::filesystem::remove("test-seq.fasta"));
    }
    SUBCASE("Empty file") {
        std::ofstream out;
        out.open("test-empty.fasta");
        REQUIRE(out);
        out << std::endl;
        out.close();
        CHECK_THROWS_AS(read_fasta("test-empty.fasta"), std::invalid_argument);
        REQUIRE(std::filesystem::remove("test-empty.fasta"));
    }
}
// GCOVR_EXCL_STOP

}  // namespace sasi::fasta
