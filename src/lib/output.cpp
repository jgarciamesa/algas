/* Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com> */

#include <doctest.h>

#include <sasi/output.hpp>

namespace sasi::gap::output {

/**
 * @brief Write result from gap::frequency to file or stdout.
 */
void frequency(const std::vector<std::pair<size_t, size_t>>& counts,
               std::ostream& out) {
    // if no gaps, print 1 gap of length zero
    if(counts.empty()) {
        out << "Gap_length,count" << std::endl;
        out << 0 << "," << 0 << std::endl;
    } else {
        out << "Gap_length,count" << std::endl;
        for(const auto& pair : counts) {
            out << pair.first << "," << pair.second << std::endl;
        }
    }
}

/**
 * @brief Write result from gap::frameshift to file or stdout.
 */
void frameshift(const std::pair<size_t, size_t>& gaps, std::ostream& out) {
    out << "frameshifting-gaps,total-gaps" << std::endl
        << gaps.first << "," << gaps.second << std::endl;
}

/**
 * @brief Write result from gap::phase to file or stdout.
 */
void phase(const std::vector<std::vector<size_t>>& phases, std::ostream& out) {
    // write counts to stdout
    out << "phase0,phase1,phase2" << std::endl;
    for(auto file : phases) {
        out << file[0] << "," << file[1] << "," << file[2] << std::endl;
    }
}

/**
 * @brief Write result from gap::position to file or stdout.
 */
void position(const std::vector<size_t>& positions, std::ostream& out) {
    // if no gaps, print 1 gap of length zero
    auto any_gaps = std::find_if(begin(positions), end(positions),
                                 [](size_t s) { return s > 0; });
    out << "position,count" << std::endl;
    if(any_gaps == std::end(positions)) {
        out << 0 << "," << 0 << std::endl;
    } else {
        // print all gaps
        for(size_t i = 1; i < positions.size(); i++) {
            out << i << "," << positions[i] << std::endl;
        }
    }
}

}  // namespace sasi::gap::output

namespace sasi::seq::output {

/**
 * @brief Write result from seq::ambiguous to file or stdout.
 */
void ambiguous(const size_t count, std::ostream& out) {
    out << "ambiguous_nucleotides" << std::endl;
    out << count << std::endl;
}

/**
 * @brief Write result from seq::frameshift to file or stdout.
 */
void frameshift(const std::pair<size_t, size_t> count, std::ostream& out) {
    out << "frameshifts,total" << std::endl;
    out << count.first << "," << count.second << std::endl;
}

/**
 * @brief Write result from seq::stop_codons to file or stdout.
 */
void stop_codons(const std::vector<std::string>& count, std::ostream& out) {
    for(const auto& line : count) {
        out << line << std::endl;
    }
}

void subst(const std::vector<std::size_t>& count, std::ostream& out) {
    out << "phase0,phase1,phase2" << std::endl
        << count[0] << ',' << count[1] << ',' << count[2] << std::endl;
}

TEST_CASE("output") {
    auto test = [](const std::vector<std::string>& expected) {
        std::ifstream in;
        in.open("test.txt");
        REQUIRE(in);
        std::string line;
        for(const auto& item : expected) {  // NOLINT
            getline(in, line);
            CHECK_EQ(line, item);
        }
        REQUIRE(std::filesystem::remove("test.txt"));
    };
    SUBCASE("gap frequency") {
        std::vector<std::pair<size_t, size_t>> counts{
            {1, 5}, {2, 4}, {3, 3}, {4, 2}, {6, 1}};
        std::vector<std::string> expected{
            "Gap_length,count", "1,5", "2,4", "3,3", "4,2", "6,1"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::gap::output::frequency(counts, outfile);
        test(expected);
    }
    SUBCASE("gap frequency - no counts") {
        std::vector<std::pair<size_t, size_t>> counts{};
        std::vector<std::string> expected{"Gap_length,count", "0,0"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::gap::output::frequency(counts, outfile);
        test(expected);
    }
    SUBCASE("gap frameshift") {
        std::pair<size_t, size_t> gaps{10, 30};
        std::vector<std::string> expected{"frameshifting-gaps,total-gaps",
                                          "10,30"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::gap::output::frameshift(gaps, outfile);
        test(expected);
    }
    SUBCASE("gap phase - single file") {
        std::vector<std::vector<size_t>> gaps{{30, 20, 10}};
        std::vector<std::string> expected{{"phase0,phase1,phase2"},
                                          {"30,20,10"}};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::gap::output::phase(gaps, outfile);
        test(expected);
    }
    SUBCASE("gap phase - multiple files") {
        std::vector<std::vector<size_t>> gaps{{30, 20, 10}, {50, 20, 30}};
        std::vector<std::string> expected{
            {"phase0,phase1,phase2"}, {"30,20,10"}, {"50,20,30"}};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::gap::output::phase(gaps, outfile);
        test(expected);
    }
    SUBCASE("gap position") {
        std::vector<size_t> positions{
            0, 2, 1, 1, 2, 3, 1, 2, 3, 1, 4, 3, 4, 2, 2, 2, 2, 2, 4, 1, 2,
            3, 2, 3, 0, 2, 2, 2, 3, 1, 1, 1, 3, 2, 2, 0, 2, 2, 2, 2, 3, 3,
            2, 2, 3, 5, 3, 3, 3, 5, 2, 1, 1, 3, 1, 1, 2, 1, 1, 3, 3, 1, 4,
            4, 2, 4, 1, 5, 2, 2, 2, 3, 1, 1, 5, 3, 4, 0, 0, 1, 1, 2, 6, 2,
            3, 3, 2, 1, 5, 1, 5, 4, 0, 0, 2, 2, 0, 4, 2, 0, 0};
        std::vector<std::string> expected{"position,count",
                                          "1,2",
                                          "2,1",
                                          "3,1",
                                          "4,2",
                                          "5,3",
                                          "6,1",
                                          "7,2",
                                          "8,3",
                                          "9,1",
                                          "10,4",
                                          "11,3",
                                          "12,4",
                                          "13,2",
                                          "14,2",
                                          "15,2",
                                          "16,2",
                                          "17,2",
                                          "18,4",
                                          "19,1",
                                          "20,2",
                                          "21,3",
                                          "22,2",
                                          "23,3",
                                          "24,0",
                                          "25,2",
                                          "26,2",
                                          "27,2",
                                          "28,3",
                                          "29,1",
                                          "30,1",
                                          "31,1",
                                          "32,3",
                                          "33,2",
                                          "34,2",
                                          "35,0",
                                          "36,2",
                                          "37,2",
                                          "38,2",
                                          "39,2",
                                          "40,3",
                                          "41,3",
                                          "42,2",
                                          "43,2",
                                          "44,3",
                                          "45,5",
                                          "46,3",
                                          "47,3",
                                          "48,3",
                                          "49,5",
                                          "50,2",
                                          "51,1",
                                          "52,1",
                                          "53,3",
                                          "54,1",
                                          "55,1",
                                          "56,2",
                                          "57,1",
                                          "58,1",
                                          "59,3",
                                          "60,3",
                                          "61,1",
                                          "62,4",
                                          "63,4",
                                          "64,2",
                                          "65,4",
                                          "66,1",
                                          "67,5",
                                          "68,2",
                                          "69,2",
                                          "70,2",
                                          "71,3",
                                          "72,1",
                                          "73,1",
                                          "74,5",
                                          "75,3",
                                          "76,4",
                                          "77,0",
                                          "78,0",
                                          "79,1",
                                          "80,1",
                                          "81,2",
                                          "82,6",
                                          "83,2",
                                          "84,3",
                                          "85,3",
                                          "86,2",
                                          "87,1",
                                          "88,5",
                                          "89,1",
                                          "90,5",
                                          "91,4",
                                          "92,0",
                                          "93,0",
                                          "94,2",
                                          "95,2",
                                          "96,0",
                                          "97,4",
                                          "98,2",
                                          "99,0",
                                          "100,0"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::gap::output::position(positions, outfile);
        test(expected);
    }
    SUBCASE("gap position - all zero") {
        std::vector<size_t> positions(101, 0);
        std::vector<std::string> expected{"position,count", "0,0"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::gap::output::position(positions, outfile);
        test(expected);
    }
    SUBCASE("sequence ambiguous") {
        std::size_t count = 10;
        std::vector<std::string> expected{"ambiguous_nucleotides", "10"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::seq::output::ambiguous(count, outfile);
        test(expected);
    }
    SUBCASE("sequence frameshift") {
        std::pair<size_t, size_t> count{9, 25};
        std::vector<std::string> expected{"frameshifts,total", "9,25"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::seq::output::frameshift(count, outfile);
        test(expected);
    }
    SUBCASE("sequence stop") {
        std::vector<std::string> count{"stop_codons", "21"};
        std::vector<std::string> expected{"stop_codons", "21"};
        std::ofstream outfile;
        outfile.open("test.txt");
        REQUIRE(outfile);
        sasi::seq::output::stop_codons(count, outfile);
        test(expected);
    }
}
}  // namespace sasi::seq::output
