//
// Created by Ilya Bryukhanov on 04.12.2017.
//
#include "InputParser.h"

namespace Parser {
    bool IsEmpty(std::ifstream &pFile) {
        return pFile.peek() == std::ifstream::traits_type::eof();
    }
    bool IsNumber(const std::string &s) {
        std::string::const_iterator it = s.begin();
        if (*it == '-') ++it;
        while (it != s.end() && std::isdigit(*it) && (*it) == '.') ++it;
        return !s.empty() && it == s.end();
    }
    std::vector<std::string> Split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        Split(s, delim, std::back_inserter(elems));
        return elems;
    }

    void InputParser(int argc, char **argv) {
        if (argc < 2) {
            throw std::logic_error("Too few arguments\n");
        }
        std::string filename = argv[1];
        std::ifstream input(filename);
        if (IsEmpty(input)) {
            throw std::logic_error("No file\n");
        }
        std::string line;
        // get first line
        // parameters of builder
        getline(input, line);

        std::string keyword;
        NTBuilder::LatticeType type;
        double x[3], y[3], z[3];
        double angle;
        double a;

        std::vector<std::string> data = Split(line, ' ');
        if (data.size() != 12) {
            throw std::logic_error("Incorrect number of parameters\n");
        }
        keyword = data[0];
        if (keyword == "FCC" || keyword == "fcc") {
            type = NTBuilder::LatticeType::FCC;
        } else if (keyword == "BCC" || keyword == "bcc") {
            type = NTBuilder::LatticeType::BCC;
        } else if (keyword == "HCP" || keyword == "hcp") {
            type = NTBuilder::LatticeType::HCP;
        } else {
            throw std::logic_error("Bad lattice type initialization\n");
        }
        for (int i = 1; i < 12; ++i) {
            if (IsNumber(data[i])) {
                throw std::logic_error("Incorrect data specification\n");
            }
        }
        x[0] = stof(data[1]);
        x[1] = stof(data[2]);
        x[2] = stof(data[3]);
        y[0] = stof(data[4]);
        y[1] = stof(data[5]);
        y[2] = stof(data[6]);
        z[0] = stof(data[7]);
        z[1] = stof(data[8]);
        z[2] = stof(data[9]);
        angle = stof(data[10]);
        a = stof(data[11]);
        NTBuilder::NanotubeBuilder builder(type, {x[0], x[1], x[2]}, {y[0], y[1], y[2]}, {z[0], z[1], z[2]}, angle, a);
        // get second line
        // parameters of what to build
        getline(input, line);
        data = Split(line, ' ');
        double length, r_in, r_out;
        NTBuilder::atoms tube_atoms;
        if (data[0] == "deform_nanotube") {
            length = stof(data[1]);
            r_in = stof(data[2]);
            r_out = stof(data[3]);
            tube_atoms = builder.BuildDeformNanotube(length, r_in, r_out);
        } else if (data[0] == "lattice_nanotube") {
            length = stof(data[1]);
            r_in = stof(data[2]);
            r_out = stof(data[3]);
            tube_atoms = builder.BuildNanotube(length, r_in, r_out);
        } else {
            throw std::logic_error("Incorrect argument for the nanotube construction\n");
        }

        // get third line
        // output parameters
        getline(input, line);
        std::string output_filename = line;
        std::ofstream out(output_filename);
        NTBuilder::OutputLammpsFormat(out, tube_atoms);
    }
}
