//
// Created by Ilya Bryukhanov on 04.12.2017.
//

#ifndef NANOTUBEBUILDER_INPUTPARSER_H
#define NANOTUBEBUILDER_INPUTPARSER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "NanotubeBuilder.h"

namespace Parser {
    void InputParser(int argc, char **argv);

    bool IsEmpty(std::ifstream &pFile);

    bool IsNumber(const std::string &s);

    template<typename Out>
    void Split(const std::string &s, char delim, Out result) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            *(result++) = item;
        }
    }

    std::vector<std::string> Split(const std::string &s, char delim);

}
#endif //NANOTUBEBUILDER_INPUTPARSER_H
