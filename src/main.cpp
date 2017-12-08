#include <iostream>
#include "NanotubeBuilder.h"
#include "InputParser.h"
using namespace std;
using namespace NTBuilder;


int main(int argc, char **argv) {
    try {
        /*NanotubeBuilder builder(LatticeType::FCC, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, 75.0, 4.05);
        atoms tube_atoms = builder.BuildDeformNanotube(200.0, 85.0, 90.0);
        string tube_file = "tube.data";
        ofstream out(tube_file);
        OutputLammpsFormat(out, tube_atoms);
        */
        Parser::InputParser(argc, argv);
    }
    catch (logic_error& e) {
        cout << "Exception caught: " << e.what() << '\n';
        exit(1);
    }
    return 0;
}