#include <iostream>
#include "NanotubeBuilder.h"
using namespace std;
using namespace NTBuilder;


int main() {
    try {

        NanotubeBuilder builder(LatticeType::FCC, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, 45.0, 3.85);
        atoms tube_atoms = builder.BuildNanotube(150.0, 20.0, 45.0);
        atoms pyr_atoms = builder.BuildPyramid(150.0, 2.0*PI*20.0, 2.0*PI*45.0, 25.0);
        string tube_file = "tube.data";
        ofstream out(tube_file);
        OutputLammpsFormat(out, tube_atoms);

        string pyr_file = "pyramid.data";
        ofstream out2(pyr_file);

        OutputLammpsFormat(out2, pyr_atoms);
    }
    catch (logic_error& e) {
        cout << "Exception caught: " << e.what() << '\n';
        exit(1);
    }
    return 0;
}