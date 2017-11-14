#include <iostream>
#include "NanotubeBuilder.h"
using namespace std;
using namespace NTBuilder;


int main() {
    try {

        NanotubeBuilder builder(LatticeType::FCC, {1, 1, -2}, {-1, 1, 0}, {1, 1, 1}, 25.0, 4.05);
        atoms tube_atoms = builder.BuildNanotube(150.0, 20.0, 45.0);
        atoms def_tube_atoms = builder.BuildDeformNanotube(300.0, 40.0, 70.0);
        atoms pyr_atoms = builder.BuildPyramid(150.0, 2*PI*20.0, 2*PI*45.0, 25.0);
        atoms def_pyr_atoms = builder.BuildDeformPyramid(150.0, 2*PI*20.0, 2*PI*45.0, 25.0);

        string tube_file = "tube.data";
        ofstream out(tube_file);
        OutputLammpsFormat(out, tube_atoms);

        string def_tube_file = "def_tube.data";
        ofstream out2(def_tube_file);
        OutputLammpsFormat(out2, def_tube_atoms);

        string pyr_file = "pyr.data";
        ofstream out3(pyr_file);
        OutputLammpsFormat(out3, pyr_atoms);

        string def_pyr_file = "def_pyr.data";
        ofstream out4(def_pyr_file);
        OutputLammpsFormat(out4, def_pyr_atoms);
    }
    catch (logic_error& e) {
        cout << "Exception caught: " << e.what() << '\n';
        exit(1);
    }
    return 0;
}