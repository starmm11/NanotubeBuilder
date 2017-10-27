#include <iostream>
#include "NanotubeBuilder.h"
using namespace std;
using namespace NTBuilder;


int main() {
    try {
        NanotubeBuilder builder(LatticeType::FCC, {1, -1, 0}, {1, 1, -2}, {1, 1, 1}, 4.05);
        atoms plane_atoms = builder.BuildSquarePlane(90.0, 100.0, 100.0, 2.0);
        ofstream out("plane.data");
        OutputLammpsFormat(out, plane_atoms);
    }
    catch (logic_error& e) {
        cout << "Exception caught: " << e.what() << '\n';
        exit(1);
    }
    return 0;
}