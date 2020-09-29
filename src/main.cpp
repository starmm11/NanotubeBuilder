#include <iostream>
#include "NanoBuilder.h"
#include "InputParser.h"
using namespace std;
using namespace NTBuilder;


int main(int argc, char **argv) {
    try {
        //double lattice = 2.8665; // Fe
        double lattice_al = 4.05; // Al
        double lattice_cu = 3.60; //Cu
        double lattice_ni = 3.52; //Ni
        double lattice_mos2 = 3.19; //MoS2
/*      int px = 44;
        int py = 59;
        int pz = 7;
*/
        //LatticeType type = LatticeType::MOS2;
        coord x = {1, 0, 0};
        coord y = {0, 1, 0};
        coord z = {0, 0, 1};
//        double lx = 1*lattice*lattice_vectors.at(type)[0][0]*getLength(x);
//        double ly = 1*lattice*lattice_vectors.at(type)[1][1]*getLength(y);
//        double lz = 1*lattice*lattice_vectors.at(type)[2][2]*getLength(z);
        double length = 700;
        double r_in = 200;
        double r_out = 210;
        //string box_file = "box.data";
        //ofstream out(box_file);
        //atoms nanotube_atoms = BuildCellBox(type, lattice, x, y, z, 1, 1, 1);
        atoms nanotube_atoms = NTBuilder::BuildNanotube(LatticeType::FCC, lattice_al, x, y, z, length, r_in, r_out);
        ofstream out("0_Al_700_200_210_001.data");
        //OutputLammpsFormat(out, nanotube_atoms, 2, lx, ly, lz);
        OutputLammpsFormat(out, nanotube_atoms, 1);

//        atoms nanotube_atoms = NTBuilder::BuildNanotube(LatticeType::MOS2, lattice_mos2, x, y, z, length, r_in, r_out);
//        ofstream out("0_MoS2_400_250_4layers.data");
//        //OutputLammpsFormat(out, nanotube_atoms, 2, lx, ly, lz);
//        OutputLammpsFormat(out, nanotube_atoms, 2);
//        out.close();
    }
    catch (logic_error& e) {
        cout << "Exception caught: " << e.what() << '\n';
        exit(1);
    }
    return 0;
}