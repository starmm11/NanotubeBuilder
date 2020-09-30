# NanotubeBuilder
This program allows to construct a nanotube rolled from an slab of crystalline material which is isotropic along
the radial directions.  

An example of code using the library is below (src/main.cpp):

#include <iostream>
#include "NanoBuilder.h"
  
int main() {
    double lattice_al = 4.05; // aluminum lattice parameter
    coord x = {1, 0, 0};
    coord y = {0, 1, 0};
    coord z = {0, 0, 1};
    double length = 700; // length
    double r_in = 200; // inner radius
    double r_out = 210; // outer radius
    atoms nanotube_atoms = NTBuilder::BuildNanotube(LatticeType::FCC, lattice_al, x, y, z, length, r_in, r_out);
    ofstream out("0_Al_700_200_210_001.data");
    OutputLammpsFormat(out, nanotube_atoms, 1);
}


Refererence:

Bryukhanov I. A., Gorodtsov V. A., Lisovenko D. S. 
Chiral Fe nanotubes with both negative Poisson’s ratio and Poynting’s effect. Atomistic simulation 
Journal of Physics: Condensed Matter 2019, 47, 475304.
https://doi.org/10.1088/1361-648X/ab3a04
