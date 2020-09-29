//
// Created by Ilya Bryukhanov on 22.10.2017.
//

#ifndef NANOTUBEBUILDER_NANOTUBEBUILDER_H
#define NANOTUBEBUILDER_NANOTUBEBUILDER_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <chrono>

#include "LatticeOperations.h"
#include "PMath.h"

namespace NTBuilder {

    /** \class NanotubeBuilder
     *  Class for Building nanotubes and some other stuff
     */

    /**
    * Build box given by the crystal dirs and number of cells along each direction
    */
    atoms BuildCellBox(LatticeType type, double a_lattice,
        r_crystal_dir x, r_crystal_dir y, r_crystal_dir z, int px, int py, int pz);

    atoms BuildYperBox(LatticeType type, double a_lattice,
                   r_crystal_dir x, r_crystal_dir y, r_crystal_dir z, double length_x, int py, double length_z);

    atoms BuildDeformPyramid(
        LatticeType type, double a,
        r_crystal_dir x, r_crystal_dir y, r_crystal_dir z,
        double length_x, int py, double def, double length_z);


    atoms BuildNanotube(LatticeType type, double a,
            r_crystal_dir x, r_crystal_dir y, r_crystal_dir z,
            double length, double r_in, double r_out);


    void  OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms, int n_types);

    void  OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms, int n_types, double lx, double ly, double lz);

    void  OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms, int n_types, double lx);

}


#endif //NANOTUBEBUILDER_NANOTUBEBUILDER_H
