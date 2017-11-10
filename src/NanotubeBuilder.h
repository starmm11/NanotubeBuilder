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
class NanotubeBuilder {
public:
    NanotubeBuilder(LatticeType type, r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal, double angle, double al) :
            type_(type), normal_(normal), chiral_angle_(angle), a_(al) {

        if (!CheckPairwiseOrthogonality(x, y, normal)) {
            throw std::logic_error("Chosen directions are not pairwise orthogonal\n");
        }
        if (!CheckRightHanded(x, y, normal)) {
            throw std::logic_error("Chosen directions are not right\n");
        }
        // rotate x and y according to the chiral angle
        double a = PI_180 * angle;
        x_ = x*cos(a) - y*sin(a);
        y_ = x*sin(a) + y*cos(a);
        rotate_row_ = Norm3VectorsRows(x_, y_, normal_);
        rotate_col_ = Norm3VectorsCols(x_, y_, normal_);
    }

    // convert atom coordinates in lattice to the box
    // rotating according to the x, y, normal directions
    coord LatticeToBox(const coord& atom);

    // convert atom coordinates in box to the lattice
    // rotating according to the x, y, normal directions
    coord BoxToLattice(const coord& atom);

    /**
     * Calculate Lattice spacing in the rotated system coordinate
     * @return Three distances in the array
     */
    coord LatticeSpacing();

    /**
     * \brief build square plane of one layer
     * \param length_x length along x direction
     * \param length_y length along y direction
     * \param length_z length along z direction ???
     * @return atom coordinates of plane
     */
    atoms BuildSquarePlane(double length_x, double length_y, double length_z);


    atoms BuildPyramid(double length_x, double length_y_in, double length_y_out, double length_z);

    /**
     * function builds nanotube by outer radius, it uses previous function
     * \param chiral_angle Turn angle
     * \param length Length of nanotube
     * \param r_in Inner radius
     * \param r_out Outer radius
     */
    atoms BuildNanotube(double length, double r_in, double r_out);

    /**
     * Builds nanotube by number of layers in it
     * \param chiral_angle Turn angle
     * \param length Length of nanotube
     * \param r_in Inner radius
     * \param n_layers Number of layers
     * @return Atoms of the nanotube
     */
    atoms BuildNanotube(double length, double r_in, int n_layers = 1);

private:
    coord BLbox(const coord& atom,
                double& xmin, double& ymin, double& zmin,
                double& xmax, double& ymax, double& zmax);

    coord LBbox(const coord& atom,
                double& xmin, double& ymin, double& zmin,
                double& xmax, double& ymax, double& zmax);
    /*! \brief Lattice type (FCC, BCC, HCP) */
    LatticeType type_;

    crystal_dir x_;
    crystal_dir y_;
    crystal_dir normal_;
    double chiral_angle_;
    double a_;

    matrix rotate_row_;

    matrix rotate_col_;
};

}


#endif //NANOTUBEBUILDER_NANOTUBEBUILDER_H
