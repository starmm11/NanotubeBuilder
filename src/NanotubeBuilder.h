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

    /**
     * \brief Constructor
     */
    NanotubeBuilder(LatticeType type, r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal, double angle, double al);


    /**
     * \brief convert atom coordinates in lattice to the box coordinates
     * @return coordinates of atom
     */
    coord LatticeToBox(const coord& atom);

    /**
     * \brief convert atom coordinates in box to the lattice coordinates
     * @return coordinates of atom
     */
    coord BoxToLattice(const coord& atom);

    /**
     * Calculate Lattice spacing in the rotated system coordinate
     * @return Three distances in the array
     */
    coord LatticeSpacing();

    /**
     * \brief build box
     * \param length_x length along x direction
     * \param length_y length along y direction
     * \param length_z length along z direction
     * @return atom coordinates of plane
     */
    atoms BuildBox(double length_x, double length_y, double length_z);

    /**
      * \brief build body along the coordinate directions
      * by multiplicate atoms in the cell along each of them
      * \param xlo min x lattice cell
      * \param xhi max x lattice cell
      * \param ylo min y lattice cell
      * \param yhi max y lattice cell
      * \param zlo min z lattice cell
      * \param zhi max z lattice cell
      * \return atoms
      */
    atoms BuildSpacingCell(int xlo, int xhi, int ylo, int yhi, int zlo, int zhi);

    /**
     * \brief build pyramid from a material without any deformation
     * \param length Length of nanotube
     * \param length_y_in shorter side
     * \param length_y_out longer side
     * \param length_z height
     * \return atoms of nanotube
     */
    atoms BuildPyramid(double length_x, double length_y_in, double length_y_out, double length_z);

    /**
     * \brief build pyramid by tension-compression bending of a box
     * median direction (ly_in + ly_out) / 2 is unstressed
     * \param length Length of nanotube
     * \param length_y_in shorter side
     * \param length_y_out longer side
     * \param length_z height
     * \return atoms of nanotube
     */
    atoms BuildDeformPyramid(double length_x, double length_y_in, double length_y_out, double length_z);

    /**
     * \brief build nanotube from a pyramid
     * which is constructed from a pyramid of material
     * without any deformation
     * \param length Length of nanotube
     * \param r_in Inner radius
     * \param r_out Outer radius
     * \return atoms of nanotube
     */
    atoms BuildNanotube(double length, double r_in, double r_out);

    /**
     * function builds nanotube from a pyramid
     * which is constructed by tension-compression bending of a box
     * \param length Length of nanotube
     * \param r_in Inner radius
     * \param r_out Outer radius
     * \return atoms of nanotube
     */
    atoms BuildDeformNanotube(double length, double r_in, double r_out);

private:
    /**
     * \brief auxiliary function that converts box to lattice coordinate
     * and refreshes boundaries in box space (min and max x, y, z values)
     */
    coord BLbox(const coord& atom,
                double& xmin, double& ymin, double& zmin,
                double& xmax, double& ymax, double& zmax);

    /**
     * \brief auxiliary function that converts lattice to box coordinate
     * and refreshes boundaries in lattice space (min and max x, y, z values)
     */
    coord LBbox(const coord& atom,
                double& xmin, double& ymin, double& zmin,
                double& xmax, double& ymax, double& zmax);

    /*! \brief Lattice type (FCC, BCC, HCP) */
    LatticeType type_;

    /*! \brief x direction */
    crystal_dir x_;
    crystal_dir y_;
    crystal_dir normal_;
    double chiral_angle_;

    /*! \brief cell size */
    double a_;

    /*! \brief matrix used for transformation from box to lattice coordinates */
    matrix rotate_row_;

    matrix rotate_col_;
};

}


#endif //NANOTUBEBUILDER_NANOTUBEBUILDER_H
