//
// Created by Ilya Bryukhanov on 01.11.2017.
//

#ifndef NANOTUBEBUILDER_LATTICEOPERATIONS_H
#define NANOTUBEBUILDER_LATTICEOPERATIONS_H

#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <array>

#include "PMath.h"

namespace NTBuilder {
    // Lattice enum class
    enum class LatticeType {
        FCC, BCC, HCP, MOS2
    };
    struct Atom {
        Atom(const std::array<double, 3>& coord, int t) : x(coord), type(t) {}
        std::array<double, 3> x;
        int type;
    };
    typedef std::map<LatticeType, std::vector<std::vector<double> > > lattice_matrices;
    typedef std::map<LatticeType, std::vector<int> > lattice_vector;
    typedef std::vector<std::vector<double> > matrix;
    typedef std::array<double, 3> crystal_dir;
    typedef const std::array<double, 3>& r_crystal_dir;
    typedef std::array<double, 3> coord;
    //typedef std::vector<std::array<double, 3> > atoms;
    typedef std::vector<Atom> atoms;

    crystal_dir operator*(const crystal_dir& dir, double num);

    crystal_dir operator+(const crystal_dir& dir1, const crystal_dir& dir2);

    crystal_dir operator-(const crystal_dir& dir1, const crystal_dir& dir2);

    int Scalar(r_crystal_dir, r_crystal_dir);

    bool CheckPairwiseOrthogonality(r_crystal_dir, r_crystal_dir, r_crystal_dir);

    bool CheckRightHanded(r_crystal_dir, r_crystal_dir, r_crystal_dir);

    double getLength(r_crystal_dir);

    matrix CreateZeroMatrix3x3();

    matrix TransposeMatrix3x3(const matrix& m);

    matrix Norm3VectorsRows(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal);

    matrix Norm3VectorsCols(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal);


    void MultiplyMatrixVector3x3(coord &, const matrix &, const coord &);

    void MultiplyMatrixMatrix3x3(matrix &, const matrix &, const matrix &);

    void MultiplyVectorNumber(coord &y, double num);

    void PrintMatrix3x3(const matrix& a, std::ostream& out);

    const lattice_vector atom_types =
            {{LatticeType::FCC, {1, 1, 1, 1}},
             {LatticeType::BCC, {1, 1}},
             {LatticeType::HCP, {1, 1, 1, 1}},
             {LatticeType::MOS2, {1, 1, 1, 1,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2}}
            };


    const lattice_matrices atom_basis =
            {{LatticeType::FCC, {{0.0, 0.0, 0.0},
                                 {0.5, 0.5, 0.0},
                                 {0.5, 0.0, 0.5},
                                 {0.0, 0.5, 0.5}}
             },
             {LatticeType::BCC, {{0.0, 0.0, 0.0},
                                 {0.5, 0.5, 0.5}}
             },
             {LatticeType::HCP, {{0.0, 0.0, 0.0},
                                 {0.5, 0.5, 0.0},
                                 {0.5, 5.0 / 6.0, 0.5},
                                 {0.0, 1.0 / 3.0, 0.5}}
             },
             {LatticeType::MOS2, {{0.99999950, 1.0 / 3.0, 0.25},
                                 {0.5, 5.0 / 6.0, 0.25},
                                 {0.5, 1.0 / 6.0, 0.75},
                                 {0.0, 2.0 / 3.0, 0.75},
                                 {0.99999950, 1.0 / 3.0, 0.855174},
                                 {0.5, 5.0 / 6.0, 0.855174},
                                 {0.5, 1.0 / 6.0, 0.144826},
                                 {0.0, 2.0 / 3.0, 0.144826},
                                 {0.5, 1.0 / 6.0, 0.355174},
                                 {0.0, 2.0 / 3.0, 0.355174},
                                 {0.99999950, 1.0 / 3.0, 0.644826},
                                 {0.5, 5.0 / 6.0, 0.644826}}
             }};

    const lattice_matrices lattice_vectors =
            {{LatticeType::FCC, {{1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {0.0, 0.0, 1.0}}
             },
             {LatticeType::BCC, {{1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {0.0, 0.0, 1.0}}
             },
             {LatticeType::HCP, {{1.0, 0.0, 0.0},
                                 {0.0, sqrt(3.0), 0.0},
                                 {0.0, 0.0, sqrt(8.0 / 3.0)}}
             },
             {LatticeType::MOS2, {{1.0, 0.0, 0.0},
                                        {0.0, sqrt(3.0), 0.0},
                                        {0.0, 0.0, 4.6638}}
             }};

    const lattice_matrices inv_lattice_vectors =
            {{LatticeType::FCC, {{1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {0.0, 0.0, 1.0}}
             },
             {LatticeType::BCC, {{1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {0.0, 0.0, 1.0}}
             },
             {LatticeType::HCP, {{1.0, 0.0, 0.0},
                                 {0.0, sqrt(1.0 / 3.0), 0.0},
                                 {0.0, 0.0, sqrt(3.0 / 8.0)}}
             },
             {LatticeType::MOS2, {{1.0, 0.0, 0.0},
                                        {0.0, sqrt(1.0 / 3.0), 0.0},
                                        {0.0, 0.0, 0.214417}}
             }};

    /**
     * \brief convert atom coordinates in lattice to the box coordinates without rotation by a chiral angle
     * @return coordinates of atom
     */
    coord LatticeToBoxNoRotation(const coord& atom, LatticeType type, double lattice_a,
                                 const matrix& unit_rows);

    /**
     * \brief convert atom coordinates in box to the lattice coordinates without rotation by a chiral angle
     * @return coordinates of atom
     */
    coord BoxToLatticeNoRotation(const coord& atom, LatticeType type, double lattice_a,
                                 const matrix& unit_cols);


    /**
     * \brief auxiliary function that converts lattice to box coordinate
     * and refreshes boundaries in lattice space (min and max x, y, z values)
     * without chiral angle and uses only initial crystal directions
     */
    coord LBboxNoRotation(const coord& atom, const matrix& unit_rows,
                          LatticeType type, double lattice_a,
                          double& xmin, double& ymin, double& zmin,
                          double& xmax, double& ymax, double& zmax);

    /**
    * \brief auxiliary function that converts box to lattice coordinate
    * and refreshes boundaries in box space (min and max x, y, z values)
    * without chiral angle and uses only initial crystal directions
    */
    coord BLboxNoRotation(const coord& atom, const matrix& unit_cols,
                          LatticeType type, double lattice_a,
                          double& xmin, double& ymin, double& zmin,
                          double& xmax, double& ymax, double& zmax);

    coord LatticeSpacing(LatticeType type, double lattice_a, r_crystal_dir x, r_crystal_dir y, r_crystal_dir z);

    /*
     * Check if orientation is 110 or 112
     */
    bool checkSpecialOrientation(r_crystal_dir x);
}

#endif //NANOTUBEBUILDER_LATTICEOPERATIONS_H
