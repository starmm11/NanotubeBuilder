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
        FCC, BCC, HCP
    };

    typedef std::map<LatticeType, std::vector<std::vector<double> > > lattice_matrices;
    typedef std::vector<std::vector<double> > matrix;
    typedef std::array<double, 3> crystal_dir;
    typedef const std::array<double, 3>& r_crystal_dir;
    typedef std::array<double, 3> coord;
    typedef std::vector<std::array<double, 3> > atoms;

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

    void OutputLammpsFormat(std::ofstream &out, const atoms &body_atoms);

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
             }};
}

#endif //NANOTUBEBUILDER_LATTICEOPERATIONS_H
