//
// Created by Ilya Bryukhanov on 22.10.2017.
//

#ifndef NANOTUBEBUILDER_NANOTUBEBUILDER_H
#define NANOTUBEBUILDER_NANOTUBEBUILDER_H

#include <stdexcept>
#include <vector>
#include <map>
#include <array>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <chrono>


const double PI_180 = 0.01745329251994329576923690768489;
const double BIG = 10000000000;

namespace NTBuilder {
    // Lattice enum class
    enum class LatticeType { FCC, BCC, HCP };
    typedef std::map<LatticeType, std::vector<std::vector<double> > > lattice_matrices;
    typedef std::vector<std::vector<double> > matrix;
    typedef std::array<int, 3> crystal_dir;
    typedef const std::array<int, 3>& r_crystal_dir;
    typedef std::array<double, 3> coord;
    typedef std::vector<std::array<double, 3> > atoms;

    int Scalar(r_crystal_dir, r_crystal_dir);
    bool CheckPairwiseOrthogonality(r_crystal_dir, r_crystal_dir, r_crystal_dir);
    bool CheckRightHanded(r_crystal_dir, r_crystal_dir, r_crystal_dir);
    double getLength(r_crystal_dir);
    matrix Norm3VectorsRows(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal);
    matrix Norm3VectorsCols(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal);

    void MultiplyMatrixVector3x3(coord&, const matrix&, const coord&);
    void MultiplyVectorNumber(coord& y, double num);

    void OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms);

    lattice_matrices CreateAtomBasis();
    lattice_matrices CreateLatticeVectors();
    lattice_matrices CreateInvertLatticeVectors();

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
                                 {0.0, 0.0, sqrt(8.0/3.0)}}
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
                                 {0.0, sqrt(3.0/8.0), 0.0},
                                 {0.0, 0.0, sqrt(1.0/3.0)}}
             }};

    // template parameters: lattice_type, lattice parameter a
    class NanotubeBuilder {
    public:
        NanotubeBuilder(LatticeType type, r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal, double al) :
                type(type), x(x), y(y), normal(normal), a(al) {

            if (!CheckPairwiseOrthogonality(x, y, normal)) {
                throw std::logic_error("Chosen directions are not pairwise orthogonal\n");
            }
            if (!CheckRightHanded(x, y, normal)) {
                throw std::logic_error("Chosen directions are not right\n");
            }
            rotate_row = Norm3VectorsRows(x, y, normal);
            rotate_col = Norm3VectorsCols(x, y, normal);
        }

        // convert atom coordinates in lattice to the box
        // rotating according to the x, y, normal directions
        coord LatticeToBox(const coord& atom);

        coord BoxToLattice(const coord& atom);

        coord BLbox(const coord& atom,
                    double& xmin, double& ymin, double& zmin,
                    double& xmax, double& ymax, double& zmax);

        coord LBbox(const coord& atom,
                    double& xmin, double& ymin, double& zmin,
                    double& xmax, double& ymax, double& zmax);

        coord LatticeSpacing();

        // create layers of the parallelogramm
        // with an angle equal to chiral_angle
        // loop over all position of atoms in the lattice coordinate
        // rotate them by new directions
        //     -----
        //    /   /
        //   /   /
        //  /   /
        // - - -  chiral angle
        // (0, 0, 0) (ly * cos(a), ly * sin(a)) (ly * cos(a) + lx, length_y * sin(a))
        atoms BuildSquarePlane(double chiral_angle,
                               double length_x, double length_y, double length_z);


        // function builds nanotube by outer radius, it uses previous function
        // input: plane normal vector, x, y, chiral angle, nanotube length (Angstrom),
        // inner radius (Angstrom), outer radius (Angstrom)
        // ouput: vector of atoms coordinates
        atoms BuildNanotube(double chiral_angle, double length, double r_in, double r_out) const;

        // function builds nanotube by number of layers in it
        // input: plane normal vector, x, y, chiral angle, nanotube length (Angstrom),
        // inner radius (Angstrom), number of layers in nanotube (Angstrom)
        // ouput: vector of atoms coordinates
        atoms BuildNanotube(double chiral_angle, double length, double r_in, int n_layers) const;

    private:
        // create atom in the box according to its coordinate in the unit cell
        LatticeType type;
        crystal_dir x;
        crystal_dir y;
        crystal_dir normal;
        double a;

        matrix rotate_row;

        matrix rotate_col;
    };

}


#endif //NANOTUBEBUILDER_NANOTUBEBUILDER_H
