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

namespace NTBuilder {
    // Lattice enum class
    enum class LatticeType { FCC, BCC, HCP };
    typedef std::map<LatticeType, std::vector<std::vector<double> > > lattice_matrices;
    typedef std::vector<std::vector<double> > matrix;
    typedef std::array<int, 3> crystal_dir;
    typedef const std::array<int, 3>& r_crystal_dir;

    template<class coord_type>
    typedef std::array<coord_type, 3> coord;

    static lattice_matrices createAtomBasis() {
        lattice_matrices m;
        m[LatticeType::FCC] = {{0.0, 0.0, 0.0},
                               {0.5, 0.5, 0.0},
                               {0.5, 0.0, 0.5},
                               {0.0, 0.5, 0.5}};
        m[LatticeType::BCC] = {{0.0, 0.0, 0.0},
                               {0.5, 0.5, 0.5}};
        m[LatticeType::HCP] = {{0.0, 0.0, 0.0},
                               {0.5, 0.5, 0.0},
                               {0.5, 5.0/6.0, 0.5},
                               {0.0, 1.0/3.0, 0.5}};
        return m;
    }

    static lattice_matrices createLatticeVectors() {
        lattice_matrices m;
        m[LatticeType::FCC] = {{0.1, 0.0, 0.0},
                               {0.0, 0.1, 0.0},
                               {0.0, 0.0, 0.1}};
        m[LatticeType::BCC] = m[LatticeType::FCC];
        m[LatticeType::HCP] = {{0.1, 0.0, 0.0},
                               {0.0, sqrt(3.0), 0.0},
                               {0.0, 0.0, sqrt(8.0/3.0)}};
        return m;
    }

    static lattice_matrices atom_basis = createAtomBasis();
    static lattice_matrices lattice_vectors = createLatticeVectors();

    int Scalar(r_crystal_dir, r_crystal_dir) const;
    bool CheckPairwiseOrthogonality(r_crystal_dir, r_crystal_dir, r_crystal_dir) const;
    double getLength(r_crystal_dir) const;
    matrix Norm3VectorsMatrix(r_crystal_dir normal, r_crystal_dir x, r_crystal_dir y) const;

    // convert atom coordinates in lattice to the box
    // rotating according to the x, y, normal directions
    template<class coord_type, LatticeType lattice_type>
    coord LatticeToBox(const coord& atom, r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal) const;

    // template parameters: lattice_type, lattice parameters a and c
    // for FCC and BCC templates will be <FCC /*BCC*/, a, a>
    template<class coord_type = double, LatticeType lattice_type, coord_type a, coord_type c>
    class NanotubeBuilder {
    public:
        // typedefs
        typedef std::array<coord_type, 3> coord;
        typedef std::vector<std::array<coord_type, 3> > atoms;

        // function builds nanotube by number of layers in it
        // input: plane normal vector, x, y, chiral angle, nanotube length (Angstrom),
        // inner radius (Angstrom), number of layers in nanotube (Angstrom)
        // ouput: vector of atoms coordinates
        atoms BuildNanotube(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal,
                            double chiral_angle, double length, double r_in, int n_layers) const;


        // create layers of the parallelogramm
        // with an angle equal to chiral_angle
        // loop over all position of atoms in the lattice coordinate
        // rotate them by new directions
        //     -----
        //    /   /
        //   /   /
        //  /   /
        // - - -  chiral angle
        // starting from (0, 0, 0)
        atoms BuildParallellogramm(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal,
                                   double chiral_angle,
                                   double length_x, double length_y, double length_z) const;

        // function builds nanotube by outer radius, it uses previous function
        // input: plane normal vector, x, y, chiral angle, nanotube length (Angstrom),
        // inner radius (Angstrom), outer radius (Angstrom)
        // ouput: vector of atoms coordinates
        atoms BuildNanotube(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal,
                            double chiral_angle, double length, double r_in, double r_out) const;


    private:
        // create atom in the box according to its coordinate in the unit cell

    };


}


#endif //NANOTUBEBUILDER_NANOTUBEBUILDER_H
