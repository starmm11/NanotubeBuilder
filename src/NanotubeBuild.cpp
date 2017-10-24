//
// Created by Ilya Bryukhanov on 23.10.2017.
//

#include "NanotubeBuilder.h"


int NTBuilder::Scalar(r_crystal_dir x, r_crystal_dir y) const {
    return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}

bool NTBuilder::CheckPairwiseOrthogonality(r_crystal_dir x, r_crystal_dir y, r_crystal_dir z) const {
    return !(NTBuilder::Scalar(x, y) || NTBuilder::Scalar(y, z) || NTBuilder::Scalar(x, z));
}

double NTBuilder::getLength(r_crystal_dir x) const {
    return sqrt(static_cast<double>(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
}

matrix NTBuilder::Norm3VectorsMatrix(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal) const {
    matrix m(3);
    for (int i = 0; i < 3; ++i) {
        m[i].resize(3);
    }
    double normal_length = NTBuilder::getLength(normal);
    double x_length = NTBuilder::getLength(x);
    double y_length = NTBuilder::getLength(y);
    for (int i = 0; i < 3; ++i) {
        m[0][i] = x[i]/x_length;
        m[1][i] = y[i]/y_length;
        m[2][i] = normal[i]/normal_length;
    }
    return m;
}

using namespace NTBuilder;

template<class coord_type, LatticeType lattice_type>
coord LatticeToBox(const coord& atom, r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal) const {
    coord atom_new;
    std::vector<std::vector<double> > primitive = lattice_vectors[lattice_type];
    coord_type x1 = primitive[0][0]*atom[0] + primitive[0][1]*atom[1] + primitive[0][2]*atom[2];
    coord_type y1 = primitive[1][0]*atom[0] + primitive[1][1]*atom[1] + primitive[1][2]*atom[2];
    coord_type z1 = primitive[2][0]*atom[0] + primitive[2][1]*atom[1] + primitive[2][2]*atom[2];

    matrix rotaterow = Norm3VectorsMatrix(x, y, normal);

    atom_new[0] = rotaterow[0][0]*x1 + rotaterow[0][1]*y1 + rotaterow[0][2]*z1;
    atom_new[1] = rotaterow[1][0]*x1 + rotaterow[1][1]*y1 + rotaterow[1][2]*z1;
    atom_new[2] = rotaterow[2][0]*x1 + rotaterow[2][1]*y1 + rotaterow[2][2]*z1;
}

template<class coord_type, LatticeType lattice_type, coord_type a>
atoms NanotubeBuilder<coord_type, lattice_type, a, a>::BuildParallellogramm(
        r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal,
        double chiral_angle,
        double length_x, double length_y, double length_z) const
{
    matrix basis = atom_basis[lattice_type];

}




template<class coord_type, LatticeType lattice_type, coord_type a>
atoms NanotubeBuilder<coord_type, lattice_type, a, a>::BuildNanotube
        (r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal,
         double chiral_angle, double length, double r_in, int n_layers) const
{

    if (lattice_type == LatticeType::HCP) {
        throw std::logic_error("HCP lattice requires two cell parameters\n");
    }
    // check orthogonality of the new directions
    if (!CheckPairwiseOrthogonality(x, y, normal)) {
        throw std::logic_error("Chosen directions are not pairwise orthogonal\n");
    }




};



template<class coord_type, LatticeType lattice_type, coord_type a>
atoms NanotubeBuilder<coord_type, lattice_type, a, a>::BuildNanotube
        (r_crystal_dir normal, r_crystal_dir x, r_crystal_dir y,
         double chiral_angle, double length, double r_in, int n_layers) const
{

    if (lattice_type == LatticeType::HCP) {
        throw std::logic_error("HCP lattice requires two cell parameters\n");
    }

    // check orthogonality of the new directions
    if (!CheckPairwiseOrthogonality(x, y, normal)) {
        throw std::logic_error("Chosen directions are not pairwise orthogonal\n");
    }


};


template<class coord_type, LatticeType lattice_type, coord_type a, coord_type c>
atoms NanotubeBuilder<coord_type, lattice_type, a, c>::BuildNanotube
        (r_crystal_dir normal, r_crystal_dir x, r_crystal_dir y,
         double chiral_angle, double length, double r_in, int n_layers) const {

    if (lattice_type == LatticeType::FCC || lattice_type == LatticeType::BCC) {
        throw std::logic_error("FCC and BCC lattice requires only one cell parameter\n");
    }

};


template<class coord_type, LatticeType lattice_type, coord_type a, coord_type c>
atoms NanotubeBuilder<coord_type, lattice_type, a, c>::BuildNanotube
        (r_crystal_dir normal, r_crystal_dir x, r_crystal_dir y,
         double chiral_angle, double length, double r_in, double r_out) const {

    if (lattice_type == LatticeType::FCC || lattice_type == LatticeType::BCC) {
        throw std::logic_error("FCC and BCC lattice requires only one cell parameter\n");
    }

};








