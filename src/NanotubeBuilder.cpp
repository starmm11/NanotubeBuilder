//
// Created by Ilya Bryukhanov on 23.10.2017.
//

#include "NanotubeBuilder.h"

namespace NTBuilder {
    int Scalar(r_crystal_dir x, r_crystal_dir y) {
        return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
    }

    bool CheckPairwiseOrthogonality(r_crystal_dir x, r_crystal_dir y, r_crystal_dir z) {
        return !(Scalar(x, y) || Scalar(y, z) || Scalar(x, z));
    }

    bool CheckRightHanded(r_crystal_dir x, r_crystal_dir y, r_crystal_dir z) {
        int xy0 = x[1] * y[2] - x[2] * y[1];
        int xy1 = x[2] * y[0] - x[0] * y[2];
        int xy2 = x[0] * y[1] - x[1] * y[0];
        return xy0 * z[0] + xy1 * z[1] + xy2 * z[2] > 0;
    }

    double getLength(r_crystal_dir x) {
        return sqrt(static_cast<double>(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
    }

    matrix Norm3VectorsRows(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal) {
        matrix m(3);
        for (int i = 0; i < 3; ++i) {
            m[i].resize(3);
        }

        double x_length = getLength(x);
        double y_length = getLength(y);
        double normal_length = getLength(normal);
        for (int i = 0; i < 3; ++i) {
            m[0][i] = x[i] / x_length;
            m[1][i] = y[i] / y_length;
            m[2][i] = normal[i] / normal_length;
        }
        return m;
    }

    matrix Norm3VectorsCols(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal) {
        matrix m(3);
        for (int i = 0; i < 3; ++i) {
            m[i].resize(3);
        }

        double x_length = getLength(x);
        double y_length = getLength(y);
        double normal_length = getLength(normal);
        for (int i = 0; i < 3; ++i) {
            m[i][0] = x[i] / x_length;
            m[i][1] = y[i] / y_length;
            m[i][2] = normal[i] / normal_length;
        }
        return m;
    }

    void MultiplyMatrixVector3x3(coord& y, const matrix& m, const coord &x) {
        y[0] = 0;
        y[1] = 0;
        y[2] = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                y[i] += m[i][j] * x[j];
            }
        }
    }

    void MultiplyVectorNumber(coord &y, double num) {
        for (int i = 0; i < 3; ++i) {
            y[i] *= num;
        }
    }

    void OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms) {
        using namespace std::chrono;
        out << "LAMMPS data file written by Nanotube Builder at " <<
              system_clock::to_time_t(system_clock::now()) << '\n';
        out << body_atoms.size() << " atoms" << '\n';
        out << " 1 atom types" << '\n';
        out << " 0.000 400 xlo xhi " << '\n';
        out << " 0.000 400 ylo yhi " << '\n';
        out << " 0.000 400 zlo zhi " << '\n';

        out << "\nAtoms \n\n";
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            out << i << " 1 "
                << body_atoms[i][0] << ' '
                << body_atoms[i][1] << ' '
                << body_atoms[i][2] << '\n';
        }
    }

    coord NanotubeBuilder::LatticeToBox(const coord &atom) {
        coord atom_new;
        coord tmp;
        MultiplyMatrixVector3x3(tmp, lattice_vectors.at(type), atom);
        MultiplyVectorNumber(tmp, a);
        MultiplyMatrixVector3x3(atom_new, rotate_row, tmp);
        return atom_new;
    }


    coord NanotubeBuilder::BoxToLattice(const coord &atom) {
        coord atom_new;
        coord tmp;
        MultiplyMatrixVector3x3(tmp, inv_lattice_vectors.at(type), atom);
        MultiplyVectorNumber(tmp, 1 / a);
        MultiplyMatrixVector3x3(atom_new, rotate_col, tmp);
        return atom_new;
    }


    coord NanotubeBuilder::BLbox(
            const coord &atom,
            double &xmin, double &ymin, double &zmin,
            double &xmax, double &ymax, double &zmax) {

        coord atom_new = BoxToLattice(atom);
        xmin = std::min(atom_new[0], xmin);
        xmax = std::max(atom_new[0], xmax);
        ymin = std::min(atom_new[1], ymin);
        ymax = std::max(atom_new[1], ymax);
        zmin = std::min(atom_new[2], zmin);
        zmax = std::max(atom_new[2], zmax);
    };


    coord NanotubeBuilder::LBbox(
            const coord &atom,
            double& xmin, double& ymin, double& zmin,
            double& xmax, double& ymax, double& zmax) {

        coord atom_new = LatticeToBox(atom);
        xmin = std::min(atom_new[0], xmin);
        xmax = std::max(atom_new[0], xmax);
        ymin = std::min(atom_new[1], ymin);
        ymax = std::max(atom_new[1], ymax);
        zmin = std::min(atom_new[2], zmin);
        zmax = std::max(atom_new[2], zmax);
    };


    coord NanotubeBuilder::LatticeSpacing() {
        double xmin, ymin, zmin;
        double xmax, ymax, zmax;
        xmin = ymin = zmin = BIG;
        xmax = ymax = zmax = -BIG;
        LBbox({0.0, 0.0, 0.0}, xmin, ymin, zmin, xmax, ymax, zmax);
        LBbox({1.0, 0.0, 0.0}, xmin, ymin, zmin, xmax, ymax, zmax);
        LBbox({0.0, 1.0, 0.0}, xmin, ymin, zmin, xmax, ymax, zmax);
        LBbox({1.0, 1.0, 0.0}, xmin, ymin, zmin, xmax, ymax, zmax);
        LBbox({0.0, 0.0, 1.0}, xmin, ymin, zmin, xmax, ymax, zmax);
        LBbox({1.0, 0.0, 1.0}, xmin, ymin, zmin, xmax, ymax, zmax);
        LBbox({0.0, 1.0, 1.0}, xmin, ymin, zmin, xmax, ymax, zmax);
        LBbox({1.0, 1.0, 1.0}, xmin, ymin, zmin, xmax, ymax, zmax);

        coord lattice;
        lattice[0] = xmax - xmin;
        lattice[1] = ymax - ymin;
        lattice[2] = zmax - zmin;
        return lattice;
    };


    atoms NanotubeBuilder::BuildSquarePlane(
                            double chiral_angle,
                            double length_x, double length_y, double length_z)
    {
        atoms square_plane_atoms;

        double a = PI_180 * chiral_angle;
        matrix box_corners = {{-0.5, -0.5, -0.5},
                              {length_x + length_y * cos(a), length_y * sin(a), length_z}
                             };

        double xmin,ymin,zmin,xmax,ymax,zmax;
        xmin = ymin = zmin = BIG;
        xmax = ymax = zmax = -BIG;

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    BLbox({box_corners[i][0],box_corners[j][1],box_corners[k][2]},
                          xmin,ymin,zmin,xmax,ymax,zmax);
                }
            }
        }

        // ilo:ihi,jlo:jhi,klo:khi = loop bounds for lattice overlap of my subbox
        int ilo,ihi,jlo,jhi,klo,khi;
        ilo = static_cast<int> (xmin) - 1;
        jlo = static_cast<int> (ymin) - 1;
        klo = static_cast<int> (zmin) - 1;
        ihi = static_cast<int> (xmax) + 1;
        jhi = static_cast<int> (ymax) + 1;
        khi = static_cast<int> (zmax) + 1;

        coord atom;
        coord atom_box;
        for (int k = klo; k < khi; k++) {
            for (int j = jlo; j < jhi; j++) {
                for (int i = ilo; i < ihi; i++) {
                    for (int m = 0; m < atom_basis.at(type).size(); m++) {
                        atom[0] = i + atom_basis.at(type)[m][0];
                        atom[1] = j + atom_basis.at(type)[m][1];
                        atom[2] = k + atom_basis.at(type)[m][2];

                        // convert from lattice coords to box coords
                        atom_box = LatticeToBox(atom);

                        if (atom_box[0] < box_corners[0][0] || atom_box[0] >= box_corners[1][0] ||
                            atom_box[1] < box_corners[0][1] || atom_box[1] >= box_corners[1][1] ||
                            atom_box[2] < box_corners[0][2] || atom_box[2] >= box_corners[1][2]) continue;

                        square_plane_atoms.push_back(atom_box);
                    }
                }
            }
        }

        return square_plane_atoms;
    }


    atoms NanotubeBuilder::BuildNanotube
            (double chiral_angle, double length, double r_in, int n_layers) const {

        // check orthogonality of the new directions




    }
};










