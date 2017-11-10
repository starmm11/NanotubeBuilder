//
// Created by Ilya Bryukhanov on 23.10.2017.
//

#include "NanotubeBuilder.h"

namespace NTBuilder {
    void OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms) {
        using namespace std::chrono;
        out << "LAMMPS data file written by Nanotube Builder at " <<
              system_clock::to_time_t(system_clock::now()) << '\n';
        out << body_atoms.size() << " atoms" << '\n';
        out << " 1 atom types" << '\n';
        out << " -200.000 200 xlo xhi " << '\n';
        out << " -200.000 200 ylo yhi " << '\n';
        out << " -200.000 200 zlo zhi " << '\n';

        out << "\nAtoms \n\n";
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            out << i+1 << " 1 "
                << body_atoms[i][0] << ' '
                << body_atoms[i][1] << ' '
                << body_atoms[i][2] << '\n';
        }
    }

    coord NanotubeBuilder::LatticeToBox(const coord &atom) {
        coord atom_new;
        coord tmp;
        MultiplyMatrixVector3x3(tmp, lattice_vectors.at(type_), atom);
        MultiplyVectorNumber(tmp, a_);
        MultiplyMatrixVector3x3(atom_new, rotate_row_, tmp);
        return atom_new;
    }


    coord NanotubeBuilder::BoxToLattice(const coord &atom) {
        coord atom_new;
        coord tmp;
        MultiplyMatrixVector3x3(tmp, inv_lattice_vectors.at(type_), atom);
        MultiplyVectorNumber(tmp, 1 / a_);
        MultiplyMatrixVector3x3(atom_new, rotate_col_, tmp);
        return atom_new;
    }

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


    atoms NanotubeBuilder::BuildSquarePlane(double length_x, double length_y, double length_z)
    {
        atoms square_plane_atoms;

        matrix box_corners = {{-0.01, -0.01, -0.01},
                              {length_x, length_y, length_z}
                             };

        double xmin,ymin,zmin,xmax,ymax,zmax;
        xmin = ymin = zmin = BIG;
        xmax = ymax = zmax = -BIG;

        // determine lattice max and min indexes corresponded to box corners
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

        // lattice coordinate of atom
        coord atom_lattice;
        // box coordinate of atom
        coord atom_box;
        for (int k = klo; k < khi; k++) {
            for (int j = jlo; j < jhi; j++) {
                for (int i = ilo; i < ihi; i++) {
                    for (int m = 0; m < atom_basis.at(type_).size(); m++) {
                        atom_lattice[0] = i + atom_basis.at(type_)[m][0];
                        atom_lattice[1] = j + atom_basis.at(type_)[m][1];
                        atom_lattice[2] = k + atom_basis.at(type_)[m][2];

                        // convert from lattice coords to box coords
                        atom_box = LatticeToBox(atom_lattice);

                        // if atom_box does not lie in our region then continue
                        if (atom_box[0] < box_corners[0][0] || atom_box[0] >= box_corners[1][0] ||
                            atom_box[1] < box_corners[0][1] || atom_box[1] >= box_corners[1][1] ||
                            atom_box[2] < box_corners[0][2] || atom_box[2] >= box_corners[1][2]) continue;

                        // if it lies than push it back to output array
                        square_plane_atoms.push_back(atom_box);
                    }
                }
            }
        }

        return square_plane_atoms;
    }

    atoms NanotubeBuilder::BuildPyramid(
            double length_x, double length_y_in, double length_y_out, double length_z)
    {

        atoms pyr_atoms = BuildSquarePlane(length_x, length_y_out, length_z);
        pyr_atoms.erase( remove_if( pyr_atoms.begin(), pyr_atoms.end(),
        [&](const coord& atom) {
            return (atom[2] - (atom[1]-length_y_out)*length_z / (length_y_in - length_y_out) > 0);
        }), pyr_atoms.end() );


        return pyr_atoms;
    }

    atoms NanotubeBuilder::BuildNanotube(double length, double r_in, double r_out) {
        // first create a plane
        std::vector<coord> pyramid = BuildPyramid(length, 2*PI*r_in, 2*PI*r_out, r_out-r_in);

        // then we need to fold a nanotube
        for (int i = 0; i < pyramid.size(); ++i) {
            double z = pyramid[i][2];
            double r_cur = r_out - z;
            double y = pyramid[i][1];
            double a_rad = y / (r_cur);
            double y_new = -r_cur*sin(a_rad);
            double z_new = r_cur*(1.0-cos(a_rad));
            pyramid[i][1] = y_new;
            pyramid[i][2] = pyramid[i][2] + z_new;
        }
        return pyramid;
    }

    atoms NanotubeBuilder::BuildNanotube(double length, double r_in, int n_layers) {
        // first create a plane
        std::vector<coord> plane = BuildSquarePlane(length, 2*PI*r_in, 2.0);

        // then we need to fold a nanotube
        for (int i = 0; i < plane.size(); ++i) {
            double y = plane[i][1];
            double a_rad = y / (r_in+0.5);
            double y_new = -r_in*sin(a_rad);
            double z_new = r_in*(1.0-cos(a_rad));
            plane[i][1] = y_new;
            plane[i][2] = z_new;
        }
        return plane;
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
};










