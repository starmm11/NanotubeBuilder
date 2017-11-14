//
// Created by Ilya Bryukhanov on 23.10.2017.
//

#include "NanotubeBuilder.h"

namespace NTBuilder {
    NanotubeBuilder::NanotubeBuilder(
            LatticeType type, r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal,
            double angle, double al) :
    type_(type), normal_(normal), chiral_angle_(angle), a_(al)
    {

        if (!CheckPairwiseOrthogonality(x, y, normal)) {
            throw std::logic_error("Chosen directions are not pairwise orthogonal\n");
        }
        if (!CheckRightHanded(x, y, normal)) {
            throw std::logic_error("Chosen directions are not right\n");
        }
        // rotate x and y according to the chiral angle
        double a = PI_180 * angle;
        matrix rotate = {{cos(a), -sin(a), 0.0},
                         {sin(a), cos(a), 0.0},
                         {0.0, 0.0, 1.0}};
        matrix unit_rows = Norm3VectorsRows(x, y, normal);
        matrix unit_cols = Norm3VectorsCols(x, y, normal);
        rotate_row_ = CreateZeroMatrix3x3();
        rotate_col_ = CreateZeroMatrix3x3();
        MultiplyMatrixMatrix3x3(rotate_row_, rotate, unit_rows);
        MultiplyMatrixMatrix3x3(rotate_col_, unit_cols, TransposeMatrix3x3(rotate));
    }

    void OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms) {
        using namespace std::chrono;
        // 0, 1: min, max values of x, y, z
        double boundary[2][3];
        for (int i = 0; i < 3; ++i) {
            boundary[0][i] = BIG;
            boundary[1][i] = -BIG;
        }
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                boundary[0][j] = std::min(body_atoms[i][j], boundary[0][j]);
                boundary[1][j] = std::max(body_atoms[i][j], boundary[1][j]);
            }
        }
        double gap = 40.0;
        for (int i = 0; i < 3; ++i) {
            boundary[0][i] -= gap;
            boundary[1][i] += gap;
        }
        out << "LAMMPS data file written by Nanotube Builder at " <<
              system_clock::to_time_t(system_clock::now()) << '\n';
        out << body_atoms.size() << " atoms" << '\n';
        out << " 1 atom types" << '\n';
        out << " " << boundary[0][0] << ' ' << boundary[1][0] << " xlo xhi " << '\n';
        out << " " << boundary[0][1] << ' ' << boundary[1][1] << " ylo yhi " << '\n';
        out << " " << boundary[0][2] << ' ' << boundary[1][2] << " zlo zhi " << '\n';

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
        MultiplyVectorNumber(tmp, 1.0 / a_);
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

    atoms NanotubeBuilder::BuildSpacingCell(int xlo, int xhi, int ylo, int yhi, int zlo, int zhi) {
        atoms spacing_atoms;
        coord atom_lattice, atom_box;
        for (int k = zlo; k <= zhi; k++) {
            for (int j = ylo; j <= yhi; j++) {
                for (int i = xlo; i <= xhi; i++) {
                    for (int m = 0; m < atom_basis.at(type_).size(); m++) {
                        atom_lattice[0] = i + atom_basis.at(type_)[m][0];
                        atom_lattice[1] = j + atom_basis.at(type_)[m][1];
                        atom_lattice[2] = k + atom_basis.at(type_)[m][2];

                        // convert from lattice coords to box coords
                        atom_box = LatticeToBox(atom_lattice);

                        // if it lies than push it back to output array
                        spacing_atoms.push_back(atom_box);
                    }
                }
            }
        }

        return spacing_atoms;
    }

    atoms NanotubeBuilder::BuildBox(double length_x, double length_y, double length_z)
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
        for (int k = klo; k <= khi; k++) {
            for (int j = jlo; j <= jhi; j++) {
                for (int i = ilo; i <= ihi; i++) {
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

        atoms pyr_atoms = BuildBox(length_x, length_y_out, length_z);
        pyr_atoms.erase( remove_if( pyr_atoms.begin(), pyr_atoms.end(),
        [&](const coord& atom) {
            return (atom[2] - (atom[1]-length_y_out)*length_z / (length_y_in - length_y_out) > 0);
        }), pyr_atoms.end() );


        return pyr_atoms;
    }

    atoms NanotubeBuilder::BuildDeformPyramid(
            double length_x, double length_y_in, double length_y_out, double length_z)
    {
        double length_y_mid = 0.5*(length_y_in + length_y_out);
        atoms pyr_atoms = BuildBox(length_x, length_y_mid, length_z);
        for (int i = 0; i < pyr_atoms.size(); ++i) {
            double z = pyr_atoms[i][2];
            double y_scale = (length_y_in-length_y_out)* z / length_z + length_y_out;
            pyr_atoms[i][1] *= y_scale/length_y_mid;
        }

        return pyr_atoms;
    }

    atoms NanotubeBuilder::BuildDeformNanotube(double length, double r_in, double r_out) {
        // first create a plane
        std::vector<coord> nanotube = BuildDeformPyramid(length, 2*PI*r_in, 2*PI*r_out, r_out-r_in);

        // fold a nanotube
        for (int i = 0; i < nanotube.size(); ++i) {
            double z = nanotube[i][2];
            double r_cur = r_out - z;
            double y = nanotube[i][1];
            double a_rad = y / (r_cur);
            double y_new = -r_cur*sin(a_rad);
            double z_new = r_cur*(1.0-cos(a_rad));
            nanotube[i][1] = y_new;
            nanotube[i][2] = nanotube[i][2] + z_new;
        }
        return nanotube;
    }

    atoms NanotubeBuilder::BuildNanotube(double length, double r_in, double r_out) {
        // first create a plane
        std::vector<coord> nanotube = BuildPyramid(length, 2*PI*r_in, 2*PI*r_out, r_out-r_in);

        // fold a nanotube
        for (int i = 0; i < nanotube.size(); ++i) {
            double z = nanotube[i][2];
            // toDO: for each z layer find maximum and minimum y position
            // toDO: and fold nanotube relating to them, not to the averaging length!!

            double r_cur = r_out - z;
            double y = nanotube[i][1];
            double a_rad = y / (r_cur);
            double y_new = -r_cur*sin(a_rad);
            double z_new = r_cur*(1.0-cos(a_rad));
            nanotube[i][1] = y_new;
            nanotube[i][2] = nanotube[i][2] + z_new;
        }
        return nanotube;
    }


    atoms NanotubeBuilder::BuildNanotube(double length, double r_in, int n_layers) {
        // first create a plane
        std::vector<coord> plane = BuildBox(length, 2 * PI * r_in, 2.0);

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










