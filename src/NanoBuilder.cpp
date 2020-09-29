//
// Created by Ilya Bryukhanov on 23.10.2017.
//

#include "NanoBuilder.h"

namespace NTBuilder {

    atoms BuildCellBox(LatticeType type, double a_lattice,
                       r_crystal_dir x, r_crystal_dir y, r_crystal_dir z, int px, int py, int pz)
    {
        if (!CheckPairwiseOrthogonality(x, y, z)) {
            throw std::logic_error("Chosen directions are not pairwise orthogonal\n");
        }
        if (!CheckRightHanded(x, y, z)) {
            throw std::logic_error("Chosen directions are not right\n");
        }
        atoms square_plane_atoms;
        bool special_orient_x = checkSpecialOrientation(x);
        bool special_orient_y = checkSpecialOrientation(y);
        bool special_orient_z = checkSpecialOrientation(z);
        double cell_lx = getLength(x);
        double cell_ly = getLength(y);
        double cell_lz = getLength(z);
        if (special_orient_x) cell_lx = cell_lx/2;
        if (special_orient_y) cell_ly = cell_ly/2;
        if (special_orient_z) cell_lz = cell_lz/2;
        double length_x = a_lattice * lattice_vectors.at(type)[0][0] * cell_lx * px;
        double length_y = a_lattice * lattice_vectors.at(type)[1][1] * cell_ly * py;
        double length_z = a_lattice * lattice_vectors.at(type)[2][2] * cell_lz * pz;


        matrix unit_rows = Norm3VectorsRows(x, y, z);
        matrix unit_cols = Norm3VectorsCols(x, y, z);

        matrix box_corners = {{-0.001, -0.001, -0.001},
                              {length_x - 0.001, length_y - 0.001, length_z - 0.001}
        };

        double xmin,ymin,zmin,xmax,ymax,zmax;
        xmin = ymin = zmin = BIG;
        xmax = ymax = zmax = -BIG;

        // determine lattice max and min indexes corresponded to box corners
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    BLboxNoRotation({box_corners[i][0],box_corners[j][1],box_corners[k][2]},
                            unit_cols, type, a_lattice,
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
        int atom_type;
        for (int k = klo; k <= khi; k++) {
            for (int j = jlo; j <= jhi; j++) {
                for (int i = ilo; i <= ihi; i++) {
                    for (unsigned m = 0; m < atom_basis.at(type).size(); m++) {
                        atom_lattice[0] = i + atom_basis.at(type)[m][0];
                        atom_lattice[1] = j + atom_basis.at(type)[m][1];
                        atom_lattice[2] = k + atom_basis.at(type)[m][2];
                        // convert from lattice coords to box coords
                        atom_box = LatticeToBoxNoRotation(atom_lattice, type, a_lattice, unit_rows);
                        atom_type = atom_types.at(type)[m];
                        // if atom_box does not lie in our region then continue
                        if (atom_box[0] < box_corners[0][0] || atom_box[0] >= box_corners[1][0] ||
                            atom_box[1] < box_corners[0][1] || atom_box[1] >= box_corners[1][1] ||
                            atom_box[2] < box_corners[0][2] || atom_box[2] >= box_corners[1][2]) continue;

                        // if it lies than push it back to output array
                        square_plane_atoms.push_back(Atom(atom_box,atom_type));
                    }
                }
            }
        }
        std::sort(square_plane_atoms.begin(), square_plane_atoms.end(),
                [](const Atom& a, const Atom& b) {
                    return a.x[2] < b.x[2];
        });

        return square_plane_atoms;
    }

    atoms BuildYperBox(LatticeType type, double a_lattice,
                       r_crystal_dir x, r_crystal_dir y, r_crystal_dir z, double length_x, int py, double length_z)
    {
        if (!CheckPairwiseOrthogonality(x, y, z)) {
            throw std::logic_error("Chosen directions are not pairwise orthogonal\n");
        }
        if (!CheckRightHanded(x, y, z)) {
            throw std::logic_error("Chosen directions are not right\n");
        }
        atoms square_plane_atoms;
        double cell_ly = getLength(y) *lattice_vectors.at(type)[1][1];
        double length_y = a_lattice * cell_ly * py;
        std::cout << "Number of cells " << py << '\n';

        matrix unit_rows = Norm3VectorsRows(x, y, z);
        matrix unit_cols = Norm3VectorsCols(x, y, z);

        matrix box_corners = {{-0.001, -0.001, -0.001},
                              {length_x - 0.001, length_y - 0.001, length_z - 0.001}
        };

        double xmin,ymin,zmin,xmax,ymax,zmax;
        xmin = ymin = zmin = BIG;
        xmax = ymax = zmax = -BIG;

        // determine lattice max and min indexes corresponded to box corners
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    BLboxNoRotation({box_corners[i][0],box_corners[j][1],box_corners[k][2]},
                                    unit_cols, type, a_lattice,
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
        int atom_type;
        for (int k = klo; k <= khi; k++) {
            for (int j = jlo; j <= jhi; j++) {
                for (int i = ilo; i <= ihi; i++) {
                    for (unsigned m = 0; m < atom_basis.at(type).size(); m++) {
                        atom_lattice[0] = i + atom_basis.at(type)[m][0];
                        atom_lattice[1] = j + atom_basis.at(type)[m][1];
                        atom_lattice[2] = k + atom_basis.at(type)[m][2];
                        atom_type = atom_types.at(type)[m];
                        // convert from lattice coords to box coords
                        atom_box = LatticeToBoxNoRotation(atom_lattice, type, a_lattice, unit_rows);

                        // if atom_box does not lie in our region then continue
                        if (atom_box[0] < box_corners[0][0] || atom_box[0] >= box_corners[1][0] ||
                            atom_box[1] < box_corners[0][1] || atom_box[1] >= box_corners[1][1] ||
                            atom_box[2] < box_corners[0][2] || atom_box[2] >= box_corners[1][2]) continue;

                        // if it lies than push it back to output array
                        square_plane_atoms.push_back(Atom(atom_box,atom_type));
                    }
                }
            }
        }
        std::sort(square_plane_atoms.begin(), square_plane_atoms.end(),
                  [](const Atom& a, const Atom& b) {
                      return a.x[2] < b.x[2];
                  });

        return square_plane_atoms;
    }

    atoms BuildDeformPyramid(
            LatticeType type, double a,
            r_crystal_dir x, r_crystal_dir y, r_crystal_dir z,
            double length_x, int py, double def, double length_z)
    {

        double ly_box = a*getLength(y)*lattice_vectors.at(type)[1][1];

        double length_y = py * ly_box;

        atoms pyr_atoms = BuildYperBox(type, a, x, y, z, length_x, py, length_z);
        std::ofstream out("box.data");
        NTBuilder::OutputLammpsFormat(out, pyr_atoms, 2, length_x+40, py * ly_box, length_z+40);
        for (unsigned i = 0; i < pyr_atoms.size(); ++i) {
            double z = pyr_atoms[i].x[2];
            double y_scale = -2*def* z / length_z + (length_y+def);
            pyr_atoms[i].x[1] *= y_scale/length_y;
        }

        return pyr_atoms;
    }

    atoms BuildNanotube(LatticeType type, double a,
                        r_crystal_dir x, r_crystal_dir y, r_crystal_dir z,
                        double length, double r_in, double r_out)
    {
        double ly_box = a*getLength(y)*lattice_vectors.at(type)[1][1];
        double r_mid = (r_in + r_out) * 0.5;
        double length_y = r_mid * 2 * PI;
        int py = static_cast<int>(length_y/ly_box);
        length_y = ly_box*py;

        double r_in_new = length_y/(2*PI) - 0.5*(r_out - r_in);
        double r_out_new = length_y/(2*PI) + 0.5*(r_out - r_in);
        // first create a plane
        atoms nanotube = BuildDeformPyramid(type, a, x, y, z, length, py, PI*(r_out - r_in), r_out-r_in);
        //std::ofstream out("def_box.data");
        //NTBuilder::OutputLammpsFormat(out, nanotube, 2);
        // fold a nanotube
        for (unsigned i = 0; i < nanotube.size(); ++i) {
            double z = nanotube[i].x[2];
            double r_cur = r_out_new - z;
            double y = nanotube[i].x[1];
            double a_rad = y / (r_cur);
            double y_new = -r_cur*sin(a_rad);
            double z_new = r_cur*(1.0-cos(a_rad));
            nanotube[i].x[1] = y_new;
            nanotube[i].x[2] = nanotube[i].x[2] + z_new;
        }
        std::cout << "Nanotube with r_in " << r_in_new << " and r_out " << r_out_new << " has been created\n";
        std::cout << "Number of atoms " << nanotube.size() << '\n';
        return nanotube;
    }

    void OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms, int n_types) {
        using namespace std::chrono;
        // 0, 1: min, max values of x, y, z
        double boundary[2][3];
        for (int i = 0; i < 3; ++i) {
            boundary[0][i] = BIG;
            boundary[1][i] = -BIG;
        }
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                boundary[0][j] = std::min(body_atoms[i].x[j], boundary[0][j]);
                boundary[1][j] = std::max(body_atoms[i].x[j], boundary[1][j]);
            }
        }
        double gap = 40.0;
        for (int i = 0; i < 3; ++i) {
            boundary[0][i] -= gap;
            boundary[1][i] += gap;
        }
        // Declaring argument for time()
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        out << "LAMMPS data file written by Nanotube Builder at " << std::ctime(&end_time) << '\n';
        out << body_atoms.size() << " atoms" << '\n';
        out << " " << n_types <<  " atom types " << '\n';
        out << " " << boundary[0][0] << ' ' << boundary[1][0] << " xlo xhi " << '\n';
        out << " " << boundary[0][1] << ' ' << boundary[1][1] << " ylo yhi " << '\n';
        out << " " << boundary[0][2] << ' ' << boundary[1][2] << " zlo zhi " << '\n';

        out << "\nAtoms \n\n";
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            out << i+1 << " " << body_atoms[i].type << " "
                << body_atoms[i].x[0] << ' '
                << body_atoms[i].x[1] << ' '
                << body_atoms[i].x[2] << '\n';
        }
    }

    void OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms, int n_types, double lx, double ly, double lz) {
        using namespace std::chrono;
        // 0, 1: min, max values of x, y, z

        // Declaring argument for time()
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        out << "LAMMPS data file written by Nanotube Builder at " << std::ctime(&end_time) << '\n';
        out << body_atoms.size() << " atoms" << '\n';
        out << " " << n_types <<  " atom types " << '\n';
        out << " " << 0.0 << ' ' << lx << " xlo xhi " << '\n';
        out << " " << 0.0 << ' ' << ly << " ylo yhi " << '\n';
        out << " " << 0.0 << ' ' << lz << " zlo zhi " << '\n';

        out << "\nAtoms \n\n";
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            out << i+1 << " " << body_atoms[i].type << " "
                << body_atoms[i].x[0] << ' '
                << body_atoms[i].x[1] << ' '
                << body_atoms[i].x[2] << '\n';
        }
    }

    void OutputLammpsFormat(std::ofstream& out, const atoms& body_atoms, int n_types, double lx) {
        using namespace std::chrono;
        // 0, 1: min, max values of x, y, z
        double boundary[2][3];
        for (int i = 0; i < 3; ++i) {
            boundary[0][i] = BIG;
            boundary[1][i] = -BIG;
        }
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                boundary[0][j] = std::min(body_atoms[i].x[j], boundary[0][j]);
                boundary[1][j] = std::max(body_atoms[i].x[j], boundary[1][j]);
            }
        }
        double gap = 40.0;
        for (int i = 0; i < 3; ++i) {
            boundary[0][i] -= gap;
            boundary[1][i] += gap;
        }
        // Declaring argument for time()
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        out << "LAMMPS data file written by Nanotube Builder at " << std::ctime(&end_time) << '\n';
        out << body_atoms.size() << " atoms" << '\n';
        out << " " << n_types <<  " atom types " << '\n';
        out << " " << 0.0 << ' ' << lx << " xlo xhi " << '\n';
        out << " " << boundary[0][1] << ' ' << boundary[1][1] << " ylo yhi " << '\n';
        out << " " << boundary[0][2] << ' ' << boundary[1][2] << " zlo zhi " << '\n';

        out << "\nAtoms \n\n";
        for (unsigned i = 0; i < body_atoms.size(); ++i) {
            out << i+1 << " " << body_atoms[i].type << " "
                << body_atoms[i].x[0] << ' '
                << body_atoms[i].x[1] << ' '
                << body_atoms[i].x[2] << '\n';
        }
    }

};










