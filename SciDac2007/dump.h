#ifndef DUMP_
#define DUMP_

#include <fstream>
#include <iomanip>
#include <filesystem>

namespace MDP
{
  /**
   * @brief Export a scalar component of an MDP field to a VTK file.
   *
   * This function writes data from an mdp_real_field defined on a 3D lattice
   * to a VTK file in the legacy STRUCTURED_POINTS format. The output can be
   * written in either ASCII or binary form and is suitable for visualization
   * in tools such as ParaView or VisIt.
   *
   * The lattice is assumed to be three-dimensional and uniformly spaced
   * with unit spacing in all directions. Each lattice site contributes one
   * scalar value to the output dataset.
   *
   * @param s
   *   Reference to the field defined on an MDP lattice. The field may contain
   *   multiple components per site.
   *
   * @param component
   *   Index of the component to be exported at each lattice site.
   *   For example, for a vector field (x, y, z), this selects which component
   *   (e.g., 0 = x, 1 = y, 2 = z) is written to the file.
   *
   * @param filename
   *   Name of the output VTK file. The file is first written to a temporary
   *   file ("tmp.vtk") and then renamed to the final filename to reduce the
   *   risk of incomplete output.
   *
   * @param ASCII
   *   If true, the file is written in ASCII format. If false, binary output
   *   is used (with endian conversion applied as required by the VTK format).
   *
   * @throws std::ios_base::failure
   *   Thrown if the temporary file cannot be opened for writing.
   *
   * @note
   *   - The function assumes the lattice has at least three dimensions.
   *   - The output uses VTK legacy format (version 2.0).
   *   - In binary mode, data is written after endian conversion to match
   *     VTK's big-endian requirement.
   *   - Existing output files with the same name are overwritten.
   *
   * @warning
   *   No explicit validation is performed on lattice dimensionality or
   *   component index bounds. Passing invalid values may result in
   *   undefined behavior.
   */
  void dump(mdp_real_field &s,
            mdp_int component = 0,
            const std::string &filename = "default.vtk",
            bool ASCII = true)
  {
    namespace fs = std::filesystem;

    const std::string tempfile = "tmp.vtk";

    fs::remove(tempfile); // remove tmp.vtk file if exists

    std::cout << "Saving file " << filename << std::endl;

    std::ofstream ofs;
    if (ASCII)
      ofs.open(tempfile, std::ios::out);
    else
      ofs.open(tempfile, std::ios::out | std::ios::binary);

    if (!ofs)
      throw std::ios_base::failure("Unable to open temporary VTK file for writing");

    const mdp_uint LX = s.lattice().size(0);
    const mdp_uint LY = s.lattice().size(1);
    const mdp_uint LZ = s.lattice().size(2);

    // VTK header
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << (ASCII ? "ASCII" : "BINARY") << "\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS " << LX << " " << LY << " " << LZ << "\n"
        << "ORIGIN 0 0 0\n"
        << "SPACING 1 1 1\n"
        << "POINT_DATA " << LX * LY * LZ << "\n"
#ifdef USE_SINGLE_PRECISION
        << "SCALARS scalar float 1\n"
#else
        << "SCALARS scalar double 1\n"
#endif
        << "LOOKUP_TABLE default\n";

    mdp_site p(s.lattice());

    for (mdp_uint k = 0; k < LZ; k++)
    {
      for (mdp_uint j = 0; j < LY; j++)
      {
        for (mdp_uint i = 0; i < LX; i++)
        {
          p.set(i, j, k);

          mdp_real fval = static_cast<mdp_real>(s(p, component));
          if (ASCII)
          {
            ofs << std::scientific << fval << "\n";
          }
          else
          {
            switch_endianess(fval);
            ofs.write(reinterpret_cast<char *>(&fval), sizeof(mdp_real));
            if (!ofs)
              error("probably out of disk space");
          }
        }
      }
    }

    ofs.close();
    fs::remove(filename);
    fs::rename(tempfile, filename);
  }
} // namespace MDP

#endif /* DUMP_ */
