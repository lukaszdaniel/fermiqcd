#ifndef DUMP_
#define DUMP_

#include <fstream>
#include <iomanip>
#include <filesystem>

namespace MDP
{
  void dump(mdp_field<mdp_real> &s,
            int site_idx = 0,
            const std::string &filename = "default.vtk",
            bool bASCII = true)
  {
    namespace fs = std::filesystem;

    const std::string tempfile = "tmp.vtk";

    fs::remove(tempfile); // remove tmp.vtk file if exists

    std::cout << "Saving file " << filename << std::endl;

    std::ofstream ofs;
    if (bASCII)
      ofs.open(tempfile, std::ios::out);
    else
      ofs.open(tempfile, std::ios::out | std::ios::binary);

    if (!ofs)
      throw std::ios_base::failure("Unable to open temporary VTK file for writing");

    const int LX = s.lattice().size(0);
    const int LY = s.lattice().size(1);
    const int LZ = s.lattice().size(2);

    // VTK header
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << (bASCII ? "ASCII" : "BINARY") << "\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS " << LX + 1 << " " << LY + 1 << " " << LZ + 1 << "\n"
        << "ORIGIN 0 0 0\n"
        << "SPACING 1 1 1\n"
        << "POINT_DATA " << (LX + 1) * (LY + 1) * (LZ + 1) << "\n"
        << "SCALARS scalar0 float 1\n"
        << "LOOKUP_TABLE default\n";

    mdp_site p(s.lattice());

    for (int k = 0; k < LZ + 1; k++)
      for (int j = 0; j < LY + 1; j++)
        for (int i = 0; i < LX + 1; i++)
        {
          p.set(i % LX, j % LY, k % LZ);

          float fval = static_cast<float>(s(p, site_idx));
          if (bASCII)
          {
            ofs << std::scientific << fval << "\n";
          }
          else
          {
            ofs.write(reinterpret_cast<char *>(&fval), sizeof(float));
          }
        }

    ofs.close();

    fs::remove(filename);
    fs::rename(tempfile, filename);
  }
} // namespace MDP

#endif /* DUMP_ */
