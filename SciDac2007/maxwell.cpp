#include "mdp.h"
#include "mdp_all.h"
#define X 0
#define Y 1
#define Z 2

void dump(mdp_field<int> &s, std::string filename = "default.vtk")
{
  char header[1024], number[1024];
  int LX = s.lattice().size(0), LY = s.lattice().size(1), LZ = s.lattice().size(2);
  snprintf(header, 1024,
           "# vtk DataFile Version 2.0\n"
           "Really cool data\n"
           "ASCII\n"
           "DATASET STRUCTURED_POINTS\n"
           "DIMENSIONS %i %i %i\n"
           "ORIGIN     0   0   0\n"
           "SPACING    1   1   1\n"
           "POINT_DATA %i\n"
           "SCALARS scalars0 float 1\n"
           "LOOKUP_TABLE default\n\n",
           LX, LY, LZ, LX * LY * LZ);

  int sfd = open("tmp.vtk", O_WRONLY);
  std::cout << "saving... " << filename << " as fd=" << sfd << std::endl;
  setFileLock(sfd);
  if (write(sfd, header, strlen(header)) == -1)
  {
    error("Error writing to file");
  }

  mdp_site p(s.lattice());
  for (int i = 0; i < LX; i++)
    for (int j = 0; j < LY; j++)
      for (int k = 0; k < LZ; k++)
      {
        p.set(i, j, k);
        snprintf(number, 1024, "%i\n", (int)s(p));
        if (write(sfd, number, strlen(number)) == -1)
        {
          error("Error writing to file");
        }
      }
  setFileUnlock(sfd);
  if (system((std::string("cp tmp.vtk ") + filename).c_str()) == -1)
  {
    error("Error copying file");
  }
  close(sfd);
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int L[] = {20, 20, 20};
  mdp_lattice cube(3, L);
  mdp_site x(cube);
  mdp_complex_field E(cube, 3);
  mdp_complex_field B(cube, 3);
  mdp_complex_field q(cube);
  mdp_complex_field j(cube, 3);

  forallsites(x)
  {
    E(x) = 0;
    B(x) = 0;
  }

  mdp.close_wormholes();
  return 0;
}
