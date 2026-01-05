#include "../Libraries/fermiqcd.h"

using namespace MDP;

void write_gauge_configuration(std::ofstream &outfp, gauge_field &U, int *L)
{
  mdp_site x(U.lattice());

  for (int x0 = 0; x0 < L[0]; x0++)
    for (int x3 = 0; x3 < L[3]; x3++)
      for (int x2 = 0; x2 < L[2]; x2++)
        for (int x1 = 0; x1 < L[1]; x1++)
        {
          x.set(x0, x1, x2, x3);

          for (int mu = 1; mu <= 4; mu++)
          {
            for (int i = 0; i < 3; i++)
              for (int j = 0; j < 3; j++)
              {
                outfp << real(U(x, mu % 4, i, j)) << " " << imag(U(x, mu % 4, i, j)) << "\n";
              }
          }
        }
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  if (argc < 3)
  {
    mdp << "Usage:\n\t" << argv[0] << " input output\n\n";
    exit(1);
  }

  if (!file_exists(argv[1]))
  {
    mdp << "Unable to open gauge configuration: " << argv[1] << "\n";
    exit(1);
  }

  std::ofstream outfp(argv[2]);
  if (!outfp)
  {
    mdp << "Unable to open output file: " << argv[2] << "\n";
    exit(1);
  }

  // read the lattice size from the gauge configuration itself
  mdp_field_file_header header = get_info(argv[1]);
  int *L = header.box;

  mdp << "Lattice size: "
      << L[0] << "x" << L[1] << "x"
      << L[2] << "x" << L[3] << "\n";

  constexpr int nc = 3; // number of colors
  mdp_lattice lattice(4, L);
  mdp_site x(lattice);

  assert(header.ndim == 4);
  assert(sizeof(mdp_complex) == 2 * sizeof(mdp_real));

  if (size_t(header.bytes_per_site) == header.ndim * nc * nc * sizeof(mdp_complex))
  {
    gauge_field U(lattice, nc);
    U.load(argv[1]);
    write_gauge_configuration(outfp, U, L);
  }

  mdp.close_wormholes();
  return 0;
}
