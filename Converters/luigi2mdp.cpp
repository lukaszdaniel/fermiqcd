#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <complex>
#include <memory>

#ifndef USE_DOUBLE_PRECISION
using Complex = std::complex<float>;
#else
using Complex = std::complex<double>;
#endif

void error(const char s[])
{
  printf("ERROR: %s\n", s);
  exit(1);
}

template <class T>
void switch_endianess_byte4(T &a)
{
  char *p = (char *)&a;
  static char q[4];
  if (sizeof(T) == 4)
  {
    q[0] = p[0];
    q[1] = p[1];
    q[2] = p[2];
    q[3] = p[3];
    p[0] = q[3];
    p[1] = q[2];
    p[2] = q[1];
    p[3] = q[0];
  }
  else
    printf("error endianess\n");
}

class _generic_field_file_header
{
public:
  char file_id[60];
  char program_version[60];
  char creation_date[60];
  unsigned int endianess;
  int ndim;
  int box_size[10];
  long bytes_per_site;
  long sites;

  _generic_field_file_header()
  {
    strcpy(file_id, "File Type: MDP FIELD\n");
  }
};

class short_field
{
public:
  std::unique_ptr<Complex[]> m_data;
  int size;
  int dim[7];

  short_field() : m_data(nullptr)
  {
  }

  ~short_field()
  {
  }

  void initialize(int x1, int x2, int x3, int a = 1, int b = 1, int c = 1,
                  int d = 1)
  {
    size = x1 * x2 * x3 * a * b * c * d;
    dim[0] = x1;
    dim[1] = x2;
    dim[2] = x3;
    dim[3] = a;
    dim[4] = b;
    dim[5] = c;
    dim[6] = d;
    m_data = std::make_unique<Complex[]>(size);
  }

  Complex &operator()(int x1, int x2, int x3,
                      int a = 0, int b = 0, int c = 0, int d = 0)
  {
    return m_data[(((((x1 * dim[1] + x2) * dim[2] + x3) * dim[3] + a) * dim[4] + b) * dim[5] + c) * dim[6] + d];
  }
};

void write_header(FILE *MDP_fp, int bytes_per_site, const char *program_version, int nx[4])
{
  _generic_field_file_header myheader;
  myheader.ndim = 4;
  for (int ii = 0; ii < 4; ii++)
    myheader.box_size[ii] = nx[ii];
  for (int ii = 4; ii < 10; ii++)
    myheader.box_size[ii] = 0;
  myheader.sites = nx[0] * nx[1] * nx[2] * nx[3];
  myheader.bytes_per_site = bytes_per_site;
  myheader.endianess = 0x87654321;
  strcpy(myheader.program_version, program_version);
  time_t time_and_date;
  time(&time_and_date);
  strcpy(myheader.creation_date, ctime(&time_and_date));
  int offset = sizeof(_generic_field_file_header) / sizeof(char);
  fwrite(&myheader, sizeof(char), offset, MDP_fp);
}

void process_gauge(FILE *LUIGI_fp, FILE *MDP_fp, int nx[4], int nc)
{
  short_field U;
  U.initialize(nx[1], nx[2], nx[3], 4, nc, nc);
  unsigned int matrix_size = nc * nc * sizeof(Complex);

  for (int x0 = 0; x0 < nx[0]; x0++)
  {
    for (int x3 = 0; x3 < nx[3]; x3++)
      for (int x2 = 0; x2 < nx[2]; x2++)
        for (int x1 = 0; x1 < nx[1]; x1++)
          for (int mu = 0; mu < 4; mu++)
            if (fread(&U(x1, x2, x3, (4 - mu) % 4, 0, 0), matrix_size, 1, LUIGI_fp) != matrix_size)
            {
              error("Error while reading from file");
            }
    fwrite(U.m_data.get(), U.size, sizeof(Complex), MDP_fp);
  }
}

int main(int argc, char **argv)
{
  printf("======================================================\n");
  printf("Program for converting LUIGI gauge configurations\n");
  printf("and propagators into MDP files\n");
  printf("Conversion of propagators not implemented yet\n");
  printf("======================================================\n");

  if (argc < 4)
  {
    printf("usage:\n\n");
    printf("luigi2mdp -gauge 16x08x08x08,NC input output\n\n");
    exit(0);
  }

  int nx[4];
  int nc;
  sscanf(argv[2], "%ix%ix%ix%i,%i", nx, nx + 1, nx + 2, nx + 3, &nc);
  long time0 = clock() / CLOCKS_PER_SEC;

  FILE *LUIGI_fp = fopen(argv[3], "r");
  if (!LUIGI_fp)
  {
    error("Cannot open input file");
  }

  FILE *MDP_fp = fopen(argv[4], "w");
  if (!MDP_fp)
  {
    fclose(LUIGI_fp);
    error("Cannot open output file");
  }

  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("opening the LUIGI file: %s (read)\n", argv[3]);
  printf("opening the MDP file: %s (write) \n", argv[4]);

  int offset = sizeof(_generic_field_file_header) / sizeof(char);
  // gauge: bytes_per_site = 4(mu)*9(SU3 matrix)*2(cplx)*8(bytes_per_double)
  //                       = 576
  // quark: bytes_per_site = 4(spin)*3(color)*2(cplx)*8(bytes_per_double)
  //                       = 192
  long bytes_per_site;
  if (strcmp(argv[1], "-gauge") == 0)
    bytes_per_site = 576;
  else
    bytes_per_site = 192;

  write_header(MDP_fp, bytes_per_site, "Converted from LUIGI", nx);
  if (strcmp(argv[1], "-gauge") == 0)
  {
    process_gauge(LUIGI_fp, MDP_fp, nx, nc);
  }
  else if (strcmp(argv[1], "-quark") == 0)
  {
    printf("I am sorry, but conversion of propagators is not implemented yet!\n");
  }

  fclose(LUIGI_fp);
  fclose(MDP_fp);

  printf("\nAll sites are OK.\n");
  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("Header size is %i bytes.\n", offset);
  printf("Output file size is %li bytes.\n", bytes_per_site * nx[0] * nx[1] * nx[2] * nx[3] + offset);
  printf("Output file name is: %s\n", argv[4]);
  printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);

  return 0;
}
