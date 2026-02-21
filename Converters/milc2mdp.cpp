#include <cstdio>
#include <cmath>
#include <cstring>
#include <complex>
#include <ctime>
#include <cstdint>
#include <memory>
#include <fstream>

using Complex = std::complex<float>;

#ifdef SHORT32
using type32 = short;
#else
using type32 = int;
#endif

void _error(const char s[])
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

/**********************************************************************/
/* Binary lattice formats                                             */
/**********************************************************************/
/* In version 5 we have two binary lattice file formats:

      serial files     Data written in coordinate natural order
      checkpoint files Data written in node dump order

      Further descriptive information is kept in a separate ASCII header
      file.  See below.

      */

/*--------------------------------------------------------------------*/
/* version 5 binary file format */

#define GAUGE_VERSION_NUMBER 20103
#define MAX_TIME_STAMP 64

/* 1. Header comes first    */

struct gauge_header_t
{
  uint32_t magic_number; /* Identifies file format */
  type32 dims[4];        /* Full lattice dimensions */
  char time_stamp[MAX_TIME_STAMP];
  /* Date and time stamp - used to
check consistency between the
ASCII header file and the
lattice file */
  type32 header_bytes; /* NOT WRITTEN TO THE FILE but
helpful for finding the data */
  type32 what_is_this;
  type32 order; /* 0 means no coordinate list is
attached and the values are in
coordinate serial order.
Nonzero means that a
coordinate list is attached,
specifying the order of values */
  type32 boo;
};
using gauge_header = struct gauge_header_t;

class _generic_field_file_header
{
public:
  char file_id[60];
  char program_version[60];
  char creation_date[60];
  uint32_t endianess;
  int32_t ndim;
  int32_t box[10];
  int32_t bytes_per_site;
  int32_t sites;

  _generic_field_file_header()
  {
    strcpy(file_id, "File Type: MDP FIELD\n");
  }
};

_generic_field_file_header get_info(const std::string &filename)
{
  _generic_field_file_header myheader;
  std::ifstream in(filename, std::ios::binary);
  if (!in)
    _error("Unable to open file");

  in.read(reinterpret_cast<char *>(&myheader), sizeof(_generic_field_file_header));

  if (!in)
    _error("Error while reading file");

  return myheader;
}

class short_field
{
public:
  std::unique_ptr<Complex[]> m_data;
  long size;
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

void write_header(std::ofstream &MDP_fp,
                  int bytes_per_site,
                  const char *program_version,
                  int nx[4],
                  const char *creation_date = nullptr)
{
  _generic_field_file_header myheader{};

  myheader.ndim = 4;

  for (int i = 0; i < 4; ++i)
    myheader.box[i] = nx[i];

  for (int i = 4; i < 10; ++i)
    myheader.box[i] = 0;

  myheader.sites = nx[0] * nx[1] * nx[2] * nx[3];

  myheader.bytes_per_site = bytes_per_site;
  myheader.endianess = 0x87654321;

  std::strcpy(myheader.program_version, program_version);

  if (creation_date)
  {
    std::strcpy(myheader.creation_date, creation_date);
  }
  else
  {
    std::time_t t = std::time(nullptr);
    std::strcpy(myheader.creation_date, std::ctime(&t));
  }

  MDP_fp.write(reinterpret_cast<const char *>(&myheader), sizeof(_generic_field_file_header));
}

void process_gauge(std::ifstream &MILC_fp, std::ofstream &MDP_fp, int nx[4], bool switch_endianess)
{
  short_field U;
  U.initialize(nx[1], nx[2], nx[3], 4, 3, 3);
  unsigned int matrix_size = 72;

  for (int x0 = 0; x0 < nx[0]; x0++)
  {
    for (int x3 = 0; x3 < nx[3]; x3++)
      for (int x2 = 0; x2 < nx[2]; x2++)
        for (int x1 = 0; x1 < nx[1]; x1++)
          for (int mu = 0; mu < 4; mu++)
          {
            char *dst = reinterpret_cast<char *>(&U(x1, x2, x3, (mu + 1) % 4, 0, 0));

            MILC_fp.read(dst, matrix_size);

            if (!MILC_fp)
            {
              _error("Error while reading from MILC file");
            }
          }

    if (switch_endianess)
    {
      for (int mu = 0; mu < U.size; mu++)
      {
        switch_endianess_byte4(*((long *)U.m_data.get() + 2 * mu));
        switch_endianess_byte4(*((long *)U.m_data.get() + 2 * mu + 1));
      }
    }

    MDP_fp.write(reinterpret_cast<const char *>(U.m_data.get()),
        static_cast<std::streamsize>(U.size * sizeof(Complex)));

    if (!MDP_fp)
    {
      _error("Error while writing to MDP file");
    }
  }
}

int main(int argc, char **argv)
{
  printf("======================================================\n");
  printf("Program for converting MILC gauge configurations\n");
  printf("and propagators into MDP files\n");
  printf("Conversion of propagators not implemented yet\n");
  printf("======================================================\n");

  if (argc < 4)
  {
    printf("usage:\n\n");
    printf("milc2mdp -gauge 16x08x08x08 input output\n\n");
    printf("milc2mdp -fermi 16x08x08x08 input output\n\n");
    exit(0);
  }

  int nx[4];
  sscanf(argv[2], "%ix%ix%ix%i", nx, nx + 1, nx + 2, nx + 3);
  long time0 = clock() / CLOCKS_PER_SEC;

  std::ifstream MILC_fp(argv[3]);
  if (!MILC_fp)
  {
    _error("Cannot open input file");
  }

  std::ofstream MDP_fp(argv[4], std::ios::binary);
  if (!MDP_fp)
  {
    _error("Cannot open output file");
  }

  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("opening the MILC file: %s (read)\n", argv[3]);
  printf("opening the MDP file: %s (write) \n", argv[4]);

  int offset = sizeof(_generic_field_file_header) / sizeof(char);
  // gauge: bytes_per_site = 4(mu)*9(SU3 matrix)*2(cplx)*4(bytes_per_float)
  //                       = 288
  // fermi: bytes_per_site = 16(spin-sq)*9(color-sq)*2(cplx)*4(bytes_per_float)
  //                       = 1152
  int32_t bytes_per_site;
  if (strcmp(argv[3], "-gauge") == 0)
    bytes_per_site = 288;
  else
    bytes_per_site = 1152;

  if (strcmp(argv[1], "-gauge") == 0)
  {
    gauge_header MILC_header;
    MILC_fp.read(reinterpret_cast<char *>(&MILC_header), sizeof(gauge_header));

    if (!MILC_fp)
    {
      _error("Error while reading from file");
    }

    printf("%x\n", MILC_header.magic_number);

    bool switch_endianess = false;
    if (MILC_header.magic_number == 0x874e0000)
    {
      printf("switching endianess\n");
      switch_endianess_byte4(MILC_header.dims[0]);
      switch_endianess_byte4(MILC_header.dims[1]);
      switch_endianess_byte4(MILC_header.dims[2]);
      switch_endianess_byte4(MILC_header.dims[3]);
      switch_endianess = true;
    }

    if ((MILC_header.dims[0] != nx[1]) &&
        (MILC_header.dims[1] != nx[2]) &&
        (MILC_header.dims[2] != nx[3]) &&
        (MILC_header.dims[3] != nx[0]))
    {
      printf("%xx%xx%xx%x => Incorrect dimensions!\n",
             MILC_header.dims[3],
             MILC_header.dims[0],
             MILC_header.dims[1],
             MILC_header.dims[2]);
      exit(1);
    }

    write_header(MDP_fp, bytes_per_site, "Converted from MILC", nx, MILC_header.time_stamp);
    process_gauge(MILC_fp, MDP_fp, nx, switch_endianess);
  }
  else if (strcmp(argv[1], "-quark") == 0)
  {
    printf("I am sorry, but conversion of propagators is not implemented yet!\n");
  }

  printf("Reading back the MDP file named %s for sanity check...\n", argv[4]);
  _generic_field_file_header header = get_info(argv[4]);

  printf("Header size is %i bytes.\n", offset);
  printf("Endianess = %x\n", header.endianess);
  printf("Dimensions = %i\n", header.ndim);
  printf("Lattice: %i x %i x %i x %i\n", header.box[0], header.box[1], header.box[2], header.box[3]);
  printf("Total sites = %i\n", header.sites);
  printf("Bytes per site = %i\n", header.bytes_per_site);
  printf("Output file size is %i bytes.\n", bytes_per_site * header.box[0] * header.box[1] * header.box[2] * header.box[3] + offset);
  printf("Output file name is: %s\n", argv[4]);
  printf("\nAll sites are OK.\n");
  printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);

  return 0;
}
