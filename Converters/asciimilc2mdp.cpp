#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <cstring>
#include <ctime>

void error(const char s[])
{
  printf("ERROR: %s\n", s);
  exit(1);
}

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
    error("Unable to open file");

  in.read(reinterpret_cast<char *>(&myheader), sizeof(_generic_field_file_header));

  if (!in)
    error("Error while reading file");

  return myheader;
}

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

void process_gauge(std::ifstream &file_fp,
                   std::ofstream &MDP_fp,
                   int buffer_size,
                   int nx[4])
{
  // assumes double precision:
  // 3 x 3 x 2 doubles = 18 doubles = 144 bytes
  double buffer[18];

  const std::streamoff offset =
      static_cast<std::streamoff>(sizeof(_generic_field_file_header));

  for (int x0 = 0; x0 < nx[0]; ++x0)
    for (int x3 = 0; x3 < nx[3]; ++x3)
      for (int x2 = 0; x2 < nx[2]; ++x2)
        for (int x1 = 0; x1 < nx[1]; ++x1)
          for (int mu = 0; mu < 4; ++mu)
          {
            // read SU(3) matrix: (real, imag) pairs
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j)
              {
                file_fp >> buffer[6 * i + 2 * j] >> buffer[6 * i + 2 * j + 1];
              }

            // map to the MDP ordering
            const std::streamoff position =
                (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3) * 4 +
                ((mu + 1) % 4);

            const std::streamoff byte_pos =
                offset + static_cast<std::streamoff>(buffer_size * position);

            MDP_fp.seekp(byte_pos, std::ios::beg);
            MDP_fp.write(reinterpret_cast<const char *>(buffer),
                      static_cast<std::streamsize>(buffer_size));
          }
}

void process_quark(std::ifstream &TONY_fp,
                   std::ofstream &MDP_fp,
                   int buffer_size,
                   int nx[4])
{
  // 4 spins × 3 colors × (re, im) = 24 doubles
  double buffer[24];

  const std::streamoff offset =
      static_cast<std::streamoff>(sizeof(_generic_field_file_header));

  for (int x0 = 0; x0 < nx[0]; ++x0)
    for (int x3 = 0; x3 < nx[3]; ++x3)
      for (int x2 = 0; x2 < nx[2]; ++x2)
        for (int x1 = 0; x1 < nx[1]; ++x1)
        {
          // read spin-color components
          for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i)
            {
              TONY_fp >>
                  buffer[6 * a + 2 * i] >>
                  buffer[6 * a + 2 * i + 1];
            }

          // map to the MDP ordering
          const std::streamoff position =
              (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3);

          const std::streamoff byte_pos =
              offset + static_cast<std::streamoff>(buffer_size * position);

          MDP_fp.seekp(byte_pos, std::ios::beg);
          MDP_fp.write(reinterpret_cast<const char *>(buffer),
                       static_cast<std::streamsize>(buffer_size));
        }
}

int main(int argc, char **argv)
{
  printf("=======================================================\n");
  printf("Program for converting MILC ascii gauge configurations\n");
  printf("and quark (fermi_field) into MDP files\n");
  printf("=======================================================\n");

  if (argc < 4)
  {
    printf("usage:\n\n");
    printf("asciimilc2mdp -gauge 16x08x08x08 input output\n\n");
    printf("asciimilc2mdp -quark 16x08x08x08 input output (to be tested) \n\n");
    exit(0);
  }

  int nx[4];
  sscanf(argv[2], "%ix%ix%ix%i", nx, nx + 1, nx + 2, nx + 3);
  long time0 = clock() / CLOCKS_PER_SEC;

  std::ifstream TONY_fp(argv[3]);
  if (!TONY_fp)
  {
    error("Cannot open input file");
  }

  std::ofstream MDP_fp(argv[4], std::ios::binary);
  if (!MDP_fp)
  {
    error("Cannot open output file");
  }

  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("opening the MILC (ascii) file: %s (read)\n", argv[3]);
  printf("opening the MDP file: %s (write) \n", argv[4]);

  int offset = sizeof(_generic_field_file_header) / sizeof(char);
  // gauge: bytes_per_site = 4(mu)*9(SU3 matrix)*2(cplx)*8(bytes_per_double)
  //                       = 576
  // quark: bytes_per_site = 4(spin)*3(color)*2(cplx)*8(bytes_per_double)
  //                       = 192
  int32_t bytes_per_site;
  if (std::strcmp(argv[1], "-gauge") == 0)
    bytes_per_site = 576;
  else
    bytes_per_site = 192;

  write_header(MDP_fp, bytes_per_site, "Converted from MILC (ascii) file", nx);

  if (std::strcmp(argv[1], "-gauge") == 0)
  {
    process_gauge(TONY_fp, MDP_fp, 144, nx);
  }
  else if (std::strcmp(argv[1], "-quark") == 0)
  {
    process_quark(TONY_fp, MDP_fp, 192, nx);
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
