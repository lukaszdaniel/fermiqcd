#include <cstdlib>
#include <cstdio>
#include <cmath>
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

void process_gauge(FILE *file_fp, FILE *MDP_fp, int buffer_size, int nx[4])
{
  double buffer[18]; // this assumes data in double precision: 144 bytes = 3 x 3 x 2 x sizeof(double)
  long position;
  int offset = sizeof(_generic_field_file_header) / sizeof(char);
  for (int x0 = 0; x0 < nx[0]; x0++)
    for (int x3 = 0; x3 < nx[3]; x3++)
      for (int x2 = 0; x2 < nx[2]; x2++)
        for (int x1 = 0; x1 < nx[1]; x1++)
          for (int mu = 0; mu < 4; mu++)
          {
            for (int i = 0; i < 3; i++)
              for (int j = 0; j < 3; j++)
              {
                if (fscanf(file_fp, "%lf%lf",
                           &(buffer[6 * i + 2 * j]),
                           &(buffer[6 * i + 2 * j + 1])))
                {
                }
              }

            // this map to the MDP ordering
            position = (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3) * 4 + ((mu + 1) % 4);
            fseek(MDP_fp, buffer_size * position + offset, SEEK_SET);
            fwrite(buffer, buffer_size, 1, MDP_fp);
          }
}

void process_quark(FILE *TONY_fp, FILE *MDP_fp, int buffer_size, int nx[4])
{
  double buffer[24];
  long position;
  int offset = sizeof(_generic_field_file_header) / sizeof(char);
  for (int x0 = 0; x0 < nx[0]; x0++)
    for (int x3 = 0; x3 < nx[3]; x3++)
      for (int x2 = 0; x2 < nx[2]; x2++)
        for (int x1 = 0; x1 < nx[1]; x1++)
        {
          for (int a = 0; a < 4; a++)
            for (int i = 0; i < 3; i++)
              if (fscanf(TONY_fp, "%lf%lf",
                         &buffer[6 * a + 2 * i],
                         &buffer[6 * a + 2 * i + 1]))
              {
              }

          // this map to the MDP ordering
          position = (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3);
          fseek(MDP_fp, buffer_size * position + offset, SEEK_SET);
          fwrite(buffer, buffer_size, 1, MDP_fp);
        }
}

int main(int argc, char **argv)
{
  printf("=======================================================\n");
  printf("Program for converting Tony Duncan gauge configurations\n");
  printf("and quark (fermi_field) into MDP files\n");
  printf("=======================================================\n");

  if (argc < 4)
  {
    printf("usage:\n\n");
    printf("duncan2mdp -gauge 16x08x08x08 input output\n\n");
    printf("duncan2mdp -quark 16x08x08x08 input output\n\n");
    exit(0);
  }

  int nx[4];
  sscanf(argv[2], "%ix%ix%ix%i", nx, nx + 1, nx + 2, nx + 3);
  long time0 = clock() / CLOCKS_PER_SEC;

  FILE *TONY_fp = fopen(argv[3], "r");
  if (!TONY_fp)
  {
    error("Cannot open input file");
  }

  FILE *MDP_fp = fopen(argv[4], "w");
  if (!MDP_fp)
  {
    fclose(TONY_fp);
    error("Cannot open output file");
  }

  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("opening the Tony Duncan file: %s (read)\n", argv[3]);
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

  write_header(MDP_fp, bytes_per_site, "Converted from Tony Duncan file", nx);
  if (strcmp(argv[1], "-gauge") == 0)
  {
    process_gauge(TONY_fp, MDP_fp, 144, nx);
  }
  else if (strcmp(argv[1], "-quark") == 0)
  {
    process_quark(TONY_fp, MDP_fp, 192, nx);
  }

  fclose(TONY_fp);
  fclose(MDP_fp);

  printf("\nAll sites are OK.\n");
  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("Header size is %i bytes.\n", offset);
  printf("Output file size is %li bytes.\n", bytes_per_site * nx[0] * nx[1] * nx[2] * nx[3] + offset);
  printf("Output file name is: %s\n", argv[4]);
  printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);

  return 0;
}
