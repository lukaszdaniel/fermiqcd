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

  int Ndim = 4;
  long position;
  long time0 = clock() / CLOCKS_PER_SEC;

  if (strcmp(argv[1], "-gauge") == 0)
  {
    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the Tony Duncan file: %s (read)\n", argv[3]);
    FILE *TONY_fp = fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp = fopen(argv[4], "w");

    _generic_field_file_header myheader;

    myheader.ndim = 4;
    for (int ii = 0; ii < 4; ii++)
      myheader.box_size[ii] = nx[ii];
    for (int ii = 4; ii < 10; ii++)
      myheader.box_size[ii] = 0;
    myheader.sites = nx[0] * nx[1] * nx[2] * nx[3];
    myheader.bytes_per_site = 288 * 2;
    myheader.endianess = 0x87654321;
    strcpy(myheader.program_version, "Converted from Tony Duncan file");
    time_t time_and_date;
    time(&time_and_date);
    strcpy(myheader.creation_date, ctime(&time_and_date));
    int offset = sizeof(_generic_field_file_header) / sizeof(char);
    fwrite(&myheader, sizeof(char), offset, MDP_fp);

    double buffer[18]; // this assumes data in double precision: 144 = 9 x 2 x sizeof(double)
    for (int x0 = 0; x0 < nx[0]; x0++)
      for (int x3 = 0; x3 < nx[3]; x3++)
        for (int x2 = 0; x2 < nx[2]; x2++)
          for (int x1 = 0; x1 < nx[1]; x1++)
            for (int mu = 0; mu < 4; mu++)
            {
              for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                  if (fscanf(TONY_fp, "%lf%lf",
                             &(buffer[6 * j + 2 * i]),
                             &(buffer[6 * j + 2 * i + 1])));
                  // printf("%e %e\n", buffer[6*j+2*i], buffer[6*j+2*i+1]);
                  buffer[6 * j + 2 * i + 1] *= -1;
                }

              // this map to the MDP ordering
              position = (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3) * 4 + ((mu + 1) % Ndim);
              fseek(MDP_fp, 144 * position + offset, SEEK_SET);
              fwrite(buffer, 144, 1, MDP_fp);
            }
    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);
  }
  else if (strcmp(argv[1], "-quark") == 0)
  {
    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the Tony Duncan file: %s (read)\n", argv[3]);
    FILE *TONY_fp = fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp = fopen(argv[4], "w");

    _generic_field_file_header myheader;

    myheader.ndim = 4;
    for (int ii = 0; ii < 4; ii++)
      myheader.box_size[ii] = nx[ii];
    for (int ii = 4; ii < 10; ii++)
      myheader.box_size[ii] = 0;
    myheader.sites = nx[0] * nx[1] * nx[2] * nx[3];
    myheader.bytes_per_site = 96;
    myheader.endianess = 0x87654321;
    strcpy(myheader.program_version, "Converted from Tony Duncan file");
    time_t time_and_date;
    time(&time_and_date);
    strcpy(myheader.creation_date, ctime(&time_and_date));
    int offset = sizeof(_generic_field_file_header) / sizeof(char);
    fwrite(&myheader, sizeof(char), offset, MDP_fp);

    float buffer[24];
    for (int x0 = 0; x0 < nx[0]; x0++)
      for (int x3 = 0; x3 < nx[3]; x3++)
        for (int x2 = 0; x2 < nx[2]; x2++)
          for (int x1 = 0; x1 < nx[1]; x1++)
          {
            for (int a = 0; a < 4; a++)
              for (int i = 0; i < 3; i++)
                if (fscanf(TONY_fp, "%f%f",
                           &buffer[6 * a + 2 * i],
                           &buffer[6 * a + 2 * i + 1]));

            // this map to the MDP ordering
            position = (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3);
            fseek(MDP_fp, 96 * position + offset, SEEK_SET);
            fwrite(buffer, 96, 1, MDP_fp);
          }

    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);
  }
  else if (strcmp(argv[1], "-gauge:d") == 0)
  {
    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the Tony Duncan file: %s (read)\n", argv[3]);
    FILE *TONY_fp = fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp = fopen(argv[4], "w");

    _generic_field_file_header myheader;

    myheader.ndim = 4;
    for (int ii = 0; ii < 4; ii++)
      myheader.box_size[ii] = nx[ii];
    for (int ii = 4; ii < 10; ii++)
      myheader.box_size[ii] = 0;
    myheader.sites = nx[0] * nx[1] * nx[2] * nx[3];
    myheader.bytes_per_site = 576;
    myheader.endianess = 0x87654321;
    strcpy(myheader.program_version, "Converted from Tony Duncan file");
    time_t time_and_date;
    time(&time_and_date);
    strcpy(myheader.creation_date, ctime(&time_and_date));
    int offset = sizeof(_generic_field_file_header) / sizeof(char);
    fwrite(&myheader, sizeof(char), offset, MDP_fp);

    double buffer[18]; // this assumes data in double precision: 144 = 9 x 2 x 8
    for (int x0 = 0; x0 < nx[0]; x0++)
      for (int x3 = 0; x3 < nx[3]; x3++)
        for (int x2 = 0; x2 < nx[2]; x2++)
          for (int x1 = 0; x1 < nx[1]; x1++)
            for (int mu = 0; mu < 4; mu++)
            {
              for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                  if (fscanf(TONY_fp, "%lf%lf",
                             &(buffer[6 * j + 2 * i]),
                             &(buffer[6 * j + 2 * i + 1])));
                  buffer[6 * j + 2 * i + 1] *= -1;
                }
              // this map to the MDP ordering
              position = (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3) * 4 + ((mu + 1) % Ndim);
              fseek(MDP_fp, 144 * position + offset, SEEK_SET);
              fwrite(buffer, 144, 1, MDP_fp);
            }
    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);
  }
  else if (strcmp(argv[1], "-quark:d") == 0)
  {
    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the Tony Duncan file: %s (read)\n", argv[3]);
    FILE *TONY_fp = fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp = fopen(argv[4], "w");

    _generic_field_file_header myheader;

    myheader.ndim = 4;
    for (int ii = 0; ii < 4; ii++)
      myheader.box_size[ii] = nx[ii];
    for (int ii = 4; ii < 10; ii++)
      myheader.box_size[ii] = 0;
    myheader.sites = nx[0] * nx[1] * nx[2] * nx[3];
    myheader.bytes_per_site = 192;
    myheader.endianess = 0x87654321;
    strcpy(myheader.program_version, "Converted from Tony Duncan file");
    time_t time_and_date;
    time(&time_and_date);
    strcpy(myheader.creation_date, ctime(&time_and_date));
    int offset = sizeof(_generic_field_file_header) / sizeof(char);
    fwrite(&myheader, sizeof(char), offset, MDP_fp);

    double buffer[24];
    for (int x0 = 0; x0 < nx[0]; x0++)
      for (int x3 = 0; x3 < nx[3]; x3++)
        for (int x2 = 0; x2 < nx[2]; x2++)
          for (int x1 = 0; x1 < nx[1]; x1++)
          {
            for (int a = 0; a < 4; a++)
              for (int i = 0; i < 3; i++)
                if (fscanf(TONY_fp, "%lf%lf",
                           &buffer[6 * a + 2 * i],
                           &buffer[6 * a + 2 * i + 1]));

            // this map to the MDP ordering
            position = (((x0 * nx[1] + x1) * nx[2] + x2) * nx[3] + x3);
            fseek(MDP_fp, 192 * position + offset, SEEK_SET);
            fwrite(buffer, 192, 1, MDP_fp);
          }

    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);
  }
}
