// ukqcd2mdp from the FermiQCD utilities by Massimo di Pierro
// last change: 27 Mar 2003 JMF

// read UKQCD binary gauge fields or propagators and save as mdp files

// compile with something like: g++ -o ukqcd2mdp ukqcd2mdp.C

// run program with no arguments for help

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <complex>
#include <ctime>
#include <memory>

using namespace std;

#define Complex std::complex<float>
#define Nspace nx[1] * nx[2] * nx[3]
#define Ndim 4

void error(const char s[])
{
  printf("ERROR: %s\n", s);
  exit(1);
}

int nx[4];

///////////////////////////////
// UKQCD byte-swapping routines
///////////////////////////////

// union lets you store values in the same memory that can be accessed
// as float or char[4] (or, double or char[8] for double precision)

void block_swap(float *buffer, long length)
{
  int i;
  union swapper
  {
    float float_number;
    char pos[4];
  } a, b;

  for (i = 0; i < length; i++)
  {
    a.float_number = *buffer;
    b.pos[0] = a.pos[3];
    b.pos[1] = a.pos[2];
    b.pos[2] = a.pos[1];
    b.pos[3] = a.pos[0];
    *buffer = b.float_number;
    buffer++;
  }
}

void block_swap_double(double *buffer, long length)
{
  int i;
  union swapper
  {
    double double_number;
    char pos[8];
  } a, b;

  for (i = 0; i < length; i++)
  {
    a.double_number = *buffer;
    b.pos[0] = a.pos[7];
    b.pos[1] = a.pos[6];
    b.pos[2] = a.pos[5];
    b.pos[3] = a.pos[4];
    b.pos[4] = a.pos[3];
    b.pos[5] = a.pos[2];
    b.pos[6] = a.pos[1];
    b.pos[7] = a.pos[0];
    *buffer = b.double_number;
    buffer++;
  }
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

  void initialize(int x1, int x2, int x3, int a = 1, int b = 1, int c = 1, int d = 1)
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

// read a single time slice of a gauge field
void read_t_gauge(short_field &U, char fileprefix[],
                  char precision, char swap, int time)
{

  char filename[200];
  std::unique_ptr<unsigned char[]> buffer;      // array to store the raw data
  long file_length = Ndim * 12 * sizeof(float); // file length in bytes (one time slice)
  long bytes_read;                              // monitor how many bytes we read
  int x1, x2, x3, mu, a, b, muf;
  // int c;
  long buf_index; // indices for the temporary buffer
  // long site_index;
  FILE *fp; // pointer to input file

  // Check input makes sense -- ie that precision is correct
  if ((precision != 'F') && (precision != 'D'))
    error("Invalid precision");
  // Check input -- whether byte swap is correct
  if ((swap != 'Y') && (swap != 'N'))
    error("Invalid swapping parameter");

  // for each time have Nspace sites, Ndim directions and 6 complex
  // numbers stored.
  // file_length=number of bytes in file
  if (precision == 'F')
    file_length = Nspace * Ndim * 12 * sizeof(float);
  if (precision == 'D')
    file_length = Nspace * Ndim * 12 * sizeof(double);

  // CARE HERE!
  snprintf(filename, 200, "%sT%02d", fileprefix, time);

  printf("Opening file: %s\n", filename);

  // Open, read, check and close the file! Allocate storage space in
  // buffer[].
  buffer = std::make_unique<unsigned char[]>(file_length);
  if (buffer == 0x0)
    error("Unable to allocate memory");
  fp = fopen(filename, "rb");
  if (fp == (FILE *)NULL)
    error("Unable to open file");
  bytes_read = fread(buffer.get(), sizeof(unsigned char), file_length, fp);
  if (bytes_read != file_length)
    error("Read wrong number of bytes");
  fclose(fp);

  // Do block swapping if necessary
  if (swap == 'Y')
  {
    if (precision == 'F')
    {
      block_swap((float *)buffer.get(), file_length / sizeof(float));
    }
    else
    {
      block_swap_double((double *)buffer.get(), file_length / sizeof(double));
    }
  }

  // for checking...
  // if (precision == 'F')
  // {
  //   printf("buffer[0] = %g\n", ((float *)buffer.get())[0]);
  // }
  // else
  // {
  //   printf("buffer[0] = %g\n", ((double *)buffer.get())[0]);
  // }

  // Copy the whole thing into the lattice structure by looping over the
  // lattice indices in the right order while continuously stepping along
  // the input buffer

  // Start of buffers
  buf_index = 0;
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Attention: the most internal loop is x1 according to Balint Joo
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for (x3 = 0; x3 < nx[3]; x3++)
    for (x2 = 0; x2 < nx[2]; x2++)
      for (x1 = 0; x1 < nx[1]; x1++)
      {
        for (muf = 1; muf < 5; muf++)
        { // ukqcd indices 1, 2, 3, 4
          if (muf == 4)
            mu = 0;
          else
            mu = muf; // mdp indices   1, 2, 3, 0
          for (a = 0; a < (3); a++)
            for (b = 0; b < (3 - 1); b++)
            {
              // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              // Attention: Balint says a are columns and b are rows
              // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (precision == 'F')
              {
                U(x1, x2, x3, mu, b, a) = Complex(((float *)buffer.get())[buf_index],
                                                  ((float *)buffer.get())[buf_index + 1]);
              }
              else
              {
                U(x1, x2, x3, mu, b, a) = Complex(((double *)buffer.get())[buf_index],
                                                  ((double *)buffer.get())[buf_index + 1]);
              }
              buf_index += 2;
            }
          // Impose special unitarity to reconstruct full SU(3) matrix
          U(x1, x2, x3, mu, 2, 0) = conj(U(x1, x2, x3, mu, 0, 1) * U(x1, x2, x3, mu, 1, 2) -
                                         U(x1, x2, x3, mu, 1, 1) * U(x1, x2, x3, mu, 0, 2));
          U(x1, x2, x3, mu, 2, 1) = conj(U(x1, x2, x3, mu, 1, 0) * U(x1, x2, x3, mu, 0, 2) -
                                         U(x1, x2, x3, mu, 0, 0) * U(x1, x2, x3, mu, 1, 2));
          U(x1, x2, x3, mu, 2, 2) = conj(U(x1, x2, x3, mu, 0, 0) * U(x1, x2, x3, mu, 1, 1) -
                                         U(x1, x2, x3, mu, 1, 0) * U(x1, x2, x3, mu, 0, 1));
        }
      }

  buffer.reset(); // free the buffer

  // for checking...
  // printf("U(0,0,0,0,0,1) = %g + I %g\n", real(U(0,0,0,0,0,1)),
  //	 real(U(0,0,0,0,0,1)));
}

/* Read all source spin/colour combinations for one time slice of a
   propagator.

   Filename prefix should contain everything up to source spin, and source
   colour so it should be something like:

   Q60U<sweep number><source fuzz><source kappa><sink fuzz><sink_kappa>

   And we will add onto it 00Ttt
                           01Ttt
         02Ttt
         10Ttt

         .
         .
         .
         32Ttt probably via some loop
*/
void read_t_prop(short_field &S, char fileprefix[],
                 char precision, char swap, int time)
{

  char filename[200];
  FILE *fp;
  std::unique_ptr<unsigned char[]> buffer;
  long bytes_to_read;
  long bytes_read;
  int x1, x2, x3;
  // int x0;
  // int comp;
  long buffer_index;

  printf("Reading propagator on timeslice %d ...\n", time);

  // Check input makes sense -- ie that precision is correct
  if ((precision != 'F') && (precision != 'D'))
    error("Invalid precision");
  // Check input -- whether byte swap is correct
  if ((swap != 'Y') && (swap != 'N'))
    error("Invalid swapping parameter");

  // File contains only sink spin and sink colour
  if (precision == 'F')
    bytes_to_read = Nspace * 3 * 4 * 2 * sizeof(float);
  else
    bytes_to_read = Nspace * 3 * 4 * 2 * sizeof(double);

  // Allocate space to buffer[]
  buffer = std::make_unique<unsigned char[]>(bytes_to_read);
  if (buffer == 0x0)
    error("Out of memory");

  for (int source_spin = 0; source_spin < 4; source_spin++)
  {
    for (int source_colour = 0; source_colour < 3; source_colour++)
    {

      // Construct filename
      snprintf(filename, 200, "%s%01d%01dT%02d", fileprefix, source_spin,
               source_colour, time);
      printf("Opening file: %s\n", filename);

      // Open, read, check and close file
      fp = fopen(filename, "rb");
      if (fp == (FILE *)NULL)
        error("Unable to open file");
      bytes_read = fread(buffer.get(), sizeof(unsigned char), bytes_to_read, fp);
      if (bytes_read != bytes_to_read)
        error("Wrong number of bytes read");
      fclose(fp);

      // Do block swapping if necessary
      if (swap == 'Y')
      {
        if (precision == 'F')
        {
          block_swap((float *)buffer.get(), bytes_to_read / sizeof(float));
        }
        else
        {
          block_swap_double((double *)buffer.get(), bytes_to_read / sizeof(double));
        }
      }

      buffer_index = 0;
      for (x3 = 0; x3 < nx[3]; x3++)
        for (x2 = 0; x2 < nx[2]; x2++)
          for (x1 = 0; x1 < nx[1]; x1++)
          {

            for (int sink_spin = 0; sink_spin < 4; sink_spin++)
            {
              for (int sink_colour = 0; sink_colour < 3; sink_colour++)
              {

                if (precision == 'F')
                {
                  S(x1, x2, x3, sink_spin, source_spin, sink_colour, source_colour) =
                      Complex(((float *)buffer.get())[buffer_index],
                              ((float *)buffer.get())[buffer_index + 1]);
                }
                else
                {
                  S(x1, x2, x3, sink_spin, source_spin, sink_colour, source_colour) =
                      Complex(((double *)buffer.get())[buffer_index],
                              ((double *)buffer.get())[buffer_index + 1]);
                }
                buffer_index += 2;
              }   // sink_colour;
            }     // sink spin
          }       // spatial loops
    }             // source colour
  }               // source spin
  buffer.reset(); // now free the buffer
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

  printf("=============================================================\n");
  printf("ukqcd2mdp: Convert UKQCD gauge configurations and propagators\n");
  printf("to MDP format.\n");
  printf("=============================================================\n");
  printf("!!! CONVERSION NOT THOROUGHLY TESTED YET !!!\n");

  if (argc < 4)
  {
    printf("Sample usage:\n\n");
    printf("ukqcd2mdp -gauge 16x8x8x8FN input_prefix output\n");
    printf("ukqcd2mdp -fermi 48x24x24x24DY input_prefix output\n\n");
    printf("UKQCD gauge filenames: [input_prefix]Txx\n");
    printf("UKQCD props filenames: [input_prefix]yyTxx\n");
    printf("where xx is the timeslice and yy is the spin source\n\n");
    printf("F=float, D=double, Y=swap, N=no swap\n\n");
    exit(0);
  };

  sscanf(argv[2], "%ix%ix%ix%i", nx, nx + 1, nx + 2, nx + 3);

  char PRECISION = argv[2][strlen(argv[2]) - 2];
  char SWAP = argv[2][strlen(argv[2]) - 1];

  long time0 = clock() / CLOCKS_PER_SEC;

  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("opening the MDP file: %s (write) \n", argv[4]);
  FILE *MDP_fp = fopen(argv[4], "w");

  _generic_field_file_header myheader;
  short_field U; // stores one timeslice

  myheader.ndim = 4;
  for (int ii = 0; ii < 4; ii++)
    myheader.box_size[ii] = nx[ii];
  for (int ii = 4; ii < 10; ii++)
    myheader.box_size[ii] = 0;
  myheader.sites = nx[0] * nx[1] * nx[2] * nx[3];
  // gauge: bytes_per_site = 4(mu)*9(SU3 matrix)*2(cplx)*4(bytes_per_float)
  //                       = 288
  // fermi: bytes_per_site = 16(spin-sq)*9(color-sq)*2(cplx)*4(bytes_per_float)
  //                       = 1152
  if (strcmp(argv[1], "-gauge") == 0)
    myheader.bytes_per_site = 288;
  if (strcmp(argv[1], "-fermi") == 0)
    myheader.bytes_per_site = 1152;
  myheader.endianess = 0x87654321;
  strcpy(myheader.program_version, "Converted from UKQCD");
  time_t time_and_date;
  time(&time_and_date);
  strcpy(myheader.creation_date, ctime(&time_and_date));
  int offset = sizeof(_generic_field_file_header) / sizeof(char);
  fwrite(&myheader, sizeof(char), offset, MDP_fp);

  if (strcmp(argv[1], "-gauge") == 0)
  {
    // U indices: space, space, space, mu, color, color
    U.initialize(nx[1], nx[2], nx[3], 4, 3, 3);
    for (int x0 = 0; x0 < nx[0]; x0++)
    {
      read_t_gauge(U, argv[3], PRECISION, SWAP, x0);
      fseek(MDP_fp, myheader.bytes_per_site * Nspace * x0 + offset, SEEK_SET);
      fwrite(U.m_data.get(), sizeof(Complex), U.size, MDP_fp);
    }
  }
  if (strcmp(argv[1], "-fermi") == 0)
  {
    U.initialize(nx[1], nx[2], nx[3], 4, 4, 3, 3);
    for (int x0 = 0; x0 < nx[0]; x0++)
    {
      read_t_prop(U, argv[3], PRECISION, SWAP, x0);
      fseek(MDP_fp, myheader.bytes_per_site * Nspace * x0 + offset, SEEK_SET);
      fwrite(U.m_data.get(), sizeof(Complex), U.size, MDP_fp);
    }
  }

  fclose(MDP_fp);

  printf("\nAll sites seem OK.\n");
  printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
  printf("Header size is %i bytes.\n", offset);
  printf("Output file size is %li bytes.\n",
         myheader.bytes_per_site * myheader.sites + offset);
  printf("Output file name is: %s\n", argv[4]);
  printf("Done in %li secs.\n", clock() / CLOCKS_PER_SEC - time0);

  return 0; // success
}
