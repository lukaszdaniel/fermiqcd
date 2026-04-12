#include "fermiqcd.h"

using namespace MDP;

int main()
{
  std::cout << "Size of various types and objects (in bytes):\n\n";
  std::cout << "sizeof(bool): " << sizeof(bool)
            << "\nsizeof(char): " << sizeof(char)
            << "\nsizeof(mdp_ssint): " << sizeof(mdp_ssint)
            << "\nsizeof(mdp_ssuint): " << sizeof(mdp_ssuint)
            << "\nsizeof(mdp_sint): " << sizeof(mdp_sint)
            << "\nsizeof(mdp_suint): " << sizeof(mdp_suint)
            << "\nsizeof(mdp_int): " << sizeof(mdp_int)
            << "\nsizeof(mdp_uint): " << sizeof(mdp_uint)
            << "\nsizeof(int): " << sizeof(int)
            << "\nsizeof(unsigned): " << sizeof(unsigned int)
            << "\nsizeof(mdp_real): " << sizeof(mdp_real)
            << "\nsizeof(double): " << sizeof(double)
            << "\nsizeof(size_t): " << sizeof(size_t)
            << "\nsizeof(mdp_complex): " << sizeof(mdp_complex)
            << "\nsizeof(mdp_field_file_header): " << sizeof(mdp_field_file_header)
            << "\nsizeof(mdp_int_field): " << sizeof(mdp_int_field)
            << "\nsizeof(mdp_real_field): " << sizeof(mdp_real_field)
            << "\nsizeof(mdp_complex_field): " << sizeof(mdp_complex_field)
            << "\nsizeof(mdp_nmatrix_field): " << sizeof(mdp_nmatrix_field)
            << "\nsizeof(mdp_matrix_field): " << sizeof(mdp_matrix_field)
            << "\nsizeof(mdp_nvector_field): " << sizeof(mdp_nvector_field)
            << "\nsizeof(mdp_vector_field): " << sizeof(mdp_vector_field)
            << "\nsizeof(gauge_field): " << sizeof(gauge_field)
            << "\nsizeof(em_field): " << sizeof(em_field)
            << "\nsizeof(dwfermi_field): " << sizeof(dwfermi_field)
            << "\nsizeof(fermi_field): " << sizeof(fermi_field)
            << "\nsizeof(fermi_propagator): " << sizeof(fermi_propagator)
            << "\nsizeof(sdwf_field): " << sizeof(sdwf_field)
            << "\nsizeof(staggered_field): " << sizeof(staggered_field)
            << "\nsizeof(staggered_propagator): " << sizeof(staggered_propagator)
            << "\nsizeof(mdp_lattice): " << sizeof(mdp_lattice)
            << '\n';
  return 0;
}