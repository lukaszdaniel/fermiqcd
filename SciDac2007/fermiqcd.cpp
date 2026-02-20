#include <memory>
#include "fermiqcd.h"

using namespace MDP;

class FermiQCD
{
public:
  std::unique_ptr<mdp_lattice> m_plattice;
  std::unique_ptr<gauge_field> m_pU;
  std::string m_prefix;
  unsigned int m_counter;
  std::ofstream m_os;

  std::string make_prefix()
  {
    time_t time0;
    time(&time0);
    std::string s = ctime(&time0);
    while (1)
    {
      int i = s.find(" ");
      if (i >= 0)
        s = s.replace(i, 1, "_");
      else
        return s.substr(0, s.size() - 1) + "/";
    }
    return s.substr(0, s.size() - 1) + "/";
  }

  FermiQCD()
  {
    m_plattice = nullptr;
    m_pU = nullptr;
    m_prefix = make_prefix();
    // if (isMainProcess())
    //   system((std::string("mkdir ") + m_prefix).c_str());
    if (isMainProcess())
      m_os.open((m_prefix + "README.log").c_str());
    m_os << "prefix=" << m_prefix << std::endl;
    m_os << "initialization completed\n";
    m_os << "running on npocs=" << Nproc << " parallel processes\n";
    m_os << std::endl;
    m_counter = 0;
  }

  template <class T>
  FermiQCD &operator<<(const T &obj)
  {
    std::cout << ME << ":" << obj;
    if (isMainProcess())
    {
      m_os << obj;
      m_os.flush();
    }
    mdp << obj;
    return *this;
  }

  ~FermiQCD()
  {
  }

  void init_cold(int LT, int LX, int LY, int LZ, int nc_)
  {
    const Box L = {LT, LX, LY, LZ};
    m_os << "making an lattice TxXxYxZ=" << LT << "x" << LX << "x" << LY << "x" << LZ << std::endl;
    m_plattice = std::make_unique<mdp_lattice>(L);
    m_os << "making a cold gauge configuration U with nc=" << nc_ << std::endl;
    m_pU = std::make_unique<gauge_field>(*m_plattice, nc_);
    std::cout << ME << std::endl;
    m_pU->update();
    // set_cold(*m_pU);
    // m_prefix = m_prefix + std::string("C") + std::to_string((int)mdp.time());
    m_counter = 0;
  }

  void init_hot(int LT, int LX, int LY, int LZ, int nc_)
  {
    const Box L = {LT, LX, LY, LZ};
    m_os << "making an lattice TxXxYxZ=" << LT << "x" << LX << "x" << LY << "x" << LZ << std::endl;
    m_plattice = std::make_unique<mdp_lattice>(L);
    m_os << "making a hot gauge configuration U with nc=" << nc_ << std::endl;
    m_pU = std::make_unique<gauge_field>(*m_plattice, nc_);
    set_hot(*m_pU);
    m_prefix = std::string("C") + std::to_string((int)mdp.time());
    m_counter = 0;
  }

  void init_load(std::string filename)
  {
    mdp_field_file_header header = get_info(filename);
    const Box L = {header.box[0], header.box[1], header.box[2], header.box[3]};
    // (4-8)*4*(1-4-9-25-36-49-64-81-100)
    int precision = 4;
    int nc = 1;
    switch (header.bytes_per_site)
    {
    case 4 * 4 * 1:
      precision = 4;
      nc = 1;
      break;
    case 8 * 4 * 1:
      precision = 8;
      nc = 1;
      break;
    case 4 * 4 * 4:
      precision = 4;
      nc = 2;
      break;
    case 8 * 4 * 4:
      precision = 8;
      nc = 2;
      break;
    case 4 * 4 * 9:
      precision = 4;
      nc = 3;
      break;
    case 8 * 4 * 9:
      precision = 8;
      nc = 3;
      break;
    case 4 * 4 * 16:
      precision = 4;
      nc = 4;
      break;
    case 8 * 4 * 16:
      precision = 8;
      nc = 4;
      break;
    case 4 * 4 * 25:
      precision = 4;
      nc = 5;
      break;
    case 8 * 4 * 25:
      precision = 8;
      nc = 5;
      break;
    case 4 * 4 * 36:
      precision = 4;
      nc = 6;
      break;
    case 8 * 4 * 36:
      precision = 8;
      nc = 6;
      break;
    case 4 * 4 * 49:
      precision = 4;
      nc = 7;
      break;
    case 8 * 4 * 49:
      precision = 8;
      nc = 7;
      break;
    case 4 * 4 * 64:
      precision = 4;
      nc = 8;
      break;
    case 8 * 4 * 64:
      precision = 8;
      nc = 8;
      break;
    case 4 * 4 * 81:
      precision = 4;
      nc = 9;
      break;
    case 8 * 4 * 81:
      precision = 8;
      nc = 9;
      break;
    case 4 * 4 * 100:
      precision = 4;
      nc = 10;
      break;
    case 8 * 4 * 100:
      precision = 8;
      nc = 10;
      break;
    }
    m_os << "making an lattice TxXxYxZ=" << L[0] << "x" << L[1] << "x" << L[2] << "x" << L[3] << std::endl;
    m_plattice = std::make_unique<mdp_lattice>(L);
    m_os << "loading gauge configuration " << filename << " into U with nc=" << nc << std::endl;
    m_pU = std::make_unique<gauge_field>(*m_plattice, nc);
    if (precision == 4)
      m_pU->load_as_float(filename);
    if (precision == 8)
      m_pU->load_as_double(filename);
    m_prefix = filename;
    m_counter = 0;
  }

  void wilson_heatbath(float beta, int steps = 1)
  {
    coefficients gauge;
    gauge["beta"] = beta;
    for (int step = 0; step < steps; step++)
    {
      (*this) << "WilsonGaugeAction::heatbath(beta=" << beta << ") step=" << step << "\n";
      WilsonGaugeAction::heatbath(*m_pU, gauge, 1);
    }
    (*this) << "average_plaquette=" << average_plaquette(*m_pU);
    m_counter += steps;
  }

  void ape_smearing(float alpha = 0.7, int iterations = 20, int cooling_steps = 10)
  {
    (*this) << "starting Ape smering...\n";
    ApeSmearing::smear(*m_pU, alpha, iterations, cooling_steps);
    (*this) << "done!\n";
  }

  void coulomb_gauge_fix(int iterations = 100)
  {
    (*this) << "starting coulomb gauge fixing...\n";
    float precision = 1e-6;
    float boost = 1;
    /*gaugefixing_stats stats = */ GaugeFixing::fix(*m_pU, GaugeFixing::Coulomb, iterations, precision, boost, false);
    (*this) << "done!\n";
  }

  void coulomb_z3_gauge_fix(int iterations = 100)
  {
    (*this) << "starting coulomb gauge fixing...\n";
    float precision = 1e-6;
    float boost = 1;
    /*gaugefixing_stats stats = */ GaugeFixing::fix(*m_pU, GaugeFixing::Coulomb, iterations, precision, boost, true);
    (*this) << "done!\n";
  }

  void landau_gauge_fix(int iterations = 100)
  {
    (*this) << "starting coulomb gauge fixing...\n";
    float precision = 1e-6;
    float boost = 1;
    /*gaugefixing_stats stats = */ GaugeFixing::fix(*m_pU, GaugeFixing::Landau, iterations, precision, boost, false);
    (*this) << "done!\n";
  }

  void compute_topological_charge()
  {
    mdp_real_scalar_field Q(*m_plattice);
    topological_charge(Q, *m_pU);
    Q.save_vtk(m_prefix + "topological_charge");
  }

  void compute_partitioning()
  {
    mdp_site x(*m_plattice);
    mdp_real_scalar_field Q(*m_plattice);
    mdp_real_scalar_field p(*m_plattice);
    forallsites(x)
        p(x) = on_which_process(*m_plattice, x(0), x(1), x(2), x(3));
    Q.save_vtk(m_prefix + "partitioning");
  }
};

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  FermiQCD fermiqcd;
  fermiqcd.init_cold(10, 10, 4, 4, 3);
  fermiqcd.compute_partitioning();
  fermiqcd.compute_topological_charge();

  mdp.close_wormholes();
  return 0;
}
