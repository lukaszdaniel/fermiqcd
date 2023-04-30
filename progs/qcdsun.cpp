#include <iomanip>
#include <iostream>
#include <memory>
#include "fermiqcd.h"
#include "my_parameters.h"
#include "my_logger.h"

using namespace std;
using namespace MDP;

int main(int argc, char **argv)
{

	if (argc < 2)
	{
		std::cerr << "usage:\t" << argv[0] << " input_file" << std::endl;
		exit(1);
	}

	mdp.open_wormholes(argc, argv);

	parameter par(argv[1]);

	string tag = par.s("tag");
	string logFileName = "log." + tag;
	string cfgFileName = "cfg." + tag;
	string plaqFileName = "plaq." + tag;
	string polyFileName = "poly." + tag;
	string polyCorFileName = "corpoly." + tag;
	string plaqCorFileName = "corplaq." + tag;
	string startConfigFile;

	mdp_int dim = 3;
	mdp_suint ncolors = 3;
	mdp_suint seed = time(0);
	mdp_uint ntherm = 0;
	mdp_uint nsweeps = 0;
	mdp_suint unit_freq = 50;
	mdp_suint relax_freq = 1;
	mdp_suint poly_meas_freq = 0;
	mdp_suint poly_cor_meas_freq = 0;
	// unsigned short int plaq_cor_meas_freq = 0;

	// unsigned int ns = par.i("nsmear");
	// double epsilon = par.d("epsilon");
	// double norm = par.d("norm");
	// unsigned int nsiter = par.i("nsiter");

	bool field_loaded = false;
	// bool z3symmetry = true;
	bool local_random = false;

	if (par.defined_i("colors"))
		ncolors = par.i("colors");
	if (par.defined_i("ndim"))
		dim = par.i("ndim");
	if (par.defined_i("relax-freq"))
		relax_freq = par.i("relax-freq");
	if (par.defined_i("unit-freq"))
		unit_freq = par.i("unit-freq");
	if (par.defined_i("nterm"))
		ntherm = par.i("nterm");
	if (par.defined_i("nsweep"))
		nsweeps = par.i("nsweep");
	if (par.defined_i("poly-meas-freq"))
		poly_meas_freq = par.i("poly-meas-freq");
	if (par.defined_i("poly-cor-meas-freq"))
		poly_cor_meas_freq = par.i("poly-cor-meas-freq");
	// if (par.defined_i("plaq-cor-meas-freq"))
	// 	plaq_cor_meas_freq = par.i("plaq-cor-meas-freq");
	if (par.defined_i("seed"))
		seed = par.i("seed");
	if (not par.defined_d("beta"))
	{
		std::cerr << "Beta undefined. Exiting..." << std::endl;
		exit(1);
	}
	if (!(par.defined_i("L_time") && par.defined_i("L_space")))
	{
		std::cerr << "L_space or L_time undefined. Exiting..." << std::endl;
		exit(1);
	}

	mdp_uint Ls = par.i("L_space");

	mdp_uint Lt = par.i("L_time");
	mdp_real stats;

	mdp_int mybox[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	mybox[0] = Lt;
	for (mdp_suint i = 1; i < dim; i++)
		mybox[i] = Ls;
	mdp_lattice lattice(dim, mybox, default_partitioning0, torus_topology, seed, local_random);
	gauge_field U(lattice, ncolors);
	coefficients coeff;
	coeff["beta"] = par.d("beta");
	FILE *stdlog = nullptr;
	stdlog = fopen(logFileName.c_str(), "w");

	if (par.defined_s("start-config-file"))
	{
		startConfigFile = par.s("start-config-file");
		U.load(startConfigFile);
		field_loaded = true;
	}
	else
		set_hot(U); // sets all U matrices by random from SU(Nc) group

	openlog(logFileName.c_str());
	header(argc, argv, "");
	par.fprint(stdlog);

	//	THERMALISATION
	if (!field_loaded)
	{
		fprintf(stdlog, "Thermalisation: %d sweeps.\n", ntherm);
		fflush(stdlog);
		for (mdp_uint i = 1; i <= ntherm; i++)
		{
			WilsonGaugeAction::heatbath(U, coeff);
			relaxation(U, relax_freq);
			if ((unit_freq > 0) and (i % unit_freq == 0))
			{
				stats = unitarize(U);
				fprintf(stdlog, "%d %.12lg\n", i, stats);
			}
			fflush(stdlog);
			U.save(cfgFileName);
		}
	}
	else
		fprintf(stdlog, "Field has been loaded. Skipping thermalisation.\n");

	ofstream plaqfile(plaqFileName);
	plaqfile << setprecision(16);
	ofstream polyfile(polyFileName);
	polyfile << setprecision(16);
	ofstream polycorfile(polyCorFileName);
	polycorfile << setprecision(16);
	// ofstream plaqcorfile(plaqCorFileName);	plaqcorfile << setprecision(16);

	mdp_real Pt, Ps;
	mdp_real S;
	const mdp_real D = 1.0 * U.ndim() * (U.ndim() - 1) / 2;
	const mdp_real as = 1.0 * (U.ndim() - 2) * (U.ndim() - 1) / 2;
	const mdp_real at = 1.0 * U.ndim() - 1;

	// const mdp_real c1 = -1.0 / 12; // Symanzik
	// const mdp_real c1 = -0.331;	   // Iwasaki
	// const mdp_real c1 = -1.4088;   // DBW2
	// const mdp_real c0 = 1 - 8 * c1;
	mdp_complex polyakov;

	//	SWEEPING
	fprintf(stdlog, "Production: %d sweeps.\n", nsweeps);
	fflush(stdlog);
	for (mdp_uint i = 1; i <= nsweeps; i++)
	{
		WilsonGaugeAction::heatbath(U, coeff);
		relaxation(U, relax_freq);

		Pt = TimePlaquette(U);
		Ps = SpacePlaquette(U);
		S = (D - (as * Ps + at * Pt) / U.nc());

		plaqfile << Pt << " " << Ps << " " << S << "\n";

		if (poly_meas_freq > 0 && i % poly_meas_freq == 0)
		{
			polyakov = PolyakovLoop(U);
			polyfile << real(polyakov) << " " << imag(polyakov) << " " << Pt << " " << Ps << "\n";
		}

		if ((poly_cor_meas_freq > 0) && (i % poly_cor_meas_freq == 0))
		{

			std::unique_ptr<double[]> cor = std::make_unique<double[]>(U.lattice().size(1));

			mdp_matrix a(U.ndim() - 1, U.lattice().size(1));
			a = PolyCor(U);

			mdp_complex average = 0;

			for (mdp_int j = 0; j < U.lattice().size(1); j++)
			{
				average += a(0, j);
			}
			average /= U.lattice().size(1);

			polycorfile << "#\t" << U.lattice().size(1) << "\n";
			polycorfile << "# PAR_0\t" << real(average) << "\n";
			polycorfile << "# PAR_1\t" << imag(average) << "\n";

			for (mdp_int k = 0; k < U.lattice().size(1); k++)
			{
				cor[k] = 0;
				for (mdp_int mu = 1; mu < U.ndim(); mu++)
				{
					for (mdp_int j = 0; j < U.lattice().size(1); j++)
					{
						cor[k] += real(a(mu - 1, j) * conj(a(mu - 1, (k + j) % U.lattice().size(1))));
					}
				}
				cor[k] /= (2.0 * U.lattice().size(1));
				polycorfile << cor[k] << "\n";
			}
			polycorfile << flush;
		}

		if (unit_freq > 0 && i % unit_freq == 0)
		{

			stats = unitarize(U);
			fprintf(stdlog, "%d %.12lg\n", i, stats);
			plaqfile << flush;
			polyfile << flush;
			fflush(stdlog);
			U.save(cfgFileName);
		}
	}

	U.save(cfgFileName);

	mdp.close_wormholes();

	//  closelog();
	plaqfile.close();
	polyfile.close();
	polycorfile.close();

	return EXIT_SUCCESS;
}
