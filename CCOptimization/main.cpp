/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: February 2022 */

/** \file main.cpp 
 * \brief Simple CLI for generating classical ground-state configurations of 'stealthy' collective coordinate potentials in 1,2,3 dimensions. */


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <nlopt.h>

/* in $(cores) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <etc.h>
/* in $(pair) */
#include <StructureFactor.h>
#include <PairCorrelation.h>
/* in $(pot) */
#include <StructureOptimization.h>
#include <CollectiveCoordinatePotential.h>
#include <MD_System.h>

size_t Verbosity = 3;

#include "RepulsiveCCPotential.h"
#include "MultiStealthyPotential.h"
#include "MultiShellS0Potential.h"

/** \brief Perform the energy optimization multiple times for a given potential energy.
Obtain ground-state configurations from multi initial configurations.
@param[in]	pConfig	A pointer of a point configuration.
@param[in]	pPotential	A pointer of a Potential object.
@param[in]	gen	A random number generator.
@param[in]	Prefix	A prefix of output files.
			It generates three distinct ConfigurationPacks;
			*_InitConfig	Initial configurations
			*				Configurations after minimizing energy
			*_Success		Configurations whose energy is lower than 1e-11 (a subset of *)
@param[in]	SampleNumber	The number of initial configurations
@param[in]	TimeLimit		When running time becomes larger than this, this function automatically ends.
@param[in]	Algorithm	Name of alogrithm for optimizations. The followings are allowed; LBFGS, LocalGradientDescent, ConjugateGradient, SteepestDescent, MINOP
@return 0 = exit without any problem. */
int CollectiveCoordinateMultiRun(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, std::string Prefix, size_t SampleNumber, time_t TimeLimit, const std::string & AlgorithmName, double E_max, size_t max_num);


/** \brief Perform a molecular dynamic simulation under NVT environment for a given potential energy.
By default, it adjust time step optimally.
@param[in]	pConfig	A pointer of a point configuration.
@param[in]	pPotential	A pointer of a Potential object.
@param[in]	gen	A random number generator.
@param[in]	TimeStep Initial time step for the simulation.
@param[in]	Temperature Temperature of the system.
@param[in]	Prefix	A prefix of output files.
@param[in]	SampleNumber	The number of target configurations.
@param[in]	StepPerSample	The number of time steps per each configuration.
@param[in]	AllowRestore	True if you want to load the previous setup (saved as *.MDDump).
@param[in]	TimeLimit		When running time becomes larger than this, this function automatically ends.
@param[in]	EquilibrateSamples	The number of sampling for equilibration.
@param[in]	MDAutoTimeStep	True if you want to adjust time step during the equilibration.
@return 0 = exit without any problem. */
int CollectiveCoordinateMD(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, double TimeStep, double Temperature, std::string Prefix, size_t SampleNumber, size_t StepPerSample, bool AllowRestore, time_t TimeLimit, size_t EquilibrateSamples, bool MDAutoTimeStep);


/// @brief Get radius of the exclusion region in Fourier space corresponding to the chi value
/// @param chi 	stealthiness parpameter
/// @param rho 	number density
/// @param d 	dimensions
/// @return 
double getK(double chi, double rho, DimensionType d);

/// @brief Get stealthiness parameter chi from K and rho.
/// @param K 
/// @param rho 
/// @param d 
/// @return 
double get_chi(double K, double rho, DimensionType d);

class ReadConfigPack
{
public:
	double Rescale;
	ConfigurationPack p;
	ReadConfigPack(std::istream & ifile, std::ostream & ofile, double Rescale)
	{
		this->Rescale=Rescale;
		ofile<<"Input Prefix:";
		std::string prefix;
		ifile>>prefix;
		p.Open(prefix);
	}
	ReadConfigPack(std::string prefix, double Rescale)
	{
		this->Rescale=Rescale;
		p.Open(prefix);
	}
	Configuration operator() (size_t i)
	{
		Configuration result=p.GetConfig((long long)i);
		result.Rescale(Rescale);
		return result;
	}
};


//masterDir = std::string("/tigress/jaeukk");


// compute ground states of stealthy hyperuniform potentials + soft-core repulsions.
int GetCCO(int argc, char ** argv){
	char tempstring[1000] = {};
	std::string tempstring2;
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;

	size_t timelimit_in_hour = 144;
	size_t beg_idx = 0;
	int seed = 0;
	std::string vtilde = "flat";
	if (argc > 1){
		timelimit_in_hour = (time_t)std::atoi(argv[1]);
		if (argc > 2){
			seed = std::atoi(argv[2]);
		}
		if (argc > 3){
			beg_idx = std::atoi(argv[3]);
		}
		if (argc > 4){
			Verbosity = (size_t) (std::atoi(argv[4]));
		}
		if (argc > 5){
			vtilde = std::string(argv[5]);
		}
		// if (argc > 2){
		// 	tolerance = std::stof(argv[2]);
		// }
		// if (argc > 3){
		// 	max_steps = std::atoi(argv[3]);
		// }
		// if (argc > 7){
		// 	algorithm = std::string(argv[7]);
		// }
	}
	//TODO: Change to use \chi as input parameters
	RandomGenerator rngGod(seed), rng(0);
	ofile << "Time limit for simulations is "<< timelimit_in_hour << " hours\n";
	ofile << "Random seed is "<< seed <<std::endl;
	ofile << "Verbosity is "<< Verbosity <<std::endl;
	ofile << "Shape of potential is "<< vtilde <<std::endl;

	nlopt_srand(999);
	auto start = std::time(nullptr);
	ProgramStart = start;
	TimeLimit = ProgramStart + timelimit_in_hour*3600 - 5*60;//default time limit of 144 hours - 5 mins (for saving progress)
	
	{/* a piece of code to compute the stealthy (hyperuniform) packing at unit number density ... */
		double val, chi = 0, L = 10, sigma, phi = 0.1, angle_deg = 90.0;
		std::vector<double> param_vtilde;	// parameter used in vtilde functions
		std::vector<std::tuple<double,double,double>> shells; // (K1a, deltaa, S0_n) per shell
		DimensionType dim = 2;
		size_t num = L * L*L, numConfig = 300, num_in_load = 0, num_in_Init = 0, num_in_Success = 0;

		int num_threads = 4;
		ofile << "Dimension = "; ifile >> dim;
		size_t M = 1;
		ofile << "M (number of stealthy shells) = "; ifile >> M;
		for (size_t n = 0; n < M; n++) {
			double K1n, deltan, S0n;
			ofile << "K1a[" << n << "] = "; ifile >> K1n;
			ofile << "delta_a[" << n << "] = "; ifile >> deltan;
			ofile << "S0[" << n << "] = "; ifile >> S0n;
			shells.emplace_back(K1n, deltan, S0n);
		}
		ofile << "val = ";	ifile >> val;
		ofile << "sigma = "; ifile >> sigma;
		ofile << "phi = "; ifile >> phi;
		ofile << "box angle (degrees, 90=square, ignored for d!=2) = "; ifile >> angle_deg;
		ofile << "# threads = "; ifile >> num_threads;
		ofile << "num particle = "; ifile >>num;
		ofile << "num configs = "; ifile >> numConfig;

		// For d=2 rhombic box: volume = L^2 * sin(angle), so L = sqrt(N/sin(angle))
		double angle_rad = angle_deg * M_PI / 180.0;
		if (dim == 2)
			L = sqrt((double)num / std::sin(angle_rad));
		else
			L = pow((double)num, 1.0/(double)dim);
		double a = pow(phi/(HyperSphere_Volume(dim, 1.0)), 1.0/(double)dim);
		ofile << "box angle = " << angle_deg << " degrees" << (dim != 2 ? " (ignored, d!=2)" : "") << "\n";

		char mode[100] = {};
		std::string loadname, savename;
		std::function<const Configuration(size_t i)> GetInitConfigs = nullptr;
		ofile << "Initial conditions (random/input)= ";
		ifile >> tempstring;

		if (strcmp (tempstring, "random") == 0){
		/* Random initial conditions */
			GetInitConfigs = [&rngGod, &dim, &num, &L, &angle_rad](size_t i) ->Configuration {
				Configuration pConfig = [&]() {
					if (dim == 2) {
						// Rhombic box: a1=(L,0), a2=(L cos θ, L sin θ)
						std::vector<GeometryVector> basis(2, GeometryVector(2));
						basis[0].x[0] = L;
						basis[1].x[0] = L * std::cos(angle_rad);
						basis[1].x[1] = L * std::sin(angle_rad);
						return Configuration(2, &basis[0], 0.1);
					} else {
						Configuration c = GetUnitCubicBox(dim, 0.1);
						c.Rescale(L);
						return c;
					}
				}();
				for (size_t i = 0; i < num; i++)
					pConfig.Insert("a", rngGod);
				return pConfig;
			};
		}
		else if (strcmp (tempstring, "input") == 0){
		/* Designated initial conditions */
			ofile << "Load a ConfigPack file named as \n";
			ifile >> loadname;
			ReadConfigPack c(loadname, 1.0);
			GetInitConfigs = c;

			num_in_load = c.p.NumConfig();
			ofile << loadname <<".ConfigPack contains " << num_in_load << " realizations\n";
			ofile << "starts from the " << beg_idx << "-th configuration\n";
		}
		else{
			ofile << tempstring << " is undefined \n";
			return 1;
		}
		
		ofile << "Save as "; ifile >> savename;
		ofile << std::endl;
		ofile << "NumParticle = "<<num<<std::endl;
		ofile << "directory = "<< savename<<std::endl;
		ofile << "a = "<< a<<std::endl;		
		ofile << "mode (MD / ground) = "; ifile >> mode;


		/* Define potential */
		ofile << "Define potential \n";
		MultiShellS0Potential * potential = nullptr;
		{
			// Convert shells from encoded (Ka, deltaa) units to absolute k-space (k1, k2, S0_n)
			std::vector<std::tuple<double,double,double>> shells_abs; // (k1, k2, S0_n)
			double k2_max = 0;
			for (size_t n = 0; n < M; n++) {
				double k1n = std::get<0>(shells[n]) / a;
				double k2n = k1n + std::get<1>(shells[n]) / a;
				double S0n = std::get<2>(shells[n]);
				shells_abs.emplace_back(k1n, k2n, S0n);
				k2_max = std::max(k2_max, k2n);
				ofile << "Shell " << n << ": k1=" << k1n << ", k2=" << k2n << ", S0=" << S0n << "\n";
			}

			Configuration Config = GetInitConfigs(beg_idx);
			potential = new MultiShellS0Potential(dim, val, sigma);
			ofile << "With repulsion, val = " << val << ", sigma = " << sigma << "\n";

			potential->ParallelNumThread = num_threads;
			potential->CCPotential1->ParallelNumThread = num_threads;
			potential->CCPotential2->ParallelNumThread = num_threads;

			// vtilde function applies to stealthy shells (S0=0); equiluminous shells always use V=1
			std::function<double(double)> v;
			if (vtilde.compare("flat") == 0) {
				ofile << "vtilde function = flat (by default)\n";
				v = [](double k) -> double { return 1.0; };
			}
			else if (vtilde.compare("overlap") == 0) {
				ofile << "vtilde function = overlap function\n";
				param_vtilde.emplace_back(k2_max);
				v = [&dim, &param_vtilde](double k) -> double { return alpha(dim, k/(param_vtilde[0])); };
			}
			else if (vtilde.compare("power-law") == 0) {
				double m;
				ofile << "vtilde function = power-law function\n";
				ofile << "exponent? = "; ifile >> m;
				param_vtilde.emplace_back(k2_max);
				param_vtilde.emplace_back(m);
				v = [&param_vtilde](double k) -> double { return std::pow(1. - k/param_vtilde[0], param_vtilde[1]); };
			}

			// Route each k-point to CCPotential1 (stealthy) or CCPotential2 (equiluminous)
			std::vector<GeometryVector> ks_temp = GetKs(Config, k2_max, k2_max, 1);
			for (auto k = ks_temp.begin(); k != ks_temp.end(); k++) {
				double km2 = k->Modulus2();
				for (auto& s : shells_abs) {
					double k1n = std::get<0>(s), k2n = std::get<1>(s), S0n = std::get<2>(s);
					if (km2 > k1n*k1n && km2 <= k2n*k2n) {
						if (S0n == 0.0) {
							std::vector<double> vals = {v(std::sqrt(km2))};
							potential->CCPotential1->AddConstraint(*k, vals);
							if (Verbosity > 3)
								ofile << std::sqrt(km2) << "\t" << vals[0] << "\n";
						} else {
							std::vector<double> vals = {1.0, S0n};
							potential->CCPotential2->AddConstraint(*k, vals);
						}
						chi++;
						break;
					}
				}
			}
			chi /= dim * (num - 1);
			ofile << "chi = " << chi << std::endl;
		}



		if (strcmp (mode, "MD") == 0)
		{
			/* NVT MD simulations */
			double MDTimeStep = 0.01, MD_Temperature = -1.0;
			size_t MDStepPerSample = 100000; 
			size_t numEquilSamples = 500;
			bool MDAutoTimeStep = false;
			bool MDAllowRestore = false;
			{
				ConfigurationPack save_(savename);				
				num_in_Init = save_.NumConfig();
			}

			/* Input parameters */
			for(;;){
				ifile >> tempstring2;
				std::transform(tempstring2.begin(), tempstring2.end(), tempstring2.begin(),::tolower); // make the input in the lower case
				if (tempstring2.compare("run") == 0){
					/* start MD simulation */
					break;
				}
				else if (tempstring2.compare("timestep") == 0){
					ofile << "\ndt = ";
					ifile >> MDTimeStep ;
				}
				else if (tempstring2.compare("samplesteps") == 0){
					ofile << "\nSample config. per \'n\' dt   ";
					ofile << "n = ";
					ifile >> MDStepPerSample ; 
				}
				else if (tempstring2.compare("temperature") == 0){
					ofile << "\nMD temperature [in the unit of v0] = ";
					ifile >> MD_Temperature ; 
					ofile << MD_Temperature <<"\n\n";
				}
				else if (tempstring2.compare("numequilibration") == 0){
					ofile << "\nThe number sampling steps to equilibrate samples = ?";
					ifile >> numEquilSamples;
				}
				else if (tempstring2.compare("autotimestep")== 0){
					ofile << "\nTime step will be determined automatically\n";
					MDAutoTimeStep = true;
				}
				else if (tempstring2.compare("restore")== 0){
					ofile << "\nSave and load the progress of MD simulations.\n";
					MDAllowRestore = true;
				}
				else {
					ofile << tempstring2 << " is an undefined command\n";
				}
			}

			if (MD_Temperature < 0.0){
				ofile << "\nTemperature is undefined;\n";
				if (dim == 1){
					MD_Temperature = 2e-4;
				}
				else if (dim == 2){
					MD_Temperature = 2e-6;
				}
				else if (dim == 3){
					MD_Temperature = 1e-6;
				}
				ofile << "Take T to be a default value, "<< MD_Temperature << "\n";
			}
			else{
				ofile << "Input T = "<< MD_Temperature << "\n";
			}
			if (strcmp (tempstring, "random") == 0){
				/* Start from a random initial condition and basic properties */
				
				Configuration pConfig = GetInitConfigs(0);
				
				size_t MDEquilibrateSamples = numEquilSamples;
				//bool MDAutoTimeStep = false;
				CollectiveCoordinateMD(&pConfig, potential, rngGod, MDTimeStep, MD_Temperature, savename, numConfig, MDStepPerSample, MDAllowRestore, TimeLimit, MDEquilibrateSamples, MDAutoTimeStep);
			}
			else if	(strcmp (tempstring, "input") == 0){
				/* For continue from the latest configuration in the previous run. */
				ofile << "MD Temperature = "<< MD_Temperature << std::endl;
			
				Configuration pConfig = GetInitConfigs(num_in_load - 1);
				//bool MDAllowRestore = true;
				size_t MDEquilibrateSamples = numEquilSamples;
				size_t numConfig_comp = (numConfig > num_in_Init)? numConfig - num_in_Init : 0 ;
				
				//bool MDAutoTimeStep = false;
				CollectiveCoordinateMD(&pConfig, potential, rngGod, MDTimeStep, MD_Temperature, savename, numConfig_comp, MDStepPerSample, MDAllowRestore, TimeLimit, MDEquilibrateSamples, MDAutoTimeStep);
			//TODO
			}
			
			
			
		}
		else if (strcmp (mode, "ground") == 0)
		{
			/* Quench to zero temperature */
			/* Parameters for minimizations */
			size_t max_steps = 10000;
			double tolerance = 1e-14;
			std::string algorithm = "LBFGS";

			ofile << "--------------------------------------- \n";
			ofile << "\t\t Input parameters for minimizations \n";
			ofile << "--------------------------------------- \n";
			for(;;){
				ifile >> tempstring2;
				std::transform(tempstring2.begin(), tempstring2.end(), tempstring2.begin(),::tolower); // make the input in the lower case
				if (tempstring2.compare("run") == 0){
					/* start minimizations simulation */
					break;
				}
				else if (tempstring2.compare("tolerance") == 0){
					ofile << "\n energy tolerance for ground states [in the unit of v0] = ";
					ifile >> tolerance ;
				}
				else if (tempstring2.compare("maxsteps") == 0){
					ofile << "\n Limit number of evaluations = ";
					ifile >> max_steps ;
				}
				else if (tempstring2.compare("algorithm") == 0){
					ofile << "\nMinimization algorithm = ";
					ifile >> algorithm ; 
				}
				else {
					ofile << tempstring2 << " is an undefined command\n";
				}
			}
			ofile << "\nNumerical optimizer is "<< algorithm <<std::endl;
			ofile << "Energy tolerance is "<< tolerance <<std::endl;
			ofile << "Max. steps of evaluations is "<< max_steps <<std::endl;


			if (strcmp (tempstring, "random") == 0){
				/* random initial conditions 
					=> use MultiRun */

				Configuration pConfig = GetInitConfigs(0);
				pConfig.PrepareIterateThroughNeighbors(sigma);

				{
					ConfigurationPack save_(savename+"_Success");				
					num_in_Success = save_.NumConfig();
				}
				size_t numConfig_comp = (numConfig > num_in_Success )? numConfig - num_in_Success : 0;

				ofile << "We already have " << num_in_Success << " Configurations in the savefile \n";
				ofile << "just compute "<< numConfig_comp <<" more Configurations \n";

				ofile << "max_eval = " << max_steps << std::endl;
				CollectiveCoordinateMultiRun(&pConfig, potential, rngGod, savename, numConfig_comp, TimeLimit, algorithm, tolerance, max_steps);
				delete potential;

			}
			else if	(strcmp (tempstring, "input") == 0){
				/* ReadConfigPack */
				size_t idx_beg = 0;
				{
					ConfigurationPack load_(savename+"_InitConfig");				
					num_in_Init = load_.NumConfig();
					ConfigurationPack success_(savename+"_InitConfig");				
					num_in_Success = success_.NumConfig();
				}
				size_t numConfig_upper = std::min(num_in_load, beg_idx+numConfig);

				ofile << "We start from the " << beg_idx << "-th Configuration from loadfile\n";  
				ofile << "Among them, we already have used " << num_in_Init << " Configurations  \n";
				num_in_Init += beg_idx;
				ofile << "compute from the "<< num_in_Init << "th Config. to the " << numConfig_upper <<"th Config.\n";
				size_t id = num_in_Init;
				for (size_t i = num_in_Success; i < numConfig; i++){
					if ( id < num_in_load ){
						Configuration pConfig = GetInitConfigs(id);
						pConfig.PrepareIterateThroughNeighbors(sigma);
						CollectiveCoordinateMultiRun(&pConfig, potential, rngGod, savename, 1, TimeLimit, algorithm, tolerance, max_steps);
					}
					id ++ ;
				}
				delete potential;				
			}

		}
		else{
			ofile << "mode "<< mode << " is undefined!\n";
			return 1;
		}
		


		


		// Some basic analyses
		/*Compute the nearest distance of each config. */
		if (false){
			ConfigurationPack ConfigSet;
			if (strcmp (mode, "MD") == 0){
				ConfigSet.Open((savename ).c_str());
			}
			else if (strcmp (mode, "ground") == 0){
				ConfigSet.Open((savename + "_Success").c_str());
			}
			double rmin = L;
			{
				for (int i = 0; i < ConfigSet.NumConfig(); i++) {
					Configuration x = ConfigSet.GetConfig(i);
					double rm = L;
					for (size_t j = 0; j < x.NumParticle(); j++) {
						double temp = x.NearestParticleDistance(j);
						rm = (temp < rm) ? temp : rm;
					}
					ofile << i << ":\t" << rm << std::endl;
				}
			}

			auto getC = [&rmin, &ConfigSet](size_t i) ->Configuration {
				Configuration c = ConfigSet.GetConfig(i);
				double R = rmin;
				for (size_t j = 0; j < c.NumParticle(); j++) {
					double temp = c.NearestParticleDistance(j);
					R = (temp < R) ? temp : R;
				}
				rmin = (R < rmin) ? R : rmin;
				return c;
			};

			/*Compute g2 function*/
			std::vector<GeometryVector> g2;
			numConfig = ConfigSet.NumConfig();
			IsotropicTwoPairCorrelation(getC, numConfig, std::min(0.2 * L, 5.0) , g2);
			WriteFunction(g2, (savename + "_g2").c_str());
			phi = 1.0*HyperSphere_Volume(dim, rmin / 2.0);
			ofile << "max packing fraction is " << phi << std::endl;

			/*Compute the structure factsor*/
			std::vector<GeometryVector> Sk;
			IsotropicStructureFactor(getC, numConfig, 0.2, 1, Sk);
			WriteFunction(Sk, (savename + "_Sk").c_str());
		}		
	}

	return 0;
}

int MC_CCO(int argc, char ** argv){
	/* common input */
	char tempstring[1000] = {};
	std::string tempstring2;
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;

	size_t timelimit_in_hour = 144;
	size_t beg_idx = 0;
	int seed = 0;
	std::string vtilde = "flat";
	if (argc > 1){
		timelimit_in_hour = (time_t)std::atoi(argv[1]);
		if (argc > 2){
			seed = std::atoi(argv[2]);
		}
		if (argc > 3){
			beg_idx = std::atoi(argv[3]);
		}
		if (argc > 4){
			Verbosity = (size_t) (std::atoi(argv[4]));
		}
		if (argc > 5){
			vtilde = std::string(argv[5]);
		}
	}
	RandomGenerator rngGod(seed), rng(0);
	ofile << "Time limit for simulations is "<< timelimit_in_hour << " hours\n";
	ofile << "Random seed is "<< seed <<std::endl;
	ofile << "Verbosity is "<< Verbosity <<std::endl;
	ofile << "Shape of potential is "<< vtilde <<std::endl;
	auto start = std::time(nullptr);
	ProgramStart = start;
	TimeLimit = ProgramStart + timelimit_in_hour*3600 - 5*60;//default time limit of 144 hours - 5 mins (for saving progress)

	/* parameters for potential */
	double K1, chi_target, val, L, sigma, phi = 0.1, S0=0.0;
	std::vector<double> param_vtilde;	// parameter used in vtilde functions
	DimensionType dim = 2;
	size_t num; 

	int num_threads = 4;
	ofile << "Dimension = "; ifile >> dim;
	ofile << "K1 = ";	ifile >> K1;
	ofile << "chi_target = ";	ifile >> chi_target;
	//ofile << "K2a = ";	ifile >> K2;
	ofile << "S0 = ";	ifile >> S0;
	ofile << "val = ";	ifile >> val;
	ofile << "sigma = "; ifile >> sigma;
	ofile << "# threads = "; ifile >> num_threads;
	ofile << "num particle = "; ifile >>num; 
	L = pow(num, 1.0/(double)dim);	

	/* Define initial conditions */
	size_t numConfig = 300, num_in_load = 0, num_in_Init = 0;
	char mode[100] = {};
	std::string loadname, savename;		
	std::function<const Configuration(size_t i)> GetInitConfigs = nullptr;

	{
		ofile << "num configs = "; ifile >> numConfig;
		ofile << "NumParticle = "<<num<<std::endl;

		ofile << "Initial conditions (random/input)= ";
		ifile >> tempstring;
		ofile << "Save as "; ifile >> savename;
		ofile << std::endl;
		ofile << "directory = "<< savename<<std::endl;

		if (strcmp (tempstring, "random") == 0){
		/* Randomg initial conditions */
			GetInitConfigs = [&rngGod, &dim, &num, &L](size_t i) ->Configuration {
				Configuration pConfig = GetUnitCubicBox(dim, 0.1);
				pConfig.Rescale(L);
				for (size_t i = 0; i < num; i++)
					pConfig.Insert("a", rngGod);
				return pConfig;
			};
		}
		else if (strcmp (tempstring, "input") == 0){
		/* Designated initial conditions */
			ofile << "Load a ConfigPack file named as \n";
			ifile >> loadname;
			ReadConfigPack c(loadname, 1.0);
			GetInitConfigs = c;
			
			num_in_load = c.p.NumConfig();
			ofile << loadname <<".ConfigPack contains " << num_in_load << " realizations\n";
			ofile << "starts from the " << beg_idx << "-th configuration\n";
		}
		else{
			ofile << tempstring << " is undefined \n";
			return 1;
		}
	}

	/* Define potential */
	//forMCPotential * pPot = nullptr;
	ShiftedCCPotential_ForMC * pCCPot = nullptr;
	{
		Configuration c = GetInitConfigs(beg_idx);
		double rho = c.NumParticle() / c.PeriodicVolume();
		/* Find a list of constrained wavevectors */
		double K2 = 1.1 * getK(chi_target + get_chi(K1, rho, dim), rho, dim);
		double K1_squared = K1*K1;
		std::vector<GeometryVector> ks_temp = GetKs(c, K2, K2, 1.0);
		std::sort(ks_temp.begin(), ks_temp.end(), [](const GeometryVector & left, const GeometryVector & right) ->bool {return (left.Modulus2()<right.Modulus2()) ; });
		ks_temp.erase(std::remove_if(ks_temp.begin(), ks_temp.end(), [&K1_squared](const GeometryVector & k) { return k.Modulus2() < K1_squared; }), ks_temp.end());
		ofile << "test = " << ks_temp.size() << std::endl;
		/* Adjust the number of constraints to be closest to the target chi value. */
		{
			size_t expected_num_constraints = (size_t) (dim * c.NumParticle() * chi_target );
			double K_search_squared = ks_temp[expected_num_constraints-1].Modulus2();
			size_t i1 = expected_num_constraints-1, i2 = i1;  // i1 < expected_num_constraints < i2.
			while (i1 > 0 && ks_temp[i1].Modulus2() == K_search_squared) i1--;
			while (i2 < ks_temp.size() && ks_temp[i2].Modulus2() == K_search_squared) i2++;
			ofile << "i1, target, i2 = " << i1 << "\t" << expected_num_constraints <<"\t" << i2 << std::endl;
			if (expected_num_constraints - i1 < i2 - expected_num_constraints){
				// i1 is the closest
				ks_temp.erase(ks_temp.begin() + i1 + 1, ks_temp.end());
			}
			else{
				// i2 is the closest
				ks_temp.erase(ks_temp.begin() + i2 + 1, ks_temp.end());
			}			
			K2 = std::sqrt(ks_temp.end()->Modulus2());
		}
		double chi_exact = ks_temp.size() / (double)(dim * (num - 1));
		ofile << "chi_exact = " << chi_exact << ",\tchi_target = " << chi_target << std::endl;

		if (S0 == 0.0)
		{
			ShiftedCCPotential pot_temp(dim);
			pot_temp.ParallelNumThread = num_threads;
			std::function<double (double k)> v;
			std::vector<double> param_vtilde;
			std::vector<double> vals;	vals.emplace_back(0);
			{
				if (vtilde.compare("flat") == 0){
					ofile << "vtilde function = flat (by default)\n";
					v = [](double k) -> double { return 1.0; };
				}
				else if (vtilde.compare("overlap") == 0){
					ofile << "vtilde function = overlap function\n";
					param_vtilde.emplace_back(K2);
					v = [&dim, &param_vtilde](double k) -> double { return alpha(dim, k/(param_vtilde[0]));};
				}
				else if (vtilde.compare("power-law") == 0){
					double m;
					ofile << "vtilde function = power-law function \n";
					ofile << "exponent? = ";	ifile >> m;	
					param_vtilde.emplace_back(K2);
					param_vtilde.emplace_back(m);
					v = [&param_vtilde](double k) -> double { return std::pow(1. - k/param_vtilde[0], param_vtilde[1]);};
				}
			
				for (auto k = ks_temp.begin(); k != ks_temp.end(); k++)
				{
					vals[0] = v(std::sqrt(k->Modulus2()));
					pot_temp.AddConstraint(*k, vals);

					if (Verbosity > 3){
						ofile << std::sqrt(k->Modulus2()) << "\t" << vals[0] <<"\n";
					}
					
				}
			}

			pCCPot = new ShiftedCCPotential_ForMC(pot_temp);	
		}
		else
		{
			//TODO: implement equiluminous case.
			return 0;
		}

	}

	/* Annealing */
	double T_init = 1e-2, T_fin = 1e-14, cooling_rate = 0.95;
	double tolerance = 1e-20;
	double stepsize = 1e-1, reduce = 0.9, incr = 1.01, stepsize_init;
	size_t equi_steps = 10000 , adjust_steps = 100, sampling_steps = 0;
	size_t samples_at_T_fin = 1;
	double acc_min = 0.25, acc_max = 0.35;
	//RandomGenerator rng (seed);

	// ofile << "stepsize = "; ifile >> stepsize ;
	// ofile << "steps for adjusting step size = "; ifile >> adjust_steps;
	// ofile << "steps for equilibration = "; ifile >> equi_steps;

	/* Input parameters */
	size_t max_inputs = 100;
	for(size_t id=0; id < max_inputs; id ++ ){
//		ofile << id << "--------------------------------------- \n";
		ifile >> tempstring2;
		std::transform(tempstring2.begin(), tempstring2.end(), tempstring2.begin(),::tolower); // make the input in the lower case
		if (tempstring2.compare("run") == 0){
			/* start MC simulation */
			break;
		}
		else if (tempstring2.compare("stepsize") == 0){
			ofile << "\n initial step size = ";
			ifile >> stepsize ;
		}
		else if (tempstring2.compare("change_stepsize") == 0){
			ofile << "\nMC increase/decrease step size by factors = ";
			ifile >> incr;
			ifile >> reduce;
		}
		else if (tempstring2.compare("equi_steps") == 0){
			ofile << "\nSample config. per \'n\' MC cycles   ";
			ofile << "n = ";
			ifile >> equi_steps; 
		}
		else if (tempstring2.compare("adjust") == 0){
			ofile << "\nMC cycles for adjusting step sizes = ?";
			ifile >> adjust_steps;
		}
		else if (tempstring2.compare("t_init") == 0){
			ofile << "\n Initial MC temperature [in the unit of v0] = ";
			ifile >> T_init ; 
			ofile << T_init <<"\n\n";
		}
		else if (tempstring2.compare("t_fin") == 0){
			ofile << "\n Final MC temperature [in the unit of v0] = ";
			ifile >> T_fin ; 
			ofile << T_fin <<"\n\n";
		}
		else if (tempstring2.compare("coolingrate")== 0){
			ofile << "\nCooling rate = \n";
			ifile >> cooling_rate;
			ofile << cooling_rate << "\n\n";
		}
		else if (tempstring2.compare("acceptance")== 0){
			ofile << "\ntarget acceptance ratio (min/max) = \n";
			ifile >> acc_min;	ifile >> acc_max;
			ofile << acc_min <<"\t" << acc_max << "\n\n";
		}
		else if (tempstring2.compare("samples_tfin")== 0){
			ofile << "\nSamples at the final temperature = \n";
			ifile >> samples_at_T_fin;
			ofile << samples_at_T_fin << "\n\n";
			ofile << "sampling_steps = ";
			ifile >> sampling_steps;
			ofile << sampling_steps << "\n\n";  
		}	
		else {
			ofile << tempstring2 << " is an undefined command\n";
			break;
		}
	}


	ConfigurationPack AfterRelaxPack(savename);
	ConfigurationPack SuccessPack(savename + "_Success");
	
	for (size_t i = 0; i < numConfig; i++ ){
		Configuration result = GetInitConfigs(i);
		pCCPot -> SetConfiguration(result);
		double E_curr = pCCPot -> Energy();
		double T_curr = T_init; 
		if (i > 0){
			/* load stepsize in the previous cooling */
			stepsize = stepsize_init;
		}
		/* cooling */
		ofile << "config " << i <<": cooling start, E_init = " << E_curr << "\n";
		while (T_curr > T_fin && E_curr > tolerance){
			ofile << "@T=" << T_curr << ",\t";
			/* adjust displacement */
			double ACCEPTANCE = 0;
			stepsize *= reduce;
			while (ACCEPTANCE > acc_max || ACCEPTANCE < acc_min)
			{	
				ACCEPTANCE = 0.0;
				for(size_t j = 0; j < adjust_steps; j++){
				
					for (size_t k = 0; k < result.NumParticle(); k++){
						GeometryVector dx (dim);
						for (size_t l = 0;  l < dim ; l++){
							dx.x[l] = stepsize*(rngGod.RandomDouble()-0.5);
						}
						GeometryVector prevCart = result.GetCartesianCoordinates(k);
						GeometryVector afterCart = prevCart + dx ;

						double dE = pCCPot -> TryMove( prevCart, afterCart);
						bool accept = (dE < 0.0) || (rngGod.RandomDouble() < std::exp(-dE/T_curr));
						//ofile << dE << " : " << std::exp(-dE/T_curr) << "\n";
						if (accept){
							pCCPot->AcceptMove();
							/* update coordinates */
							GeometryVector afterRelative = result.CartesianCoord2RelativeCoord(afterCart);
							result.MoveParticle(k, afterRelative);
							E_curr += dE;

							ACCEPTANCE ++ ;
						}
					}

				}
				ofile << "-";
				double E_actual = pCCPot -> Energy();
				ACCEPTANCE /= (double)adjust_steps * result.NumParticle();
				if (Verbosity > 3){
					ofile << "stepsize = " << stepsize << ", ";
					ofile << "err_E = " << std::abs(E_curr-E_actual) << ", ";
					ofile << "acceptance ratio = " << ACCEPTANCE << std::endl;
				}
				E_curr = E_actual;

				if (ACCEPTANCE < acc_min){
					stepsize *= reduce;
				}
				if (ACCEPTANCE > acc_max){
					stepsize *= incr;
				}
			}
			if (T_curr == T_init){
				stepsize_init = stepsize;
			}

			if (Verbosity > 2){
				/* determined parameters */
				ofile << "(stepsize=" << stepsize << ", ";
				ofile << "curent E=" << E_curr << ", ";
				ofile << "acceptance ratio=" << ACCEPTANCE << ") ";
			}

			/* equilibration */
			double E1 = 0.0, E2 =0.0, count  = 0;
			for(size_t j = 0; j < equi_steps; j++){
				
				for (size_t k = 0; k < result.NumParticle(); k++){
					GeometryVector dx (dim);
					for (size_t l = 0;  l < dim ; l++){
						dx.x[l] = stepsize*(rngGod.RandomDouble()-0.5);
					}
					GeometryVector prevCart = result.GetCartesianCoordinates(k);
					GeometryVector afterCart = prevCart + dx ;

					double dE = pCCPot -> TryMove( prevCart, afterCart);
					bool accept = (dE < 0.0) || (rngGod.RandomDouble() < std::exp(-dE/T_curr));
					if (accept){
						pCCPot->AcceptMove();
						/* update coordinates */
						GeometryVector afterRelative = result.CartesianCoord2RelativeCoord(afterCart);
						result.MoveParticle(k, afterRelative);
						E_curr += dE;

						E1 += E_curr;
						E2 += E_curr*E_curr;
						count++;
						if ( E_curr < tolerance ){
							break ;
						}
					}
				}

				if ( j % 500 == 99){
					ofile << "*";
					double E_actual = pCCPot -> Energy();
					if (Verbosity > 3){
						ofile << std::abs(E_curr-E_actual) ;
					}
					E_curr = E_actual;
				}
			}

			{
				/* intermeidate statistics */
				double Emean = E1 / count; 
				double Evar = std::sqrt(E2 / count - Emean*Emean);
				ofile << "\t E=" << Emean << " (" << Evar << ", " << count << ")" << std::endl;
			}

			T_curr *= cooling_rate;
		}


		/* sampling at T = T_fin*/
		ofile << "sampling at the final temperature, T=" << T_curr << "\n";

		for (size_t j = 0; j < samples_at_T_fin ; j++){
			ofile << ".";
			for (size_t l = 0 ; l < sampling_steps; l++){

				for (size_t k = 0; k < result.NumParticle(); k++){
					GeometryVector dx (dim);
					for (size_t l = 0;  l < dim ; l++){
						dx.x[l] = stepsize*(rngGod.RandomDouble()-0.5);
					}
					GeometryVector prevCart = result.GetCartesianCoordinates(k);
					GeometryVector afterCart = prevCart + dx ;

					double dE = pCCPot -> TryMove( prevCart, afterCart);
					bool accept = (dE < 0.0) || (rngGod.RandomDouble() < std::exp(-dE/T_curr));
					if (accept){
						pCCPot->AcceptMove();
						/* update coordinates */
						GeometryVector afterRelative = result.CartesianCoord2RelativeCoord(afterCart);
						result.MoveParticle(k, afterRelative);
						E_curr += dE;
					}
				}
			}

			ofile << "E_relax="<< E_curr << "\n";
			AfterRelaxPack.AddConfig(result);
			if (E_curr < tolerance ){
				ofile << "Add to success pack\n";
				SuccessPack.AddConfig(result);
			}

		}
		
	}


	delete pCCPot;

	return 0;
}

// int CollectiveCoordinateMD(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, double TimeStep, double Temperature, std::string Prefix, size_t SampleNumber, size_t StepPerSample, bool AllowRestore, time_t TimeLimit, size_t EquilibrateSamples, bool MDAutoTimeStep);

/**
 * @brief CLI of multi-stealthy hyperuniform packings
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int GetMultiCCO(int argc, char ** argv){

	char tempstring[1000] = {};
	std::string tempstring2;
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;

	size_t timelimit_in_hour = 144;
	size_t beg_idx = 0;
	int seed = 0;
	std::string vtilde = "flat";
	if (argc > 2){
		timelimit_in_hour = (time_t)std::atoi(argv[2]);
		if (argc > 3){
			seed = std::atoi(argv[3]);
		}
		if (argc > 4){
			beg_idx = std::atoi(argv[4]);
		}
		if (argc > 5){
			Verbosity = (size_t) (std::atoi(argv[5]));
		}
		if (argc > 6){
			vtilde = std::string(argv[6]);
		}
	}
	RandomGenerator rngGod(seed), rng(0);
	ofile << "Time limit for simulations is "<< timelimit_in_hour << " hours\n";
	ofile << "Random seed is "<< seed <<std::endl;
	ofile << "Verbosity is "<< Verbosity <<std::endl;
	ofile << "Shape of potential is "<< vtilde <<std::endl;
	
	nlopt_srand(999);
	auto start = std::time(nullptr);
	ProgramStart = start;
	TimeLimit = ProgramStart + timelimit_in_hour*3600 - 5*60;//default time limit of 144 hours - 5 mins (for saving progress)
	
	/* 1. define sizes of a packing */
	DimensionType dim = 2;
	double val, L = 10, phi = 0.1, rho;
	std::vector<double> radii;
	std::vector<size_t> sizes;
	size_t num_species = 1, num_radii = 1;
	size_t N = L * L*L;
	std::string InitCondition;
	std::function<const SpherePacking(size_t i)> GetInitConfigs = nullptr;
	std::string loadname, savename;
	{
		ofile << "Define an initial packing (random / input)\n";
		ifile >> InitCondition;
		if (strcmp(InitCondition.c_str(), "random") == 0){
			/* random initial conditions */

			ofile << "Dimension = "; 		ifile >> dim;
			ofile << "Number of different radii = ";	ifile >> num_radii; 
			radii.reserve(num_radii);
			sizes.reserve(num_radii);
			ofile << "Input the number of particles and the relative radius \n";
			N = 0;
			for (int i=0; i < num_radii; i++){
				ofile << "number: "; ifile >> tempstring; 
				sizes.emplace_back(atoi(tempstring));
				N += atoi(tempstring);

				ofile << "relative radius: "; ifile >> tempstring;
				radii.emplace_back(atof(tempstring));
			}
			ofile << "Target packing fraction = "; 	ifile >> phi;
			
			/* scale particle radii */
			ofile << "-- Particle number and radius --\n";
			{

				double curr_phi = 0;
				for (int i=0; i < num_radii; i++){
					curr_phi += sizes[i] * HyperSphere_Volume(dim, radii[i]);
				}
				curr_phi /= N;

				double factor = std::pow( phi/curr_phi , 1.0/dim);
				for (int i=0; i < num_radii; i++){
					radii[i] *= factor;
					ofile << sizes[i] << "\t" << radii[i] <<"\n";
				}
				L = pow(N, 1./dim);				
			}		

			/* a function to generate a random configuration of spheres. */
			GetInitConfigs = [&rngGod, &dim, &sizes, &radii, L](size_t idx) ->SpherePacking {
				SpherePacking pConfig(GetUnitCubicBox(dim, 0.1),0);
				pConfig.Rescale(L);
				for (size_t i = 0; i < sizes.size(); i++){					
					for (size_t j = 0; j < sizes[i]; j++){
						GeometryVector x;
						RandomVector(dim, rngGod, x);
						pConfig.Insert(radii[i], x);
					}
				}
				return pConfig;
			};
		}
		else if (strcmp(InitCondition.c_str(), "input") == 0){
			ofile << "Load a txt file for a SpherePacking (include a %d symbol to identify a configuration index) ";
			ifile >> loadname;

			GetInitConfigs = [&loadname](size_t idx) -> SpherePacking{
				char name[300];
				std::sprintf(name, loadname.c_str(), idx);
				SpherePacking temp; 
				ReadConfiguration(temp, name);
				return temp;
			};

			// ReadConfigPack c(loadname, 1.0);
			// GetInitConfigs = c;
			
			// num_in_load = c.p.NumConfig();
			// ofile << loadname <<".ConfigPack contains " << num_in_load << " realizations\n";
			// ofile << "starts from the " << beg_idx << "-th configuration\n";
		}
		else{
			ofile << tempstring << " is undefined \n";
			return 1;
		}
	}
	ofile << "Input the relative strength of contact potential =  ";	ifile >> val; 	ofile << val << std::endl;
	
	/* 2. Define stealthy potential */
	ofile << "Define a stealthy potential \n";
	MultiStealthyPotential * potential = new MultiStealthyPotential(dim, val);
	{
		ofile << "number of different collective potentials = ";
		ifile >> num_species;	ofile << " (" << num_species <<" )\n";

		SpherePacking c = GetInitConfigs(beg_idx);
		N = c.NumParticle();
		rho = N / c.PeriodicVolume();
		potential -> SetPacking(c);

		potential -> SetNumSpecies(num_species);
		int N_curr = 0;
		for (size_t s = 0; s < num_species; s++){
			size_t N_s;
			
			ofile << "species [" << s << "]: particle number = ";
			ifile >> N_s;

			if(N_s + N_curr > N){
				ofile << "The input is larger than the unspecified particle number. All unspecified particles are taken to be a new species.\n";
				N_s = (N-N_curr);
			}

			for (size_t i = 0; i < N_s ; i++){
				potential->SetSpeciesOfParticle(N_curr + i, s);
				ofile << N_curr + i <<"\n";
			}
			N_curr += N_s;

			double lower, upper, S0, chi;
			std::vector<double> param_vtilde;	// parameter used in vtilde functions

			ofile << "S(0) = ";	ifile >> S0;	ofile << S0 << std::endl;
			ofile << "How to define the constrained region? (chi / wavenumber) = ";	ifile >> tempstring;	ofile << tempstring << std::endl;
			ofile << "lower bound = ";	ifile >> lower;	ofile << lower <<"\n";
			ofile << "upper bound = ";	ifile >> upper;	ofile << upper <<"\n";

			/* convert to the wavenumbers */
			if (strcmp(tempstring, "chi")==0){
				ofile << "Note that the values of chi here are normalized by the total number of particles.\n";
				lower = getK(lower, rho, dim);
				upper = getK(upper, rho, dim);
			}

			/* define functional form. */
			std::function<std::vector<double> (double k)> v;
			if (S0 == 0.0)
			{	
				ofile << "Stealthy:\n";
				if (vtilde.compare("flat") == 0){
					ofile << "vtilde function = flat (by default)\n";
					v = [](double k) -> std::vector<double> { return std::vector<double>({1.0}); };
				}
				else if (vtilde.compare("overlap") == 0){
					ofile << "vtilde function = overlap function\n";
					param_vtilde.emplace_back(upper);
					v = [&dim, &param_vtilde](double k) -> std::vector<double> { return std::vector<double>({alpha(dim, k/(param_vtilde[0]))}); };
				}
				else if (vtilde.compare("power-law") == 0){
					double m;
					ofile << "vtilde function = power-law function \n";
					ofile << "exponent? = ";	ifile >> m;	
					param_vtilde.emplace_back(upper);
					param_vtilde.emplace_back(m);
					v = [&param_vtilde](double k) -> std::vector<double> { return std::vector<double>({std::pow(1. - k/param_vtilde[0], param_vtilde[1])}); };
				}
			}
			else{
				ofile << "equiluminous patterns with S0 = " << S0 << "\n";
				v = [S0](double k) -> std::vector<double> { return std::vector<double>(1.0, S0); }; 
			}

			double k1_modulus = lower*lower;
			std::vector<GeometryVector> ks_temp = GetKs(c, upper, upper, 1);
			for (auto k = ks_temp.begin(); k != ks_temp.end(); k++) {
				if (k->Modulus2() > k1_modulus) {
					potential->AddConstraint(s,*k, v(std::sqrt(k->Modulus2())));
					chi++;
				}					
			}
			chi /= dim * (N_s - 1);				
			ofile << "chi (of a subconfig) = " << chi << std::endl;			
		}
	}
		
	size_t N_config = 300, N_load = 0, N_init = 0, N_success = 0;
	int num_threads = 4;
	
	ofile << "# threads = "; ifile >> num_threads;	ofile << "\t "<< num_threads << "\n";
	potential->ParallelNumThread = num_threads;
	ofile << "Save as "; ifile >> savename; ofile << "\t "<< savename << "\n";
	ofile << "target num configs = "; ifile >> N_config; ofile << "\t "<< N_config << "\n";

	//ofile << std::endl;
	//ofile << "NumParticle = "<<num<<std::endl;
	//potential->CCPotential->ParallelNumThread = num_threads;
	
	{/* a piece of code to compute the stealthy (hyperuniform) packing at unit number density ... */
		
		char mode[100] = {};
		ofile << "directory = "<< savename<<std::endl;
		ofile << "mode (MD / ground) = "; ifile >> mode; ofile << std::endl;
		
		if (strcmp (mode, "MD") == 0)
		{
			/* NVT MD simulations */
			double MDTimeStep = 0.01, MD_Temperature = -1.0;
			size_t MDStepPerSample = 100000; 
			size_t numEquilSamples = 500;
			bool MDAutoTimeStep = false;
			bool MDAllowRestore = false;
			{
				ConfigurationPack save_(savename);				
				N_init = save_.NumConfig();
			}

			/* Input parameters */
			for(;;){
				ifile >> tempstring2;
				std::transform(tempstring2.begin(), tempstring2.end(), tempstring2.begin(),::tolower); // make the input in the lower case
				if (tempstring2.compare("run") == 0){
					/* start MD simulation */
					break;
				}
				else if (tempstring2.compare("timestep") == 0){
					ofile << "\ndt = ";
					ifile >> MDTimeStep ;
				}
				else if (tempstring2.compare("samplesteps") == 0){
					ofile << "\nSample config. per \'n\' dt   ";
					ofile << "n = ";
					ifile >> MDStepPerSample ; 
				}
				else if (tempstring2.compare("temperature") == 0){
					ofile << "\nMD temperature [in the unit of v0] = ";
					ifile >> MD_Temperature ; 
					ofile << MD_Temperature <<"\n\n";
				}
				else if (tempstring2.compare("numequilibration") == 0){
					ofile << "\nThe number sampling steps to equilibrate samples = ?";
					ifile >> numEquilSamples;
				}
				else if (tempstring2.compare("autotimestep")== 0){
					ofile << "\nTime step will be determined automatically\n";
					MDAutoTimeStep = true;
				}
				else if (tempstring2.compare("restore")== 0){
					ofile << "\nSave and load the progress of MD simulations.\n";
					MDAllowRestore = true;
				}
				else {
					ofile << tempstring2 << " is an undefined command\n";
				}
			}

			if (MD_Temperature < 0.0){
				ofile << "\nTemperature is undefined;\n";
				if (dim == 1){
					MD_Temperature = 2e-4;
				}
				else if (dim == 2){
					MD_Temperature = 2e-6;
				}
				else if (dim == 3){
					MD_Temperature = 1e-6;
				}
				ofile << "Take T to be a default value, "<< MD_Temperature << "\n";
			}
			else{
				ofile << "Input T = "<< MD_Temperature << "\n";
			}
			if (strcmp (tempstring, "random") == 0){
				/* Start from a random initial condition and basic properties */
				
				SpherePacking pck = GetInitConfigs(0);
				potential -> SetPacking(pck);
				Configuration pConfig(pck, "a");
				size_t MDEquilibrateSamples = numEquilSamples;
				//bool MDAutoTimeStep = false;
				CollectiveCoordinateMD(&pConfig, potential, rngGod, MDTimeStep, MD_Temperature, savename, N_config, MDStepPerSample, MDAllowRestore, TimeLimit, MDEquilibrateSamples, MDAutoTimeStep);
			}
			else if	(strcmp (tempstring, "input") == 0){
				/* For continue from the latest configuration in the previous run. */
				ofile << "MD Temperature = "<< MD_Temperature << std::endl;

				SpherePacking pck = GetInitConfigs(N_load - 1);
				potential -> SetPacking(pck);
				Configuration pConfig(pck, "a");

				//bool MDAllowRestore = true;
				size_t MDEquilibrateSamples = numEquilSamples;
				size_t numConfig_comp = (N_config > N_init)? N_config - N_init : 0 ;
				
				//bool MDAutoTimeStep = false;
				CollectiveCoordinateMD(&pConfig, potential, rngGod, MDTimeStep, MD_Temperature, savename, numConfig_comp, MDStepPerSample, MDAllowRestore, TimeLimit, MDEquilibrateSamples, MDAutoTimeStep);
			//TODO
			}
			
			
			
		}
		else if (strcmp (mode, "ground") == 0)
		{
			/* Quench to zero temperature */
			/* Parameters for minimizations */
			size_t max_steps = 10000;
			double tolerance = 1e-14;
			bool split_species = true;
			std::string algorithm = "LBFGS";

			ofile << "--------------------------------------- \n";
			ofile << "\t\t Input parameters for minimizations \n";
			ofile << "--------------------------------------- \n";
			for(;;){
				ifile >> tempstring2;
				std::transform(tempstring2.begin(), tempstring2.end(), tempstring2.begin(),::tolower); // make the input in the lower case
				if (tempstring2.compare("run") == 0){
					/* start minimizations simulation */
					break;
				}
				else if (tempstring2.compare("tolerance") == 0){
					ofile << "\n energy tolerance for ground states [in the unit of v0] = ";
					ifile >> tolerance ;
				}
				else if (tempstring2.compare("maxsteps") == 0){
					ofile << "\n Limit number of evaluations = ";
					ifile >> max_steps ;
				}
				else if (tempstring2.compare("algorithm") == 0){
					ofile << "\nMinimization algorithm = ";
					ifile >> algorithm ; 
				}
				else if (tempstring2.compare("combined") == 0){
					ofile << "\nWrite the configurations of all species\n";
					split_species = false;
				}
				else {
					ofile << tempstring2 << " is an undefined command\n";
				}
			}
			ofile << "\nNumerical optimizer is "<< algorithm <<std::endl;
			ofile << "Energy tolerance is "<< tolerance <<std::endl;
			ofile << "Max. steps of evaluations is "<< max_steps <<std::endl;


			if (strcmp (InitCondition.c_str(), "random") == 0){
				/* random initial conditions 
					=> use MultiRun */

				SpherePacking pck = GetInitConfigs(0);
				potential -> SetPacking(pck);
				Configuration pConfig(pck, "a");

				pConfig.PrepareIterateThroughNeighbors(potential->GetRc());

				{
					ConfigurationPack save_(savename+"_Success");				
					N_success = save_.NumConfig();
				}
				size_t numConfig_comp = (N_config > N_success )? N_config - N_success : 0;

				ofile << "We already have " << N_success << " Configurations in the savefile \n";
				ofile << "just compute "<< numConfig_comp <<" more Configurations \n";

				ofile << "max_eval = " << max_steps << std::endl;
				CollectiveCoordinateMultiRun(&pConfig, potential, rngGod, savename, numConfig_comp, TimeLimit, algorithm, tolerance, max_steps);

				// separate configurations 
				if (split_species){
					std::vector<ConfigurationPack> cp;
					for(size_t s = 0; s  < potential->species_constraints.size(); s++){
						char name[300];
						sprintf(name, "%s_s%d_Success", savename.c_str(), s);
						cp.emplace_back(name);
					}
					ConfigurationPack success(savename+"_Success");

					for (size_t i = N_success; i < success.NumConfig(); i++){
						Configuration c = success.GetConfig(i);
						Configuration empty_cell(c);
						empty_cell.RemoveParticles();
						for(size_t s = 0; s  < potential->species_constraints.size(); s++){
							Configuration c_spec(empty_cell);
							for (auto & id: potential->species_constraints[s].idx){
								c_spec.Insert("a", c.GetRelativeCoordinates(id));
							}
							cp[s].AddConfig(c_spec);						
						}
					}
				}
				delete potential;

			}
			else if	(strcmp (InitCondition.c_str(), "input") == 0){
				/* ReadConfigPack */
				size_t idx_beg = 0;
				{
					ConfigurationPack load_(savename+"_InitConfig");				
					N_init = load_.NumConfig();
					ConfigurationPack success_(savename+"_InitConfig");				
					N_success = success_.NumConfig();
				}
				size_t numConfig_upper = std::min(N_load, beg_idx+N_config);

				ofile << "We start from the " << beg_idx << "-th Configuration from loadfile\n";  
				ofile << "Among them, we already have used " << N_init << " Configurations  \n";
				N_init += beg_idx;
				ofile << "compute from the "<< N_init << "th Config. to the " << numConfig_upper <<"th Config.\n";
				size_t id = N_init;
				for (size_t i = N_success; i < N_config; i++){
					if ( id < N_load ){
						SpherePacking pck = GetInitConfigs(id);
						Configuration pConfig(pck, "a");
						pConfig.PrepareIterateThroughNeighbors(potential->GetRc());
						CollectiveCoordinateMultiRun(&pConfig, potential, rngGod, savename, 1, TimeLimit, algorithm, tolerance, max_steps);
					}
					id ++ ;
				}
				delete potential;				
			}

		}
		else{
			ofile << "mode "<< mode << " is undefined!\n";
			return 1;
		}
		// Some basic analyses
		/*Compute the nearest distance of each config. */
		if (false){
			ConfigurationPack ConfigSet;
			if (strcmp (mode, "MD") == 0){
				ConfigSet.Open((savename ).c_str());
			}
			else if (strcmp (mode, "ground") == 0){
				ConfigSet.Open((savename + "_Success").c_str());
			}
			double rmin = L;
			{
				for (int i = 0; i < ConfigSet.NumConfig(); i++) {
					Configuration x = ConfigSet.GetConfig(i);
					double rm = L;
					for (size_t j = 0; j < x.NumParticle(); j++) {
						double temp = x.NearestParticleDistance(j);
						rm = (temp < rm) ? temp : rm;
					}
					ofile << i << ":\t" << rm << std::endl;
				}
			}

			auto getC = [&rmin, &ConfigSet](size_t i) ->Configuration {
				Configuration c = ConfigSet.GetConfig(i);
				double R = rmin;
				for (size_t j = 0; j < c.NumParticle(); j++) {
					double temp = c.NearestParticleDistance(j);
					R = (temp < R) ? temp : R;
				}
				rmin = (R < rmin) ? R : rmin;
				return c;
			};

			/*Compute g2 function*/
			std::vector<GeometryVector> g2;
			N_config = ConfigSet.NumConfig();
			IsotropicTwoPairCorrelation(getC, N_config, std::min(0.2 * L, 5.0) , g2);
			WriteFunction(g2, (savename + "_g2").c_str());
			phi = 1.0*HyperSphere_Volume(dim, rmin / 2.0);
			ofile << "max packing fraction is " << phi << std::endl;

			/*Compute the structure factsor*/
			std::vector<GeometryVector> Sk;
			IsotropicStructureFactor(getC, N_config, 0.2, 1, Sk);
			WriteFunction(Sk, (savename + "_Sk").c_str());
		}		
	}



	return 0;
}



int main(int argc, char ** argv){
	char tempstring[1000] = {};
	char tempstring2[300] = {};
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	RandomGenerator rngGod(0), rng(0);

	if (argc > 1 && strcmp(argv[1], "multi") == 0){
		GetMultiCCO(argc, argv);
	}
	else
	{
		GetCCO(argc, argv);
	}
	return 0;
	{
		MC_CCO(argc, argv);
	}
	return 0;
	{
		// generate stitched point configurations from an ensemble.

		size_t num_selected = 0, num_generated = 0;
		std::string loadname, savename;

		ofile << "load an ensemble from \n";
		ifile >> loadname;
		ofile << "save as \n";
		ifile >> savename;
		ConfigurationPack loaded_ensemble(loadname);
		ConfigurationPack result(savename);
		result.Clear();

		ofile << "num. of selected configurations from " + loadname <<"\n";
		ifile >> num_selected;
		ofile << "num. of generated configurations, saved as " + savename <<"\n";
		ifile >> num_generated;

		int option = 0;
		ofile << "how to (0:identical, 1:distinct) \n";
		ifile >> option;

		double L = 0.0, l = 0.0;
		size_t num_configs = 0;
		{
			Configuration c = loaded_ensemble.GetConfig(0);
			l = c.GetBasisVector(0).x[0];
			L = l * num_selected;
			num_configs = loaded_ensemble.NumConfig();
			ofile << "Using " << num_configs << "configurations in an ensemble, we generate " << num_generated << "configurations of N = "<< c.NumParticle() * num_selected << " and L = " << L << "\n";
 		}
		Configuration c = GetUnitCubicBox(1);
		c.Resize(L);

		for (size_t i = 0 ; i < num_generated; i++){
			Configuration finC(c);
			size_t id = 0;
			if (option == 0){
				id = static_cast<size_t> (std::floor(std::floor(num_configs * rng.RandomDouble())));
			}

			for (size_t j = 0; j < num_selected; j++){
				if (option == 1){
					size_t id = static_cast<size_t> (std::floor(std::floor(num_configs * rng.RandomDouble())));
				}
				Configuration patch = loaded_ensemble.GetConfig(id);

				for (size_t k = 0; k < patch.NumParticle(); k++){
					
					finC.Insert("a", 
						finC.CartesianCoord2RelativeCoord( 
							patch.GetCartesianCoordinates(k) + GeometryVector(static_cast<double> (j*l))
						) 
					);

				}

			}
			result.AddConfig(finC);
		}


	}

	return 0;
}


int CollectiveCoordinateMultiRun(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, std::string Prefix, size_t SampleNumber, time_t TimeLimit, const std::string & AlgorithmName, double E_max, size_t max_num)
{
	DimensionType dim=pConfig->GetDimension();
	size_t Num=pConfig->NumParticle();
	ConfigurationPack AfterRelaxPack(Prefix);
	ConfigurationPack SuccessPack(Prefix+"_Success");
	ConfigurationPack InitConfigPack(Prefix+"_InitConfig");
	if(AfterRelaxPack.NumConfig()>0)
		std::cout<<"Found "<<AfterRelaxPack.NumConfig()<<" configurations, continue MultiRun.\n";

	std::cout<<"MultiRun start."<<'\n';
	for(size_t i=0; i<SampleNumber; i++)
	{
		std::cout<<"at time "<<std::time(nullptr)-ProgramStart<<",\tgenerating config"<<i<<"\n";
		Configuration result(*pConfig);
		pPotential->SetConfiguration(*pConfig);
		double E_init=pPotential->Energy();

		InitConfigPack.AddConfig(result);
		if(AlgorithmName=="LBFGS")
			RelaxStructure_NLOPT(result, *pPotential, 0.0, 0, 0.0, max_num);
		else if(AlgorithmName=="LocalGradientDescent")
			RelaxStructure_LocalGradientDescent(result, *pPotential, 0.0, 0, 0.0, max_num);
		else if(AlgorithmName=="ConjugateGradient")
			RelaxStructure_ConjugateGradient(result, *pPotential, 0.0, 0, 0.0, max_num);
		else if(AlgorithmName=="SteepestDescent")
			RelaxStructure_SteepestDescent(result, *pPotential, 0.0, 0, 0.0, max_num);
		else if(AlgorithmName=="MINOP")
			RelaxStructure_MINOP_withoutPertubation(result, *pPotential, 0.0, 0, 0.0, max_num);
		else {
			std::cout << "Algorithm is undefined!!"<< std::endl;
			return -1;
		}
		AfterRelaxPack.AddConfig(result);
		pPotential->SetConfiguration(result);
		double E=pPotential->Energy();
		std::cout<<" \tE_init="<<E_init<<" \n";
		std::cout<<" \tE_relax="<<E<<" \n";
		if(E<E_max)
		{
			std::cout<<"Add to success pack\n";
			SuccessPack.AddConfig(result);
		}
		if(std::time(nullptr) > TimeLimit)
		{
			//std::fstream ofile2("continue.txt", std::fstream::out);
			//ofile2<<"Not Completed. Please run again\n";
			return 0;
		}

		for(size_t n=0; n<Num; n++)
		{
			GeometryVector temp(dim);
			for(DimensionType j=0; j<dim; j++)
				temp.x[j]=gen.RandomDouble();
			pConfig->MoveParticle(n, temp);
		}
	}
	return 0;
}


int CollectiveCoordinateMD(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, double TimeStep, double Temperature, std::string Prefix, size_t SampleNumber, size_t StepPerSample, bool AllowRestore, time_t TimeLimit, size_t EquilibrateSamples, bool MDAutoTimeStep)
{
	ConfigurationPack BeforeRelaxPack(Prefix);

	size_t OneTenthStepPerSample = StepPerSample/10;
	if(OneTenthStepPerSample==0)
		OneTenthStepPerSample=1;

	DimensionType dim=pConfig->GetDimension();
	size_t Num=pConfig->NumParticle();
	size_t dimTensor=dim*Num;
	Potential * pPot=pPotential;
	//unsigned long tempNumThread=1;//use this when doing MD because there is parallelization in MD code, so no need to parallelize in the potential
	ParticleMolecularDynamics * psystem = NULL;

	bool Restart = false;
	//signed char stage=0;
	signed int stage=0;
	long long step=0;
	if(AllowRestore)
	{
		std::fstream ifile( Prefix+std::string(".MDDump"), std::fstream::in | std::fstream::binary);
		if(ifile.good())
		{
			ifile.read( (char*)(&stage), sizeof(stage) );
			ifile.read( (char*)(&step), sizeof(step) );
			psystem = new ParticleMolecularDynamics(ifile);
			pPot->SetConfiguration(psystem->Position);
			Restart=true;
			std::cout << "Continue from " << Prefix+std::string(".MDDump") <<"\n";
			
			double Ek = psystem->GetKineticEnergy();
			double Ep = pPot->Energy();
			std::cout << "stage "<< stage << ",\tstep " << step <<": \tE_relax="<< Ep <<" \t E_k=" << Ek <<"\n";
			step ++;
		}
	}
	if(Restart == false)
	{
		BeforeRelaxPack.Clear();
		/* equilibration starts from a ground state */
		RelaxStructure_NLOPT(*pConfig, *pPotential, 0.0, 0, 0.0, 100000);
		pPot->SetConfiguration(* pConfig);
		double E=pPot->Energy();

		psystem = new ParticleMolecularDynamics(*pConfig, TimeStep, 1.0); 
	}

	size_t NumExistConfig=0;
	std::cout<<"CCMD start. Temperature="<<Temperature<<'\n';
	
	if(stage==0)
	{
		/* equilibration starts from a ground state. */
		stage++;
		step=0;
	}
	if(stage==1)
	{
		//stage 1: equilibration after adjusting time step
		std::cout << "------------------------\n";
		std::cout << "	Equilibration step    \n";
		std::cout << "------------------------\n";	
		for(long long i=step; i<EquilibrateSamples; i++)
		{
			if (std::time(nullptr) > TimeLimit || std::time(nullptr) > ::TimeLimit)
			{
				std::fstream ofile( Prefix+std::string(".MDDump"), std::fstream::out | std::fstream::binary);
				ofile.write( (char*)(&stage), sizeof(stage) );
				ofile.write( (char*)(&i), sizeof(i) );
				psystem->WriteBinary(ofile);
				delete psystem;
				std::fstream ofile2( Prefix+std::string("_continue.txt"), std::fstream::out);
				ofile2<<"Not Completed. Please run again\n";
				return 0;
			}
			double c0=psystem->Position.GetCartesianCoordinates(0).x[0];

			for(size_t ii=0; ii<2; ii++)
			{
				psystem->SetRandomSpeed(Temperature, gen);
				//std::swap(tempNumThread, pPot->ParallelNumThread);
				if(MDAutoTimeStep)
//					psystem->Evolve_AutoTimeStep(StepPerSample/2, *pPot, 0.0001/SampleNumber);
					psystem->Evolve_AutoTimeStep(StepPerSample/2, *pPot, 0.0001);
				else
					psystem->Evolve(StepPerSample/2, *pPot);
				//std::swap(tempNumThread, pPot->ParallelNumThread);
			}

			pPot->SetConfiguration(psystem->Position);
			std::cout<<"1:"<<i<<"/"<<(EquilibrateSamples)<<", x0="<<psystem->Position.GetCartesianCoordinates(0).x[0]<<", Ep="<<pPot->Energy()<<", Ek="<<psystem->GetKineticEnergy()<<", dt="<<psystem->TimeStep<<'\n';

		}
		stage++;
		step=0;
	}

	//stage 2: sample
	std::cout << "------------------------\n";
	std::cout << "		Sampling step     \n";
	std::cout << "------------------------\n";
	long long final_step=0;
	for(long long i=step; i<SampleNumber; i++)
	{
		if (std::time(nullptr) > TimeLimit || std::time(nullptr) > ::TimeLimit)
		{
			std::fstream ofile( Prefix + std::string(".MDDump"), std::fstream::out | std::fstream::binary);
			ofile.write( (char*)(&stage), sizeof(stage) );
			ofile.write( (char*)(&i), sizeof(i) );
			psystem->WriteBinary(ofile);
			delete psystem;
			std::fstream ofile2("continue.txt", std::fstream::out);
			ofile2<<"Not Completed. Please run again\n";
			return 0;
		}
		std::cout<<"at time"<<std::time(nullptr)-ProgramStart;
		//std::swap(tempNumThread, pPot->ParallelNumThread);
		psystem->AndersonEvolve(StepPerSample/2, *pPot, Temperature, 0.01, gen);
		psystem->Evolve(StepPerSample/2, *pPot);
		//std::swap(tempNumThread, pPot->ParallelNumThread);
		Configuration result(psystem->Position);
		
		pPot->SetConfiguration(result);
		double Ek = psystem->GetKineticEnergy(), Ep = pPot->Energy();
		std::cout<<", 2:"<<i<<"/"<<(SampleNumber)<<" \tE_relax="<< Ep <<" \t";
		std::cout<<"E_k="<<Ek<<", dt="<<psystem->TimeStep<<'\n';

		//AfterRelaxPack.AddConfig(result);
		//if(QuenchAfterMD)
		BeforeRelaxPack.AddConfig(result);

		std::cout.flush();
		final_step = i;
	}

	if (AllowRestore){
		std::cout << "Store the last configuration and simulation parameters for future use.\n";
		std::fstream ofile( Prefix+std::string(".MDDump"), std::fstream::out | std::fstream::binary);
		ofile.write( (char*)(&stage), sizeof(stage) );
		ofile.write( (char*)(&final_step), sizeof(final_step) );
		psystem->WriteBinary(ofile);
	}
	//::Output("FinalConfiguration", psystem->Position);
	delete psystem;

	return 0;
}

double getK(double chi, double rho, DimensionType d)
{
	return 2*pi *std::pow(2*d*rho*chi/HyperSphere_Volume(d,1.0), 1./d);
}

double get_chi(double K, double rho, DimensionType d)
{
	return HyperSphere_Volume(d, K) /(2*d*std::pow(2*pi,d)) /rho;
}
