/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: March 2022 */

/** \file main.cpp 
 * \brief A simple CLI for computing various pair statistics from point configurations. 
 * This is a simplified version of the PairStatisticsCLI code written by Ge.*/

#include "Computation.h"
//#include "ContactStatistics.h"

/* header fields in $(cores) */
#include <PeriodicCellList.h>
#include <etc.h>
#include <algorithm>

size_t Verbosity = 2;

class ReadConfigPackWeights
{
public:
	double Rescale;
	size_t StartIdx;
	ConfigurationPack p;
	RandomGenerator rng;
	std::string prefix;
	size_t column = 0;	
	
	std::function<std::complex<double> (std::complex<double> input)> filter; 

	bool sequential_sampling = true;
	bool weight_shuffling = false;
	ReadConfigPackWeights(std::istream & ifile, std::ostream & ofile, double Rescale, size_t StartIdx = 0)
	{
		this->Rescale=Rescale;
		this->StartIdx = StartIdx;
		this->filter = [](std::complex<double> x)-> std::complex<double> { return x;};
		ofile<<"Input Prefix (ConfigPack):";
		std::string temp;
		ifile>>temp;
		p.Open(temp);

		ofile<<"Input Prefix (Weights):";
		ifile>> this->prefix;
		ofile<<"column to read:";
		ifile>> column;
	}

	void replace_filter(std::function<const std::complex<double> (std::complex<double>)> new_filter){
		this->filter = new_filter;
	}

	PeriodicCellList<std::complex<double>> operator() (size_t i)
	{
		Configuration result;
		size_t idx; 
		if (sequential_sampling){
			idx = i+this->StartIdx;
		}
		else{
			idx = (int)(std::floor(p.NumConfig() * rng.RandomDouble()));
		}
		result=p.GetConfig(idx);
		result.Rescale(Rescale);
		PeriodicCellList<std::complex<double>> out (result, std::complex<double>(0.,0.));
		out.RemoveParticles();

		{
			char name[300] = {};//std::string name;
			sprintf(name, this->prefix.c_str(), idx);
			std::fstream IFILE(name, std::fstream::in);

			if (IFILE.good() == false)
			{

				std::cerr << "Error in ReadFunction : input stream is not good! " << name <<"\n";
				return MarkedPointPattern();
			}

			std::string buffer;	std::istringstream iss;
			size_t id = 0;
			while (std::getline(IFILE, buffer) && (id < result.NumParticle())) {
				if (buffer[0] != '#') {/* ignore comments or field names. */
					iss.str(buffer);
					iss.clear();	//THIS IS VERY CRUCIAL!!!
					std::complex<double> temp;
					for (size_t i = 0; i < column; i++) {
						iss >> temp;
						if (i == column - 1){
							temp = this->filter(temp);
							out.Insert(temp, result.GetRelativeCoordinates(id));
						}
					}
					//	if (IFILE.eof() == false)
						//result.push_back(temp);
					id++;
				}
			}

			IFILE.close();
		}

		if (weight_shuffling)
		{
			std::vector<std::complex<double>> weights;
			for (size_t i = 0; i < out.NumParticle(); ++i)
			{
				weights.push_back(out.GetCharacteristics(i));
			}
			size_t N = weights.size();

			/* shuffle the weights */
			for (int i = weights.size() - 1; i > 0; --i) {
				int j = rng.RandomInt()%N; // random index in [0, i]
				std::swap(weights[i], weights[j]);
			}

			PeriodicCellList<std::complex<double>> temp (out);
			out.RemoveParticles();
			for (size_t i = 0; i < temp.NumParticle(); ++i)
			{
				out.Insert(weights[i], temp.GetRelativeCoordinates(i));
			}
		}
		return out;
	}
};

//take every configuration from GetConfigsFunction, then rescale it to unit number density, then output
class RescaleToUnitDensity
{
public:
	std::function<const Configuration(size_t i)> GetConfigsFunction;
	RescaleToUnitDensity(std::function<const Configuration(size_t i)> GetConfigsFunction) : GetConfigsFunction(GetConfigsFunction)
	{}
	Configuration operator() (size_t i)
	{
		Configuration result = GetConfigsFunction(i);
		result.Resize(result.NumParticle());
		return result;
	}
};

int PairStatisticsCLI();

int main(int argc, char ** argv){
	int signal;
	signal = PairStatisticsCLI();
	return signal;
}

int PairStatisticsCLI(){
	char tempstring[1000];
	/* An extern variable defined in ${core}/etc.h */
	Verbosity = 3;

	std::istream & ifile=std::cin;
	std::ostream & ofile=std::cout;

	double Rescale=1.0;
	size_t NumConfig =0;
	size_t StartIdx = 0;
	size_t NumThreads = 1;
	std::string OutputPrefix;
	std::string GraceTitle;
	std::function<const PeriodicCellList<std::complex<double>> (size_t i)>  GetConfigsFunction = nullptr;	
	std::function<std::complex<double> (std::complex<double>)> new_filter;
	bool doReplaceFilter = false;
	bool weight_shuffling = false;

	bool AlsoWriteGrace = false;
	bool RandomSampling = false;

	std::vector<std::string> available_computations{
		"g2","spectraldensity", "Variance"};
	std::vector<std::string> todo_computations{
		// "HpComputation",
		// "HvComputation",
		// "VoronoiVolumeComputation",
		// "VoronoiNumSidesComputation",
		// "CoveringRadiusDistributionComputation",
		// "PackingRadiusDistributionComputation",
		// "PercolationVolumeFractionDistribution",
		// "PercolationP1Calculation",
		// "PercolationP1Calculation_v2",
		// "VoronoiVolumeCorrelationComputation",
		// "Psi6CorrelationComputation",
		// "LocalQ6DistributionComputation",
		// "AnalyticWindowNumberVarianceComputation",
		// "WindowNumberDistributionComputation",
		// "AverageClusterSizeComputation",
		// "NeighborLinkTortuosityComputation",
		// "M2Computation",
		// "DiscretizationVolumeFractionComputation",
		// "FirstPassageTimeComputation"
	};

	std::vector<Computation *> vpComputations;
	for(;;)
	{
		ifile>>tempstring;
		/* Basic parameters */
		if(strcmp(tempstring, "Exit")==0)
		{
			for(auto iter=vpComputations.begin(); iter!=vpComputations.end(); iter++)
				delete *iter;
			return 0;
		}
		else if(strcmp(tempstring, "Available")==0){
			for(auto iter=available_computations.begin(); iter!=available_computations.end(); iter++)
				ofile << *iter <<";\n";
		}
		else if(strcmp(tempstring, "ToDo")==0){
			for(auto iter=todo_computations.begin(); iter!=todo_computations.end(); iter++)
				ofile << *iter <<";\n";
		}
		else if(strcmp(tempstring, "GraceTitle")==0)
		{
			ifile>>GraceTitle;
		}
		else if(strcmp(tempstring, "OutputPrefix")==0)
		{
			ifile>>OutputPrefix;
			ofile << "Output Prefix = " << OutputPrefix << std::endl;
		}
		else if(strcmp(tempstring, "AlsoWriteEPS")==0)
		{
			ifile>>AlsoWriteEPS;
		}
		else if(strcmp(tempstring, "AlsoWriteGrace")==0)
		{
			ifile>>AlsoWriteGrace;
		}
		else if(strcmp(tempstring, "Rescale")==0)
		{
			std::cerr << "Warning : Rescale is deprecated\n";
			ifile>>Rescale;
		}
		else if(strcmp(tempstring, "NumConfig")==0)
		{
			ifile>>NumConfig;
		}
		else if (strcmp(tempstring, "StartIndex") == 0)
		{
			ifile>> StartIdx;
			ofile<< "Start Index = "<< StartIdx <<std::endl;
		}
		else if (strcmp(tempstring, "RandomSampling") == 0)
		{
			ifile >> RandomSampling;
			if (RandomSampling){
				ofile << "Turn on random sampling in GetConfigsFunction" << std::endl;
				ofile << "This option is effective only for ReadConfigPack" << std::endl;
				//TODO: It is applied only to the ReadConfigPack case. Please apply it ot other cases

			}
			else{
				ofile << "Turn off random sampling in GetConfigsFunction" << std::endl;				
			}
		}

		/* temporary filter functions */
		else if(strcmp(tempstring, "filter5")==0)
		{
			if (GetConfigsFunction == nullptr){
				doReplaceFilter=true;
				new_filter = [](std::complex<double> x) -> std::complex<double> {
					if ((4.9<x.real()) && (x.real()<5.1)){
						return std::complex<double>(1.,0.0);							
					}
					else{
						return std::complex<double> (0.,0.);
					}
				};
			}
			else{
				std::cerr << "Input this command before GetConfigsFunction !!\n";
			}
		}
		else if(strcmp(tempstring, "filter6")==0)
		{
			if (GetConfigsFunction == nullptr){
				doReplaceFilter=true;
				new_filter = [](std::complex<double> x) -> std::complex<double> {
					if ((5.9<x.real()) && (x.real()<6.1)){
						return std::complex<double>(1.,0.0);							
					}
					else{
						return std::complex<double> (0.,0.);
					}
				};
			}
			else{
				std::cerr << "Input this command before GetConfigsFunction !!\n";
			}
		}
		else if(strcmp(tempstring, "filter7")==0)
		{
			if (GetConfigsFunction == nullptr){
				doReplaceFilter=true;
				new_filter = [](std::complex<double> x) -> std::complex<double> {
					if ((6.9<x.real()) && (x.real()<7.1)){
						return std::complex<double>(1.,0.0);							
					}
					else{
						return std::complex<double> (0.,0.);
					}
				};
			}
			else{
				std::cerr << "Input this command before GetConfigsFunction !!\n";
			}
		}
		
		/* shuffling weights */
		else if (strcmp(tempstring, "ShuffleWeights") == 0)
		{
			ifile >> weight_shuffling;
			if (weight_shuffling)
			{
				ofile << "Weight shuffling is enabled.\n";
			}
			else
			{
				ofile << "Weight shuffling is disabled.\n";
			}
		}

		/* define computations */
		else if(strcmp(tempstring, "g2Computation")==0)
		{
			vpComputations.push_back( new g2Computation(ifile, ofile) );
		}
		else if (strcmp(tempstring, "SkComputation") == 0)
		{
			vpComputations.push_back(new SkComputation(ifile, ofile));
		}
		else if (strcmp(tempstring, "LocalVariance") == 0)
		{
			vpComputations.push_back( new LocalVariance(ifile, ofile) );
		}

		else if (strcmp(tempstring, "GetConfigsFunction") == 0)
		{

			ofile<<"Centers:";
			ifile>>tempstring;
			if (strcmp(tempstring, "ReadConfigPack") == 0)
			{
				ReadConfigPackWeights c(ifile, ofile, Rescale, StartIdx);
				if (RandomSampling){
					c.sequential_sampling = false;
					NumConfig=c.p.NumConfig();
					std::cout<<"The config pack contains "<<NumConfig<<" configurations, set NumConfig to this value.\n";
				}
				//automatic sets NumConfig since Configuration Pack contains this information
				else{
					c.sequential_sampling = true;
					std::cout<<"The config pack contains "<<c.p.NumConfig()<<" configurations.";
					if (c.p.NumConfig() < StartIdx){
						GetConfigsFunction = nullptr;
						NumConfig = 0;
						std::cout << "StartIdx exceeds the size of ConfigPack\n";
					}
					else{
						size_t num = c.p.NumConfig()-StartIdx;
	
						if (NumConfig == 0){
							std::cout << "Set NumConfig to this value because the former is undefined.\n";
							NumConfig= num;
						}
						else if (NumConfig > num){
							std::cout << "Set NumConfig to this value because the former is greater than the latter.\n";
							NumConfig= num;
						}
					}
				}
				if(doReplaceFilter){
					ofile << "\n\nFilter will transform the particle weights\n";
					c.replace_filter(new_filter);
				}
				if(weight_shuffling){
					c.weight_shuffling = weight_shuffling;
				}

				GetConfigsFunction = c;
			}
		}
		// else if (strcmp(tempstring, "PrintTheta") == 0){
		// 	if (GetConfigsFunction == nullptr)
		// 		std::cerr << "Specify GetConfigsFunction before printing!\n";
		// 	else
		// 	{
		// 		for (int index = 0; index < NumConfig; index++) {
		// 			std::string file = OutputPrefix + std::string("__") + std::to_string(index) + ".txt";
		// 			std::fstream OFILE(file.c_str(), std::fstream::out);

		// 			MarkedPointPattern c = GetConfigsFunction(index);
		// 			ofile << "---------------------------\n";
		// 			ofile << "\t" << index << "\n";
		// 			ofile << "---------------------------\n";
		// 			for (int i=0; i<c.NumParticle(); i++){
		// 				ofile << i << ": ";
		// 				for (size_t j = 0; j < c.GetDimension(); j++){
		// 					ofile << c.GetCartesianCoordinates(i).x[j] << " ";
		// 				}
		// 				ofile << c.GetCharacteristics(i) << "\n";
		// 			}

		// 			OFILE.close();
		// 		}
		// 	}
		// }
		else if (strcmp(tempstring, "PrintConfigs") == 0) { //print ConfigurationPack => readable txt files
			if (GetConfigsFunction == nullptr)
				std::cerr << "Specify GetConfigsFunction before printing!\n";
			else 
			{
				
				for (int index = 0; index < NumConfig; index++) {
					std::string file = OutputPrefix + std::string("__") + std::to_string(index);
					//sprintf(file, "%s__%d", OutputPrefix, index);
					//WriteConfiguration(GetConfigsFunction(index), file);
					MarkedPointPattern c = GetConfigsFunction(index);
					ofile << "---------------------------\n";
					ofile << "\t" << index << "\n";
					ofile << "---------------------------\n";
					for (int i=0; i<c.NumParticle(); i++){
						ofile << i << ": ";
						for (size_t j = 0; j < c.GetDimension(); j++){
							ofile << c.GetCartesianCoordinates(i).x[j] << " ";
						}
						ofile << c.GetCharacteristics(i) << "\n";
					}
				}
			}
		}

		// else if (strcmp(tempstring, "RescaleToUnitDensity") == 0)
		// {
		// 	if (GetConfigsFunction == nullptr)
		// 		std::cerr << "Specify GetConfigsFunction before specifying rescaling!\n";
		// 	else
		// 	{
		// 		RescaleToUnitDensity c(GetConfigsFunction);
		// 		GetConfigsFunction = c;
		// 	}
		// }
		else if(strcmp(tempstring, "NumThreads")==0){
			ifile >> NumThreads;
		}

		else if(strcmp(tempstring, "Calculation")==0)
		{
			if(GetConfigsFunction == nullptr)
				std::cerr<<"Specify GetConfigsFunction before compute!\n";
			else
			{
				for(auto iter=vpComputations.begin(); iter!=vpComputations.end(); iter++)
				{
					(*iter)->SetNumThreads(NumThreads);
					if (NumConfig > 0){
						(*iter)->Compute(GetConfigsFunction, NumConfig);
						(*iter)->Write(OutputPrefix);
						if(AlsoWriteGrace)
							(*iter)->Plot(OutputPrefix, GraceTitle);
					}
					else{
						ofile << "Skip calculations \n";
					}
				}
			}
		}
		
		
		else
		{
			if (vpComputations.empty())
				std::cout << "Unrecognized command!\n";
			else
				vpComputations.back()->ProcessAdditionalOption(tempstring, std::cin, std::cout);

			std::cin.clear();
		}
		ifile.ignore(1000,'\n');
		if(ifile.eof()) return 2;
	}
	return 1;
}

