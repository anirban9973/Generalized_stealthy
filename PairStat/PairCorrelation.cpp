#include "PairCorrelation.h"
#include <omp.h>

//adaptive bins calculated after calculating partial pair distances, use less memory
void IsotropicTwoPairCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double MaxDistance, std::vector<GeometryVector> & Result, size_t SampleDistanceSize, double ResolutionPreference)
{
	if(!(MaxDistance>0.0))
		return;
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"computing g_2";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	DimensionType dim;

	size_t NumConfigsSampled=0;
	size_t TotalAtomCount=0;
	double AverageDensity=0;
	size_t StartParticle=0;
	size_t NumParticle=1;
	for(size_t j=0; j<NumConfigs; j++)
	{
		Configuration Config = GetConfigsFunction(j);
		Config.SetCellSize(0.35*MaxDistance);
		NumParticle=Config.NumParticle();
		if(j==0)
			dim=Config.GetDimension();
		//pre-allocate
		Config.IterateThroughNeighbors(GeometryVector(Config.GetDimension()), MaxDistance, 
			[](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
		{
		});
		for(size_t i=0; i<NumParticle; i++)
		{
			//Configuration::particle * a=Config.GetParticle(i);
			Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(i), MaxDistance, 
				[&distances, &i](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
			{
				if(SourceAtom>i || (SourceAtom==i && LatticeShift.Modulus2()!=0) )
					distances.push_back(std::sqrt(shift.Modulus2()));
			});
			if(distances.size()>SampleDistanceSize)
			{
				StartParticle=i+1;
				break;
			}
		}
		if(distances.size()>SampleDistanceSize)
			break;

		NumConfigsSampled++;
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();
		pd++;
	}

	
	//if(Verbosity>2)
	//	std::cout<<"generating bins...";

	//find the lowest pair distance
	if(distances.empty())
	{
		std::cerr<<"Warning in IsotropicTwoPairCorrelation : Pair distance not found in specified cutoff!\n";
		Result.clear();
		return;
	}

	//prepare for the iteration
	Result.clear();
	//iteration
	HistogramGenerator HGen;
	HGen.GenerateBins(distances, static_cast<double>(NumConfigs)/((double)(NumConfigsSampled)+(double)(StartParticle)/NumParticle), ResolutionPreference);

	//discard empty bins at the end. Such bins are useful when we don't have an upper bound of data (e.g. nearest neighbor distance binning)
	//here, they are useless
	while(HGen.bins.size()>1 && HGen.bins.back().Count==0)
		HGen.bins.pop_back();

	//the bin counts
	for(size_t j=NumConfigsSampled; j<NumConfigs; j++)
	{
		Configuration Config = GetConfigsFunction(j);
		Config.SetCellSize(0.35*MaxDistance);
		Config.RefineCellList();
		size_t NumParticle=Config.NumParticle();
		//pre-allocate
		Config.PrepareIterateThroughNeighbors(MaxDistance);

#pragma omp parallel
		{
#pragma omp for schedule(guided)
			for(size_t i=StartParticle; i<NumParticle; i++)
			{
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(i), MaxDistance, 
					[&HGen, &i](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
				{
					if(SourceAtom>i || (SourceAtom==i && LatticeShift.Modulus2()!=0) )
					{
						double dist=std::sqrt(shift.Modulus2());
						HGen.Report(dist);
					}
				});
			}
		}
		StartParticle=0;
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();
		pd++;
	}
	AverageDensity/=NumConfigs;

	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
	{
		double BinEnd;
		if(iter->Start>=MaxDistance)
			continue;
		else if(iter!=HGen.bins.end()-1)
			BinEnd=std::min(iter[1].Start, MaxDistance);
		else
			BinEnd=MaxDistance;
		
		GeometryVector temp(4);
		temp.x[0]=(BinEnd+iter->Start)/2;
		temp.x[1]=(double)iter->Count/( (HyperSphere_Volume(dim, BinEnd)-HyperSphere_Volume(dim, iter[0].Start))*0.5*AverageDensity*TotalAtomCount );//g_2(r)
		temp.x[2]=(temp.x[0]-iter[0].Start);//\delta r
		if(iter->Count!=0)
			temp.x[3]=temp.x[1]/std::sqrt(static_cast<double>(iter->Count));//\delta g_2(r)
		else
			temp.x[3]=0.0;

		//if bin is too small and inaccurate, discard
		//if(Result.empty()==false)
		//	if(temp.x[2]<1e-3*Result.back().x[2] && iter->Count< 3)
		//		continue;
		if(temp.x[0]>0)
			Result.push_back(temp);
	}

	if(Verbosity>2)
		std::cout<<"done!\n";

}

void IsotropicTwoPairCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double MaxDistance, std::vector<GeometryVector> & Result, HistogramGenerator & HGen)
{
	if (!(MaxDistance>0))
		return;
	//get the pair distances list
	//if(Verbosity>2)
	//	std::cout<<"Generate g_2, sampling pair distances...";
	if (Verbosity>2)
		std::cout << "computing g_2";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	DimensionType dim;

	size_t TotalAtomCount = 0;
	double AverageDensity = 0;
	//prepare for the iteration
	Result.clear();

	//the bin counts
	//if(Verbosity>2)
	//	std::cout<<"counting pair distances...";
	for (size_t j = 0; j<NumConfigs; j++)
	{
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		if (j == 0)
			dim = Config.GetDimension();
		Config.SetCellSize(0.35*MaxDistance);
		Config.RefineCellList();
		size_t NumParticle = Config.NumParticle();
		//pre-allocate
		Config.IterateThroughNeighbors(GeometryVector(Config.GetDimension()), MaxDistance,
			[](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
		{
		});

		//Configuration::particle * p0=Config.GetParticle(0);
#pragma omp parallel
		{
#pragma omp for schedule(guided)
				for (size_t i = 0; i<NumParticle; i++)
			{
				//Configuration::particle * a=Config.GetParticle(i);
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(i), MaxDistance,
					[&HGen, &i](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
				{
					if (SourceAtom>i || (SourceAtom == i && LatticeShift.Modulus2() != 0))
					{
						double dist = std::sqrt(shift.Modulus2());
						HGen.Report(dist);

						//debug temp
						//if(dist<3e-4)
						//	std::cerr<<"Particle "<<SourceAtom-p0<<" and "<<a-p0<<" overlap !\n";
					}
				});
			}
		}
		TotalAtomCount += Config.NumParticle();
		AverageDensity += Config.NumParticle() / Config.PeriodicVolume();
		pd++;
	}
	AverageDensity /= NumConfigs;

	for (auto iter = HGen.bins.begin(); iter != HGen.bins.end(); iter++)
	{
		double BinEnd;
		if (iter->Start >= MaxDistance)
			continue;
		else if (iter != HGen.bins.end() - 1)
			BinEnd = std::min(iter[1].Start, MaxDistance);
		else
			BinEnd = MaxDistance;

		GeometryVector temp(4);
		temp.x[0] = (BinEnd + iter->Start) / 2;
		temp.x[1] = iter->Count / ((HyperSphere_Volume(dim, BinEnd) - HyperSphere_Volume(dim, iter[0].Start))*0.5*AverageDensity*TotalAtomCount);//g_2(r)
		temp.x[2] = (temp.x[0] - iter[0].Start);//\delta r
		if (iter->Count != 0)
			temp.x[3] = temp.x[1] / std::sqrt(static_cast<double>(iter->Count));//\delta g_2(r)
		else
			temp.x[3] = 0.0;

		//if bin is too small and inaccurate, discard
		//if(Result.empty()==false)
		//	if(temp.x[2]<1e-3*Result.back().x[2] && iter->Count< 3)
		//		continue;
		if (temp.x[0]>0)
			Result.push_back(temp);
	}

	if (Verbosity>2)
		std::cout << "done!\n";
}


void TwoPairCorrelation_2DAnisotropic(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double MaxDistance, std::vector<std::vector<GeometryVector> > & Result, size_t NumDirectionalBins, std::vector<std::vector<GeometryVector> > * pError, size_t SampleDistanceSize, double ResolutionPreference)
{
	if(!(MaxDistance>0))
		return;
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"Generate g_2, sampling pair distances...";

	std::vector<double> distances;

	DimensionType dim;

	size_t NumConfigsSampled=0;
	size_t TotalAtomCount=0;
	double AverageDensity=0;
	for(size_t j=0; j<NumConfigs; j++)
	{
		if(Verbosity>3 || (Verbosity>2&&j%100==0) )
			std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		if(j==0)
		{
			dim=Config.GetDimension();

			if(dim!=2)
			{
				std::cerr<<"Error in TwoPairCorrelation_2DAnisotropic : Configuration is not two-dimensional!\n";
				std::cerr<<"This function only support two dimensional configurations!\n";
				assert(false);
			}
		}
		size_t NumParticle=Config.NumParticle();
		//pre-allocate
		Config.IterateThroughNeighbors(GeometryVector(Config.GetDimension()), MaxDistance, 
			[](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
		{
		});
		for(size_t i=0; i<NumParticle; i++)
		{
			//Configuration::particle * a=Config.GetParticle(i);
			Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(i), MaxDistance, 
				[&distances, &i](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
			{
				if(SourceAtom>i || (SourceAtom==i && LatticeShift.Modulus2()!=0) )
					distances.push_back(std::sqrt(shift.Modulus2()));
			});
		}
		NumConfigsSampled++;
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();

		if(distances.size()>SampleDistanceSize)
			break;
	}

	if(Verbosity>2)
		std::cout<<"generating bins...";

	//find the lowest pair distance
	if(distances.empty())
	{
		std::cerr<<"Warning in IsotropicTwoPairCorrelation : Pair distance not found in specified cutoff!\n";
		Result.clear();
		return;
	}

	//prepare for the iteration
	Result.clear();
	//iteration
	HistogramGenerator HGen;
	HGen.GenerateBins(distances, static_cast<double>(NumConfigs)/NumConfigsSampled/NumDirectionalBins, ResolutionPreference);

	//discard empty bins at the end. Such bins are useful when we don't have an upper bound of data (e.g. nearest neighbor distance binning)
	//here, they are useless
	while(HGen.bins.size()>1 && HGen.bins.back().Count==0)
		HGen.bins.pop_back();

	//clear the counts, we need directional g2(r), so the previous isotropic counts are useless
	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
		iter->Count=0;

	//Generate directional bins
	std::vector<HistogramGenerator> DirectionalHGen(NumDirectionalBins, HGen);

	//Re-do the bin counts
	if(Verbosity>2)
		std::cout<<"counting pair distances...";
	TotalAtomCount=0;
	AverageDensity=0.0;
	for(size_t j=0; j<NumConfigs; j++)
	{
		if(Verbosity>3 || (Verbosity>2&&j%100==0) )
			std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		size_t NumParticle=Config.NumParticle();
		//pre-allocate
		Config.IterateThroughNeighbors(GeometryVector(Config.GetDimension()), MaxDistance, 
			[](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
		{
		});


#pragma omp parallel
		{
#pragma omp for schedule(guided)
			for(size_t i=0; i<NumParticle; i++)
			{
				//Configuration::particle * a=Config.GetParticle(i);
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(i), MaxDistance, 
					[&DirectionalHGen, &i, &NumDirectionalBins](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
				{
					if(SourceAtom>i || (SourceAtom==i && LatticeShift.Modulus2()!=0) )
					{
						double dist=std::sqrt(shift.Modulus2());

						//angle/pi, within range [0, 1]
						double angle= 0.5 - std::atan(shift.x[0]/shift.x[1])/::pi;

						//debug temp
						//std::cout<<angle<<'\n';

						//turns 1 into 0 for binning
						while(! (angle < 1) )
							angle-=1;

						size_t AngleBin=std::floor(angle*NumDirectionalBins);
						DirectionalHGen[AngleBin].Report(dist);
					}
				});
			}
		}
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();
	}
	AverageDensity/=NumConfigs;

	Result.resize(NumDirectionalBins, std::vector<GeometryVector>() );
	if(pError!=nullptr)
		pError->resize(NumDirectionalBins, std::vector<GeometryVector>() );
	for(size_t i=0; i<NumDirectionalBins; i++)
	{
		for(auto iter=DirectionalHGen[i].bins.begin(); iter!=DirectionalHGen[i].bins.end(); iter++)
		{
			double BinEnd;
			if(iter->Start>=MaxDistance)
				continue;
			else if(iter!=DirectionalHGen[i].bins.end()-1)
				BinEnd=std::min(iter[1].Start, MaxDistance);
			else
				BinEnd=MaxDistance;

			GeometryVector temp(4);
			temp.x[0]=(BinEnd+iter->Start)/2;//r
			temp.x[1]=( (double)(i)/NumDirectionalBins + 0.5 )* ::pi;//theta
			temp.x[2]=iter->Count/( (HyperSphere_Volume(dim, BinEnd)-HyperSphere_Volume(dim, iter->Start))*0.5*AverageDensity*TotalAtomCount )*NumDirectionalBins;//g_2(r)
			Result[i].push_back(temp);

			//debug temp
			//std::cout<<i<<" \t"<<iter-DirectionalHGen[i].bins.begin()<<" \t"<<iter->Count<<" \t"<<BinEnd<<" \t"<<iter->Start<<" \t"<<temp.x[2]<<'\n';

			if(pError != nullptr)
			{
				temp.x[0]=(temp.x[0]-iter[0].Start);//\delta r
				temp.x[1]= ::pi/(2*NumDirectionalBins);//\delta theta
				if(iter->Count!=0)
					temp.x[2]=temp.x[2]/std::sqrt(static_cast<double>(iter->Count));//\delta g_2(r)
				else
					temp.x[2]=0.0;
				pError->at(i).push_back(temp);
			}
		}
	}
	if(Verbosity>2)
		std::cout<<"done!\n";
}



void NearestNeighborDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleDistanceSize, double ResolutionPreference)
{
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"Computing H_p(r)";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	size_t NumConfigsSampled=0;
	size_t TotalAtomCount=0;
	double AverageDensity=0;
	for(size_t j=0; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);

		const size_t NumParticle = Config.NumParticle();
		double TypicalLength=std::pow(Config.PeriodicVolume()/NumParticle, 1.0/Config.GetDimension())*1.5;
		for(size_t particle=0; particle<NumParticle; particle++)
		{
			std::vector<GeometryVector> neighbors;
			double l=TypicalLength;
			while(neighbors.size()<2)
			{
				neighbors.clear();
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(particle), l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
				{
					neighbors.push_back(shift);
				});
				l*=2;
			}

			std::partial_sort(neighbors.begin(), neighbors.begin()+2, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2()<right.Modulus2();});
			distances.push_back(std::sqrt(neighbors[1].Modulus2()));
		}
		NumConfigsSampled++;
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();

		if(distances.size()>SampleDistanceSize)
			break;
	}

	//if(Verbosity>2)
	//	std::cout<<"generating bins...";
	//prepare for the iteration
	Result.clear();
	//iteration
	HistogramGenerator HGen;
	HGen.GenerateBins(distances, static_cast<double>(NumConfigs)/NumConfigsSampled, ResolutionPreference);

	//discard empty bins at the end. Such bins are useful when we don't have an upper bound of data (e.g. nearest neighbor distance binning)
	//here, they are useless
//	while(HGen.bins.size()>1 && HGen.bins.back().Count==0)
//		HGen.bins.pop_back();

	//Re-do the bin counts
	//if(Verbosity>2)
	//	std::cout<<"counting pair distances...";
	for(size_t j=NumConfigsSampled; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		const size_t NumParticle=Config.NumParticle();
		//Configuration::particle * p0=Config.GetParticle(0);

		for(size_t particle=0; particle<NumParticle; particle++)
		{
			std::vector<GeometryVector> neighbors;
			double TypicalLength=std::pow(Config.PeriodicVolume()/NumParticle, 1.0/Config.GetDimension())*1.5;
			double l=TypicalLength;
			while(neighbors.size()<2)
			{
				neighbors.clear();
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(particle), l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
				{
					neighbors.push_back(shift);
				});
				l*=2;
			}

			std::partial_sort(neighbors.begin(), neighbors.begin()+2, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2()<right.Modulus2();});
			HGen.Report(std::sqrt(neighbors[1].Modulus2()));
		}
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();
	}
	AverageDensity/=NumConfigs;

	size_t TotalCount=0;
	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
		TotalCount+=iter->Count;

	bool IgnoreBin = true;
	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
	{
		GeometryVector temp(4);
		if(iter!=HGen.bins.end()-1)
			temp.x[0]=((iter[0].Start)+(iter[1].Start))/2.0;//r 
		else
			continue;

		if(iter==HGen.bins.begin())
		{
			if(iter->Start==0.0 && (iter[2].Start-iter[1].Start)>iter[1].Start && iter->Count==0)
				continue;//there shouldn't be a bin at here
		}
		temp.x[2]=(temp.x[0]-iter[0].Start);//\delta r

		temp.x[1]=static_cast<double>(iter->Count)/TotalCount/(2*temp.x[2]);//p(r)

		if(iter->Count!=0)
			temp.x[3]=temp.x[1]/std::sqrt(static_cast<double>(iter->Count));//\delta p(r)
		else
			temp.x[3]=0.0;

		if (temp.x[0]>=0.0){
			IgnoreBin &= (temp.x[1] == 0.0);	/* Once the second column is non-zero, IgnoreBin becomes false.*/
			if (! IgnoreBin ){
				Result.push_back(temp);
			}
		}
	}

	if(Verbosity>2)
		std::cout<<"done!\n";
}

void HvDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleDistanceSize, double ResolutionPreference, double OverSampling)
{
	RandomGenerator gen(81479058);
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"Computing H_v(r)";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	DimensionType dim;

	size_t NumConfigsSampled=0;
	size_t TotalAtomCount=0;
	double AverageDensity=0;
	for(size_t j=0; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		if(j==0)
			dim=Config.GetDimension();
		const size_t NumParticle=Config.NumParticle();
		double TypicalLength=std::pow(Config.PeriodicVolume()/NumParticle, 1.0/Config.GetDimension())*1.5;
		for (size_t j = 0; j<NumParticle*OverSampling; j++)
		{
			std::vector<GeometryVector> neighbors;
			double l=TypicalLength;
			GeometryVector temp(dim);
			for(int i=0; i<dim; i++)
				temp.x[i]=gen.RandomDouble();
			while(neighbors.size()<2)
			{
				neighbors.clear();
				Config.IterateThroughNeighbors(temp, l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
				{
					neighbors.push_back(shift);
				});
				l*=2;
			}

			std::partial_sort(neighbors.begin(), neighbors.begin()+1, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2()<right.Modulus2();});
			distances.push_back(std::sqrt(neighbors[0].Modulus2()));
		}
		NumConfigsSampled++;
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();

		if(distances.size()>SampleDistanceSize)
			break;
	}

	//if(Verbosity>2)
	//	std::cout<<"generating bins...";
	//prepare for the iteration
	Result.clear();
	//iteration
	HistogramGenerator HGen;
	HGen.GenerateBins(distances, static_cast<double>(NumConfigs)/NumConfigsSampled, ResolutionPreference);

	//Re-do the bin counts
	//if(Verbosity>2)
	//	std::cout<<"counting pair distances...";
	for(size_t j=NumConfigsSampled; j<NumConfigs; j++)
	{
			pd++;
			//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
			//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
			Configuration Config = GetConfigsFunction(j);
			const size_t NumParticle=Config.NumParticle();
			//Configuration::particle * p0=Config.GetParticle(0);

			for(size_t j=0; j<NumParticle; j++)
			{
				std::vector<GeometryVector> neighbors;
				double TypicalLength=std::pow(Config.PeriodicVolume()/NumParticle, 1.0/Config.GetDimension())*1.5;
			double l=TypicalLength;
			while(neighbors.size()<2)
			{
				neighbors.clear();
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(j), l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
				{
					neighbors.push_back(shift);
				});
				l*=2;
			}

			std::partial_sort(neighbors.begin(), neighbors.begin()+2, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2()<right.Modulus2();});
			HGen.Report(std::sqrt(neighbors[1].Modulus2()));
		}
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();
	}
	AverageDensity/=NumConfigs;

	size_t TotalCount=0;
	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
		TotalCount+=iter->Count;

	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
	{
		GeometryVector temp(4);
		if(iter!=HGen.bins.end()-1)
			temp.x[0]=((iter[0].Start)+(iter[1].Start))/2.0;//r 
		else
			continue;

		if(iter==HGen.bins.begin())
		{
			if(iter->Start==0.0 && (iter[2].Start-iter[1].Start)>iter[1].Start && iter->Count==0)
				continue;//there shouldn't be a bin at here
		}
		temp.x[2]=(temp.x[0]-iter[0].Start);//\delta r

		temp.x[1]=static_cast<double>(iter->Count)/TotalCount/(2*temp.x[2]);//p(r)

		if(iter->Count!=0)
			temp.x[3]=temp.x[1]/std::sqrt(static_cast<double>(iter->Count));//\delta p(r)
		else
			temp.x[3]=0.0;
		Result.push_back(temp);
	}

	if(Verbosity>2)
		std::cout<<"done!\n";
}



//return the minimum pair distance in the configuration
double MinDistance(const Configuration & a, size_t * pi, size_t * pj)
{
	DimensionType dim=a.GetDimension();
	double UnitSphereVolume= ::HyperSphere_Volume(dim, 1.0);
	double result=std::pow(a.PeriodicVolume()/a.NumParticle()/UnitSphereVolume, (double)(1.0)/dim)*2;
	for(size_t i=0; i<a.NumParticle(); i++)
	{
		GeometryVector loc=a.GetRelativeCoordinates(i);
		a.IterateThroughNeighbors(loc, result, [&](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
		{
			if(i==SourceAtom)
				return;
			double dis=std::sqrt(shift.Modulus2());
			if(dis<result)
			{
				result=dis;
				if(pi!=nullptr)
					(*pi)=i;
				if(pj!=nullptr)
					(*pj)=SourceAtom;
			}
		});
	}
	return result;
}


double MeanNearestNeighborDistance(const Configuration & a)
{
	double sum=0.0;
	for(size_t j=0; j<a.NumParticle(); j++)
	{
		std::vector<GeometryVector> neighbors;
		double TypicalLength=std::pow(a.PeriodicVolume()/a.NumParticle(), 1.0/a.GetDimension())*1.5;
		double l=TypicalLength;
		while(neighbors.size()<2)
		{
			neighbors.clear();
			a.IterateThroughNeighbors(a.GetRelativeCoordinates(j), l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l*=2;
		}

		std::partial_sort(neighbors.begin(), neighbors.begin()+2, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2()<right.Modulus2();});
		sum+=std::sqrt(neighbors[1].Modulus2());
	}
	return sum/a.NumParticle();
}

size_t HistogramToPDF(const HistogramGenerator & HGen, std::vector<GeometryVector> & result)
{
	result.clear();
	size_t TotalCount = 0;
	for (auto iter = HGen.bins.begin(); iter != HGen.bins.end(); iter++)
		TotalCount += iter->Count;

	for (auto iter = HGen.bins.begin(); iter != HGen.bins.end(); iter++)
	{
		GeometryVector temp(4);
		if (iter != HGen.bins.end() - 1)
			temp.x[0] = ((iter[0].Start) + (iter[1].Start)) / 2.0;//r 
		else
			continue;

		if (iter == HGen.bins.begin())
		{
			if (iter->Start == 0.0 && (iter[2].Start - iter[1].Start)>iter[1].Start && iter->Count == 0)
				continue;//there shouldn't be a bin at here
		}
		temp.x[2] = (temp.x[0] - iter[0].Start);//\delta r

		temp.x[1] = static_cast<double>(iter->Count) / TotalCount / (2 * temp.x[2]);//p(r)

		if (iter->Count != 0)
			temp.x[3] = temp.x[1] / std::sqrt(static_cast<double>(iter->Count));//\delta p(r)
		else
			temp.x[3] = 0.0;
		result.push_back(temp);
	}
	return TotalCount;
}
size_t HistogramToCDF(const HistogramGenerator & HGen, std::vector<GeometryVector> & result)
{
	result.clear();
	size_t TotalCount = 0;
	for (auto iter = HGen.bins.begin(); iter != HGen.bins.end(); iter++)
		TotalCount += iter->Count;

	//debug temp
	//std::cout << "TotalCount=" << TotalCount;
	size_t AlreadyCount = 0;

	for (auto iter = HGen.bins.begin(); iter != HGen.bins.end(); iter++)
	{
		GeometryVector temp(4);
		if (iter != HGen.bins.end() - 1)
			temp.x[0] = ((iter[0].Start) + (iter[1].Start)) / 2.0;//r 
		else
			continue;

		if (iter == HGen.bins.begin())
		{
			if (iter->Start == 0.0 && (iter[2].Start - iter[1].Start)>iter[1].Start && iter->Count == 0)
				continue;//there shouldn't be a bin at here
		}
		temp.x[2] = (temp.x[0] - iter[0].Start);//\delta r
		temp.x[1] = (AlreadyCount + 0.5*iter->Count) / TotalCount;//p(r)
		temp.x[3] = (0.5*iter->Count) / TotalCount;//\delta p(r)
		result.push_back(temp);
		AlreadyCount += iter->Count;
	}
	return TotalCount;
}


