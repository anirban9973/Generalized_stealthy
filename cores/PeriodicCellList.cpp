#include "PeriodicCellList.h"
#include <exception>



void Output(std::string prefix, const Configuration & List)
{
	prefix+=std::string(".pos");
	std::fstream out(prefix.c_str(), std::fstream::out);
	out.precision(17);

	Output(out, List);
}

void Output(std::ostream & out, const Configuration & List)
{
	const Configuration * list = & List;
	//determine Rc of output
	//GeometryVector temp(list->GetDimension());
	//for( ::DimensionType i=0; i<list->GetDimension(); i++)
	//{
	//	GeometryVector t=list->GetBasisVector(i);
	//	if(temp.Dot(t)>0)
	//		temp.AddFrom(t);
	//	else
	//		temp.MinusFrom(t);
	//}
	//double Rc=2*std::sqrt(temp.Modulus2())/std::pow(list->NumParticle(), 1.0/list->GetDimension());
	double Rc=std::pow(100*list->PeriodicVolume()/list->NumParticle(), 1.0/list->GetDimension());
	//////////////////////////////////////
	GeometryVector origin(list->GetDimension());
	size_t count=0;
	list->IterateThroughNeighbors(origin, Rc, [&count](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void {count++;} );
	out<<count<<'\n';
	//the second line is a comment
	out<<'\n';
	//.xyz part, used by visualization sofwares
	list->IterateThroughNeighbors(origin, Rc, [&out, &list](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void 
	{
		out<<list->GetCharacteristics(SourceAtom).name<<" \t"; 
		shift.OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	} );


	out<<"\n\n****************************************\n\n";

	out<<"Primitive vectors\n";
	for(::DimensionType i=0; i<list->GetDimension(); i++)
	{
		out<<"a(";
		out<<i;
		out<<") = \t";
		list->GetBasisVector(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}

	out<<"\n\nVolume =  "<<list->PeriodicVolume()<<"\n\n";
	out<<"Reciprocal vectors\n";
	for(::DimensionType i=0; i<list->GetDimension(); i++)
	{
		out<<"b(";
		out<<i;
		out<<") = \t";
		list->GetReciprocalBasisVector(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}
	out<<"\n\nBasis Vectors:\nAtom    Lattice Coordinates                Cartesian Coordinates\n\n";
	for(size_t i=0; i<list->NumParticle(); i++)
	{
		//Configuration::particle * now=list->GetParticle(i);
		out<<list->GetCharacteristics(i).name<<" \t";
		list->GetRelativeCoordinates(i).OutputCoordinate(out, list->GetDimension());
		list->GetCartesianCoordinates(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}
}
Configuration ReadPos(std::string prefix)
{
	std::string name=prefix+std::string(".pos");
	std::fstream file(name.c_str(), std::fstream::in);
	return ReadPos(file);
}
Configuration ReadPos(std::istream & ifile)
{
	const std::string flag1("Primitive");
	const std::string flag2("Basis");

	std::string current;

	size_t NumParticle=0;

	::DimensionType dimension=0;
	do
	{
		std::getline(ifile, current);
		if(ifile.eof() || ifile.fail())
		{
			return Configuration(0, nullptr, 0.0, false);
		}
	}
	while(current.find(flag1)!=0);
	for(;;)
	{
		std::getline(ifile, current);
		if(current.size()<3)
			break;
		if(current[2]>='0' && current[2]<='9')
			dimension++;
		else
			break;
	}

	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag2)!=0);
	std::getline(ifile, current);
	std::getline(ifile, current);
	while(ifile.eof()==false)
	{
		char junk[1000];
		ifile.getline(junk, 1000);

		if(ifile.eof()==false)
			NumParticle++;
	}
	ifile.clear();
	ifile.seekg(0, std::ios::beg);


	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag1)!=0);

	std::vector<GeometryVector> base;
	for(int i=0; i<dimension; i++)
	{
		char a;
		do
		{
			a=ifile.get();
		}
		while(a!='=');
		base.push_back(GeometryVector(dimension));
		base.back().InputCoordinate(ifile, dimension);
	}
	double volume= ::Volume(&base[0], dimension);
	double cellsize = std::pow(volume/(NumParticle/ParticlePerCell), static_cast<double>(1)/dimension);
	Configuration result(dimension, &base[0], cellsize);

	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag2)!=0);
	std::getline(ifile, current);
	std::getline(ifile, current);
	while(ifile.eof()==false)
	{
		char name[128];
		ifile >> name;

		GeometryVector loc(dimension);
		loc.InputCoordinate(ifile, dimension);

		GeometryVector junk(dimension);
		junk.InputCoordinate(ifile, dimension);

		if(ifile.eof()==false)
			result.Insert(name, loc);
	}
	return result;
}
Configuration ReadHoomdXml(std::string prefix)
{
	prefix+=std::string(".xml");
	std::fstream ifile(prefix.c_str());
	return ReadHoomdXml(ifile);
}
Configuration ReadHoomdXml(std::istream & ifile)
{
	std::string temp;
	std::string flag1("<configuration");
	std::string flag2("<type");

	DimensionType dim;
	size_t NumParticle;
	double volume=1;

	std::vector<GeometryVector> basis;
	std::vector<double> coords;

	do
	{
		std::getline(ifile, temp);
	}
	while(temp.find(flag1)!=0);

	{
		std::stringstream stream;
		stream<<temp;

		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>NumParticle;//this is actually timestep. junk it
		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>dim;
		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>NumParticle;
	}

	std::getline(ifile, temp);
	{
		std::stringstream stream;
		stream<<temp;
		double length;

		for(DimensionType d=0; d<dim; d++)
		{
			//read a number
			for(;;)
			{
				char c=stream.peek();
				if(c>='0' && c<='9')
					break;
				if(c=='.')
					break;
				if(stream.eof())
					break;
				stream>>c;
			}
			stream>>length;
			GeometryVector tempVec(dim);
			tempVec.x[d]=length;
			basis.push_back(tempVec);
			volume=volume*length;
		}
	}

	double cellsize = std::pow(volume/(NumParticle/ParticlePerCell), static_cast<double>(1)/dim);
	Configuration result(dim, &basis[0], cellsize);
	std::getline(ifile, temp);
	for(size_t i=0; i<dim*NumParticle; i++)
	{
		double tempD;
		ifile>>tempD;
		coords.push_back(tempD);
	}


	do
	{
		std::getline(ifile, temp);
	}
	while(temp.find(flag2)!=0);
	for(size_t i=0; i<NumParticle; i++)
	{
		std::getline(ifile, temp);
		assert(temp.length()<3);
		GeometryVector rel(dim);
		for(DimensionType j=0; j<dim; j++)
			rel.x[j]=coords[i*dim+j]/basis[j].x[j];

		result.Insert(temp.c_str(), rel);
	}

	return result;
}
Configuration ReadCoordinate(std::istream & file)
{
	return ReadCoordinate(file, 1.0);
}
Configuration ReadCoordinate(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;

	//determine dimension
	std::getline(file, temp);
	std::stringstream str;
	str<<temp;
	for(;;)
	{
		double dtemp=-1.0;
		str>>dtemp;
		if(dtemp!=-1.0)
			dim++;
		else
			break;
	}

	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
			file>>rel.x[j];
		result.Insert("A", rel);
	}

	return result;
}
Configuration ReadCartesianCoordinate(std::istream & file)
{
	return ReadCartesianCoordinate(file, 1.0);
}
Configuration ReadCartesianCoordinate(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;

	//determine dimension
	std::getline(file, temp);
	std::stringstream str;
	str<<temp;
	for(;;)
	{
		double dtemp=-1.0;
		str>>dtemp;
		if(dtemp!=-1.0)
			dim++;
		else
			break;
	}

	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
			file>>rel.x[j];
		rel=result.CartesianCoord2RelativeCoord(rel);
		result.Insert("A", rel);
	}

	return result;
}

Configuration ReadGro(std::istream & file)
{
	GeometryVector basis[ ::MaxDimension];
	DimensionType dim=3;
	size_t NumP=0;
	std::string temp;

	//determine NumP
	std::getline(file, temp);
	file>>NumP;
	for(size_t i=0; i<NumP; i++)
	{
		std::getline(file, temp);
	}
	std::getline(file, temp);
	//determine box size
	double Volume=1.0;
	for(int i=0; i<3; i++)
	{
		double l;
		file>>l;
		basis[i]=GeometryVector(3);
		basis[i].x[i]=l;
		Volume*=l;
	}
	file.clear();
	file.seekg(0, std::ios::beg);
	std::getline(file, temp);
	std::getline(file, temp);

	double CellSize=Volume/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		//Only read atoms named 'O'!
		file>>temp;
		file>>temp;
		if(temp.at(0)=='O')
		{
			file>>temp;
			GeometryVector rel(dim);
			for(size_t j=0; j<dim; j++)
				file>>rel.x[j];
			rel=result.CartesianCoord2RelativeCoord(rel);
			result.Insert("A", rel);
		}
		else
			std::getline(file, temp);
	}

	return result;
}
Configuration ReadStealthOutput(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;
	std::vector<double> OneOverSize;//the side length of the input file

	//determine dimension
	std::getline(file, temp);
	{
		std::stringstream str;
		str<<temp;
		str>>temp;
		str>>dim;
	}

	//determine size
	std::getline(file, temp);
	{
		std::stringstream str;
		str<<temp;
		str>>temp;
		for(DimensionType i=0; i<dim; i++)
		{
			double tt;
			str>>tt;
			OneOverSize.push_back(1.0/tt);
		}
	}

	std::getline(file, temp);
	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);
	std::getline(file, temp);
	std::getline(file, temp);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
		{
			file>>rel.x[j];
			rel.x[j]*=OneOverSize[j];
		}
		result.Insert("A", rel);
	}

	return result;
}
Configuration ReadStealthOutput(std::istream & file)
{
	return ReadStealthOutput(file, 1.0);
}

SpherePacking ReadTJOutput(std::istream & inFile)
{
	std::string input;
	std::getline(inFile,input);
	std::stringstream ss1(input);
	int dim, N;
	ss1 >> dim;
	std::string buffer;
	ss1 >> buffer; //"HS"
	std::string packType;
	ss1 >> packType; //"poly" or "mono"
	getline(inFile,input); //Line 2 is junk.
	inFile >> N;

	double radMono;
	if (packType == "mono") 
	{//This line is the diameter value for monodisperse packings
		inFile >> radMono;
		radMono /= 2.0; //Convert from diameter to radius
	}
	else if (packType != "poly")
	{
		std::cerr << "Error in ReadTJOutput : Unsupported packType!\n";
		return SpherePacking();
	}

	// lattice vectors next
	GeometryVector basis[ ::MaxDimension];
	for (int i=0; i<dim; i++) 
	{
		basis[i]=GeometryVector(dim);
		for (int j=0; j<dim; j++) 
		{
			double val;
			inFile >> val;
			basis[i].x[j]=val;
		}
	}
	SpherePacking packing(dim, basis, ::MaxDistance);

	// then more trash
	getline(inFile,input);
	getline(inFile,input);

	// then coordinates and radii
	for (int i=0; i<N; i++) 
	{
		double radii;
		GeometryVector car(dim);
		for (int d=0; d<dim; d++) 
		{
			inFile >> car.x[d];
		}
		if (packType == "poly") 
		{
			inFile >> radii;
		}
		else 
		{
			radii = radMono;
		}
		packing.Insert(radii, packing.CartesianCoord2RelativeCoord(car));
	}
	packing.SetCellSize(packing.GetMaxRadius()*1.05);

	return packing; //Got a packing!

}
SpherePacking ReadTJOutput(const std::string & prefix)
{
	std::string name = prefix + std::string(".dat");
	std::fstream ifile(name.c_str(), std::fstream::in);
	return ReadTJOutput(ifile);
}
void WriteTJOutput(std::fstream & ofile, const SpherePacking & pk)
{
	ofile.precision(17);
	DimensionType dim = pk.GetDimension();
	size_t num = pk.NumParticle();
	if (dim == 0 || num == 0)
	{
		std::cerr << "Error in WriteTJOutput : invalid packing!\n";
		return;
	}
	//determine mono or poly
	bool mono = true;
	for(size_t i=1; i<num; i++)
		if (pk.GetCharacteristics(i) != pk.GetCharacteristics(0))
		{
			mono = false;
			break;
		}

	ofile << dim << "	hs	";
	if (mono)
		ofile << "mono";
	else
		ofile << "poly";
	ofile << "\n" << num << " 1\n" << num << '\n';
	if (mono)
		ofile << 2 * pk.GetCharacteristics(0) << '\n';

	for (DimensionType i = 0; i < dim; i++)
	{
		for (DimensionType j = 0; j < dim; j++)
			ofile << pk.GetBasisVector(i).x[j] << " \t";
		ofile << '\n';
	}
	ofile << "T T T";
	for (size_t n = 0; n < num; n++)
	{
		ofile << "\n";
		for (DimensionType j = 0; j < dim; j++)
			ofile << pk.GetCartesianCoordinates(n).x[j] << " \t";
		if (!mono)
			ofile << pk.GetCharacteristics(n);
	}
}
void WriteTJOutput(const std::string & prefix, const SpherePacking & pk)
{
	std::string name = prefix + std::string(".dat");
	std::fstream ifile(name.c_str(), std::fstream::out);
	return WriteTJOutput(ifile, pk);
}


Configuration MultiplicateStructure(const Configuration & src, size_t NumberInEachSide)
{
	double dNumberInEachSide=static_cast<double>(NumberInEachSide);
	DimensionType dim=src.GetDimension();
	std::vector<GeometryVector> bases;
	for(DimensionType i=0; i<dim; i++)
		bases.push_back(src.GetBasisVector(i)*dNumberInEachSide);
	Configuration result(src);
	result.ChangeBasisVector(&bases[0]);
	result.RemoveParticles();

	for(size_t i=0; i<src.NumParticle(); i++)
	{
		GeometryVector OriginalRelative=(src.GetRelativeCoordinates(i))*(1.0/dNumberInEachSide);
		const char * name=src.GetCharacteristics(i).name;

		//d-dimensional loop
		std::vector<size_t> indexes(dim, 0);
		while(indexes.back()!=NumberInEachSide)
		{
			GeometryVector a(OriginalRelative);
			for(DimensionType j=0; j<dim; j++)
				a.x[j]+=indexes[j]/dNumberInEachSide;
			result.Insert(name, a);

			//loop end
			indexes[0]++;
			for(DimensionType j=0; j<dim-1; j++)
			{
				if(indexes[j]==NumberInEachSide)
				{
					indexes[j]=0;
					indexes[j+1]++;
				}
			}
		}
	}
	return result;
}
Configuration MultiplicateStructure(const Configuration & src, const std::vector<size_t> & NumberInEachSide)
{
	assert(NumberInEachSide.size()==src.GetDimension());
	std::vector<double> dNumberInEachSide;
	for(auto iter=NumberInEachSide.begin(); iter!=NumberInEachSide.end(); iter++)
		dNumberInEachSide.push_back((double)(*iter));

	DimensionType dim=src.GetDimension();
	std::vector<GeometryVector> bases;
	for(DimensionType i=0; i<dim; i++)
		bases.push_back(src.GetBasisVector(i)*dNumberInEachSide[i]);
	Configuration result(src);
	result.ChangeBasisVector(&bases[0]);
	result.RemoveParticles();

	for(size_t i=0; i<src.NumParticle(); i++)
	{
		GeometryVector OriginalRelative=(src.GetRelativeCoordinates(i));
		for(int i=0; i< src.GetDimension(); i++)
			OriginalRelative.x[i]/=dNumberInEachSide[i];

		const char * name=src.GetCharacteristics(i).name;

		//d-dimensional loop
		std::vector<size_t> indexes(dim, 0);
		while(indexes.back()!=NumberInEachSide[dim-1])
		{
			GeometryVector a(OriginalRelative);
			for(DimensionType j=0; j<dim; j++)
				a.x[j]+=indexes[j]/dNumberInEachSide[j];
			result.Insert(name, a);

			//loop end
			indexes[0]++;
			for(DimensionType j=0; j<dim-1; j++)
			{
				if(indexes[j]==NumberInEachSide[j])
				{
					indexes[j]=0;
					indexes[j+1]++;
				}
			}
		}
	}
	return result;
}
Configuration GetUnitCubicBox(DimensionType d, double CellSize)
{
	std::vector<GeometryVector> vbas;
	for(DimensionType i=0; i<d; i++)
	{
		GeometryVector t(d);
		t.x[i]=1.0;
		vbas.push_back(t);
	}
	Configuration result(d, &vbas[0], CellSize);

	return result;
}



void ConfigurationPack::Open(const std::string & prefix)
{
	this->IndexName = prefix + ".ConfigPackIndex";
	this->PackName = prefix + ".ConfigPack";
	std::fstream ifile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	if (ifile.good())
	{
		ifile.read((char*)(&NConfig), sizeof(NConfig));
		if (ifile.eof())
			NConfig = 0;
	}
	else
		NConfig = 0;
}
long long ConfigurationPack::NumConfig(void)
{
	return NConfig;
}
void ConfigurationPack::AddConfig(const Configuration & c)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary);
	if (indexFile.good() == false)
	{
		indexFile.clear();
		indexFile.close();
		indexFile.clear();
		indexFile.open(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	}
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::ate);
	if (packFile.good() == false)
	{
		packFile.clear();
		packFile.close();
		packFile.clear();
		packFile.open(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	}
	this->NConfig++;
	indexFile.seekp(0, std::fstream::beg);
	indexFile.write((char*)(&NConfig), sizeof(NConfig));
	long long loc = packFile.tellp();
	indexFile.seekp(0, std::fstream::end);
	indexFile.write((char*)(&loc), sizeof(loc));
	c.WriteBinary(packFile);
	indexFile.close();
	if (indexFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Index file failing!\n";
	packFile.close();
	if (packFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Pack file failing!\n";
}
Configuration ConfigurationPack::GetConfig(long long i)
{
	if (i >= this->NConfig)
	{
		std::cerr << "Error in ConfigurationPack::GetConfig : Index is larger than number of configuration!\n";
		exit(1);
	}
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::binary);
	long long loc;
	indexFile.seekg(sizeof(long long)*(i + 1));
	indexFile.read((char*)(&loc), sizeof(loc));
	packFile.seekg(std::streampos(loc));
	return packFile.good() ? Configuration(packFile) : Configuration();
}
void ConfigurationPack::Clear(void)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	this->NConfig = 0;
}

// Implementations of member functions of the SpherePacking class
// for First-Passage-Time simulations
int SpherePacking::IsInsideParticle(const GeometryVector & RelativeCoordinate) const
{
	int idx = -1;
	bool inside = false;
	const std::vector<double> * pSphereRadii = &this->ParticleCharacteristics;

	this->IterateThroughNeighbors(RelativeCoordinate, this->GetMaxRadius(), 
	[&pSphereRadii, &inside, &idx, this] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void{
		//GeometryVector pos_Cartesian = this->RelativeCoord2CartesianCoord(shift);
		
//		double r = pSphereRadii->at(Sourceparticle)
//		double R = std::sqrt(shift.Modulus2());

		if(pSphereRadii->at(Sourceparticle)*pSphereRadii->at(Sourceparticle) > shift.Modulus2() ){
			idx = Sourceparticle;
			inside = true;
		}
	});
	return idx;
}

void SpherePacking::FindNearestNeighbors (const GeometryVector & RelativeCoordinate, std::vector<GeometryVector> & list, GeometryVector & x0_Relative, GeometryVector & n2host) const{
	//TODO: fix this!
	size_t num = 1;
	for (size_t i=0; i < this->Dimension; i++){
		num *= 2;
	}
	//std::min((size_t)3, this->NumParticle());
	double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension());
	std::vector<size_t> neighbors;
	std::vector<GeometryVector> shifts;
	// search neighbors
	while (neighbors.size() < num)
	{
		neighbors.clear();
		shifts.clear();
		(*this).IterateThroughNeighbors(RelativeCoordinate, TypicalLength, [&shifts, &neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
		{	
			neighbors.push_back(Sourceparticle);
			shifts.push_back(shift);
		});
		TypicalLength *= 1.5;
	}

	list.resize(0);
	//std::vector<GeometryVector> origins, ns;
	double smallest_distance = 0.0;

	//search for the closest points to each neighbor
	for (size_t i = 0; i < neighbors.size(); i++){
		size_t prt_idx = neighbors[i];
		GeometryVector pos_wrt_center = -1.0* shifts[i];	

		std::vector<GeometryVector> list_i;
		GeometryVector x0, n;
		double r = std::sqrt(pos_wrt_center.Modulus2());
		list_i.emplace_back(std::abs(r-this->GetCharacteristics(prt_idx)), 0.0, prt_idx);
		list_i.emplace_back(std::abs(r+this->GetCharacteristics(prt_idx)), 0.0, prt_idx);

		// find the smallest distance to neighboring particles.
		if (i==0)
			smallest_distance = 1.0 + list_i[0].x[0];
		
		if ( 
			(list_i[0].x[0]<list_i[1].x[0])
		&& (list_i[0].x[0] < smallest_distance) 
		){
			smallest_distance = list_i[0].x[0];
			n2host = 1.0/r*pos_wrt_center ;
			x0_Relative = (1.0 - this->GetCharacteristics(prt_idx)/r) * shifts[i];	//In Cartesian coordinates.
		}
		else if(
			(list_i[0].x[0] > list_i[1].x[0])
		&& (list_i[1].x[0] < smallest_distance)
		){
			smallest_distance = list_i[1].x[0];
			n2host = -1.0/r*pos_wrt_center ;
			x0_Relative = (1.0 + this->GetCharacteristics(prt_idx)/r) * shifts[i];
			// In Cartesian coordinates.
		}
		list.insert(list.end(), list_i.begin(), list_i.end());
	}

	// sort the list of distances...
	std::sort(list.begin(), list.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.x[0] < right.x[0];} );
	list.resize(3);

	// up to here, x0_Relative = the Cartesian coordinates of interface point closest to a test point with respect to this test point. 
	x0_Relative = this->CartesianCoord2RelativeCoord(x0_Relative) + RelativeCoordinate;
	this->RelativeCoordToMinimumImage(x0_Relative);
}


/* ConfigurationPack -> SpherePack
void SpherePackingPack::Open(const std::string &prefix) {
	this->IndexName = prefix + ".SpherePackIndex";
	this->PackName = prefix + ".SpherePack";
	std::fstream ifile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	if (ifile.good())
	{
		ifile.read((char*)(&NConfig), sizeof(NConfig));
		if (ifile.eof())
			NConfig = 0;
	}
	else
		NConfig = 0;
}
long long SpherePackingPack::NumConfig(void)
{
	return NConfig;
}
void SpherePackingPack::AddConfig(const SpherePacking & c)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary);
	if (indexFile.good() == false)
	{
		indexFile.clear();
		indexFile.close();
		indexFile.clear();
		indexFile.open(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	}
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::ate);
	if (packFile.good() == false)
	{
		packFile.clear();
		packFile.close();
		packFile.clear();
		packFile.open(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	}
	this->NConfig++;
	indexFile.seekp(0, std::fstream::beg);
	indexFile.write((char*)(&NConfig), sizeof(NConfig));
	long long loc = packFile.tellp();
	indexFile.seekp(0, std::fstream::end);
	indexFile.write((char*)(&loc), sizeof(loc));
	c.WriteBinary(packFile);
	indexFile.close();
	if (indexFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Index file failing!\n";
	packFile.close();
	if (packFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Pack file failing!\n";
}
SpherePacking SpherePackingPack::GetConfig(long long i)
{
	if (i >= this->NConfig)
	{
		std::cerr << "Error in ConfigurationPack::GetConfig : Index is larger than number of configuration!\n";
		exit(1);
	}
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::binary);
	long long loc;
	indexFile.seekg(sizeof(long long)*(i + 1));
	indexFile.read((char*)(&loc), sizeof(loc));
	packFile.seekg(std::streampos(loc));
	return packFile.good() ? SpherePacking(packFile) : SpherePacking();
}
void SpherePackingPack::Clear(void)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	this->NConfig = 0;
}
*/
//polydisperse sphere packing
SpherePacking::SpherePacking(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList) : PeriodicCellList<double>(Dimension, BasisVectors, CellSize, UseSortedList)
{
	//std::cerr<<"SpherePacking CellSize="<<CellSize<<'\n';
	MaxRadius=0;
}

void SpherePacking::Insert(double radius, const GeometryVector & RelativeCoordinate)
{
	UpdateMaxRadius();
	this->ParticleCartesians.push_back(GeometryVector(this->Dimension));
	this->ParticleRelatives.push_back(GeometryVector(this->Dimension));
	this->ParticleCharacteristics.push_back(radius);
	for(DimensionType i=0; i<this->Dimension; i++)
	{
		assert( RelativeCoordinate.x[i]==RelativeCoordinate.x[i] );

		//NEED TO subtract the floor TWICE
		this->ParticleRelatives.back().x[i]=RelativeCoordinate.x[i]-std::floor(RelativeCoordinate.x[i]);
		this->ParticleRelatives.back().x[i]-=std::floor(this->ParticleRelatives.back().x[i]);
	}
	this->UpdateCartesianCoordinates(this->NumParticle()-1);
	this->Cells[this->GetIndex(this->NumParticle()-1)].push_back(this->NumParticle()-1);

	if(this->MaxRadius<radius)
		this->MaxRadius=radius;
}
//relative coord and halfsize 
//Radius : cartesian radius of the upcoming sphere
bool SpherePacking::CheckVoxelOverlap(const GeometryVector & cr, double halfsize, double Radius) const
{
	UpdateMaxRadius();
	bool Overlap=false;
	const GeometryVector * ba=this->BasisVector;
	DimensionType dim=this->Dimension;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius, [&Overlap, &halfsize, &ba, dim, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
	{
		//additional processing of shift
		GeometryVector ss(shift);
		for(DimensionType d=0; d<dim; d++)
			if(ss.Dot(ba[d])>0.0)
				ss.AddFrom( halfsize*ba[d]);
			else
				ss.MinusFrom( halfsize*ba[d]);

		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);
		if(ss.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	return Overlap;
}

//returns 0 if voxel is completely empty
//2 if voxel is fully occupied
//1 otherwise
int SpherePacking::CheckVoxelOverlap2(const GeometryVector & cr, double halfsize, double Radius) const
{		
	UpdateMaxRadius();
	bool Overlap=false;
	bool PartialOverlap=false;
	const GeometryVector * ba=this->BasisVector;
	DimensionType dim=this->Dimension;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	double maxVoxelDiagonal=this->MaxDiagonalLength()*halfsize;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius+maxVoxelDiagonal, [&maxVoxelDiagonal, &PartialOverlap, &Overlap, &halfsize, &ba, dim, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
	{
		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);

		double temp=OverlapRadius+maxVoxelDiagonal;
		if(shift.Modulus2()<temp*temp)
			PartialOverlap=true;
		//additional processing of shift
		GeometryVector ss(shift);
		for(DimensionType d=0; d<dim; d++)
			if(ss.Dot(ba[d])>0.0)
				ss.AddFrom( halfsize*ba[d]);
			else
				ss.MinusFrom( halfsize*ba[d]);

		if(ss.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	if(Overlap)
		return 2;
	else if(PartialOverlap)
		return 1;
	else
		return 0;
}


//relative coord and halfsize 
//Radius : cartesian radius of the upcoming sphere
bool SpherePacking::CheckOverlap(const GeometryVector & cr, double Radius) const
{
	UpdateMaxRadius();
	bool Overlap=false;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius, [&Overlap, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
	{
		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);
		if(shift.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	return Overlap;
}


//pixelizing the configuration to MeshSide pixes in each side
void DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> & occu)
{
	DimensionType dim = c.GetDimension();
	size_t MeshTotal = std::pow(MeshSide, dim);
	struct VicinityMesh
	{
		signed long DeltaIndex[::MaxDimension] = {};
	};

	//get the list of mesh points occupied by a sphere
	std::vector<VicinityMesh> mesh;

	GeometryVector newbas[::MaxDimension] = {};
	for (int i = 0; i<dim; i++)
		newbas[i] = c.GetBasisVector(i)*(1.0 / MeshSide);
	Configuration * pt = new Configuration(dim, newbas, ::MaxDistance, false);
	pt->Insert("A", GeometryVector(dim));
	pt->IterateThroughNeighbors(GeometryVector(dim), radius, [&](const GeometryVector & s, const GeometryVector & l, const signed long * p, size_t source) ->void
	{
		VicinityMesh temp;
		for (int i = 0; i<dim; i++)
			temp.DeltaIndex[i] = p[i];
		mesh.push_back(temp);
	});
	delete pt;

	occu.assign(MeshTotal, 0);
	signed long m = MeshSide;
	int num = c.NumParticle();
	for (int i = 0; i<num; i++)
	{
		GeometryVector rel = c.GetRelativeCoordinates(i);
		for (auto iter = mesh.begin(); iter != mesh.end(); ++iter)
		{
			size_t index = 0;
			for (int j = 0; j<dim; j++)
			{
				index *= MeshSide;
				signed long temp = std::floor(rel.x[j] * MeshSide);
				temp += iter->DeltaIndex[j];
				while (temp >= m)
					temp -= m;
				while (temp<0)
					temp += m;
				index += temp;
			}
			occu[index] = 1;
		}
	}
}
//calculate volume fraction if each particle is replaced with a hypersphere of certain radius
//do this by pixelizing the configuration to MeshSide pixes in each side
double Volume(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> * pbuffer)
{
	size_t MeshTotal = std::pow(MeshSide, c.GetDimension());
	size_t InsideCount=0;
	if (pbuffer == nullptr)
	{
		std::vector<VoxelType> occu;
		DigitizeConfiguration(c, radius, MeshSide, occu);
		for (auto iter = occu.begin(); iter != occu.end(); ++iter)
			InsideCount += *iter;
	}
	else
	{
		std::vector<VoxelType> & occu = *pbuffer;
		DigitizeConfiguration(c, radius, MeshSide, occu);
		for (auto iter = occu.begin(); iter != occu.end(); ++iter)
			InsideCount += *iter;
	}
	return c.PeriodicVolume()*InsideCount/MeshTotal;
}


#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <set>
//return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected.
double MinConnectDiameter(const Configuration & c, size_t a, size_t b)
{
	double clength = 0.0;
	double mlength = 2 * std::pow(c.PeriodicVolume() / c.NumParticle(), 1.0 / c.GetDimension());
	std::set<size_t> cluster;
	struct toSearchStruct
	{
		double dist;
		size_t i;
		toSearchStruct(double dist, size_t i) : dist(dist), i(i)
		{}
	};
	typedef boost::multi_index_container<
		toSearchStruct,
		boost::multi_index::indexed_by<
		boost::multi_index::ordered_unique<boost::multi_index::member<toSearchStruct, size_t, &toSearchStruct::i> >,
		boost::multi_index::ordered_non_unique<boost::multi_index::member<toSearchStruct, double, &toSearchStruct::dist> >
		>
	> toSearch_Set;
	toSearch_Set toAdd;
	toSearch_Set::nth_index<0>::type & iview = toAdd.get<0>();
	toSearch_Set::nth_index<1>::type & dview = toAdd.get<1>();

	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
	{
		if (cluster.find(SourceAtom) == cluster.end())
		{
			double dist = std::sqrt(shift.Modulus2());
			auto iter = iview.find(SourceAtom);
			if (iter == iview.end())
				toAdd.insert(toSearchStruct(dist, SourceAtom));
			else if (iter->dist > dist)
			{
				toSearchStruct temp(dist, SourceAtom);
				toAdd.replace(iter, temp);
			}
		}
	};
	cluster.insert(a);
	GeometryVector rel = c.GetRelativeCoordinates(a);
	c.IterateThroughNeighbors(rel, mlength, IterateFunction);
	for (;;)
	{
		if (toAdd.empty())
		{
			mlength *= 2;
			for (auto iter = cluster.begin(); iter != cluster.end(); iter++)
			{
				GeometryVector rel = c.GetRelativeCoordinates(*iter);
				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
			}
			//std::cout<<"Warning in MinConnectDiameter : resizing mlength. toAdd size="<<toAdd.size()<<'\n';
		}
		else
		{
			auto ibegin = dview.begin();
			toSearchStruct temp = *ibegin;
			dview.erase(ibegin);
			clength = std::max(clength, temp.dist);
			size_t np = temp.i;
			if (np == b)
				return clength;
			if (cluster.find(np) == cluster.end())
			{
				cluster.insert(np);
				GeometryVector rel = c.GetRelativeCoordinates(np);
				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
				//std::cerr<<"Explore particle "<<np<<'\n';
			}
		}
	}

	//shouldn't reach this line
	std::cerr << "Error in MinConnectDiameter\n";
	return 0.0;
}
//double oldMinConnectDiameter(const Configuration & c, size_t a, size_t b)
//{
//	double clength=0.0;
//	double mlength=1.5*std::pow(c.PeriodicVolume()/c.NumParticle(), 1.0/c.GetDimension());
//	std::set<size_t> cluster;
//	std::multimap<double, size_t> toAdd;
//	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
//	{
//		if(cluster.find(SourceAtom)==cluster.end())
//		{
//			double dist=std::sqrt(shift.Modulus2());
//			toAdd.insert(std::make_pair(dist, SourceAtom));
//		}
//	};
//	cluster.insert(a);
//	GeometryVector rel=c.GetRelativeCoordinates(a);
//	c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//	for(;;)
//	{
//		if(toAdd.empty())
//		{
//			mlength*=2;
//			for(auto iter=cluster.begin(); iter!=cluster.end(); iter++)
//			{
//				GeometryVector rel=c.GetRelativeCoordinates(*iter);
//				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//			}
//			//std::cout<<"Warning in MinConnectDiameter : resizing mlength. \n";
//		}
//		else
//		{
//			auto ibegin=toAdd.begin();
//			std::pair<double, size_t> temp=*ibegin;
//			toAdd.erase(ibegin);
//			clength=std::max(clength, temp.first);
//			size_t np=temp.second;
//			if(np==b)
//				return clength;
//			if(cluster.find(np)==cluster.end())
//			{
//				cluster.insert(np);
//				GeometryVector rel=c.GetRelativeCoordinates(np);
//				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//			}
//		}
//	}
//}

//return the minimum D such that if each particle is replaced with a sphere with diameter D, 
//particle n will connect to its periodic image in direction dir.
double PercolationDiameter(const Configuration & c, size_t n, DimensionType dir)
{
	std::vector<size_t> m(c.GetDimension(), 1);
	m[dir] = 2;
	Configuration t = MultiplicateStructure(c, m);
	return MinConnectDiameter(t, 2 * n, 2 * n + 1);
}

//Replace each particle in c with a hypersphere with radius,
//what's the cluster size InitParticle is in?
//FoundPeriodicImages: whether or not found more than 1 periodic images of a particle, 
//which means the cluster size is actually infinite rather than the returned value.
#include <list>
size_t ClusterSize(const Configuration & c, size_t InitParticle, double radius, bool & FoundPeriodicImages)
{
	FoundPeriodicImages = false;
	std::map<size_t, GeometryVector> Found;
	struct toSearchStruct
	{
		size_t i;
		GeometryVector LatticeShift;
		toSearchStruct(size_t i, const GeometryVector & LatticeShift) : i(i), LatticeShift(LatticeShift)
		{}
	};
	std::list<toSearchStruct> toSearch;
	GeometryVector alreadyLatticeShift = GeometryVector(c.GetDimension());
	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
	{
		auto iter = Found.find(SourceAtom);
		GeometryVector l = alreadyLatticeShift + LatticeShift;
		if (iter != Found.end())
		{
			if (!FoundPeriodicImages)
			{
				//this particle has been found before. Test if they belong to the same periodic cell
				GeometryVector deltaLatticeShift = l - iter->second;
				if (deltaLatticeShift.Modulus2() > LengthPrecision*LengthPrecision)
					FoundPeriodicImages = true;
			}
		}
		else
		{
			toSearch.push_back(toSearchStruct(SourceAtom, l));
			Found.insert(std::make_pair(SourceAtom, l));
		}
	};
	Found.insert(std::make_pair(InitParticle, GeometryVector(c.GetDimension())));
	c.IterateThroughNeighbors(c.GetRelativeCoordinates(InitParticle), radius, IterateFunction);
	while (!toSearch.empty())
	{
		toSearchStruct s = toSearch.front();
		toSearch.pop_front();
		alreadyLatticeShift = s.LatticeShift;
		c.IterateThroughNeighbors(c.GetRelativeCoordinates(s.i), radius, IterateFunction);
	}
	return Found.size();
}


void WriteConfiguration(const Configuration & conf, const std::string & OutputName)
{
	// dimension 0	0	0
	// basis vectors ...
	// Cartessian coordinates
	std::string name = OutputName + std::string(".txt");
	std::fstream output(name.c_str(), std::fstream::out);
	
	output << std::setprecision(17);
	output << conf.GetDimension() << "\n";
	for (int i = 0; i < conf.GetDimension(); i++) {
		output <<conf.GetBasisVector(i);
	}
	for (size_t i = 0; i < conf.NumParticle(); i++) {
		output << conf.GetCartesianCoordinates(i);
	}
	output.close();
}

void WriteConfiguration(const SpherePacking & conf, const std::string & OutputName)
{
	//dimension	0	0	0
	//basisvector1
	// ...
	//Relative coordinates1	radius 1
	//
	//....

	std::string name = OutputName + std::string(".txt");
	std::fstream ofile(name.c_str(), std::fstream::out);
	ofile << std::setprecision(10);
	ofile << conf.GetDimension() << "\n";
	for (int i = 0; i < conf.GetDimension(); i++) {
		ofile << conf.GetBasisVector(i);
	}
	//coordinates
	for (size_t i = 0; i < conf.NumParticle(); i++) {
		GeometryVector temp(conf.GetCartesianCoordinates(i));
		ofile << temp.x[0]<<"\t"<<temp.x[1]<<"\t"<<temp.x[2]<<"\t"<<temp.x[3]<<"\t"<<conf.GetCharacteristics(i)<<"\n";
	}
	ofile.close();
}

template<typename T> void WriteConfiguration(const PeriodicCellList<T> & conf, const std::string & OutputName)
{

	std::string name = OutputName + std::string(".txt");
	std::fstream ofile(name.c_str(), std::fstream::out);
	ofile << std::setprecision(10);
	ofile << conf.GetDimension() << "\n";
	for (int i = 0; i < conf.GetDimension(); i++) {
		ofile << conf.GetBasisVector(i);
	}
	//coordinates
	for (size_t i = 0; i < conf.NumParticle(); i++) {
		GeometryVector temp(conf.GetCartesianCoordinates(i));
		ofile << temp.x[0] << "\t" << temp.x[1] << "\t" << temp.x[2] << "\t" << temp.x[3] << "\t" << conf.GetCharacteristics(i) << "\n";
	}
	ofile.close();
}
//assume # columns = 4
void ReadConfiguration(Configuration & conf, const std::string & InputName) {
	
	std::string name = InputName + std::string(".txt");
	std::fstream ifile(name.c_str(), std::fstream::in);
	
	if (ifile.good() == false)
	{
		std::cerr << "Error in ReadFunction : input stream is not good!\n" << InputName << "\n";
		return;
	}

	
	DimensionType d = 0;
	GeometryVector * basis;
	// set dimension
	ifile >> d;
	basis = new GeometryVector[d];

	// set basis vectors
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			ifile >> basis[i].x[j];
		}
		ifile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		basis[i].SetDimension(d);
	}
	//initialize conf
	conf = Configuration(d, basis, 1.0);
	while (ifile.eof() == false)
	{
		char dump[255] = {};		
		//input particle positions
		GeometryVector temp(static_cast<int> (d));
		for (size_t i = 0; i < ::MaxDimension; i++) {
			if (i < d) {
				ifile >> temp.x[i];
			}
			else {
				ifile.getline(dump, 255); i = ::MaxDimension;
			}
		}

		if (ifile.eof() == false)
			conf.Insert("a", conf.CartesianCoord2RelativeCoord(temp));
	}

	//std::cout << "the number of particles = " << conf.NumParticle() << "\n";

}

//assume # of columns is 5
void ReadConfiguration(SpherePacking & conf, const std::string & InputName) {
	std::string name = InputName + std::string(".txt");
	std::fstream ifile(name.c_str(), std::fstream::in);

	if (ifile.good() == false)
	{
		std::cerr << "Error in ReadFunction : input stream is not good!\n"<<InputName<<"\n";
		return;
	}
	
	DimensionType d = 0;
	GeometryVector * basis;
	// set dimension
	ifile >> d;
	basis = new GeometryVector[d];

	// set basis vectors
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			ifile >> basis[i].x[j];
		}
		ifile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		basis[i].SetDimension(d);
	}
	//initialize conf
	conf = SpherePacking(d, basis, 1);
	while (ifile.eof() == false)
	{
		char dump[255] = {};
		//input particle positions
		GeometryVector temp(static_cast<int> (d));
		double radius = 0;
		for (size_t i = 0; i < ::MaxDimension; i++) {
			if (i < d) {
				ifile >> temp.x[i];
			}
			else {
				ifile >> dump;
			}
		}
		ifile >> radius;

		if (ifile.eof() == false)
			conf.Insert(radius, conf.CartesianCoord2RelativeCoord(temp));
	}
	
	
	/*std::vector<GeometryVector> temp;
	ReadFunction(temp, InputName);

	DimensionType d = (DimensionType)temp[0].x[0];
	GeometryVector *basis;
	basis = new GeometryVector[d];
	for (int i = 0; i < d; i++) {
		basis[i] = temp[i + 1];
		basis[i].SetDimension(d);
	}

	conf = SpherePacking(d, basis, 1.0);
	int numPrt = (temp.size() - (1 + d)) / 2;
	for (int i = 0; i < numPrt; i++) {
		conf.Insert(temp[2 * i + d + 2].x[0], temp[2 * i + d + 1]);
	}

	std::cout << "the number of particles = " << conf.NumParticle() << "\n";*/

}


std::vector<GeometryVector> GetKs(const PeriodicCellList<Empty> & tempList, double CircularKMax, double LinearKMax, double SampleProbability, bool use_iteratethrough, double CircularKMin)
{
	//Exclude probablematic cases.
	if (CircularKMax <= 0.0) {
		std::cerr << "GetKs():: CircularKMax should be positive." << std::endl;
		return std::vector<GeometryVector>();
	}
	if (LinearKMax <= 0.0) {
		std::cerr << "GetKs():: Linear KMax should be positive." << std::endl;
		return std::vector<GeometryVector>();
	}

	RandomGenerator gen(98765);
	std::vector<GeometryVector> results;
	DimensionType dim = tempList.GetDimension();
	if (dim == 0) {
		std::cerr << "Dimension of PeriodicCellList must be positive!\n";
		exit(1);
	}
	else if (dim > 1)
	{
		std::vector<GeometryVector> bs;
		for (DimensionType i = 0; i<dim; i++)
			bs.push_back(tempList.GetReciprocalBasisVector(i));
		
		//a periodic list of reciprocal lattice
		PeriodicCellList<Empty> reciprocal(dim, &bs[0], std::sqrt(bs[0].Modulus2()));
		reciprocal.Insert(Empty(), GeometryVector(dim));
		std::vector<GeometryVector> ks;
		std::vector<GeometryVector> preKs;
		double KMin2 = CircularKMin*CircularKMin;

		if (use_iteratethrough)
		{
			reciprocal.IterateThroughNeighbors(GeometryVector(dim), CircularKMax, [&ks, &dim, &KMin2](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void{
			bool Trivial = true;
			for (DimensionType d = 0; d<dim; d++)
			{
				if (PeriodicShift[d]>0)
				{
					Trivial = false;
					break;
				}
				else if (PeriodicShift[d]<0)
				{
					Trivial = true;
					break;
				}
			}
			if (!Trivial && (shift.Modulus2() > KMin2))
				ks.push_back(shift);
			});
		}
		else{
			/* Get wavevectors without IterateThroughNeighbors */
			size_t num_ks = (size_t) round(0.5*HyperSphere_Volume(dim,CircularKMax) / reciprocal.PeriodicVolume());
			ks.reserve(num_ks);

			std::vector<int> limits;
			std::vector<int> PeriodicShift;
			size_t iter = 1;
			for (size_t i =0; i < dim; i++){
				int val = int(std::ceil( CircularKMax* std::sqrt( tempList.GetBasisVector(i).Modulus2()) /(2.*::pi) ));
				limits.push_back(val);
				PeriodicShift.push_back(-val);
				iter *= (2*val + 1);
			}

			/* iterate over all reciprocal lattice points */
			double K2 = CircularKMax*CircularKMax;
			for (size_t i=0; i < iter; i++){
				bool Trivial = false;
				
				for (DimensionType d = 0; d<dim; d++)
				{
					if (PeriodicShift[d]>0)
					{
						Trivial = false;
						break;
					}
					else if (PeriodicShift[d]<0)
					{
						Trivial = true;
						break;
					}
				}
				if (!Trivial){
					/* Ignore inverted vectors */
					GeometryVector k(dim);
					for (DimensionType d = 0; d<dim; d++)
					{
						k = k + PeriodicShift[d] * bs[d];
					} 
					double k2 = k.Modulus2();
					if ((k2 < K2) && (k2 > KMin2))
						ks.emplace_back(k);
				}

				/* Change PeriodicShift */
				PeriodicShift[0] ++;
				for (DimensionType d = 0; d<dim; d++)
				{	
					if (PeriodicShift[d] > limits[d]){
						PeriodicShift[d] = - limits[d];
						PeriodicShift[d+1] ++ ; 
						/* Note that PeriodicShift[dim-1] cannot exceed limits[dim-1] */ 
					}
				} 
			}	
		}
		std::sort(ks.begin(), ks.end(), [](const GeometryVector & left, const GeometryVector & right) ->bool {return (left.Modulus2()<right.Modulus2()) && (left.x[0]<right.x[0]); });
		//ks contains K points in that circle

		if (LinearKMax == CircularKMax)
			return ks;

		std::vector<GeometryVector> tempBase;
		for (DimensionType i = 0; i<dim; i++)
		{
			GeometryVector t(dim);
			t.x[i] = 2;
			tempBase.push_back(t);
		}
		PeriodicCellList<Empty> KDirection(dim, &tempBase[0], std::sqrt(tempBase[0].Modulus2())*std::pow(1.0 / ks.size(), 1.0 / (dim))*10.0);//use this list to discover k points of the same direction
		for (auto iter = ks.begin(); iter != ks.end(); iter++)
		{
			double Length = std::sqrt(iter->Modulus2());
			if (Length == 0)
				continue;
			GeometryVector dir(dim);
			for (DimensionType i = 0; i<dim; i++)
				dir.x[i] = iter->x[i] / Length / 2.0;

			bool NeighborFound = false;
			KDirection.IterateThroughNeighbors(dir, 1e-10, [&NeighborFound](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void{NeighborFound = true; });
			if (NeighborFound == false)
			{
				preKs.push_back(*iter);
				KDirection.Insert(Empty(), dir);
			}
		}
		//preKs contains K points in that circle with different directions

		ks.clear();
		for (auto iter = preKs.begin(); iter != preKs.end(); iter++)
		{
			for (size_t i = 1;; i++)
			{
				GeometryVector temp = static_cast<double>(i)*(*iter);
				if (temp.Modulus2()>LinearKMax*LinearKMax)
					break;
				else
					if(SampleProbability>=1.0 || gen.RandomDouble()<SampleProbability)
						results.push_back(temp);
			}
		}
		preKs.clear();
	}
	else
	{
		/* 1D case */
		//std::cout<< "GetKs(...):: 1D case. \n";
		GeometryVector b = tempList.GetReciprocalBasisVector(0);
		if (CircularKMin <= 0.0){
			for (double i = 1.0;; i += 1.0)
			{
				GeometryVector k = i*b;
				if (k.x[0] > LinearKMax)
					break;
				if (SampleProbability >= 1.0 || gen.RandomDouble()<SampleProbability)
					results.push_back(k);
			}
		}
		else{
			/* search wavevectors between CircularKMin and CircularKMax */
			//std::cout<< "GetKs(...):: CircularKMin is positive. \n";
			size_t num_inside_circular = 0;
			for (double i = 1.0;; i += 1.0)
			{
				GeometryVector k = i*b;
				if (k.x[0] > CircularKMax)
					break;
				if (k.x[0] > CircularKMin){
					if (SampleProbability >= 1.0 || gen.RandomDouble()<SampleProbability)
						results.push_back(k);
					num_inside_circular ++;
				}
			}
			
			if (results.size() > 0){
				/* consider integer multiples in results. */
				bool OutsideRange = false;
				for (double i = 2.0;; i+= 1.0){
					for (size_t j = 0; j<num_inside_circular; j++){
						GeometryVector k = i * results[j];
						if (k.x[0] > LinearKMax){
							OutsideRange = true;
							break;
						}
						if (SampleProbability >= 1.0 || gen.RandomDouble()<SampleProbability)
							results.push_back(k);
					}
					if (OutsideRange)
						break;
				}
			}
			else{
				std::cerr << "GetKs(...):: There is no wavevector bewteen CircularKMin and CircularKMax " <<std::endl;
			}

		}
	}
	return results;
}


#ifdef VORONOI_INCLUDED
#include "Voronoi.h"
#include <gsl/gsl_sf_legendre.h>
double LegenderSphericalPlm(int l, int m, double x)
{
	if (m >= 0)
		return gsl_sf_legendre_sphPlm(l, m, x);
	else
	{
		int mm = (-1)*m;
		double result = gsl_sf_legendre_sphPlm(l, mm, x);
		if (mm % 2 == 1)
			result *= (-1.0);
		for (int t = l - mm + 1; t <= l + mm; t++)
			result /= t;
		return result;
	}
}




#endif
