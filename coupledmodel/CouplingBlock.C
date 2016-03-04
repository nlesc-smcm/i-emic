#include "CouplingBlock.H"
#include "GlobalDefinitions.H"
#include "Utils.H"

//------------------------------------------------------------------
CouplingBlock::CouplingBlock()
	:
	initialized_(false)
{}

//------------------------------------------------------------------
// constructor
CouplingBlock::CouplingBlock(std::string const &coupling,
							 std::vector<double> const &values,
							 std::vector<int> const &row_ind,
							 std::vector<int> const &col_ind)
	:
	coupling_   (coupling),
	values_     (values),
	row_ind_    (row_ind),
	col_ind_    (col_ind),
	initialized_(false)
{
	unique_row_ind_ = find_unique(row_ind);
	unique_col_ind_ = find_unique(col_ind);
}

//------------------------------------------------------------------
// default constructor
void CouplingBlock::initialize()
{
}

//------------------------------------------------------------------
std::vector<int> CouplingBlock::find_unique(std::vector<int> const &in)
{
	std::vector<int> tmp(values_.size(), 0);
	for (auto &i: in)
		tmp[i] = 1;

	std::vector<int> out;
	for (int i = 0; i != (int) tmp.size(); ++i)
		if (tmp[i])
			out.push_back(i);
	
	return out;		
}

//------------------------------------------------------------------
void CouplingBlock::initializeOA(Epetra_BlockMap const &map)
{
	// --> here we should use the UNIQUE col indices, not this array
	indexMap_  = Utils::CreateSubMap(map, col_ind_);
	gatherMap_ = Utils::AllGather(*indexMap_);
	restrVec_  = Teuchos::rcp(new Epetra_Vector(*indexMap_));
	restrImp_  = Teuchos::rcp(new Epetra_Import(*indexMap_, map));
	gatherImp_ = Teuchos::rcp(new Epetra_Import(*gatherMap_, *indexMap_));
	gathered_  = Teuchos::rcp(new Epetra_MultiVector(*gatherMap_, 1));

	initialized_ = true;
}

//------------------------------------------------------------------
void CouplingBlock::initializeAO()
{
	// I don't believe we need any initializing here
	initialized_ = true;
}

//------------------------------------------------------------------
void CouplingBlock::setValues(std::vector<double> const  &values)
{
	values_ = values;
}

//------------------------------------------------------------------
void CouplingBlock::setColInd(std::vector<int> const &col_ind)
{
	col_ind_ = col_ind;
}

//------------------------------------------------------------------
void CouplingBlock::setRowInd(std::vector<int> const &row_ind)
{
	row_ind_ = row_ind;
}

//------------------------------------------------------------------
void CouplingBlock::applyMatrix(SuperVector const &v, SuperVector &out)
{
	if (coupling_.compare("OA") == 0)
		applyOA(v, out);
	else if (coupling_.compare("AO") == 0)
		applyAO(v, out);
	else
	{
		ERROR("Invalid coupling mode!", __FILE__, __LINE__);
	}
}

//------------------------------------------------------------------
void CouplingBlock::applyOA(SuperVector const &in, SuperVector &out)
{
	TIMER_START("CouplingBlock: applyOA");
	if (!initialized_) // initialize
		initializeOA(in.getOceanVector()->Map());
	
	// obtain restricted vector
	restrVec_->Import(*in.getOceanVector(), *restrImp_, Insert);

	// gather vector
	gathered_->Import(*restrVec_, *gatherImp_, Insert);

	// get the ocean values
	std::vector<double> oceanValues(col_ind_.size(), 0.0);
	gathered_->ExtractCopy(&oceanValues[0], col_ind_.size());

	// compute the new atmosphere values (sparse matrix-vector mult)
	out.zeroAtmos(); // clear
	std::shared_ptr<std::vector<double> > result = out.getAtmosVector();
	
	for (size_t i = 0; i != values_.size(); ++i)
		(*result)[row_ind_[i]] += values_[i] * oceanValues[i]; 

	// the above is wrong
	// instead of oceanValues[i] we nee oceanValues[restricted_col_ind_[i]]
	
	TIMER_STOP("CouplingBlock: applyOA");
}

//------------------------------------------------------------------
void CouplingBlock::applyAO(SuperVector const &in, SuperVector &out)
{
	TIMER_START("CouplingBlock: applyAO");
	
	if (!initialized_) // initialize
		initializeAO();
	
	out.zeroOcean();

	std::vector<double> result(values_.size(),0.0);
	std::shared_ptr<std::vector<double> > atmosValues =	in.getAtmosVector();
	
	// we need a restricted_row_ind_[i] in result[i]
	for (size_t i = 0; i != values_.size(); ++i)
		result[i] += values_[i] * (*atmosValues)[col_ind_[i]];

	out.getOceanVector()->ReplaceGlobalValues(
		atmosValues->size(),
		&result[0], &row_ind_[0]);
	
	TIMER_STOP("CouplingBlock: applyAO");
}

//------------------------------------------------------------------
void CouplingBlock::info()
{
	INFO("Coupling block: " << coupling_);
	INFO(" |     #values: " << values_.size());
}
