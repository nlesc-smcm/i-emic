#include "CouplingBlock.H"
#include "GlobalDefinitions.H"
#include "Utils.H"

//------------------------------------------------------------------
// constructor
CouplingBlock::CouplingBlock(std::string const &coupling,
							 std::vector<double> const &values,
							 std::vector<int> const &row_ind,
							 std::vector<int> const &col_ind)
	:
	coupling_  (coupling),
	values_    (values),
	row_ind_   (row_ind),
	col_ind_   (col_ind)
{
}

void CouplingBlock::initialize()
{
}

//------------------------------------------------------------------
void CouplingBlock::initializeOA(Epetra_BlockMap const &map)
{
	// --> here we should use the UNIQUE col indices, not this array
	indexMap_ = Utils::CreateSubMap(map, col_ind_);
	restrVec_ = Teuchos::rcp(new Epetra_Vector(*indexMap_));
	restrImp_ = Teuchos::rcp(new Epetra_Import(*indexMap_, map));


	initialized_ = true;
}

//------------------------------------------------------------------
void CouplingBlock::initializeAO()
{
	
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
	TIMER_START("CouplingBlock: applyOA restrict vector");
	restrVec_->Import(*in.getOceanVector(), *restrImp_, Insert);
	TIMER_STOP("CouplingBlock: applyOA restrict vector");
	
	// gather the restricted vector
	TIMER_START("CouplingBlock: applyOA gather");
	gathered_ = Utils::AllGather(*restrVec_);
	TIMER_STOP("CouplingBlock: applyOA gather");

	// get the ocean values
	std::vector<double> oceanValues(col_ind_.size(), 0.0);
	gathered_->ExtractCopy(&oceanValues[0], col_ind_.size());

	// compute the new atmosphere values (sparse matrix-vector mult)
	out.zeroAtmos(); // clear
	std::shared_ptr<std::vector<double> > result = out.getAtmosVector();
	
	for (size_t i = 0; i != values_.size(); ++i)
		(*result)[row_ind_[i]] += values_[i] * oceanValues[i]; 

	// this is wrong
	// instead of fullSol[i] we nee fullSol[restricted_col_ind_[i]]
	TIMER_STOP("CouplingBlock: applyOA");
}

//------------------------------------------------------------------
void CouplingBlock::applyAO(SuperVector const &v, SuperVector &out)
{
	
}

//------------------------------------------------------------------
void CouplingBlock::info()
{
	INFO("Coupling block: " << coupling_);
	INFO(" |     #values: " << values_.size());
}
