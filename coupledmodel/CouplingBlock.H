#include <vector>
#include <string>

class CouplingBlock
{
	// Storage for sparse coupling block
	std::vector<double> val_;
	std::vector<int>    col_ind_;
	std::vector<int>    row_ind_;

	// string defining the coupling type
	std::string coupling_;

public:
	// Constructor
	CouplingBlock(std::string const &coupling);

	// Destructor
	~CouplingBlock(){};

	// Fill our storage
	void setValues(std::vector<double> const  &val)  {val_     = val    ;}
	void setColInd(std::vector<int> const &col_ind)  {col_ind_ = col_ind;}
	void setRowInd(std::vector<int> const &row_ind)  {row_ind_ = row_ind;}

	// Apply the coupling block
	void applyMatrix(SuperVector const &v, SuperVector &out);

private:
	void initialize();
	void applyOA(SuperVector const &v, SuperVector &out);
	void applyAO(SuperVector const &v, SuperVector &out);
};
