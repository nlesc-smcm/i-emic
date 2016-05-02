#ifndef MYPARAMETERLIST_H
#define MYPARAMETERLIST_H

//------------------------------------------------------------------
// Simple parameterList for testing purposes
class MyParameterList
{
	std::map<std::string, double>  list_;

public:

	// constructor
	MyParameterList() {}
	
	MyParameterList(std::map<std::string, double> const &other)
		:
		list_(other)
		{}

	template<typename T>
	T get(std::string const &key, T value)
		{
			if (list_.find(key) != list_.end())				
				return list_[key];
			else
				return value;
		}

	void set(std::string const &key, double value)
		{
			list_[key] = value;
		}	
};

#endif
