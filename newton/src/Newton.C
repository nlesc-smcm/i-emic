#include "GlobalDefinitions.H"
#include "Newton.H"

template<typename Model, typename Vector>
Newton<Model, Vector>::Newton(Model model, char *parameterFile)
	:
	model_(model),
	isInitialized(false)
{
	INFO("Entering Newton constructor");
	INFO("Leaving Newton constructor");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::initialize()
{
	INFO("Entering Newton::initialize()");
	//
	isInitialized = true;
	INFO("Leaving Newton::initialize()");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::run()
{
	INFO("Entering Newton::run()");
	INFO("Leaving Newton::run()");
}

template<typename Model, typename Vector>
void Newton<Model, Vector>::step()
{
	INFO("Entering Newton::step()");
	if (!isInitialized) initialize();
	
	INFO("Leaving Newton::step()");
}
