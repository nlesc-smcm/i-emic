#include "Newton.H"

Newton::Newton(ModelType model, char *parameterFile)
	:
	model_(model),
	isInitialized(false),
{
	INFO("Entering Newton constructor");
	INFO("Leaving Newton constructor");
}

void Newton::initialize()
{
	INFO("Entering Newton::initialize()");
	//
	isInitialized = true;
	INFO("Leaving Newton::initialize()");
}

void Newton::run()
{
	INFO("Entering Newton::run()");
	INFO("Leaving Newton::run()");
}

void Newton::step()
{
	INFO("Entering Newton::step()");
	if (!isInitialized) initialize();
	
	INFO("Leaving Newton::step()");
}

