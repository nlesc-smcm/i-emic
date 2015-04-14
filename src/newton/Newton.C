#include "Newton.H"

Newton::Newton(ModelType model, Statetype state_, char *parameterFile)
	:
	model_(model),
	state_(state)
{
	INFO("Entering Newton constructor");
	INFO("Leaving Newton constructor");
}

void Newton::solve()
{
	INFO("Entering Newton::solve()");
	INFO("Leaving Newton::solve()");
}
