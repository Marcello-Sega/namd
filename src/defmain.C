extern "C" {
#include "converse.h"
}

#include "ProcessorPrivate.h"

extern void _initCharm(int, char**);

void charm_init(int argc, char **argv)
{
  ProcessorPrivateInit();
  _initCharm(argc, argv);
}

int main(int argc, char **argv)
{
  ConverseInit(argc, argv, charm_init,0,0);
}

