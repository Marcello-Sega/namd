#include "Priorities.h"

int Priorities::numBits = 16; // 32 causes crash

int Priorities::comm_urgent = 0;
int Priorities::comm_high = 1;
int Priorities::comm_default = 2;
int Priorities::comp_nonlocal_base = 3;
int Priorities::comp_nonlocal_range = 100;
int Priorities::comp_default = Priorities::comp_nonlocal_base;
int Priorities::comm_low =
	Priorities::comp_nonlocal_base + Priorities::comp_nonlocal_range;
int Priorities::comp_local_large = Priorities::comm_low + 1;
int Priorities::comp_local_base = Priorities::comp_local_large + 1;
int Priorities::comp_local_range = 100;
int Priorities::comp_synchronizing =
	Priorities::comp_local_base + Priorities::comp_local_range;
