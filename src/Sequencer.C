/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include "Sequencer.h"

Sequencer::thread_initForce(void)
{
	patch->positionsReady();
	Cth_Suspend();
}

Sequencer::thread_run(int number_of_cycles, int steps_per_cycle)
{
	int step, cycle;
	for ( cycle = 0; cycle < number_of_cycles; ++cycle )
	{
		for ( step = 0; step < steps_per_cycle; ++step )
		{
			patch->addForceToMomentum(0.5*timestep);
			patch->addVelocityToPosition(timestep);
			patch->positionsReady();
			Cth_Suspend();
			patch->addForceToMomentum(0.5*timestep);
		}
	}
}
