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

#ifndef SEQUENCER_H
#define SEQUENCER_H

class Sequencer
{
public:
	Sequencer(HomePatch *p) { patch = p };
	void awaken(void)
	{
		Cth_Awaken(thread);
	}
	void initForce(void)
	{
		thread = Cth_Create(this, thread_initForce());
	}
	void run(int number_of_cycles, int steps_per_cycle)
	{
		thread = Cth_Create(this,
			thread_run(number_of_cycles, steps_per_cycle));
	}

private:
	HomePatch *patch;
	Cthread thread;
	virtual void thread_initForce(void);
	virtual void thread_run(int number_of_cycles, int steps_per_cycle);
};

#endif // SEQUENCER_H

