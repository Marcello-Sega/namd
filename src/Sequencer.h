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

#include "converse.h"

extern class HomePatch;

class Sequencer
{
public:
	Sequencer(HomePatch *p) : patch(p) { };
	void run(int numberOfCycles);             // spawn thread, etc.
	void awaken(void){ CthAwaken(thread); };

protected:
	virtual void threadRun(void);  // subclasses redefine this method
	int numberOfCycles;            // stores argument to run()
	const HomePatch *patch;        // access methods in patch

private:
	CthThread thread;
	friend void Sequencer_threadRun(Sequencer*);
};

#endif // SEQUENCER_H

