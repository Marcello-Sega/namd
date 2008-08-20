/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef SEQUENCER_H
#define SEQUENCER_H

#include "converse.h"
#include "PatchTypes.h"

class HomePatch;
class SimParameters;
class SubmitReduction;
class CollectionMgr;
class ControllerBroadcasts;
class LdbCoordinator;
class Random;

class Sequencer
{
    friend class HomePatch;
public:
    Sequencer(HomePatch *p);
    virtual ~Sequencer(void);
    void run(void);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); }
    void suspend(void) { CthSuspend(); }

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void integrate(); // Verlet integrator
    void minimize(); // CG minimizer

    void runComputeObjects(int migration = 1, int pairlists = 0);
    int pairlistsAreValid;
    int pairlistsAge;

    void submitReductions(int);
    void submitHalfstep(int);
    void submitMinimizeReductions(int);
    void submitCollections(int step, int zeroVel = 0);

    void submitMomentum(int step);
    void correctMomentum(int step, BigReal drifttime);

    void saveForce(const int ftag = Results::normal);
    void addForceToMomentum(BigReal, const int ftag = Results::normal,
						const int useSaved = 0, const int pressure = 0);
    void addVelocityToPosition(BigReal);
    
    // For TIP4
    void redistrib_tip4p_forces(const int, const int);

    void addRotDragToPosition(BigReal);
    void addMovDragToPosition(BigReal);

    void newMinimizeDirection(BigReal);
    void newMinimizePosition(BigReal);
    void quenchVelocities();

    void rattle1(BigReal,int);
    void rattle2(BigReal,int);

    void maximumMove(BigReal);
    void minimizationQuenchVelocity(void);

    void reloadCharges();

    void rescaleVelocities(int);
    int rescaleVelocities_numTemps;
    void reassignVelocities(BigReal,int);
    void reinitVelocities(void);
    void rescaleVelocitiesByFactor(BigReal);
    void tcoupleVelocities(BigReal,int);
    void berendsenPressure(int);
      int berendsenPressure_count;
      int checkpoint_berendsenPressure_count;
    void langevinPiston(int);
      int slowFreq;
    void langevinVelocities(BigReal);
    void langevinVelocitiesBBK1(BigReal);
    void langevinVelocitiesBBK2(BigReal);

    void cycleBarrier(int,int);
    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    HomePatch *const patch;		// access methods in patch
    SubmitReduction *reduction;
    SubmitReduction *pressureProfileReduction;

    CollectionMgr *const collection;
    ControllerBroadcasts * broadcast;

    int ldbSteps;
    void rebalanceLoad(int timestep);


private:
    CthThread thread;
    static void threadRun(Sequencer*);

    LdbCoordinator *ldbCoordinator;
};

#endif
