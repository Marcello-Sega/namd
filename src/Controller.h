/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "converse.h"
#include "Node.h"
#include "common.h"
#include <fstream.h>

class ControllerBroadcasts;
class NamdState;
class SimParameters;
class RequireReduction;
class CollectionMaster;
class Random;

class Controller
{
public:
    Controller(NamdState *s);
    virtual ~Controller(void);
    void run(void);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); };

protected:
    friend class ScriptTcl;
    int scriptSeq;
    virtual void algorithm(void);	// subclasses redefine this method

    void receivePressure(int);
      Tensor pressure_normal;
      Tensor pressure_nbond;
      Tensor pressure_slow;
      Tensor groupPressure_normal;
      Tensor groupPressure_nbond;
      Tensor groupPressure_slow;
      Tensor controlPressure_normal;
      Tensor controlPressure_nbond;
      Tensor controlPressure_slow;
      int nbondFreq;
      int slowFreq;
    void printEnergies(int);
      int computeChecksum;
      int numDegFreedom;
      BigReal electEnergy;
      BigReal electEnergySlow;
      BigReal ljEnergy;
      BigReal kineticEnergy;
      BigReal temperature;
      Tensor pressure;
      Tensor groupPressure;
      int controlNumDegFreedom;
      Tensor controlPressure;
    void enqueueCollections(int);
    void rescaleVelocities(int);
      BigReal rescaleVelocities_sumTemps;
      int rescaleVelocities_numTemps;
    void reassignVelocities(int);
    void tcoupleVelocities(int);
    void berendsenPressure(int);
    void langevinPiston1(int);
    void langevinPiston2(int);
      Tensor langevinPiston_strainRate;

    int ldbSteps;
    void rebalanceLoad(int);
    void cycleBarrier(int,int);

    // void suspend(void) { CthSuspend(); };
    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    NamdState *const state;		// access data in state
    RequireReduction *reduction;
    CollectionMaster *const collection;
    ControllerBroadcasts * broadcast;
    ofstream xstFile;
    void outputExtendedSystem(int step);
    void writeExtendedSystemLabels(ofstream &file);
    void writeExtendedSystemData(int step, ofstream &file);

private:
    CthThread thread;
    static void threadRun(Controller*);

    double startCTime;
    double startWTime;
    double startBenchTime;

};

#endif // CONTROLLER_H

