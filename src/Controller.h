//-*-c++-*-
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
    virtual void algorithm(void);	// subclasses redefine this method

    void receivePressure(int);
      Vector pressure_normal;
      Vector pressure_nbond;
      Vector pressure_slow;
      Vector groupPressure_normal;
      Vector groupPressure_nbond;
      Vector groupPressure_slow;
      Vector controlPressure_normal;
      Vector controlPressure_nbond;
      Vector controlPressure_slow;
      int nbondFreq;
      int slowFreq;
    void printEnergies(int);
      int computeChecksum;
      int numDegFreedom;
      BigReal kineticEnergy;
      BigReal temperature;
      Vector pressure;
      Vector groupPressure;
      int controlNumDegFreedom;
      Vector controlPressure;
    void enqueueCollections(int);
    void rescaleVelocities(int);
      BigReal rescaleVelocities_sumTemps;
      int rescaleVelocities_numTemps;
    void reassignVelocities(int);
    void tcoupleVelocities(int);
    void berendsenPressure(int);
    void langevinPiston1(int);
    void langevinPiston2(int);
      Vector langevinPiston_strainRate;

    // void suspend(void) { CthSuspend(); };
    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    NamdState *const state;		// access data in state
    RequireReduction *reduction;
    CollectionMaster *const collection;
    ControllerBroadcasts * broadcast;
    ofstream xstFile;

private:
    CthThread thread;
    static void threadRun(Controller*);

    double startCTime;
    double startWTime;
    double startBenchTime;

};

#endif // CONTROLLER_H

