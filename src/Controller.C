/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "Controller.h"
#include "ReductionMgr.h"
#include "CollectionMaster.h"
#include "Output.h"
#include "strlib.h"
#include "BroadcastObject.h"
#include "NamdState.h"
#include "ScriptTcl.h"
#include "Broadcasts.h"
#include "LdbCoordinator.h"
#include "Thread.h"
#include <math.h>
#include <signal.h>
#include "NamdOneTools.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "Random.h"
#include "imd.h"
#include "IMDOutput.h"
#include "InfoStream.h"
#include "BackEnd.h"

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
extern "C" void CApplicationDepositNode0Data(char *);
#endif

#ifndef cbrt
  // cbrt() not in math.h on goneril
  #define cbrt(x)  pow(x,(double)(1.0/3.0))
#endif

//#define DEBUG_PRESSURE
#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#define XXXBIGREAL 1.0e32

Controller::Controller(NamdState *s) :
	computeChecksum(0), marginViolations(0), pairlistWarnings(0),
	simParams(Node::Object()->simParameters),
	state(s),
	collection(CollectionMaster::Object()),
        startCTime(0),
        firstCTime(CmiTimer()),
        startWTime(0),
        firstWTime(CmiWallTimer()),
        startBenchTime(0),
	ldbSteps(0)

{
    broadcast = new ControllerBroadcasts;
    reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_BASIC);
    if (simParams->pressureProfileOn) {
      pressureProfileReduction = 
        ReductionMgr::Object()->willRequire(REDUCTIONS_PPROFILE);
      pressureProfileAverage = new BigReal[3*simParams->pressureProfileSlabs];
      memset(pressureProfileAverage, 0, 
          3*simParams->pressureProfileSlabs*sizeof(BigReal));
    } else {
      pressureProfileReduction = NULL;
      pressureProfileAverage = NULL;
    }
    random = new Random(simParams->randomSeed);
    random->split(0,PatchMap::Object()->numPatches()+1);

    rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
    berendsenPressure_avg = 0; berendsenPressure_count = 0;
    // strainRate tensor is symmetric to avoid rotation
    langevinPiston_strainRate =
	Tensor::symmetric(simParams->strainRate,simParams->strainRate2);
    if ( ! simParams->useFlexibleCell ) {
      BigReal avg = trace(langevinPiston_strainRate) / 3.;
      langevinPiston_strainRate = Tensor::identity(avg);
    } else if ( simParams->useConstantRatio ) {
#define AVGXY(T) T.xy = T.yx = 0; T.xx = T.yy = 0.5 * ( T.xx + T.yy );\
		 T.xz = T.zx = T.yz = T.zy = 0.5 * ( T.xz + T.yz );
      AVGXY(langevinPiston_strainRate);
#undef AVGXY
    }
    smooth2_avg = XXXBIGREAL;
    temp_avg = 0;
    pressure_avg = 0;
    groupPressure_avg = 0;
    avg_count = 0;
    pressure_tavg = 0;
    groupPressure_tavg = 0;
    tavg_count = 0;
    checkpoint_stored = 0;
}

Controller::~Controller(void)
{
    delete broadcast;
    delete reduction;
    delete pressureProfileReduction;
    delete [] pressureProfileAverage;
    delete random;
}

void Controller::threadRun(Controller* arg)
{
    arg->algorithm();
}

void Controller::run(void)
{
    // create a Thread and invoke it
    DebugM(4, "Starting thread in controller on this=" << this << "\n");
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),CTRL_STK_SZ);
    CthSetStrategyDefault(thread);
    awaken();
}


void Controller::algorithm(void)
{
  int scriptTask;
  int scriptSeq = 0;
  BackEnd::awaken();
  while ( (scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) != SCRIPT_END) {
    switch ( scriptTask ) {
      case SCRIPT_OUTPUT:
        enqueueCollections(FILE_OUTPUT);
        outputExtendedSystem(FILE_OUTPUT);
        break;
      case SCRIPT_MEASURE:
        enqueueCollections(EVAL_MEASURE);
        break;
      case SCRIPT_REINITVELS:
        iout << "REINITIALIZING VELOCITIES AT STEP " << simParams->firstTimestep
          << " TO " << simParams->initialTemp << " KELVIN.\n" << endi;
        break;
      case SCRIPT_RESCALEVELS:
        iout << "RESCALING VELOCITIES AT STEP " << simParams->firstTimestep
          << " BY " << simParams->scriptArg1 << "\n" << endi;
        break;
      case SCRIPT_CHECKPOINT:
        iout << "CHECKPOINTING POSITIONS AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        checkpoint_stored = 1;
        checkpoint_lattice = state->lattice;
        checkpoint_langevinPiston_strainRate = langevinPiston_strainRate;
        checkpoint_berendsenPressure_avg = berendsenPressure_avg;
        checkpoint_berendsenPressure_count = berendsenPressure_count;
        break;
      case SCRIPT_REVERT:
        iout << "REVERTING POSITIONS AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        if ( ! checkpoint_stored )
          NAMD_die("Unable to revert, checkpoint was never called!");
        state->lattice = checkpoint_lattice;
        langevinPiston_strainRate = checkpoint_langevinPiston_strainRate;
        berendsenPressure_avg = checkpoint_berendsenPressure_avg;
        berendsenPressure_count = checkpoint_berendsenPressure_count;
        break;
      case SCRIPT_MINIMIZE:
        minimize();
        break;
      case SCRIPT_RUN:
        integrate();
        break;
    }
    BackEnd::awaken();
  }
  enqueueCollections(END_OF_RUN);
  outputExtendedSystem(END_OF_RUN);
  terminate();
}


extern int eventEndOfTimeStep;

// Handle SIGINT so that restart files get written completely.
static int gotsigint = 0;
static void my_sigint_handler(int sig) {
  if (sig == SIGINT) gotsigint = 1;
}
extern "C" {
  typedef void (*namd_sighandler_t)(int);
}

void Controller::integrate() {

    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;

    nbondFreq = simParams->nonbondedFrequency;
    const int dofull = ( simParams->fullDirectOn ||
			simParams->FMAOn || simParams->PMEOn );
    if (dofull)
      slowFreq = simParams->fullElectFrequency;
    else
      slowFreq = simParams->nonbondedFrequency;

    reassignVelocities(step);  // only for full-step velecities
    receivePressure(step);
    printFepMessage(step);
    printDynamicsEnergies(step);
    outputFepEnergy(step);
    traceUserEvent(eventEndOfTimeStep);
    outputExtendedSystem(step);
    rebalanceLoad(step);

    // Handling SIGINT doesn't seem to be working on Lemieux, and it
    // sometimes causes the net-xxx versions of NAMD to segfault on exit, 
    // so disable it for now.
    // namd_sighandler_t oldhandler = signal(SIGINT, 
    //  (namd_sighandler_t)my_sigint_handler);
    for ( ++step ; step <= numberOfSteps; ++step )
    {

        rescaleVelocities(step);
	tcoupleVelocities(step);
	berendsenPressure(step);
	langevinPiston1(step);
	enqueueCollections(step);  // after lattice scaling!
	receivePressure(step);
	langevinPiston2(step);
        reassignVelocities(step);
        printDynamicsEnergies(step);
        outputFepEnergy(step);
        traceUserEvent(eventEndOfTimeStep);
  // if (gotsigint) {
  //   iout << iINFO << "Received SIGINT; shutting down.\n" << endi;
  //   NAMD_quit();
  // }
        outputExtendedSystem(step);
#if CYCLE_BARRIER
        cycleBarrier(!((step+1) % stepsPerCycle),step);
#elif  PME_BARRIER
        cycleBarrier(dofull && !(step%slowFreq),step);
#endif

        rebalanceLoad(step);

#if  PME_BARRIER
        cycleBarrier(dofull && !((step+1)%slowFreq),step);   // step before PME
#endif
    }
    // signal(SIGINT, oldhandler);
}


#define CALCULATE \
  printMinimizeEnergies(step); \
  rebalanceLoad(step); \
  if ( step == numberOfSteps ) return; \
  else ++step;

#define MOVETO(X) \
  if ( step == numberOfSteps ) { \
    if ( 0 ) { iout << "LINE MINIMIZER: RETURNING TO " << mid.x << " FROM " << last.x << "\n" << endi; } \
    if ( newDir || (mid.x-last.x) ) { \
      broadcast->minimizeCoefficient.publish(minSeq++,mid.x-last.x); \
    } else { \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
    } \
    enqueueCollections(step); \
    CALCULATE \
  } else if ( (X)-last.x ) { \
    if ( 0 ) { iout << "LINE MINIMIZER: MOVING FROM " << last.x << " TO " << (X) << "\n" << endi; } \
    broadcast->minimizeCoefficient.publish(minSeq++,(X)-last.x); \
    newDir = 0; \
    last.x = (X); \
    enqueueCollections(step); \
    CALCULATE \
    last.u = min_energy; \
    last.dudx = -1. * min_f_dot_v; \
  }

struct minpoint {
  BigReal x, u, dudx;
  minpoint() : x(0), u(0), dudx(0) { ; }
};

void Controller::minimize() {
  // iout << "Controller::minimize() called.\n" << endi;

  const int numberOfSteps = simParams->N;
  int step = simParams->firstTimestep;
  slowFreq = nbondFreq = 1;
  BigReal tinystep = simParams->minTinyStep;  // 1.0e-6
  BigReal babystep = simParams->minBabyStep;  // 1.0e-2
  BigReal linegoal = simParams->minLineGoal;  // 1.0e-4
  BigReal initstep = tinystep;
  const BigReal goldenRatio = 0.5 * ( sqrt(5.0) - 1.0 );

  CALCULATE

  int minSeq = 0;
  int atStart = 2;
  BigReal old_f_dot_f = min_f_dot_f;
  broadcast->minimizeCoefficient.publish(minSeq++,0.);
  broadcast->minimizeCoefficient.publish(minSeq++,0.); // v = f
  int newDir = 1;
  min_f_dot_v = min_f_dot_f;
  min_v_dot_v = min_f_dot_f;
  while ( 1 ) {
    // line minimization
    // bracket minimum on line
    minpoint lo,hi,mid,last;
    BigReal x = 0;
    lo.x = x;
    lo.u = min_energy;
    lo.dudx = -1. * min_f_dot_v;
    int noGradients = min_huge_count;
    mid = lo;
    last = mid;
    BigReal tol = fabs( linegoal * min_f_dot_v );
    if ( initstep > babystep ) initstep = babystep;
    if ( initstep < 1.0e-300 ) initstep = 1.0e-300;
    iout << "INITIAL STEP: " << initstep << "\n" << endi;
    iout << "GRADIENT TOLERANCE: " << tol << "\n" << endi;
    x = initstep;
    x *= sqrt( min_f_dot_f / min_v_dot_v ); MOVETO(x)
    // bracket minimum on line
    initstep *= 0.25;
    BigReal maxinitstep = initstep * 16.0;
    while ( last.u < mid.u ) {
      initstep *= 2.0;
      lo = mid; mid = last;
      // when bracketed, need to know if midpoint gradient is valid
      noGradients = min_huge_count;
      x *= 2.0; MOVETO(x)
    }
    if ( initstep > maxinitstep ) initstep = maxinitstep;
    hi = last;
    iout << "BRACKET: " << (hi.x-lo.x) << " " << ((hi.u>lo.u?hi.u:lo.u)-mid.u) << " " << lo.dudx << " " << mid.dudx << " " << hi.dudx << " \n" << endi;
    // converge on minimum on line
    int itcnt;
    for ( itcnt = 10; fabs(last.dudx) > tol && itcnt > 0 ; --itcnt ) {
      // select new position
      if ( noGradients ) {
       if ( ( mid.x - lo.x ) > ( hi.x - mid.x ) ) {  // subdivide left side
	x = (1.0 - goldenRatio) * lo.x + goldenRatio * mid.x;
	MOVETO(x)
	if ( last.u <= mid.u ) {
	  hi = mid; mid = last; noGradients = min_huge_count;
	} else {
	  lo = last;
	}
       } else {  // subdivide right side
	x = (1.0 - goldenRatio) * hi.x + goldenRatio * mid.x;
	MOVETO(x)
	if ( last.u <= mid.u ) {
	  lo = mid; mid = last; noGradients = min_huge_count;
	} else {
	  hi = last;
	}
       }
      } else {
       if ( mid.dudx > 0. ) {  // subdivide left side
        BigReal altxhi = 0.1 * lo.x + 0.9 * mid.x;
        BigReal altxlo = 0.9 * lo.x + 0.1 * mid.x;
        x = mid.dudx*(mid.x*mid.x-lo.x*lo.x) + 2*mid.x*(lo.u-mid.u);
        x /= 2*(mid.dudx*(mid.x-lo.x)+(lo.u-mid.u));
        if ( x > altxhi ) x = altxhi;
        if ( x < altxlo ) x = altxlo;
        if ( x-last.x == 0 ) break;
        MOVETO(x)
        if ( last.u <= mid.u ) {
	  hi = mid; mid = last; noGradients = min_huge_count;
	} else {
	  lo = last;
	}
       } else {  // subdivide right side
        BigReal altxlo = 0.1 * hi.x + 0.9 * mid.x;
        BigReal altxhi = 0.9 * hi.x + 0.1 * mid.x;
        x = mid.dudx*(mid.x*mid.x-hi.x*hi.x) + 2*mid.x*(hi.u-mid.u);
        x /= 2*(mid.dudx*(mid.x-hi.x)+(hi.u-mid.u));
        if ( x < altxlo ) x = altxlo;
        if ( x > altxhi ) x = altxhi;
        if ( x-last.x == 0 ) break;
        MOVETO(x)
        if ( last.u <= mid.u ) {
	  lo = mid; mid = last; noGradients = min_huge_count;
	} else {
	  hi = last;
	}
       }
      }
      iout << "BRACKET: " << (hi.x-lo.x) << " " << ((hi.u>lo.u?hi.u:lo.u)-mid.u) << " " << lo.dudx << " " << mid.dudx << " " << hi.dudx << " \n" << endi;
    }
    // new direction
    broadcast->minimizeCoefficient.publish(minSeq++,0.);
    BigReal c = min_f_dot_f / old_f_dot_f;
    c = ( c > 1.5 ? 1.5 : c );
    if ( atStart ) { c = 0; --atStart; }
    if ( c*c*min_v_dot_v > 100*min_f_dot_f ) { c = 0; }
    if ( c == 0 ) {
      iout << "RESTARTING CONJUGATE GRADIENT ALGORITHM\n" << endi;
    } else {
      iout << "NEW SEARCH DIRECTION\n" << endi;
    }
    broadcast->minimizeCoefficient.publish(minSeq++,c); // v = c*v+f
    newDir = 1;
    old_f_dot_f = min_f_dot_f;
    min_f_dot_v = c * min_f_dot_v + min_f_dot_f;
    min_v_dot_v = c*c*min_v_dot_v + 2*c*min_f_dot_v + min_f_dot_f;
  }

}

#undef MOVETO
#undef CALCULATE

void Controller::berendsenPressure(int step)
{
  if ( simParams->berendsenPressureOn ) {
   berendsenPressure_count += 1;
   berendsenPressure_avg += controlPressure;
   const int freq = simParams->berendsenPressureFreq;
   if ( ! (berendsenPressure_count % freq) ) {
    Tensor factor = berendsenPressure_avg / berendsenPressure_count;
    berendsenPressure_avg = 0;
    berendsenPressure_count = 0;
    // We only use on-diagonal terms (for now)
    factor = Tensor::diagonal(diagonal(factor));
    factor -= Tensor::identity(simParams->berendsenPressureTarget);
    factor *= ( ( simParams->berendsenPressureCompressibility / 3.0 ) *
       simParams->dt * freq / simParams->berendsenPressureRelaxationTime );
    factor += Tensor::identity(1.0);
#define LIMIT_SCALING(VAR,MIN,MAX,FLAG) {\
         if ( VAR < (MIN) ) { VAR = (MIN); FLAG = 1; } \
         if ( VAR > (MAX) ) { VAR = (MAX); FLAG = 1; } }
    int limited = 0;
    LIMIT_SCALING(factor.xx,1./1.03,1.03,limited)
    LIMIT_SCALING(factor.yy,1./1.03,1.03,limited)
    LIMIT_SCALING(factor.zz,1./1.03,1.03,limited)
#undef LIMIT_SCALING
    if ( limited ) {
      iout << iERROR << "Step " << step <<
	" cell rescaling factor limited.\n" << endi;
    }
    broadcast->positionRescaleFactor.publish(step,factor);
    state->lattice.rescale(factor);
   }
  } else {
    berendsenPressure_avg = 0;
    berendsenPressure_count = 0;
  }
}

void Controller::langevinPiston1(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMAN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

#ifdef DEBUG_PRESSURE
    iout << iINFO << "entering langevinPiston1, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( (step-1) % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell ) {
        // We only use on-diagonal terms (for now)
        if ( simParams->useConstantRatio ) {
	  BigReal r = f2 * random->gaussian();
	  strainRate.xx += r;
	  strainRate.yy += r;
	  strainRate.zz += f2 * random->gaussian();
        } else {
	  strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
        }
      } else {
	strainRate += f2 * Tensor::identity(random->gaussian());
      }
#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
    
    Tensor ptarget;
    ptarget.zz = simParams->langevinPistonTarget;
    ptarget.xx = ptarget.yy = simParams->langevinPistonTarget - 
        simParams->surfaceTensionTarget / state->lattice.c().z;

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );

#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( (step-1-slowFreq/2) % slowFreq ) )
    {
      // We only use on-diagonal terms (for now)
      Tensor factor;
      if ( !simParams->useConstantArea ) {
        factor.xx = exp( dt_long * strainRate.xx );
        factor.yy = exp( dt_long * strainRate.yy );
      } else {
        factor.xx = factor.yy = 1;
      }
      factor.zz = exp( dt_long * strainRate.zz );
      broadcast->positionRescaleFactor.publish(step,factor);
      state->lattice.rescale(factor);
#ifdef DEBUG_PRESSURE
      iout << iINFO << "rescaling by: " << factor << "\n";
#endif
    }

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }

  }
}

void Controller::langevinPiston2(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMAN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
   
    Tensor ptarget;
    ptarget.zz = simParams->langevinPistonTarget;
    ptarget.xx = ptarget.yy = simParams->langevinPistonTarget -
        simParams->surfaceTensionTarget / state->lattice.c().z;

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );
 
#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( step % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell ) {
        // We only use on-diagonal terms (for now)
        if ( simParams->useConstantRatio ) {
	  BigReal r = f2 * random->gaussian();
	  strainRate.xx += r;
	  strainRate.yy += r;
	  strainRate.zz += f2 * random->gaussian();
        } else {
	  strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
        }
      } else {
	strainRate += f2 * Tensor::identity(random->gaussian());
      }
#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

#ifdef DEBUG_PRESSURE
    iout << iINFO << "exiting langevinPiston2, strain rate: " << strainRate << "\n";
#endif
  }
}

void Controller::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 ) {
    rescaleVelocities_sumTemps += temperature;  ++rescaleVelocities_numTemps;
    if ( rescaleVelocities_numTemps == rescaleFreq ) {
      BigReal avgTemp = rescaleVelocities_sumTemps / rescaleVelocities_numTemps;
      const BigReal rescaleTemp = simParams->rescaleTemp;
      BigReal factor = sqrt(rescaleTemp/avgTemp);
      broadcast->velocityRescaleFactor.publish(step,factor);
      iout << "RESCALING VELOCITIES AT STEP " << step
           << " FROM AVERAGE TEMPERATURE OF " << avgTemp
           << " TO " << rescaleTemp << " KELVIN.\n" << endi;
      rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
    }
  }
}

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
void Controller::printFepMessage(int step)
{
  if (simParams->fepOn) {
    const BigReal lambda = simParams->lambda;
    const BigReal lambda2 = simParams->lambda2;
    const int fepEquilSteps = simParams->fepEquilSteps;
    iout << "FEP: RESETTING FOR NEW FEP WINDOW "
         << "LAMBDA SET TO " << lambda << " LAMBDA2 " << lambda2
         << "\nFEP: WINDOW TO HAVE " << fepEquilSteps
         << " STEPS OF EQUILIBRATION PRIOR TO FEP DATA COLLECTION.\n" << endi;
  }
} 
//fepe

void Controller::reassignVelocities(int step)
{
  const int reassignFreq = simParams->reassignFreq;
  if ( ( reassignFreq > 0 ) && ! ( step % reassignFreq ) ) {
    BigReal newTemp = simParams->reassignTemp;
    newTemp += ( step / reassignFreq ) * simParams->reassignIncr;
    if ( simParams->reassignIncr > 0.0 ) {
      if ( newTemp > simParams->reassignHold && simParams->reassignHold > 0.0 )
        newTemp = simParams->reassignHold;
    } else {
      if ( newTemp < simParams->reassignHold )
        newTemp = simParams->reassignHold;
    }
    iout << "REASSIGNING VELOCITIES AT STEP " << step
         << " TO " << newTemp << " KELVIN.\n" << endi;
  }
}

void Controller::tcoupleVelocities(int step)
{
  if ( simParams->tCoupleOn )
  {
    const BigReal tCoupleTemp = simParams->tCoupleTemp;
    BigReal coefficient = 1.;
    if ( temperature > 0. ) coefficient = tCoupleTemp/temperature - 1.;
    broadcast->tcoupleCoefficient.publish(step,coefficient);
  }
}

static char *FORMAT(BigReal X)
{
  static char tmp_string[25];
  const double maxnum = 99999999.9999;
  if ( X > maxnum ) X = maxnum;
  if ( X < -maxnum ) X = -maxnum;
  sprintf(tmp_string," %14.4f",X); 
  return tmp_string;
}

static char *FORMAT(const char *X)
{
  static char tmp_string[25];
  sprintf(tmp_string," %14s",X); 
  return tmp_string;
}

static char *ETITLE(int X)
{
  static char tmp_string[21];
  sprintf(tmp_string,"ENERGY: %7d",X); 
  return  tmp_string;
}

void Controller::receivePressure(int step, int minimize)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    reduction->require();
    if (pressureProfileReduction) 
      pressureProfileReduction->require();


    Tensor virial;
    Tensor virial_normal;
    Tensor virial_nbond;
    Tensor virial_slow;
#ifdef ALTVIRIAL
    Tensor altVirial_normal;
    Tensor altVirial_nbond;
    Tensor altVirial_slow;
#endif
    Tensor intVirial;
    Tensor intVirial_normal;
    Tensor intVirial_nbond;
    Tensor intVirial_slow;
    Vector extForce_normal;
    Vector extForce_nbond;
    Vector extForce_slow;
    BigReal volume;

    int numAtoms = molecule->numAtoms;
    numDegFreedom = 3 * numAtoms;
    int numGroupDegFreedom = 3 * molecule->numHydrogenGroups;
    int numFixedAtoms =
	( simParameters->fixedAtomsOn ? molecule->numFixedAtoms : 0 );
    int numFixedGroups = ( numFixedAtoms ? molecule->numFixedGroups : 0 );
    if ( numFixedAtoms ) numDegFreedom -= 3 * numFixedAtoms;
    if ( numFixedGroups ) numGroupDegFreedom -= 3 * numFixedGroups;
    if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
      numDegFreedom -= 3;
      numGroupDegFreedom -= 3;
    }
    int numRigidBonds = molecule->numRigidBonds;
    int numFixedRigidBonds =
	( simParameters->fixedAtomsOn ? molecule->numFixedRigidBonds : 0 );
    numDegFreedom -= ( numRigidBonds - numFixedRigidBonds );

    kineticEnergyHalfstep = reduction->item(REDUCTION_HALFSTEP_KINETIC_ENERGY);
    kineticEnergyCentered = reduction->item(REDUCTION_CENTERED_KINETIC_ENERGY);

    BigReal groupKineticEnergyHalfstep = kineticEnergyHalfstep -
	reduction->item(REDUCTION_INT_HALFSTEP_KINETIC_ENERGY);
    BigReal groupKineticEnergyCentered = kineticEnergyCentered -
	reduction->item(REDUCTION_INT_CENTERED_KINETIC_ENERGY);

    BigReal atomTempHalfstep = 2.0 * kineticEnergyHalfstep
					/ ( numDegFreedom * BOLTZMAN );
    BigReal atomTempCentered = 2.0 * kineticEnergyCentered
					/ ( numDegFreedom * BOLTZMAN );
    BigReal groupTempHalfstep = 2.0 * groupKineticEnergyHalfstep
					/ ( numGroupDegFreedom * BOLTZMAN );
    BigReal groupTempCentered = 2.0 * groupKineticEnergyCentered
					/ ( numGroupDegFreedom * BOLTZMAN );

    /*  test code for comparing different temperatures
    iout << "TEMPTEST: " << step << " " << 
	atomTempHalfstep << " " <<
	atomTempCentered << " " <<
	groupTempHalfstep << " " <<
	groupTempCentered << "\n" << endi;
    */

    GET_TENSOR(virial_normal,reduction,REDUCTION_VIRIAL_NORMAL);
    GET_TENSOR(virial_nbond,reduction,REDUCTION_VIRIAL_NBOND);
    GET_TENSOR(virial_slow,reduction,REDUCTION_VIRIAL_SLOW);

#ifdef ALTVIRIAL
    GET_TENSOR(altVirial_normal,reduction,REDUCTION_ALT_VIRIAL_NORMAL);
    GET_TENSOR(altVirial_nbond,reduction,REDUCTION_ALT_VIRIAL_NBOND);
    GET_TENSOR(altVirial_slow,reduction,REDUCTION_ALT_VIRIAL_SLOW);
#endif

    GET_TENSOR(intVirial_normal,reduction,REDUCTION_INT_VIRIAL_NORMAL);
    GET_TENSOR(intVirial_nbond,reduction,REDUCTION_INT_VIRIAL_NBOND);
    GET_TENSOR(intVirial_slow,reduction,REDUCTION_INT_VIRIAL_SLOW);

    GET_VECTOR(extForce_normal,reduction,REDUCTION_EXT_FORCE_NORMAL);
    GET_VECTOR(extForce_nbond,reduction,REDUCTION_EXT_FORCE_NBOND);
    GET_VECTOR(extForce_slow,reduction,REDUCTION_EXT_FORCE_SLOW);
    Vector extPosition = lattice.origin();
    virial_normal -= outer(extForce_normal,extPosition);
    virial_nbond -= outer(extForce_nbond,extPosition);
    virial_slow -= outer(extForce_slow,extPosition);

    kineticEnergy = kineticEnergyCentered;
    temperature = 2.0 * kineticEnergyCentered / ( numDegFreedom * BOLTZMAN );

    if ( (volume=lattice.volume()) != 0. )
    {
      // kinetic energy component included in virials
      pressure_normal = virial_normal / volume;
      groupPressure_normal = ( virial_normal - intVirial_normal ) / volume;

      if ( minimize || ! ( step % nbondFreq ) )
      {
        pressure_nbond = virial_nbond / volume;
        groupPressure_nbond = ( virial_nbond - intVirial_nbond ) / volume;
      }

      if ( minimize || ! ( step % slowFreq ) )
      {
        pressure_slow = virial_slow / volume;
        groupPressure_slow = ( virial_slow - intVirial_slow ) / volume;
      }

/*
      iout << "VIRIALS: " << virial_normal << " " << virial_nbond << " " <<
	virial_slow << " " << ( virial_normal - intVirial_normal ) << " " <<
	( virial_nbond - intVirial_nbond ) << " " <<
	( virial_slow - intVirial_slow ) << "\n";
*/

      pressure = pressure_normal + pressure_nbond + pressure_slow; 
      groupPressure = groupPressure_normal + groupPressure_nbond +
						groupPressure_slow;
    }
    else
    {
      pressure = Tensor();
      groupPressure = Tensor();
    }

    if ( simParameters->useGroupPressure )
    {
      controlPressure_normal = groupPressure_normal;
      controlPressure_nbond = groupPressure_nbond;
      controlPressure_slow = groupPressure_slow;
      controlPressure = groupPressure;
      controlNumDegFreedom = molecule->numHydrogenGroups - numFixedGroups;
      if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
        controlNumDegFreedom -= 1;
      }
    }
    else
    {
      controlPressure_normal = pressure_normal;
      controlPressure_nbond = pressure_nbond;
      controlPressure_slow = pressure_slow;
      controlPressure = pressure;
      controlNumDegFreedom = numDegFreedom / 3;
    }

    if ( simParameters->useFlexibleCell ) {
      // use symmetric pressure to control rotation
      // controlPressure_normal = symmetric(controlPressure_normal);
      // controlPressure_nbond = symmetric(controlPressure_nbond);
      // controlPressure_slow = symmetric(controlPressure_slow);
      // controlPressure = symmetric(controlPressure);
      // only use on-diagonal components for now
      controlPressure_normal = Tensor::diagonal(diagonal(controlPressure_normal));
      controlPressure_nbond = Tensor::diagonal(diagonal(controlPressure_nbond));
      controlPressure_slow = Tensor::diagonal(diagonal(controlPressure_slow));
      controlPressure = Tensor::diagonal(diagonal(controlPressure));
      if ( simParameters->useConstantRatio ) {
#define AVGXY(T) T.xy = T.yx = 0; T.xx = T.yy = 0.5 * ( T.xx + T.yy );\
		 T.xz = T.zx = T.yz = T.zy = 0.5 * ( T.xz + T.yz );
        AVGXY(controlPressure_normal);
        AVGXY(controlPressure_nbond);
        AVGXY(controlPressure_slow);
        AVGXY(controlPressure);
#undef AVGXY
      }
    } else {
      controlPressure_normal =
		Tensor::identity(trace(controlPressure_normal)/3.);
      controlPressure_nbond =
		Tensor::identity(trace(controlPressure_nbond)/3.);
      controlPressure_slow =
		Tensor::identity(trace(controlPressure_slow)/3.);
      controlPressure =
		Tensor::identity(trace(controlPressure)/3.);
    }

#ifdef DEBUG_PRESSURE
    iout << iINFO << "Control pressure = " << controlPressure <<
      " = " << controlPressure_normal << " + " <<
      controlPressure_nbond << " + " << controlPressure_slow << "\n" << endi;
#endif

  if (pressureProfileReduction) {
    int i;

    // accumulate the pressure profile computed for this step into the average.
    for (i=0; i<3*simParameters->pressureProfileSlabs; i++)
      pressureProfileAverage[i] += pressureProfileReduction->item(i);
    if (!(step % simParameters->pressureProfileFreq)) {
      // convert NAMD internal virial to pressure in units of bar by 
      // multiplying by PRESSUREFACTOR and dividing by the volume of one slab.
      BigReal scalefac = PRESSUREFACTOR * 
        simParameters->pressureProfileSlabs / lattice.volume();

      // output pressure profile for this step
      iout << "PRESSUREPROFILE: " << step << " ";
      for (i=0; i<3*simParameters->pressureProfileSlabs; i++) 
        iout << pressureProfileReduction->item(i)*scalefac << " ";
      iout << "\n" << endi; 

      if (step != simParameters->firstTimestep) {
        scalefac /= simParameters->pressureProfileFreq;
      }
      // output pressure profile averaged over the last Freq steps.
      iout << "PRESSUREPROFILEAVG: " << step << " ";
      for (i=0; i<3*simParameters->pressureProfileSlabs; i++) 
        iout << pressureProfileAverage[i]*scalefac << " ";
      iout << "\n" << endi; 
      // Clear the average for the next block
      memset(pressureProfileAverage, 0, 
          3*simParameters->pressureProfileSlabs*sizeof(BigReal));
    }
  }
}

void Controller::compareChecksums(int step, int forgiving) {
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;

    // Some basic correctness checking
    BigReal checksum;

    checksum = reduction->item(REDUCTION_ATOM_CHECKSUM);
    if ( ((int)checksum) != molecule->numAtoms )
      NAMD_bug("Bad global atom count!\n");

    checksum = reduction->item(REDUCTION_COMPUTE_CHECKSUM);
    if ( ((int)checksum) != computeChecksum ) {
      if ( computeChecksum )
        NAMD_bug("Bad global compute count!\n");
      else
        computeChecksum = ((int)checksum);
    }

    checksum = reduction->item(REDUCTION_BOND_CHECKSUM);
    if ( checksum && (((int)checksum) != molecule->numCalcBonds) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcBonds) )
        iout << iWARN << "Bad global bond count!\n" << endi;
      else NAMD_bug("Bad global bond count!\n");
    }

    checksum = reduction->item(REDUCTION_ANGLE_CHECKSUM);
    if ( checksum && (((int)checksum) != molecule->numCalcAngles) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcAngles) )
        iout << iWARN << "Bad global angle count!\n" << endi;
      else NAMD_bug("Bad global angle count!\n");
    }

    checksum = reduction->item(REDUCTION_DIHEDRAL_CHECKSUM);
    if ( checksum && (((int)checksum) != molecule->numCalcDihedrals) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcDihedrals) )
        iout << iWARN << "Bad global dihedral count!\n" << endi;
      else NAMD_bug("Bad global dihedral count!\n");
    }

    checksum = reduction->item(REDUCTION_IMPROPER_CHECKSUM);
    if ( checksum && (((int)checksum) != molecule->numCalcImpropers) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcImpropers) )
        iout << iWARN << "Bad global improper count!\n" << endi;
      else NAMD_bug("Bad global improper count!\n");
    }

    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) > molecule->numCalcExclusions &&
         ( ! simParams->mollyOn || step % slowFreq ) )
      iout << iWARN << "Not all atoms have unique coordinates.\n" << endi;

    checksum = reduction->item(REDUCTION_MARGIN_VIOLATIONS);
    if ( ((int)checksum) && ! marginViolations ) {
      iout << iERROR << "Margin is too small for " << ((int)checksum) <<
        " atoms during timestep " << step << ".\n" << iERROR <<
      "Incorrect nonbonded forces and energies may be calculated!\n" << endi;
    }
    marginViolations += (int)checksum;

    checksum = reduction->item(REDUCTION_PAIRLIST_WARNINGS);
    if ( simParams->outputPairlists && ((int)checksum) && ! pairlistWarnings ) {
      iout << iINFO <<
        "Pairlistdist is too small for " << ((int)checksum) <<
        " computes during timestep " << step << ".\n" << endi;
    }
    if ( simParams->outputPairlists )  pairlistWarnings += (int)checksum;

    checksum = reduction->item(REDUCTION_STRAY_CHARGE_ERRORS);
    if ( checksum ) {
      if ( forgiving )
        iout << iWARN << "Stray PME grid charges ignored!\n" << endi;
      else NAMD_bug("Stray PME grid charges detected!\n");
    }
}

void Controller::printTiming(int step) {

    if ( simParams->outputTiming && ! ( step % simParams->outputTiming ) )
    {
      const double endWTime = CmiWallTimer() - firstWTime;
      const double endCTime = CmiTimer() - firstCTime;

      const double elapsedW = 
	(endWTime - startWTime) / simParams->outputTiming;
      const double elapsedC = 
	(endCTime - startCTime) / simParams->outputTiming;

      const double remainingW = elapsedW * (simParams->N - step);
      const double remainingW_hours = remainingW / 3600;

      startWTime = endWTime;
      startCTime = endCTime;

      if ( step >= (simParams->firstTimestep + simParams->outputTiming) ) {
	CmiPrintf("TIMING: %d  CPU: %g, %g/step  Wall: %g, %g/step"
		  ", %g hours remaining, %d kB of memory in use.\n",
		  step, endCTime, elapsedC, endWTime, elapsedW,
		  remainingW_hours, (memusage()/1024));
#if 0
        iout << "TIMING: " << step
             << "  CPU: " << endCTime << ", " << elapsedC << "/step"
             << "  Wall: " << endWTime << ", " << elapsedW << "/step"
             << ", " << remainingW_hours << " hours remaining"
             << ", " << (memusage()/1024) << " kB of memory in use"
             << ".\n" << endi;
#endif
      }
    }
}

void Controller::printMinimizeEnergies(int step) {

    receivePressure(step,1);

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    BigReal checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) && ((int)checksum) < molecule->numCalcExclusions )
      iout << iWARN << "Bad global exclusion count, possible error!\n" << iWARN
        << "Increasing cutoff during minimization may avoid this.\n" << endi;
    compareChecksums(step,1);

    printEnergies(step,1);

/*
    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal totalEnergy;

    bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);
    electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
    ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);
    electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);

    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	electEnergy + electEnergySlow + ljEnergy + boundaryEnergy + miscEnergy;
*/

    min_energy = totalEnergy;
    min_f_dot_f = reduction->item(REDUCTION_MIN_F_DOT_F);
    min_f_dot_v = reduction->item(REDUCTION_MIN_F_DOT_V);
    min_v_dot_v = reduction->item(REDUCTION_MIN_V_DOT_V);
    min_huge_count = reduction->item(REDUCTION_MIN_HUGE_COUNT);

/*
    if ( ( step % 10 ) == 0 ) {
	iout << "ETITLE:     TS    BOND        ANGLE       "
	     << "DIHED       IMPRP       ELECT       VDW       "
	     << "BOUNDARY    MISC        TOTAL       GRADIENT\n" << endi;
    }

    // N.B.  HP's aCC compiler merges FORMAT calls in the same expression.
    //       Need separate statements because data returned in static array.

    iout << ETITLE(step);
    iout << FORMAT(bondEnergy);
    iout << FORMAT(angleEnergy);
    iout << FORMAT(dihedralEnergy);
    iout << FORMAT(improperEnergy);
    iout << FORMAT(electEnergy+electEnergySlow);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(miscEnergy);
    iout << FORMAT(totalEnergy);
    iout << FORMAT(sqrt(min_f_dot_f));
    iout << "\n" << endi;
*/
}

void Controller::printDynamicsEnergies(int step) {

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    BigReal checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) && ((int)checksum) < molecule->numCalcExclusions )
      NAMD_bug("Bad global exclusion count!\n");
    compareChecksums(step);

    printEnergies(step,0);
}

void Controller::printEnergies(int step, int minimize)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal potentialEnergy;
    BigReal flatEnergy;
    BigReal smoothEnergy;
    Vector momentum;
    Vector angularMomentum;
    BigReal volume = lattice.volume();

    bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);

    if ( minimize || ! ( step % nbondFreq ) )
    {
      electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
      ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);
//fepb
      electEnergy_f = reduction->item(REDUCTION_ELECT_ENERGY_F);
      ljEnergy_f = reduction->item(REDUCTION_LJ_ENERGY_F);
//fepe
    }

    if ( minimize || ! ( step % slowFreq ) )
    {
      electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
//fepb
      electEnergySlow_f = reduction->item(REDUCTION_ELECT_ENERGY_SLOW_F);
//fepe
    }

    momentum.x = reduction->item(REDUCTION_MOMENTUM_X);
    momentum.y = reduction->item(REDUCTION_MOMENTUM_Y);
    momentum.z = reduction->item(REDUCTION_MOMENTUM_Z);
    angularMomentum.x = reduction->item(REDUCTION_ANGULAR_MOMENTUM_X);
    angularMomentum.y = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Y);
    angularMomentum.z = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Z);

    potentialEnergy = bondEnergy + angleEnergy + dihedralEnergy +
	improperEnergy + electEnergy + electEnergySlow + ljEnergy +
	boundaryEnergy + miscEnergy;
    totalEnergy = potentialEnergy + kineticEnergy;
    flatEnergy = totalEnergy +
        (1.0/3.0)*( kineticEnergyHalfstep - kineticEnergyCentered);
    if ( !(step%slowFreq) ) {
      // only adjust based on most accurate energies
      BigReal s = (4.0/3.0)*( kineticEnergyHalfstep - kineticEnergyCentered);
      if ( smooth2_avg == XXXBIGREAL ) smooth2_avg = s;
      smooth2_avg *= 0.9375;
      smooth2_avg += 0.0625 * s;
    }
    smoothEnergy = flatEnergy + smooth2_avg -
        (4.0/3.0)*( kineticEnergyHalfstep - kineticEnergyCentered);

    if ( simParameters->outputMomenta && ! minimize &&
         ! ( step % simParameters->outputMomenta ) )
    {
      iout << "MOMENTUM: " << step 
           << " P: " << momentum
           << " L: " << angularMomentum
           << "\n" << endi;
    }

    if ( simParameters->outputPressure ) {
      pressure_tavg += pressure;
      groupPressure_tavg += groupPressure;
      tavg_count += 1;
      if ( minimize || ! ( step % simParameters->outputPressure ) ) {
        iout << "PRESSURE: " << step << " "
           << PRESSUREFACTOR * pressure << "\n"
           << "GPRESSURE: " << step << " "
           << PRESSUREFACTOR * groupPressure << "\n";
        if ( tavg_count > 1 ) iout << "PRESSAVG: " << step << " "
           << (PRESSUREFACTOR/tavg_count) * pressure_tavg << "\n"
           << "GPRESSAVG: " << step << " "
           << (PRESSUREFACTOR/tavg_count) * groupPressure_tavg << "\n" << endi;
        pressure_tavg = 0;
        groupPressure_tavg = 0;
        tavg_count = 0;
      }
    }

    if (simParameters->IMDon && !(step % simParameters->IMDfreq)) {
      IMDEnergies energies;
      energies.tstep = step;
      energies.T = temperature;
      energies.Etot = totalEnergy;
      energies.Epot = potentialEnergy;
      energies.Evdw = ljEnergy;
      energies.Eelec = electEnergy + electEnergySlow;
      energies.Ebond = bondEnergy;
      energies.Eangle = angleEnergy;
      energies.Edihe = dihedralEnergy;
      energies.Eimpr = improperEnergy;
      Node::Object()->imd->gather_energies(&energies);
    }
  
    int stepInRun = step - simParams->firstTimestep;
    if ( stepInRun % simParams->firstLdbStep == 0 ) {
     int benchPhase = stepInRun / simParams->firstLdbStep;
     if ( benchPhase > 0 && benchPhase < 7 ) {
      iout << iINFO;
      if ( benchPhase < 4 ) iout << "Initial time: ";
      else iout << "Benchmark time: ";
      iout << CkNumPes() << " CPUs ";
      {
        BigReal wallTime = CmiWallTimer() - startBenchTime;
        BigReal wallPerStep =
		(CmiWallTimer() - startBenchTime) / simParams->firstLdbStep;
	BigReal ns = simParams->dt / 1000000.0;
	BigReal days = 1.0 / (24.0 * 60.0 * 60.0);
	BigReal daysPerNano = wallPerStep * days / ns;
	iout << wallPerStep << " s/step " << daysPerNano << " days/ns ";
        iout << (memusage()/1024) << " kB memory\n" << endi;
      }
     }
     startBenchTime = CmiWallTimer();
    }

    printTiming(step);

    // callback to Tcl with whatever we can
#ifdef NAMD_TCL
#define CALLBACKDATA(LABEL,VALUE) \
		labels << (LABEL) << " "; values << (VALUE) << " ";
#define CALLBACKLIST(LABEL,VALUE) \
		labels << (LABEL) << " "; values << "{" << (VALUE) << "} ";
    if (node->getScript() && node->getScript()->doCallback()) {
      ostrstream labels, values;
      CALLBACKDATA("TS",step);
      CALLBACKDATA("BOND",bondEnergy);
      CALLBACKDATA("ANGLE",angleEnergy);
      CALLBACKDATA("DIHED",dihedralEnergy);
      CALLBACKDATA("IMPRP",improperEnergy);
      CALLBACKDATA("ELECT",electEnergy+electEnergySlow);
      CALLBACKDATA("VDW",ljEnergy);
      CALLBACKDATA("BOUNDARY",boundaryEnergy);
      CALLBACKDATA("MISC",miscEnergy);
      CALLBACKDATA("KINETIC",kineticEnergy);
      CALLBACKDATA("TOTAL",totalEnergy);
      CALLBACKDATA("TEMP",temperature);
      CALLBACKLIST("PRESSURE",pressure*PRESSUREFACTOR);
      CALLBACKLIST("GPRESSURE",groupPressure*PRESSUREFACTOR);
      CALLBACKDATA("VOLUME",lattice.volume());
      CALLBACKLIST("CELL_A",lattice.a());
      CALLBACKLIST("CELL_B",lattice.b());
      CALLBACKLIST("CELL_C",lattice.c());
      CALLBACKLIST("CELL_O",lattice.origin());
      labels << "PERIODIC"; values << "{" << lattice.a_p() << " "
		<< lattice.b_p() << " " << lattice.c_p() << "}";

      labels << '\0';  values << '\0';  // insane but makes Linux work
      char *labelstr = labels.str();
      char *valuestr = values.str();
      node->getScript()->doCallback(labelstr,valuestr);
      delete [] labelstr;
      delete [] valuestr;
    }
#undef CALLBACKDATA
#endif

    temp_avg += temperature;
    pressure_avg += trace(pressure)/3.;
    groupPressure_avg += trace(groupPressure)/3.;
    avg_count += 1;

    Vector pairVDWForce, pairElectForce;
    if ( simParameters->pairInteractionOn ) {
      GET_VECTOR(pairVDWForce,reduction,REDUCTION_PAIR_VDW_FORCE);
      GET_VECTOR(pairElectForce,reduction,REDUCTION_PAIR_ELECT_FORCE);
    }

    if ( simParams->outputPairlists && pairlistWarnings &&
				! (step % simParams->outputPairlists) ) {
      iout << iINFO << pairlistWarnings <<
        " pairlist warnings in past " << simParams->outputPairlists <<
	" steps.\n" << endi;
      pairlistWarnings = 0;
    }
    
    // NO CALCULATIONS OR REDUCTIONS BEYOND THIS POINT!!!
    if ( ! minimize &&  step % simParameters->outputEnergies ) return;
    // ONLY OUTPUT SHOULD OCCUR BELOW THIS LINE!!!

    if ( marginViolations ) {
      iout << iERROR << marginViolations <<
        " margin violations detected since previous energy output.\n" << endi;
    }
    marginViolations = 0;

    if ( (step % (10 * (minimize?1:simParameters->outputEnergies) ) ) == 0 )
    {
	iout << "ETITLE:      TS";
	iout << FORMAT("BOND");
	iout << FORMAT("ANGLE");
	iout << FORMAT("DIHED");
	iout << FORMAT("IMPRP");
        iout << "     ";
	iout << FORMAT("ELECT");
	iout << FORMAT("VDW");
	iout << FORMAT("BOUNDARY");
	iout << FORMAT("MISC");
	iout << FORMAT("KINETIC");
        iout << "     ";
	iout << FORMAT("TOTAL");
	iout << FORMAT("TEMP");
	iout << FORMAT("TOTAL2");
	iout << FORMAT("TOTAL3");
	iout << FORMAT("TEMPAVG");
	if ( volume != 0. ) {
          iout << "     ";
	  iout << FORMAT("PRESSURE");
	  iout << FORMAT("GPRESSURE");
	  iout << FORMAT("VOLUME");
	  iout << FORMAT("PRESSAVG");
	  iout << FORMAT("GPRESSAVG");
	}
	iout << "\n\n" << endi;
    }

    // N.B.  HP's aCC compiler merges FORMAT calls in the same expression.
    //       Need separate statements because data returned in static array.

    iout << ETITLE(step);
    iout << FORMAT(bondEnergy);
    iout << FORMAT(angleEnergy);
    iout << FORMAT(dihedralEnergy);
    iout << FORMAT(improperEnergy);
    iout << "     ";
    iout << FORMAT(electEnergy+electEnergySlow);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(miscEnergy);
    iout << FORMAT(kineticEnergy);
    iout << "     ";
    iout << FORMAT(totalEnergy);
    iout << FORMAT(temperature);
    iout << FORMAT(flatEnergy);
    iout << FORMAT(smoothEnergy);
    iout << FORMAT(temp_avg/avg_count);
    if ( volume != 0. )
    {
        iout << "     ";
	iout << FORMAT(trace(pressure)*PRESSUREFACTOR/3.);
	iout << FORMAT(trace(groupPressure)*PRESSUREFACTOR/3.);
	iout << FORMAT(volume);
	iout << FORMAT(pressure_avg*PRESSUREFACTOR/avg_count);
	iout << FORMAT(groupPressure_avg*PRESSUREFACTOR/avg_count);
    }
    iout << "\n\n" << endi;

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
     char webout[80];
     sprintf(webout,"%d %d %d %d",(int)totalEnergy,
	     (int)(potentialEnergy),
	     (int)kineticEnergy,(int)temperature);
     CApplicationDepositNode0Data(webout);
#endif

    if (simParameters->pairInteractionOn) {
      iout << "PAIR INTERACTION:";
      iout << " STEP: " << step;
      iout << " VDW_FORCE: ";
      iout << FORMAT(pairVDWForce.x);
      iout << FORMAT(pairVDWForce.y);
      iout << FORMAT(pairVDWForce.z);
      iout << " ELECT_FORCE: ";
      iout << FORMAT(pairElectForce.x);
      iout << FORMAT(pairElectForce.y);
      iout << FORMAT(pairElectForce.z);
      iout << "\n" << endi;
    }
    temp_avg = 0;
    pressure_avg = 0;
    groupPressure_avg = 0;
    avg_count = 0;

}

void Controller::writeExtendedSystemLabels(ofstream &file) {
  Lattice &lattice = state->lattice;
  file << "#$LABELS step";
  if ( lattice.a_p() ) file << " a_x a_y a_z";
  if ( lattice.b_p() ) file << " b_x b_y b_z";
  if ( lattice.c_p() ) file << " c_x c_y c_z";
  file << " o_x o_y o_z";
  if ( simParams->langevinPistonOn ) {
    file << " s_x s_y s_z s_u s_v s_w";
  }
  file << endl;
}

void Controller::writeExtendedSystemData(int step, ofstream &file) {
  Lattice &lattice = state->lattice;
  file << step;
    if ( lattice.a_p() ) file << " " << lattice.a().x << " " << lattice.a().y << " " << lattice.a().z;
    if ( lattice.b_p() ) file << " " << lattice.b().x << " " << lattice.b().y << " " << lattice.b().z;
    if ( lattice.c_p() ) file << " " << lattice.c().x << " " << lattice.c().y << " " << lattice.c().z;
    file << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
  if ( simParams->langevinPistonOn ) {
    Vector strainRate = diagonal(langevinPiston_strainRate);
    Vector strainRate2 = off_diagonal(langevinPiston_strainRate);
    file << " " << strainRate.x;
    file << " " << strainRate.y;
    file << " " << strainRate.z;
    file << " " << strainRate2.x;
    file << " " << strainRate2.y;
    file << " " << strainRate2.z;
  }
  file << endl;
}

void Controller::enqueueCollections(int timestep)
{
  if ( Output::coordinateNeeded(timestep) )
    collection->enqueuePositions(timestep,state->lattice);
  if ( Output::velocityNeeded(timestep) )
    collection->enqueueVelocities(timestep);
}

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
static char *FEPTITLE(int X)
{ 
  static char tmp_string[21];
  sprintf(tmp_string, "FepEnergy:%6d   ",X);
  return tmp_string;
}

void Controller::outputFepEnergy(int step) {
 if (simParams->fepOn) {
  const int stepInRun = step - simParams->firstTimestep;
  const int fepEquilSteps = simParams->fepEquilSteps;
  const BigReal lambda = simParams->lambda;
  const BigReal lambda2 = simParams->lambda2;
  if (stepInRun == 0 || stepInRun == fepEquilSteps) {
    FepNo = 0;
    exp_dE_ByRT = 0.0;
    net_dE = 0.0;
  }
  if (stepInRun == fepEquilSteps) {
    fepFile << "#" << fepEquilSteps << " STEPS OF EQUILIBRATION AT "
            << "LAMBDA " << simParams->lambda << " COMPLETED\n"
            << "#STARTING COLLECTION OF ENSEMBLE AVERAGE" << endl;
  }
  if (simParams->fepOutFreq && ((step%simParams->fepOutFreq)==0)) {
    if (!fepFile.rdbuf()->is_open()) {
      fepSum = 0.0;
      NAMD_backup_file(simParams->fepOutFile);
      fepFile.open(simParams->fepOutFile);
      iout << "OPENING FEP ENERGY OUTPUT FILE\n" << endi;
      fepFile << "#           STEP           Elec          "
              << "           vdW              dE        dE_avg        Temp        dG\n"
              << "#                    l         l+dl   "
              << "       l          l+dl      E(l+dl)-E(l)" << endl;
    }
    if (stepInRun == 0) {
      fepFile << "#RESCALED CHARGE FOR NEW FEP WINDOW "
              << "LAMBDA SET TO " << lambda << " LAMBDA2 " << lambda2 << endl;
    }
    writeFepEnergyData(step, fepFile);
    fepFile.flush();
  }
  if (step == simParams->N) {
    fepSum = fepSum + dG;
    fepFile << "#Net free energy change at end of lambda@"<<lambda <<" is " << fepSum << endl;
  }
 }
}

void Controller::writeFepEnergyData(int step, ofstream &file) {
  BigReal eeng = electEnergy+electEnergySlow;
  BigReal eeng_f = electEnergy_f + electEnergySlow_f;
  BigReal dE = eeng_f + ljEnergy_f - eeng - ljEnergy;
  BigReal RT = BOLTZMAN *temperature;
  FepNo++;
  exp_dE_ByRT += exp(-dE/RT);
  net_dE += dE;
  dG = -(RT * log(exp_dE_ByRT/FepNo));
  BigReal dE_avg = net_dE/FepNo;
  fepFile << FEPTITLE(step);
  fepFile << FORMAT(eeng);
  fepFile << FORMAT(eeng_f);
  fepFile << FORMAT(ljEnergy);
  fepFile << FORMAT(ljEnergy_f);
  fepFile << FORMAT(dE);
  fepFile << FORMAT(dE_avg);
  fepFile << FORMAT(temperature);
  fepFile << FORMAT(dG);
  fepFile << endl;
}
//fepe

void Controller::outputExtendedSystem(int step)
{

  if ( step >= 0 ) {

    // Write out eXtended System Trajectory (XST) file
    if ( simParams->xstFrequency &&
         ((step % simParams->xstFrequency) == 0) )
    {
      if ( ! xstFile.rdbuf()->is_open() )
      {
        NAMD_backup_file(simParams->xstFilename);
        xstFile.open(simParams->xstFilename);
        iout << "OPENING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
        xstFile << "# NAMD extended system trajectory file" << endl;
        writeExtendedSystemLabels(xstFile);
      }
      writeExtendedSystemData(step,xstFile);
      xstFile.flush();
    }

    // Write out eXtended System Configuration (XSC) files
    //  Output a restart file
    if ( simParams->restartFrequency &&
         ((step % simParams->restartFrequency) == 0) &&
         (step != simParams->firstTimestep) )
    {
      iout << "WRITING EXTENDED SYSTEM TO RESTART FILE AT STEP "
		<< step << "\n" << endi;
      char fname[140];
      strcpy(fname, simParams->restartFilename);
      if ( simParams->restartSave ) {
        char timestepstr[20];
        sprintf(timestepstr,".%d",step);
        strcat(fname, timestepstr);
      }
      strcat(fname, ".xsc");
      NAMD_backup_file(fname,".old");
      ofstream xscFile(fname);
      if (!xscFile) {
        char err_msg[257];
        sprintf(err_msg, "Error opening XSC restart file %s",fname);
        NAMD_err(err_msg);
      } 
      xscFile << "# NAMD extended system configuration restart file" << endl;
      writeExtendedSystemLabels(xscFile);
      writeExtendedSystemData(step,xscFile);
      if (!xscFile) {
        char err_msg[257];
        sprintf(err_msg, "Error writing XSC restart file %s",fname);
        NAMD_err(err_msg);
      } 
    }

  }

  //  Output final coordinates
  if (step == FILE_OUTPUT || step == END_OF_RUN)
  {
    iout << "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP "
		<< simParams->N << "\n" << endi;
    static char fname[140];
    strcpy(fname, simParams->outputFilename);
    strcat(fname, ".xsc");
    NAMD_backup_file(fname);
    ofstream xscFile(fname);
    if (!xscFile) {
      char err_msg[257];
      sprintf(err_msg, "Error opening XSC output file %s",fname);
      NAMD_err(err_msg);
    } 
    xscFile << "# NAMD extended system configuration output file" << endl;
    writeExtendedSystemLabels(xscFile);
    writeExtendedSystemData(simParams->N,xscFile);
    if (!xscFile) {
      char err_msg[257];
      sprintf(err_msg, "Error writing XSC output file %s",fname);
      NAMD_err(err_msg);
    } 
  }

  //  Close trajectory file
  if (step == END_OF_RUN) {
    if ( xstFile.rdbuf()->is_open() ) {
      xstFile.close();
      iout << "CLOSING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
    }
  }

}

void Controller::rebalanceLoad(int)
{
  if ( ! ldbSteps ) { 
    ldbSteps = LdbCoordinator::Object()->steps();
  }
  if ( ! --ldbSteps ) {
    startBenchTime -= CmiWallTimer();
    LdbCoordinator::Object()->rebalance(this);
    startBenchTime += CmiWallTimer();
  }
}

void Controller::cycleBarrier(int doBarrier, int step) {
#if USE_BARRIER
	if (doBarrier) {
	  broadcast->cycleBarrier.publish(step,1);
	  CkPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);
	}
#endif
}

void Controller::terminate(void) {
  BackEnd::awaken();
  CthFree(thread);
  CthSuspend();
}

