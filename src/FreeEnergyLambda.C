// written by David Hurwitz, March to May 1998.

#include <string.h>
#include <iostream.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyLambda.h"
#include "Vector.h"
#include "FreeEnergyVector.h"

// initialize static member variables
int  ALambdaControl::m_CurrStep = 0;

ALambdaControl::ALambdaControl() {
//------------------------------------------------------------------------
// initialize member variables
//------------------------------------------------------------------------
  // initialize these for the first block.
  m_StartStep = 0;
  m_LambdaKf  = 1.0;   // applied to Kf
  m_LambdaRef = 0.0;   // applied to reference (pos, dist, angle, dihe)

  // set the rest to illegal values.
  m_Task =           kUnknownTask;
  m_NumSteps =      -1;
  m_NumEquilSteps = -1;
  m_NumAccumSteps = -1;
  m_NumPrintSteps = -1;
  m_NumRepeats =    -1;
  m_StopStep =      -1;
}


void ALambdaControl::Init(ALambdaControl& PriorBlock) {
//------------------------------------------------------------------------
// initialize this object using the settings for the prior object.
//------------------------------------------------------------------------
  m_Task =          PriorBlock.m_Task;
  m_NumSteps =      PriorBlock.m_NumSteps;
  m_NumEquilSteps = PriorBlock.m_NumEquilSteps;
  m_NumAccumSteps = PriorBlock.m_NumAccumSteps;
  m_NumPrintSteps = PriorBlock.m_NumPrintSteps;
  m_NumRepeats =    PriorBlock.m_NumRepeats;
  switch (PriorBlock.m_Task) {

    case kUp:       m_LambdaKf=1.0;  m_LambdaRef=1.0; break;
    case kStepUp:   m_LambdaKf=1.0;  m_LambdaRef=1.0; break;
    case kDown:     m_LambdaKf=1.0;  m_LambdaRef=0.0; break;
    case kStepDown: m_LambdaKf=1.0;  m_LambdaRef=0.0; break;
    case kStop:     m_LambdaKf=1.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;

    case kGrow:     m_LambdaKf=1.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kStepGrow: m_LambdaKf=1.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kFade:     m_LambdaKf=0.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kStepFade: m_LambdaKf=0.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kNoGrow:   m_LambdaKf =  PriorBlock.m_LambdaKf;
                    m_LambdaRef = PriorBlock.m_LambdaRef;  break;
    default:        ASSERT(kFalse); break;  //should never get here

  }
  m_StartStep = PriorBlock.GetLastStep() + 1;
}


ALambdaControl& ALambdaControl::operator= (ALambdaControl& PmfBlock) {
//------------------------------------------------------------------------
// copy everything from PmfBlock to this block.
//------------------------------------------------------------------------
  m_NumSteps =      PmfBlock.m_NumSteps;
  m_NumEquilSteps = PmfBlock.m_NumEquilSteps;
  m_NumAccumSteps = PmfBlock.m_NumAccumSteps;
  m_NumRepeats =    PmfBlock.m_NumRepeats;
  m_NumPrintSteps = PmfBlock.m_NumPrintSteps;
  m_StartStep =     PmfBlock.m_StartStep;
  m_StopStep =      PmfBlock.m_StopStep;
  m_LambdaKf =      PmfBlock.m_LambdaKf;
  m_LambdaRef =     PmfBlock.m_LambdaRef;
  m_Task =          PmfBlock.m_Task;
  return(*this);
}


int ALambdaControl::GetNumSteps() {
//------------------------------------------------------------------------
// get the number of steps needed for this pmf or mcti block
//------------------------------------------------------------------------
  // make sure m_StopStep is calculated
  GetLastStep();
  
  return( (m_StopStep - m_StartStep) + 1 );
}


int ALambdaControl::GetLastStep() {
//------------------------------------------------------------------------
// get the last step of this task
//------------------------------------------------------------------------
  // if it's already calculated, just return it
  if (m_StopStep > 0) {
    return(m_StopStep);
  }
  // otherwise calculate it
  switch (m_Task) {
    case kStepUp:
    case kStepDown:
    case kStepGrow:
    case kStepFade:
      m_StopStep = m_StartStep +
                  (m_NumAccumSteps+m_NumEquilSteps) * (m_NumRepeats+1) - 1;
      break;
    default:
      m_StopStep = m_StartStep + m_NumSteps - 1;
      break;
  }
  // and return it
  return(m_StopStep);
}


void ALambdaControl::GetTaskStr(char* Str) {
//------------------------------------------------------------------------
// get a string that describes this task
//------------------------------------------------------------------------
  switch (m_Task) {
    case kUp:        strcpy(Str, "Up");            break;
    case kDown:      strcpy(Str, "Down");          break;
    case kStop:      strcpy(Str, "Stop");          break;
    case kGrow:      strcpy(Str, "Grow");          break;
    case kFade:      strcpy(Str, "Fade");          break;
    case kNoGrow:    strcpy(Str, "NoGrow");        break;
    case kStepUp:    strcpy(Str, "StepUp");        break;
    case kStepDown:  strcpy(Str, "StepDown");      break;
    case kStepGrow:  strcpy(Str, "StepGrow");      break;
    case kStepFade:  strcpy(Str, "StepFade");      break;
    default:         strcpy(Str, "Bug Alert!!!");  break;
  }
}


Bool_t ALambdaControl::IsTimeToPrint(double dT) {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's time
// to print out dU/dLambda.
// if it's time to print, write out the current time
//------------------------------------------------------------------------
  double  Time;
  Bool_t  RetVal=kFalse;
  char    Str[100];

  // if printing is required
  if (m_NumPrintSteps > 0) {
    // if number-of-steps from StartStep is an even multiple of NumPrintSteps
    // then it might be time to print
    if ( ((m_CurrStep-m_StartStep) % m_NumPrintSteps) == 0) {
      // for mcti blocks
      if ((m_Task==kStepUp)   || (m_Task==kStepDown) ||
          (m_Task==kStepGrow) || (m_Task==kStepFade)) {
        // it's only time to print if we're no longer equilibrating
        if ( ((m_CurrStep-m_StartStep) % (m_NumEquilSteps+m_NumAccumSteps))
             >= m_NumEquilSteps) {
          RetVal = kTrue;
        }
      }
      // for the other tasks, it's ok to print
      else {
        RetVal = kTrue;
      }
    }
  }
  if (RetVal) {
    iout << "FreeEnergy: " << endl << endi;
    iout << "FreeEnergy: ";
    iout << "Time Step = "  << m_CurrStep            <<    ",  ";
    iout << "Time = ";
    // calculate current time in femto-seconds
    Time = (double)m_CurrStep * dT;
    // write out time in either fs, ps, or ns
    if      (Time < 1000)    {iout << Time           << " fs,  ";}
    else if (Time < 1000000) {iout << Time/1000.0    << " ps,  ";}
    else                     {iout << Time/1000000.0 << " ns,  ";}
    iout << "Lambda_Kf = "  << m_LambdaKf            <<    ",  ";
    iout << "Lambda_Ref = " << m_LambdaRef           <<     "  ";
    GetTaskStr(Str);
    iout << "(" << Str << ")";
    iout << endl << endi;
    iout << "FreeEnergy: ";
    iout << "------------------------------------------------";
    iout << "-------------------";
    iout << endl << endi;
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsActive() {
//------------------------------------------------------------------------
// determine if this object is currently active
//------------------------------------------------------------------------
  if ( (m_CurrStep>=m_StartStep) && (m_CurrStep<=GetLastStep()) ) {
    return(kTrue);
  }
  else {
    return(kFalse);
  }
}


double ALambdaControl::GetLambdaKf() {
//------------------------------------------------------------------------
// calculate LambdaKf for grow, fade, stepgrow, and stepfade.
// LambdaKf=1.0, for up, down, stepup, stepdown, and stop.
// for nogrow, LambdaKf is Lambda from the config file.
//------------------------------------------------------------------------
  int     N;
  double  CurrStep  = m_CurrStep;
  double  StartStep = m_StartStep;
  double  NumSteps  = m_NumSteps;
  double  NumEquilSteps = m_NumEquilSteps;
  double  NumAccumSteps = m_NumAccumSteps;
  double  NumRepeats = m_NumRepeats;

  if (IsActive()) {
    switch (m_Task) {
      case kGrow:
        m_LambdaKf = (CurrStep-StartStep)/NumSteps;
        break;
      case kFade:
        m_LambdaKf = 1.0-(CurrStep-StartStep)/NumSteps;
        break;
      case kStepGrow:
        N = (int) ( (CurrStep-StartStep) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaKf = N/NumRepeats;
        break;
      case kStepFade:
        N = (int) ( (CurrStep-StartStep) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaKf = 1.0 - N/NumRepeats;
        break;
      case kNoGrow:
        break;              // return prior setting of m_LambdaKf
      default:
        m_LambdaKf=1.0;
    }
  }
  else {
    m_LambdaKf=1.0;
  }
  return(m_LambdaKf);
}


double ALambdaControl::GetLambdaRef() {
//------------------------------------------------------------------------
// calculate LambdaRef for up, down, stepup, and stepdown.
// for stop, LambdaRef is Lambda from the config file.
// for grow, fade, stepgrow, stepfade, and nogrow,
//   LambdaRef is LambdaT from the config file.
//------------------------------------------------------------------------
  int     N;
  double  CurrStep  = m_CurrStep;
  double  StartStep = m_StartStep;
  double  NumSteps  = m_NumSteps;
  double  NumEquilSteps = m_NumEquilSteps;
  double  NumAccumSteps = m_NumAccumSteps;
  double  NumRepeats = m_NumRepeats;

  if (IsActive()) {
    switch (m_Task) {
      case kUp:
        m_LambdaRef = (CurrStep-StartStep)/NumSteps;
        break;
      case kDown:
        m_LambdaRef = 1.0-(CurrStep-StartStep)/NumSteps;
        break;
      case kStepUp:
        N = (int) ( (CurrStep-StartStep) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaRef = N/NumRepeats;
        break;
      case kStepDown:
        N = (int) ( (CurrStep-StartStep) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaRef = 1.0 - N/NumRepeats;
      default: 
        break;             // return prior setting of m_LambdaRef
    }
  }
  else {
    m_LambdaRef=0.0;
  }
  return(m_LambdaRef);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker $		$State $
 *	$Revision $	$Date $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: FreeEnergyLambda.C,v $
 * Revision 1.2  1998/05/22 19:08:31  hurwitz
 * Do NAMD_die if there aren't enough steps to complete all pmf & mcti blocks
 *
 *
 ***************************************************************************/
