// written by David Hurwitz, March to May 1998.

#include <memory.h>
#include <iostream.h>
#include <iomanip.h>
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"
#include "FreeEnergyGroup.h"

ALambdaManager::ALambdaManager() {
//------------------------------------------------------------------------
// make room for some LambdaControl objects.
//------------------------------------------------------------------------
  m_ActiveIndex = 0;
  m_NumObjects = 0;
  m_pPmfBlocks = new ALambdaControl[kLambdaNumToStart];
  m_MaxNum = kLambdaNumToStart;
}


ALambdaManager::~ALambdaManager() {
//------------------------------------------------------------------------
// return borrowed memory to the free store.
//------------------------------------------------------------------------
  delete []m_pPmfBlocks;
}


Bool_t ALambdaManager::IsTimeToPrint(double dT) {
//------------------------------------------------------------------------
// ASSUMING that the m_ActiveIndex'th LambdaControl is active,
// decide if it's time to print out dU/dLambda
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].IsTimeToPrint(dT));
}


Bool_t ALambdaManager::GetLambdas(double& LambdaKf, double& LambdaRef) {
//------------------------------------------------------------------------
// get LambdaKf and LambdaRef from the active LambdaControl
//
// return(kTrue) if an active LambdaControl is found
// return(kFalse) if all LambdaControls have expired
//------------------------------------------------------------------------
  // don't continue if all LamdaControl's have expired
  if (m_ActiveIndex == m_NumObjects) {
    return(kFalse);
  }

  // if the m_ActiveIndex'th LambdaControl is no longer active
  if ( !(*this)[m_ActiveIndex].IsActive()) {
    // move on to the next one
    m_ActiveIndex++;
    // if there is no next object, return KFalse
    if (m_ActiveIndex == m_NumObjects) {
      return(kFalse);
    }
    // otherwise, make sure the next one's active
    else {
      ASSERT( (*this)[m_ActiveIndex].IsActive() );
    }
  }
  // return LambdaKf and LambdaRef from the active LambdaControl
  LambdaKf =  (*this)[m_ActiveIndex].GetLambdaKf();
  LambdaRef = (*this)[m_ActiveIndex].GetLambdaRef();
  return(kTrue);
}


void ALambdaManager::Clear() {
//------------------------------------------------------------------------
// leave memory allocation alone.
//------------------------------------------------------------------------
  m_NumObjects = 0;
}


int ALambdaManager::Add(ALambdaControl& PmfBlock) {
//------------------------------------------------------------------------
// add an object to the list.  if there's not enough room, make room.
// return an index to the added oject.
//------------------------------------------------------------------------
  ALambdaControl*  pPmfBlocks;

  // if there's no room for another object
  if (m_NumObjects == m_MaxNum) {
    // create an array with more space
    m_MaxNum *= kLambdaMultiplier;
    pPmfBlocks = new ALambdaControl[m_MaxNum];
    // copy from the full array to the new one
    for (int i=0; i<m_NumObjects; i++) {
      pPmfBlocks[i] = m_pPmfBlocks[i];
    }
    // return the space used for the full array
    delete []m_pPmfBlocks;
    // point to the bigger array
    m_pPmfBlocks = pPmfBlocks;
  }
  // add the object to the array
  m_pPmfBlocks[m_NumObjects] = PmfBlock;
  m_NumObjects++;
  return(m_NumObjects-1);
}


ALambdaControl& ALambdaManager::operator[] (int Index) {
//------------------------------------------------------------------------
// return an object from this group of objects.
//------------------------------------------------------------------------
  ASSERT((Index>=0) && (Index<m_NumObjects));
  return(m_pPmfBlocks[Index]);
}
