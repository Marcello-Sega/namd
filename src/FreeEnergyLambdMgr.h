//-----------------------------------------------------------------------------
// ALambdaManager contains a (potentially long) list of LambdaControl objects.
// written by David Hurwitz, March to May 1998.
//-----------------------------------------------------------------------------
#if !defined(LAMBDMGR_HPP)
  #define LAMBDMGR_HPP

const int kLambdaNumToStart = 16; // to start, there's room for this num objects.
const int kLambdaMultiplier = 4;  // each time array size is exceeded,
                                  // its size is increased by this many times.

class ALambdaManager {
private:
  ALambdaControl*  m_pPmfBlocks;  // the list of objects
  int  m_NumObjects;              // the number of objects in the list
  int  m_MaxNum;                  // the maximum number of objects allowed in
                                  // the list without allocating more memory
  int  m_ActiveIndex;             // keep track of which LambdaControl is active
  ALambdaControl m_Dummy;

public:
  ALambdaManager();
  ~ALambdaManager();
  ALambdaControl&  operator[] (int index);
  void    Clear();
  int     Add(ALambdaControl& PmfBlock);
  int     GetNumObjects() {return(m_NumObjects);}
  Bool_t  GetLambdas(double& LambdaKf, double& LambdaRef);
  Bool_t  IsTimeToPrint();
  Bool_t  IsTimeToPrint_dU_dLambda();
  Bool_t  IsTimeToClearAccumulator();
  void    PrintHeader(double dT);
  void    Print_dU_dLambda_Summary(double Sum_dU_dLambdas);
  void    IncCurrStep() {m_Dummy.IncCurrStep();}
  int     GetTotalNumSteps();
  void    Accumulate(double dU_dLambda);
  double  GetAccumulation();
  void    ZeroAccumulator();
  int     GetNum_dU_dLambda();
};

#endif
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
 * $Log: FreeEnergyLambdMgr.h,v $
 * Revision 1.3  1998/06/05 22:54:40  hurwitz
 * accumulate dU/dLambda for free energy calculation
 *
 * Revision 1.2  1998/05/22 19:08:30  hurwitz
 * Do NAMD_die if there aren't enough steps to complete all pmf & mcti blocks
 *
 *
 ***************************************************************************/
