// written by David Hurwitz, March to May 1998.

#if !defined(LAMBDA_HPP)
  #define LAMBDA_HPP

class ALambdaControl {
private:

  // don't forget to change operator= if member variables change
  int     m_NumSteps;         // for pmf block
  int     m_NumEquilSteps;    // for mcti block
  int     m_NumAccumSteps;    // "
  int     m_NumRepeats;       // "
  int     m_NumPrintSteps;    // for pmf & mcti blocks
  int     m_StartStep;        // "
  int     m_StopStep;         // "
  double  m_LambdaKf;         // "
  double  m_LambdaRef;        // "
  task_t  m_Task;             // "
  double  m_Sum_dU_dLambda;   // for accumulating dU/dLambda
  int     m_Num_dU_dLambda;   // number averaged
  double  m_MCTI_Integration; // for accumulating <dU/dLambda> * dLambda

  static int  m_CurrStep;     // for all pmf & mcti blocks

public:
  ALambdaControl();
  void    Init(ALambdaControl& PriorBlock);
  double  GetLambdaKf();
  double  GetLambdaRef();
  Bool_t  IsActive();
  Bool_t  IsTimeToPrint();
  Bool_t  IsFirstStep();
  Bool_t  IsTimeToPrint_dU_dLambda();
  Bool_t  IsTimeToClearAccumulator();
  Bool_t  IsEndOf_MCTI_Step();
  Bool_t  IsEndOf_MCTI();
  void    PrintHeader(double dT);
  void    PrintLambdaHeader(double dT);
  void    IncCurrStep() {m_CurrStep++;}
  ALambdaControl&  operator= (ALambdaControl& PmfBlock);
  void    GetTaskStr(char* Str);
  void    GetPaddedTaskStr(char* Str);
  void    Integrate_MCTI();
  void    Accumulate(double dU_dLambda);
  double  GetIntegration();
  double  GetAccumulation();
  void    ZeroAccumulator() {
    m_Sum_dU_dLambda = 0.0;
    m_Num_dU_dLambda = 0;
  }

  int    GetNumSteps();
  int    GetNumStepsSoFar()            {return(m_CurrStep-m_StartStep);}
  int    GetNumAccumStepsSoFar();
  int    GetNum_dU_dLambda()           {return(m_Num_dU_dLambda);}
  void   SetNumSteps(int Steps)        {m_NumSteps=Steps;}
  void   SetNumEquilSteps(int Steps)   {m_NumEquilSteps=Steps;}
  void   SetNumAccumSteps(int Steps)   {m_NumAccumSteps=Steps;}
  void   SetNumPrintSteps(int Steps)   {m_NumPrintSteps=Steps;}
  void   SetNumRepeats(int Repeats)    {m_NumRepeats=Repeats;}
  void   SetStartStep(int Step)        {m_StartStep=Step;}
  void   SetStopStep(int Step)         {m_StopStep=Step;}
  void   SetLambdaKf(double LambdaKf)  {m_LambdaKf=LambdaKf;}
  void   SetLambdaRef(double LambdaRef){m_LambdaRef=LambdaRef;}
  void   SetTask(task_t Task)          {m_Task=Task;}
  task_t GetTask()                     {return(m_Task);}

private:
  Bool_t IsLastStep();
  int    GetLastStep();
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
 * $Log: FreeEnergyLambda.h,v $
 * Revision 1.4  1998/09/20 16:34:58  hurwitz
 * make sure Lambda control objects start and stop on just the right step.
 * made output shorter and more readable (compile with _VERBOSE_PMF for old output)
 * : ----------------------------------------------------------------------
 *
 * Revision 1.3  1998/06/05 22:54:41  hurwitz
 * accumulate dU/dLambda for free energy calculation
 *
 * Revision 1.2  1998/05/22 19:08:31  hurwitz
 * Do NAMD_die if there aren't enough steps to complete all pmf & mcti blocks
 *
 *
 ***************************************************************************/
