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

  static int  m_CurrStep;     // for all pmf & mcti blocks

public:
  ALambdaControl();
  void    Init(ALambdaControl& PriorBlock);
  double  GetLambdaKf();
  double  GetLambdaRef();
  Bool_t  IsActive();
  Bool_t  IsTimeToPrint(double dT);
  void    IncCurrStep() {m_CurrStep++;}
  ALambdaControl&  operator= (ALambdaControl& PmfBlock);
  void    GetTaskStr(char* Str);

  int    GetNumSteps();
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
  int  GetLastStep();
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
 * Revision 1.2  1998/05/22 19:08:31  hurwitz
 * Do NAMD_die if there aren't enough steps to complete all pmf & mcti blocks
 *
 *
 ***************************************************************************/
