//-------------------------------------------------------------------------
// ARestraintManager contains a (potentially long) list of restraint POINTERS
// written by David Hurwitz, March to May 1998.
//-------------------------------------------------------------------------
#if !defined(RMGR_HPP)
  #define RMGR_HPP

typedef ARestraint* pRestr;

// to start, there's room for this number of restraint pointers
// each time array size is exceeded, its size is increased by this many times.
const int kNumToStart = 1024;
const int kMultiplier = 4;

class ComputeFreeEnergy;

class ARestraintManager {
private:
  ARestraint** m_ppRestraints; // list of restraint pointers
  int  m_NumRestraints;        // number of pointers in the list
  int  m_MaxNum;               // max num pointers without allocating more mem
  AFixedPosRestraint  m_Dummy; // for setting ARestraint statics

public:
  ARestraintManager();
  ~ARestraintManager();
  ARestraint*  operator[] (int Index);
  void  Add(ARestraint* pRestraint);
  int   GetNumRestraints() {return(m_NumRestraints);}
  void  UpdateCOMs(ComputeFreeEnergy& CFE);
  void  AddForces(ComputeFreeEnergy& CFE);
  void  PrintInfo();
  void  PrintPreInfo(int Index);
  void  SetLambdaKf(double LambdaKf)   {m_Dummy.SetLambdaKf(LambdaKf);}
  void  SetLambdaRef(double LambdaRef) {m_Dummy.SetLambdaRef(LambdaRef);}
  void  SetLambdas(double LambdaKf, double LambdaRef)
  {
     m_Dummy.SetLambdaKf(LambdaKf);
     m_Dummy.SetLambdaRef(LambdaRef);
  }
};

#endif
