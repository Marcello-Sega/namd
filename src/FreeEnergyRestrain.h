// written by David Hurwitz, March to May 1998.

#if !defined(RESTRAINT_HPP)
  #define RESTRAINT_HPP

class ComputeFreeEnergy;

//****************************************************************************
//  ARestraint:
//****************************************************************************
class ARestraint {
//-----------------------------------------------------------------------------
// ARestraint is the base class for all restraints
// Here's a tree of derived and base classes:
//             -----------------ARestraint-----------------
//           /                /           \                 \
//         /                /               \                 \
//        |                 |                |                 |
//   APosRestraint   ADistRestraint   AnAngleRestraint   ADiheRestraint
//        |                 |                |                 |
//        |                 |                |                 |
//   AFixedPosRestraint     |         AFixedAngleRestraint     |            
//   ABoundPosRestraint     |         ABoundAngleRestraint     |
//   AForcingPosRestraint   |         AForcingAngleRestraint   |
//                          |                                  |
//                   AFixedDistRestraint                 AFixedDiheRestraint
//                   ABoundDistRestraint                 ABoundDiheRestraint
//                   AForcingDistRestraint               AForcingDiheRestraint
//-----------------------------------------------------------------------------
protected:
  double    m_Kf;
  int       m_NumGroups;
  AGroup*   m_pGroups;
  AVector*  m_pCOMs;

  // lambda is the same for all forcing restraints
  static double  m_LambdaKf;
  static double  m_LambdaRef;

public:
  ARestraint();
  ~ARestraint();
  int     GetNumGroups()    {return(m_NumGroups);}
  void    SetKf(double Kf)  {m_Kf=Kf;}
  double  GetKf()           {return(m_Kf);}
  void    SetLambdaKf(double LambdaKf)   {m_LambdaKf=LambdaKf;}
  void    SetLambdaRef(double LambdaRef) {m_LambdaRef=LambdaRef;}
  double  GetLambdaKf()                  {return(m_LambdaKf);}
  double  GetLambdaRef()                 {return(m_LambdaRef);}
  void    SetGroup(AGroup& Group, int GroupIndex);
  void    SetGroups(AGroup& Group1);
  void    SetGroups(AGroup& Group1, AGroup& Group2);
  void    SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3);
  void    SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3, AGroup& Group4);
  void    UpdateCOMs(ComputeFreeEnergy& CFE);
  void    DistributeForce(int WhichGroup, AVector Force, ComputeFreeEnergy& CFE);

//---------------- only for testing------------//
#if defined(_DEBUG)                            //
  void  SetCOM(int Index, AVector& Pos);       // 
#endif                                         //
//---------------------------------------------//

  // pure virtual functions
  virtual AVector  GetGradient(int WhichGroup) = 0;
  virtual double   GetEnergy() = 0;
  virtual void     GetStr(char* Str) = 0;
  virtual void     PrintInfo() = 0;
  // virtual functions that are only meaningful in the forcing restraints
  virtual Bool_t   IsForcing()           {return(kFalse);}
  virtual double   Get_dU_dLambda()      {return(0.0);}
  // virtual functions that should only be called in the derived classes
  virtual void  SetRefPos(AVector)       {ASSERT(kFalse);}
  virtual void  SetRefDist(double)       {ASSERT(kFalse);}
  virtual void  SetRefAngle(double)      {ASSERT(kFalse);}
  virtual void  SetBound(Bound_t)        {ASSERT(kFalse);}
  virtual void  SetLowerAngle(double)    {ASSERT(kFalse);}
  virtual void  SetUpperAngle(double)    {ASSERT(kFalse);}
  virtual void  SetIntervalAngle(double) {ASSERT(kFalse);}
  virtual void  SetStartPos(AVector)     {ASSERT(kFalse);}
  virtual void  SetStopPos(AVector)      {ASSERT(kFalse);}
  virtual void  SetStartDist(double)     {ASSERT(kFalse);}
  virtual void  SetStopDist(double)      {ASSERT(kFalse);}
  virtual void  SetStartAngle(double)    {ASSERT(kFalse);}
  virtual void  SetStopAngle(double)     {ASSERT(kFalse);}
protected:
  double  GetAngle(AVector& A, AVector& B, AVector& C);
  double  GetDihe(AVector& A, AVector& B, AVector& C, AVector& D);
  void    EarlyExit(char* Str, int AtomID);
};


//****************************************************************************
//  APosRestraint, ADistRestraint, AnAngleRestraint, ADiheRestraint:
//****************************************************************************
class APosRestraint : public ARestraint {
//-------------------------------------------------------------------
// APosRestraint is a derived class of ARestraint
//-------------------------------------------------------------------
public:
  APosRestraint();
  void    PrintInfo();
protected:
  double  GetE(AVector RefPos, double LambdaKf=1.0);
  AVector GetGrad(int WhichGroup, AVector RefPos, double LambdaKf=1.0);
};

class ADistRestraint : public ARestraint {
//-------------------------------------------------------------------
// ADistRestraint is a derived class of ARestraint
//-------------------------------------------------------------------
public:
  ADistRestraint();
  void    PrintInfo();
protected:
  double  GetE(double RefDist, double LambdaKf=1.0);
  AVector GetGrad(int WhichGroup, double RefDist, double LambdaKf=1.0);
};

class AnAngleRestraint : public ARestraint {
//-------------------------------------------------------------------
// AnAngleRestraint is a derived class of ARestraint
//-------------------------------------------------------------------
public:
  AnAngleRestraint();
  void    PrintInfo();
protected:
  double  GetE(double RefAngle, double LambdaKf=1.0);
  AVector GetGrad(int WhichGroup, double RefAngle, double LambdaKf=1.0);
};

class ADiheRestraint : public ARestraint {
//-------------------------------------------------------------------
// ADiheRestraint is a derived class of ARestraint
//-------------------------------------------------------------------
public:
  ADiheRestraint();
  void    PrintInfo();
protected:
  double  GetE(double RefDihe, double Const);
  AVector GetGrad(int WhichGroup, double RefDihe, double Const);
  AVector gradU(AVector& P1P2P3, AVector& P4P5P6,
                AVector& dP1,    AVector& dP2,    AVector& dP3,
                AVector& dP4,    AVector& dP5,    AVector& dP6);
};


//****************************************************************************
//  AFixedPosRestraint, ABoundPosRestraint, AForcingPosRestraint:
//****************************************************************************
class AFixedPosRestraint : public APosRestraint {
//-------------------------------------------------------------------
// AFixedPosRestraint is a derived class of APosRestraint
//-------------------------------------------------------------------
private:
  AVector m_RefPos;
public:
  void    SetRefPos(AVector Pos) {m_RefPos=Pos;}
  AVector GetRefPos()            {return(m_RefPos);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Fixed   Position Restraint");
  }
};

class ABoundPosRestraint : public APosRestraint {
//-------------------------------------------------------------------
// ABoundPosRestraint is a derived class of APosRestraint
//-------------------------------------------------------------------
private:
  AVector m_RefPos;
  double  m_RefDist;
  Bound_t m_Bound;
public:
  void    SetRefPos(AVector Pos)  {m_RefPos=Pos;}
  void    SetRefDist(double Dist) {m_RefDist=Dist;}
  void    SetBound(Bound_t Bound) {m_Bound=Bound;}
  AVector GetRefPos()             {return(m_RefPos);}
  double  GetRefDist()            {return(m_RefDist);}
  Bound_t GetBound()              {return(m_Bound);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Bound   Position Restraint");
  }
};

class AForcingPosRestraint : public APosRestraint {
//-------------------------------------------------------------------
// AForcingPosRestraint is a derived class of APosRestraint
//-------------------------------------------------------------------
private:
  AVector m_StartPos;
  AVector m_StopPos;
public:
  void    SetStartPos(AVector Pos)       {m_StartPos=Pos;}
  void    SetStopPos(AVector Pos)        {m_StopPos=Pos;}
  AVector GetStartPos()                  {return(m_StartPos);}
  AVector GetStopPos()                   {return(m_StopPos);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  double  Get_dU_dLambda();
  Bool_t  IsForcing() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Forcing Position Restraint");
  }
};


//****************************************************************************
//  AFixedDistRestraint, ABoundDistRestraint, AForcingDistRestraint:
//****************************************************************************
class AFixedDistRestraint : public ADistRestraint {
//-------------------------------------------------------------------
// AFixedDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  double  m_RefDist;
public:
  void    SetRefDist(double Dist)  {m_RefDist=Dist;}
  double  GetRefDist()             {return(m_RefDist);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Fixed   Distance Restraint");
  }
};

class ABoundDistRestraint : public ADistRestraint {
//-------------------------------------------------------------------
// ABoundDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  double  m_RefDist;
  Bound_t m_Bound;
public:
  void    SetRefDist(double Dist)  {m_RefDist=Dist;}
  void    SetBound(Bound_t Bound)  {m_Bound=Bound;}
  double  GetRefDist()             {return(m_RefDist);}
  Bound_t GetBound()               {return(m_Bound);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Bound   Distance Restraint");
  }
};

class AForcingDistRestraint : public ADistRestraint {
//-------------------------------------------------------------------
// AForcingDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  double  m_StartDist;
  double  m_StopDist;
public:
  void    SetStartDist(double Dist)      {m_StartDist=Dist;}
  void    SetStopDist(double Dist)       {m_StopDist=Dist;}
  double  GetStartDist()                 {return(m_StartDist);}
  double  GetStopDist()                  {return(m_StopDist);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  double  Get_dU_dLambda();
  Bool_t  IsForcing() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Forcing Distance Restraint");
  }
};


//****************************************************************************
//  AFixedAngleRestraint, ABoundAngleRestraint, AForcingAngleRestraint:
//****************************************************************************
class AFixedAngleRestraint : public AnAngleRestraint {
//-------------------------------------------------------------------
// AFixedAngleRestraint is a derived class of AnAngleRestraint
//-------------------------------------------------------------------
private:
  double  m_RefAngle;     // in radians
public:
  void    SetRefAngle(double Angle)  {m_RefAngle=Angle;}
  double  GetRefAngle()              {return(m_RefAngle);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Fixed   Angle    Restraint");
  }
};

class ABoundAngleRestraint : public AnAngleRestraint {
//-------------------------------------------------------------------
// ABoundAngleRestraint is a derived class of AnAngleRestraint
//-------------------------------------------------------------------
private:
  double  m_RefAngle;     // in radians
  Bound_t m_Bound;
public:
  void    SetRefAngle(double Angle)  {m_RefAngle=Angle;}
  void    SetBound(Bound_t Bound)    {m_Bound=Bound;}
  double  GetRefAngle()              {return(m_RefAngle);}
  Bound_t GetBound()                 {return(m_Bound);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Bound   Angle    Restraint");
  }
};

class AForcingAngleRestraint : public AnAngleRestraint {
//-------------------------------------------------------------------
// AForcingAngleRestraint is a derived class of AnAngleRestraint
//-------------------------------------------------------------------
private:
  double  m_StartAngle;     // in radians
  double  m_StopAngle;      // in radians
public:
  void    SetStartAngle(double Angle)    {m_StartAngle=Angle;}
  void    SetStopAngle(double Angle)     {m_StopAngle=Angle;}
  double  GetStartAngle()                {return(m_StartAngle);}
  double  GetStopAngle()                 {return(m_StopAngle);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  double  Get_dU_dLambda();
  Bool_t  IsForcing() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Forcing Angle    Restraint");
  }
};


//****************************************************************************
//  AFixedDiheRestraint, ABoundDiheRestraint, AForcingDiheRestraint:
//****************************************************************************
class AFixedDiheRestraint : public ADiheRestraint {
//-------------------------------------------------------------------
// AFixedDiheRestraint is a derived class of ADiheRestraint
//-------------------------------------------------------------------
private:
  double  m_RefAngle;     // in radians
public:
  void    SetRefAngle(double Angle)  {m_RefAngle=Angle;}
  double  GetRefAngle()              {return(m_RefAngle);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Fixed   Dihedral Restraint");
  }
};

class ABoundDiheRestraint : public ADiheRestraint {
//-------------------------------------------------------------------
// ABoundDiheRestraint is a derived class of ADiheRestraint
//-------------------------------------------------------------------
private:
  double  m_LowerAngle;     // in radians (between 0 and 2pi)
  double  m_UpperAngle;     // in radians (between 0 and 2pi)
  double  m_IntervalAngle;  // in radians
public:
  void    SetLowerAngle(double Angle)     {m_LowerAngle=Angle;}
  void    SetUpperAngle(double Angle)     {m_UpperAngle=Angle;}
  void    SetIntervalAngle(double Angle)  {m_IntervalAngle=Angle;}
  double  GetLowerAngle()                 {return(m_LowerAngle);}
  double  GetUpperAngle()                 {return(m_UpperAngle);}
  double  GetIntervalAngle()              {return(m_IntervalAngle);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Bound   Dihedral Restraint");
  }
};

class AForcingDiheRestraint : public ADiheRestraint {
//-------------------------------------------------------------------
// AForcingDiheRestraint is a derived class of ADiheRestraint
//-------------------------------------------------------------------
private:
  double  m_StartAngle;     // in radians
  double  m_StopAngle;      // in radians
public:
  void    SetStartAngle(double Angle)    {m_StartAngle=Angle;}
  void    SetStopAngle(double Angle)     {m_StopAngle=Angle;}
  double  GetStartAngle()                {return(m_StartAngle);}
  double  GetStopAngle()                 {return(m_StopAngle);}
  double  GetEnergy();
  AVector GetGradient(int WhichGroup);
  double  Get_dU_dLambda();
  Bool_t  IsForcing() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Forcing Dihedral Restraint");
  }
};

#endif
