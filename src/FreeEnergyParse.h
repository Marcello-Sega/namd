// written by David Hurwitz, March to May 1998.

#if !defined(PARSE_HPP)
  #define PARSE_HPP

enum item_t  {kAtom, kResidue, kAtomName, kAtomNameList,
              kResidueRange, kStartGroup, kEndGroup, kUnknownItem};
enum restr_t {kPosi,      kDist,      kAngle,      kDihe,
              kPosiBound, kDistBound, kAngleBound, kDiheBound,
              kPosiPMF,   kDistPMF,   kAnglePMF,   kDihePMF, kUnknownRestr};
enum pmf_t   {kTask, kTime, kLambda, kLambdaT, kPrint, kNoPrint,
              kEquilTime, kAccumTime, kNumRepeats, kUnknownPmf};
enum TimeUnits_t  {k_fs, k_ps, k_ns, kUnknownTime};

// the "main" function for parsing
void   ReadInput(char* Str, ARestraintManager& RMgr,
                            ALambdaManager& LMgr,
                            ComputeFreeEnergy& CFE,
                            double dT);

// get an initialized restraint that's read from the input String
int    ReadRestraints(char* Str, ARestraintManager& RMgr,
                                 ComputeFreeEnergy& CFE);
ARestraint* GetRestraint(char* Str, int& NumChars, ComputeFreeEnergy& CFE);

// for reading pmf/mcti blocks
int    ReadPmfBlock(char* Str, ALambdaControl& PmfBlock, double dT);
int    ReadNextPmfSpec(char* Str, pmf_t& PmfSpec);
int    ReadTaskType(char* Str, task_t& Task);
int    ReadTimeUnits(char* Str, TimeUnits_t& Units, TimeUnits_t DefaultUnits);
double GetTime(double Val, TimeUnits_t Units);

// functions for parsing the config file
void    CheckParentheses(char* Str);
void    ProblemParsing(char* Message, char* Str, Bool_t Terminate=kTrue);
void    ToLower(char* Str);
int     ReadWhite(char* Str);
int     ReadAlphaNum(char* Str);
int     ReadAlpha(char* Str);
int     ReadDigits(char* Str);
int     ReadParentheses(char* Str);
int     IsStartGroup(char* Str);
int     IsEndGroup(char* Str);
int     IsAtomName(char* Str);
int     IsAtomNameList(char* Str);
int     IsAtom(char* Str);
int     IsAResidue(char* Str);
int     IsResidue(char* Str);
int     IsResidueRange(char* Str);
int     ReadWord(char* Str, char* Word, Bool_t ErrMsg=kFalse);
int     ReadChar(char* Str, char Char, Bool_t ErrMsg=kFalse);
int     ReadAValue(char* Str, double& Value, Bool_t ErrMsg=kFalse);
int     ReadBound(char* Str, Bound_t& Bound);
item_t  ReadNextItem(char* Str, int& NumChars);
restr_t ReadNextRestraintType(char* Str, int& NumChars);

// functions for putting the specified atoms into Group's
int  AddAtoms(AGroup& Group, char* Str, ComputeFreeEnergy& CFE);
void AddAtom(AGroup& Group, char* Atom, ComputeFreeEnergy& CFE);
void AddResidues(AGroup& Group, char* ResRange, ComputeFreeEnergy& CFE);
void AddAtomsInResidues(AGroup& Group, char* AtomNames, char* ResRange,
                        ComputeFreeEnergy& CFE);

// functions for helping with the above
void AddAtom(AGroup& Group, char* ResRange, char* AtomName,
             ComputeFreeEnergy& CFE);
void GetResRange(char* ResRange, int& ResNum1, int& ResNum2);
int  GetSegName(char* Str, char* SegName);
int  GetResNum(char* Str, int& ResNum);
int  GetAtomName(char* Str, char* AtomName);

#endif
