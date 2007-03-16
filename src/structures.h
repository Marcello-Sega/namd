/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "common.h"
#include <vector>
using namespace std;


// status elements, used for Atom status 
#define UnknownAtom      0x00
#define HydrogenAtom     0x01
#define OxygenAtom       0x02
#define HBDonorAtom      0x04
#define HBAcceptorAtom   0x08
#define HBAntecedentAtom 0x10
#define HBHydrogenAtom   0x20


typedef unsigned short Index;		//  Used for index into arrays
					//  or parameters

typedef struct atom_name_info
{
	char *resname;
	char *atomname;
	char *atomtype;
} AtomNameInfo;

typedef struct atom_constants
{
	Real mass;
	Real charge;
	Index vdw_type;
	int32 status;	         // flags telling about this atom
	int32 partner;             // connecting atom, for hydrogens
	int32 hydrogenList;	// index of atom in hydrogenGroup list
} Atom;

typedef struct bond
{
	int32 atom1;
	int32 atom2;
	Index bond_type;
} Bond;

typedef struct angle
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	Index angle_type;
} Angle;

typedef struct dihedral
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	Index dihedral_type;
} Dihedral;

typedef struct improper
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	Index improper_type;
} Improper;

typedef struct crossterm
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	int32 atom5;
	int32 atom6;
	int32 atom7;
	int32 atom8;
	Index crossterm_type;
} Crossterm;

class Exclusion
{
public:
	Exclusion(void) : modified(0) {;}
	Exclusion(int a1, int a2, int mod = 0) :
		atom1(a1), atom2(a2), modified(mod) {;}
	int32 atom1;
	int32 atom2;
	Index modified;
	int hash(void) const
	{
		return atom1 + atom2;
	}
	int operator==(const Exclusion &o) const
	{
		return atom1 == o.atom1 && atom2 == o.atom2;
	}
	int operator<(const Exclusion &o) const
	{
		return
		(
		  ( atom1 < o.atom1 ) ||
		  ( atom1 == o.atom1 && atom2 < o.atom2 )
		);
	}
};


class MOStream;
class MIStream;

enum TupleSigType {BOND=0, ANGLE, DIHEDRAL, IMPROPER, DONOR, ACCEPTOR, CROSSTERM, EXCLUSION};

class TupleSignature{
public:
    //This field indicates which class this tuple belongs to (bond or angle or ...)
    TupleSigType tupleType;
    int numOffset;
    //the relative offset from the atom index
    int *offset;

    //This type is determined by the Parameter object
    //for Exclusion, this indicates the modified field
    Index tupleParamType;

public:
    TupleSignature(){        
        offset = NULL;
    }

    TupleSignature(int n, TupleSigType t, Index paramType){
        tupleType = t;
        numOffset = n;
        offset = new int[n];
        tupleParamType = paramType;
    }
    TupleSignature(const TupleSignature& oneSig){
        tupleType = oneSig.tupleType;
        numOffset = oneSig.numOffset;
        offset = new int[numOffset];
        setOffsets(oneSig.offset);
        tupleParamType = oneSig.tupleParamType;
    }
    TupleSignature &operator=(const TupleSignature& oneSig){
        tupleType = oneSig.tupleType;
        numOffset = oneSig.numOffset;
        if(offset) delete [] offset;
        offset = new int[numOffset];
        setOffsets(oneSig.offset);
        tupleParamType = oneSig.tupleParamType;
        return *this;
    }
    ~TupleSignature(){
        if(offset) delete[] offset;
    }
    void setOffsets(int *offs){
        for(int i=0; i<numOffset; i++)
            offset[i] = offs[i];
	//sort the offset in increasing order
	//based on the input files, this offset is almost sorted increasingly
	//Therefore using insertion sort
	/*switch(numOffset){
	    case 1:
		break;
	    case 2:
		if(offset[0]>offset[1]){
		    int tmp = offset[0];
		    offset[0] = offset[1];
		    offset[1] = tmp;
		}
		break;
	    default: //insertion sort
		for(int ii=1; ii<numOffset; ii++){
		    int val = offset[ii];
		    for(int jj=ii-1; jj>=0; jj--){
			if(offset[jj]<val){
			    offset[jj+1] = offset[jj];
			    offset[jj] = val;
			}else
			    break;
		    }
		}
	}*/

    }
    void setEmpty(){
        delete [] offset;
        offset = NULL;
    }
    int isEmpty(){
        return offset==NULL;
    }
    void output(FILE *ofp){
        for(int i=0; i<numOffset; i++)
            fprintf(ofp, "%d ", offset[i]);
        fprintf(ofp, "| %d\n", tupleParamType);         
    }
    void pack(MOStream *msg);
    void unpack(MIStream *msg);
};

int operator==(const TupleSignature &sig1, const TupleSignature& sig2){
    if(sig1.tupleType!=sig2.tupleType)
        return 0;

    if(sig1.tupleParamType != sig2.tupleParamType)
        return 0;

    int equalCnt=0;
    int numOffset = sig1.numOffset;
    for(int i=0; i<numOffset; i++){
        equalCnt += (sig1.offset[i]==sig2.offset[i]);      
    }
    return equalCnt==numOffset;
}

//represents the signatures for atoms
class AtomSignature{
    //Missing crossterms (which is usually not presented in psf file)
    //Missing exclusions (which is very important!!!)
public:
    int bondCnt;
    int angleCnt;
    int dihedralCnt;
    int improperCnt;
    //TODO: support crossterm!!!!
    int crosstermCnt;

    TupleSignature *bondSigs;
    TupleSignature *angleSigs;
    TupleSignature *dihedralSigs;
    TupleSignature *improperSigs;
    //TODO: support crossterm!!!
    TupleSignature *crosstermSigs;

    AtomSignature(){
        bondCnt=angleCnt=dihedralCnt=improperCnt=crosstermCnt=0;
        bondSigs = NULL;
        angleSigs = NULL;
        dihedralSigs = NULL;
        improperSigs = NULL;
        crosstermSigs = NULL;
    }
    AtomSignature(const AtomSignature &sig){
        bondSigs = NULL;
        angleSigs = NULL;
        dihedralSigs = NULL;
        improperSigs = NULL;
        crosstermSigs = NULL;

        bondCnt = sig.bondCnt;
        if(bondCnt>0){            
            bondSigs = new TupleSignature[bondCnt];
            for(int i=0; i<bondCnt; i++)
                bondSigs[i] = sig.bondSigs[i];
        }

        angleCnt = sig.angleCnt;
        if(angleCnt>0){            
            angleSigs = new TupleSignature[angleCnt];
            for(int i=0; i<angleCnt; i++)
                angleSigs[i] = sig.angleSigs[i];
        }
        
        dihedralCnt = sig.dihedralCnt;
        if(dihedralCnt>0){            
            dihedralSigs = new TupleSignature[dihedralCnt];
            for(int i=0; i<dihedralCnt; i++)
                dihedralSigs[i] = sig.dihedralSigs[i];
        }
        
        improperCnt = sig.improperCnt;
        if(improperCnt>0){
            improperSigs = new TupleSignature[improperCnt];
            for(int i=0; i<improperCnt; i++)
                improperSigs[i] = sig.improperSigs[i];
        }      

        crosstermCnt = sig.crosstermCnt;
        if(crosstermCnt>0){
            crosstermSigs = new TupleSignature[crosstermCnt];
            for(int i=0; i<crosstermCnt; i++)
                crosstermSigs[i] = sig.crosstermSigs[i];
        }        
    }
    AtomSignature& operator=(const AtomSignature& sig){        
        bondCnt = sig.bondCnt;
        if(bondSigs) delete [] bondSigs;
        if(bondCnt>0){
            bondSigs = new TupleSignature[bondCnt];
            for(int i=0; i<bondCnt; i++)
                bondSigs[i] = sig.bondSigs[i];
        }else
            bondSigs = NULL;

        angleCnt = sig.angleCnt;
        if(angleSigs) delete [] angleSigs;
        if(angleCnt>0){
            angleSigs = new TupleSignature[angleCnt];
            for(int i=0; i<angleCnt; i++)
                angleSigs[i] = sig.angleSigs[i];
        }else
            angleSigs = NULL;
        
        dihedralCnt = sig.dihedralCnt;
        if(dihedralSigs) delete [] dihedralSigs;
        if(dihedralCnt>0){
            dihedralSigs = new TupleSignature[dihedralCnt];
            for(int i=0; i<dihedralCnt; i++)
                dihedralSigs[i] = sig.dihedralSigs[i];
        }else
            dihedralSigs = NULL;
        
        improperCnt = sig.improperCnt;
        if(improperSigs) delete [] improperSigs;
        if(improperCnt>0){
            improperSigs = new TupleSignature[improperCnt];
            for(int i=0; i<improperCnt; i++)
                improperSigs[i] = sig.improperSigs[i];
        }else
            improperSigs = NULL;      

        crosstermCnt = sig.crosstermCnt;
        if(crosstermSigs) delete [] crosstermSigs;
        if(crosstermCnt>0){
            crosstermSigs = new TupleSignature[crosstermCnt];
            for(int i=0; i<crosstermCnt; i++)
                crosstermSigs[i] = sig.crosstermSigs[i];
        }else
            crosstermSigs = NULL;

        return *this;
    }
    ~AtomSignature(){
        if(bondSigs) delete[] bondSigs;
        if(angleSigs) delete[] angleSigs;
        if(dihedralSigs) delete[] dihedralSigs;
        if(improperSigs) delete[] improperSigs;
        if(crosstermSigs) delete[] crosstermSigs;
    }    

    void removeEmptyTupleSigs();
    void pack(MOStream *msg);
    void unpack(MIStream *msg);
};

int operator==(const AtomSignature &sig1, const AtomSignature& sig2){
    if(sig1.bondCnt!=sig2.bondCnt) return 0;
    if(sig1.angleCnt!=sig2.angleCnt) return 0;
    if(sig1.dihedralCnt!=sig2.dihedralCnt) return 0;
    if(sig1.improperCnt!=sig2.improperCnt) return 0;
    if(sig1.crosstermCnt!=sig2.crosstermCnt) return 0;

#define CMPSIGS(TUPLE) \
for(int i=0; i<sig1.TUPLE##Cnt; i++){ \
    if(!(sig1.TUPLE##Sigs[i]==sig2.TUPLE##Sigs[i])) return 0; \
} \

    CMPSIGS(bond)
    CMPSIGS(angle)
    CMPSIGS(dihedral)
    CMPSIGS(improper)
    CMPSIGS(crossterm)
}

struct AtomNameIdx{
    Index resnameIdx;
    Index atomnameIdx;
    Index atomtypeIdx;
};
struct AtomCstInfo{
    Index vdw_type;
    int32 status;
    int32 partner;
    int32 hydrogenList;
};

struct ExclusionSignature{
    int fullExclCnt; //1-2, 1-3 exclusion
    int *fullOffset; //should be in increasing order
    int modExclCnt; //1-4 exclusion
    int *modOffset; //should be in increasing order

    ExclusionSignature(){
    	fullExclCnt = modExclCnt = 0;
    	fullOffset = modOffset = NULL;
    }    
    ExclusionSignature(const ExclusionSignature& sig){
        fullOffset = modOffset = NULL;
    	fullExclCnt = sig.fullExclCnt;
        if(fullExclCnt>0){
            fullOffset = new int[fullExclCnt];
            for(int i=0; i<fullExclCnt; i++)
                fullOffset[i] = sig.fullOffset[i];
        }
    	
    	modExclCnt = sig.modExclCnt;
        if(modExclCnt>0){
            modOffset = new int[modExclCnt];
            for(int i=0; i<modExclCnt; i++)
                modOffset[i] = sig.modOffset[i];
        }
    }
    ~ExclusionSignature(){
    	if(fullOffset) delete [] fullOffset;
    	if(modOffset) delete [] modOffset;
    }
    
    ExclusionSignature& operator=(const ExclusionSignature& sig){
        fullExclCnt = sig.fullExclCnt;
        if(fullOffset) delete [] fullOffset;
        if(fullExclCnt>0){
            fullOffset = new int[fullExclCnt];
            for(int i=0; i<fullExclCnt; i++)
                fullOffset[i] = sig.fullOffset[i];
        }else
            fullOffset = NULL;
    
    	modExclCnt = sig.modExclCnt;
        if(modOffset) delete [] modOffset;
        if(modExclCnt>0){
            modOffset = new int[modExclCnt];
            for(int i=0; i<modExclCnt; i++)
                modOffset[i] = sig.modOffset[i];
        }else
            modOffset = NULL;

        return *this;
    }

    //both input should be sorted in increasing order
    void setOffsets(vector<int>& fullVec, vector<int>& modVec){
    	fullExclCnt = fullVec.size();
    	modExclCnt = modVec.size();
    	if(fullExclCnt>0) {        
            fullOffset = new int[fullExclCnt];
    	    for(int i=0; i<fullExclCnt; i++)
                fullOffset[i] = fullVec[i];
        }

    	if(modExclCnt>0) {        
            modOffset = new int[modExclCnt];
    	    for(int i=0; i<modExclCnt; i++)
                modOffset[i] = modVec[i];	
        }
    }

    void removeEmptyOffset();
    //assuming offsets in the signature have been sorted increasingly
    int findOffset(int offset, int  *fullOrMod);
    void pack(MOStream *msg);
    void unpack(MIStream *msg);
};
int operator==(const ExclusionSignature& sig1, const ExclusionSignature& sig2){
    if(sig1.fullExclCnt!=sig2.fullExclCnt) return 0;
    if(sig1.modExclCnt!=sig2.modExclCnt) return 0;
    
    for(int i=0; i<sig1.fullExclCnt; i++){
	if(sig1.fullOffset[i]!=sig2.fullOffset[i]) return 0;
    }
    for(int i=0; i<sig1.modExclCnt; i++){
	if(sig1.modOffset[i]!=sig2.modOffset[i]) return 0;
    }
    return 1;
}

#endif

