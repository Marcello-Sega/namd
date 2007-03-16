#ifndef ATOMSDISINFO_H
#define ATOMSDISINFO_H

#include "structures.h"
#include "AtomsDisInfo.decl.h"

class AtomStaticInfoMsg: public CMessage_AtomStaticInfoMsg{
    public:
        int actualNumAtoms;
        Atom *atoms;
};

class AtomsCoorMsg: public CMessage_AtomsCoorMsg{
    public:
        int x;
        //Vector coors
};

class AtomsForcesMsg: public CMessage_AtomsForcesMsg{
    public:
        int x;
        //Vector forces
};


class AtomsDisInfo : public CBase_AtomsDisInfo
{
public:
    static const int ATOMDISNUM = 500;

    /* 
     * Blindly returns the array index of the atom (specified by atomId).
     * In reality, the atomId should be less than the total number of atoms 
     */ 
    static int getAtomArrayIndex(int atomId){
        if(atomId<0){
            return -1;
        }

        return atomId/ATOMDISNUM;
    }

/* TODO: Need to figure out what data structures are needed later */

    //AtomInfo's array (atoms' static info, mass/charge etc) with size ATOMDISNUM
    int atomsCnt;
    Atom *disAtoms;
    
    //Atom's force array 

    //Atom's coordinates' array (Points[ATOMNUM])

public: //methods
    AtomsDisInfo(){
	//initialize the atoms' static info array/force array/coordinates array
        disAtoms = new Atom[ATOMDISNUM];
    }

    AtomsDisInfo(CkMigrateMessage *m){}    

    ~AtomsDisInfo(){
	//free allocated memory for static info/force array/coor array
        delete disAtoms;
    }

    void recvStaticInfo(AtomStaticInfoMsg *msg);
    void recvAtomsCoor(AtomsCoorMsg *msg);
    void recvAtomsForces(AtomsForcesMsg *msg);
};

#include "AtomsDisInfo.def.h"

#endif
