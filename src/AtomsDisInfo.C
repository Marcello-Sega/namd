#include "AtomsDisInfo.h"

void AtomsDisInfo::recvStaticInfo(AtomStaticInfoMsg *msg){    
    atomsCnt = msg->actualNumAtoms;
    Atom *msgAtoms = msg->atoms;
    for(int i=0; i<atomsCnt; i++){
        disAtoms[i] = msgAtoms[i];
    }
    delete msg;
}

void AtomsDisInfo::recvAtomsCoor(AtomsCoorMsg *msg){
}

void AtomsDisInfo::recvAtomsForces(AtomsForcesMsg *msg){
}
