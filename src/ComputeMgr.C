/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ProcessorPrivate.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

#include "BOCgroup.h"
#include "ComputeMgr.decl.h"
#include "ComputeMgr.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"

#include "Node.h"
#include "ComputeMap.h"
#include "PatchMap.h"
#include "PatchMap.inl"

#include "Compute.h"
#include "ComputeNonbondedUtil.h"
#include "ComputeNonbondedSelf.h"
#include "ComputeNonbondedPair.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeCrossterms.h"
#include "ComputeBonds.h"
#include "ComputeFullDirect.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "ComputeExt.h"
#include "ComputeDPMTA.h"
#include "ComputeDPME.h"
#include "ComputeDPMEMsgs.h"
#include "ComputePme.h"
#include "ComputeEwald.h"
#include "ComputeEField.h"
/* BEGIN gf */
#include "ComputeGridForce.h"
/* END gf */
#include "ComputeStir.h"
#include "ComputeSphericalBC.h"
#include "ComputeCylindricalBC.h"
#include "ComputeTclBC.h"
#include "ComputeRestraints.h"
#include "ComputeConsForce.h"
#include "ComputeConsForceMsgs.h"
#include "WorkDistrib.h"

/* include all of the specific masters we need here */
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"
#include "FreeEnergyVector.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"

#include "GlobalMasterTest.h"
#include "GlobalMasterIMD.h"
#include "GlobalMasterTcl.h"
#include "GlobalMasterSMD.h"
#include "GlobalMasterTMD.h"
#include "GlobalMasterEasy.h"
#include "GlobalMasterMisc.h"
#include "GlobalMasterFreeEnergy.h"

ComputeMgr::ComputeMgr()
{
    CkpvAccess(BOCclass_group).computeMgr = thisgroup;
    computeGlobalObject = 0;
    computeDPMEObject = 0;
    computeEwaldObject = 0;
    computeNonbondedWorkArrays = new ComputeNonbondedWorkArrays;
}

ComputeMgr::~ComputeMgr(void)
{
    delete computeNonbondedWorkArrays;
}

void ComputeMgr::updateComputes(int ep, CkGroupID chareID)
{
    updateComputesReturnEP = ep;
    updateComputesReturnChareID = chareID;
    updateComputesCount = CkNumPes();

    if (CkMyPe())
    {
        iout << iPE << iERRORF << "updateComputes signaled on wrong Pe!\n" << endi;
        CkExit();
        return;
    }

#if CHARM_VERSION > 050402
    CkStartQD(CkIndex_ComputeMgr::updateComputes2((CkQdMsg*)0),&thishandle);
#else
    CkStartQD(CProxy_ComputeMgr::ckIdx_updateComputes2((CkQdMsg*)0),&thishandle);
#endif
}

void ComputeMgr::updateComputes2(CkQdMsg *msg)
{
    delete msg;

    CProxy_WorkDistrib wd(CkpvAccess(BOCclass_group).workDistrib);
    WorkDistrib  *workDistrib = wd.ckLocalBranch();
#if CHARM_VERSION > 050402
    workDistrib->saveComputeMapChanges(CkIndex_ComputeMgr::updateComputes3(),thisgroup);
#else
    workDistrib->saveComputeMapChanges(CProxy_ComputeMgr::ckIdx_updateComputes3(),thisgroup);
#endif
}

void ComputeMgr::updateComputes3()
{
    CProxy_ComputeMgr(thisgroup).updateLocalComputes();
}

void ComputeMgr::updateLocalComputes()
{
    ComputeMap *computeMap = ComputeMap::Object();
    CProxy_ProxyMgr pm(CkpvAccess(BOCclass_group).proxyMgr);
    ProxyMgr *proxyMgr = pm.ckLocalBranch();

    computeFlag = new int[computeMap->numComputes()];

    for (int i=0; i<computeMap->numComputes(); i++)
    {
        DebugM(3, "updateLocalComputes("<<i<<") curnode="<<computeMap->node(i)
               <<" newnode="<<computeMap->newNode(i)<<"\n");
        computeFlag[i] = 0;

        if (computeMap->newNode(i) == CkMyPe() && computeMap->node(i) != CkMyPe())
        {
            DebugM(4, "updateLocal - creating new computeID("<<i<<")\n");
            computeFlag[i] = 1;
            computeMap->setNode(i,computeMap->newNode(i));
            for (int n=0; n < computeMap->numPids(i); n++)
            {
                proxyMgr->createProxy(computeMap->pid(i,n));
            }
        }
        else if (computeMap->node(i) == CkMyPe() &&
                 (computeMap->newNode(i) != -1 && computeMap->newNode(i) != CkMyPe() ))
        {
            DebugM(4, "updateLocal - deleting computeID("<<i<<")\n");
            computeFlag[i] = -1;
            computeMap->setNode(i,computeMap->newNode(i));
        }
        else if (computeMap->newNode(i) != -1)
        {
            computeMap->setNode(i,computeMap->newNode(i));
        }
        computeMap->setNewNode(i,-1);
    }

    if (!CkMyPe())
    {
#if CHARM_VERSION > 050402
        CkStartQD(CkIndex_ComputeMgr::updateLocalComputes2((CkQdMsg*)0), &thishandle);
#else
        CkStartQD(CProxy_ComputeMgr::ckIdx_updateLocalComputes2((CkQdMsg*)0), &thishandle);
#endif
    }
}

void
ComputeMgr::updateLocalComputes2(CkQdMsg *msg)
{
    delete msg;
    CProxy_ComputeMgr(thisgroup).updateLocalComputes3();
}

void
ComputeMgr::updateLocalComputes3()
{
    ComputeMap *computeMap = ComputeMap::Object();
    CProxy_ProxyMgr pm(CkpvAccess(BOCclass_group).proxyMgr);
    ProxyMgr *proxyMgr = pm.ckLocalBranch();

    ProxyMgr::nodecount = 0;

    for (int i=0; i<computeMap->numComputes(); i++)
    {
        if (1 == computeFlag[i])
        {
            DebugM(4, "updateLocalCompute3() - create computeID(" << i << ")\n");
            createCompute(i, computeMap);
        }
        else if (-1 == computeFlag[i])
        {
            // remove this compute
            DebugM(4, "updateLocalCompute3() - delete computeID(" << i << ")\n");
            delete computeMap->compute(i);
            computeMap->registerCompute(i,NULL);
        }
    }
    delete[] computeFlag;

    proxyMgr->removeUnusedProxies();

    DebugM(4, "msg to doneUpdateLocalComputes on Pe("<<CkMyPe()<<")\n");
    ComputeMap::Object()->checkMap();
    PatchMap::Object()->checkMap();

    if (!CkMyPe())
    {
#if CHARM_VERSION > 050402
        CkStartQD(CkIndex_ComputeMgr::updateLocalComputes4((CkQdMsg*)0), &thishandle);
#else
        CkStartQD(CProxy_ComputeMgr::ckIdx_updateLocalComputes4((CkQdMsg*)0), &thishandle);
#endif
// added a new phase to build spanning tree after load balance
// was
//    CkStartQD(CProxy_ComputeMgr::ckIdx_doneUpdateLocalComputes(), &thishandle);
    }
    //CSendMsgBranch(ComputeMgr, doneUpdateLocalComputes, thisgroup, 0);
}

void
ComputeMgr::updateLocalComputes4(CkQdMsg *msg)
{
    delete msg;
    CProxy_ComputeMgr(thisgroup).updateLocalComputes5();
}

int firstphase = 1;

void
ComputeMgr::updateLocalComputes5()
{
    // we always use the centralized building of spanning tree
    // distributed building of ST called in Node.C only
    if (proxySendSpanning || proxyRecvSpanning)
        ProxyMgr::Object()->buildProxySpanningTree2();

    // this code needs to be turned on if we want to
    // shift the creation of ST to the load balancer

#if 0
    if (proxySendSpanning || proxyRecvSpanning)
    {
        if (firstphase)
            ProxyMgr::Object()->buildProxySpanningTree2();
        else
            if (CkMyPe() == 0)
                ProxyMgr::Object()->sendSpanningTrees();

        firstphase = 0;
    }
#endif

    if (!CkMyPe())
#if CHARM_VERSION > 050402
        CkStartQD(CkIndex_ComputeMgr::doneUpdateLocalComputes(), &thishandle);
#else
        CkStartQD(CProxy_ComputeMgr::ckIdx_doneUpdateLocalComputes(), &thishandle);
#endif
}

void ComputeMgr::doneUpdateLocalComputes()
{

//  if (!--updateComputesCount) {
    DebugM(4, "doneUpdateLocalComputes on Pe("<<CkMyPe()<<")\n");
    void *msg = CkAllocMsg(0,0,0);
    CkSendMsgBranch(updateComputesReturnEP,msg,0,updateComputesReturnChareID);
//  }
}

//
void
ComputeMgr::createCompute(ComputeID i, ComputeMap *map)
{
    Compute *c;
    PatchID pid2[2];
    PatchIDList pids;
    int trans2[2];
    SimParameters *simParams = Node::Object()->simParameters;

    switch ( map->type(i) )
    {
    case computeNonbondedSelfType:
        c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0].pid,
                                     computeNonbondedWorkArrays,
                                     map->partition(i),map->partition(i)+1,
                                     map->numPartitions(i)); // unknown delete
        ++numNonbondedSelf;
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeNonbondedPairType:
        pid2[0] = map->computeData[i].pids[0].pid;
        trans2[0] = map->computeData[i].pids[0].trans;
        pid2[1] = map->computeData[i].pids[1].pid;
        trans2[1] = map->computeData[i].pids[1].trans;
        c = new ComputeNonbondedPair(i,pid2,trans2,
                                     computeNonbondedWorkArrays,
                                     map->partition(i),map->partition(i)+1,
                                     map->numPartitions(i)); // unknown delete
        ++numNonbondedPair;
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeBondsType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeBonds(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeAnglesType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeAngles(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeDihedralsType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeDihedrals(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeImpropersType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeImpropers(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeCrosstermsType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeCrossterms(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfBondsType:
        c = new ComputeSelfBonds(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfAnglesType:
        c = new ComputeSelfAngles(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfDihedralsType:
        c = new ComputeSelfDihedrals(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfImpropersType:
        c = new ComputeSelfImpropers(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfCrosstermsType:
        c = new ComputeSelfCrossterms(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
#ifdef DPMTA
    case computeDPMTAType:
        c = new ComputeDPMTA(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
#endif
#ifdef DPME
    case computeDPMEType:
        c = computeDPMEObject = new ComputeDPME(i,this); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
#endif
    case computePmeType:
        c = new ComputePme(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeEwaldType:
        c = computeEwaldObject = new ComputeEwald(i,this); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeFullDirectType:
        c = new ComputeFullDirect(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeGlobalType:
        c = computeGlobalObject = new ComputeGlobal(i,this); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeStirType:
        c = new ComputeStir(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeExtType:
        c = new ComputeExt(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeEFieldType:
        c = new ComputeEField(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
        /* BEGIN gf */
    case computeGridForceType:
        c = new ComputeGridForce(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
        /* END gf */
    case computeSphericalBCType:
        c = new ComputeSphericalBC(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeCylindricalBCType:
        c = new ComputeCylindricalBC(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeTclBCType:
        c = new ComputeTclBC(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeRestraintsType:
        c = new ComputeRestraints(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeConsForceType:
        c = new ComputeConsForce(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeConsTorqueType:
        c = new ComputeConsTorque(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    default:
        NAMD_bug("Unknown compute type in ComputeMgr::createCompute().");
        break;
    }
}

void registerUserEventsForAllComputeObjs()
{
#ifdef TRACE_COMPUTE_OBJECTS
    ComputeMap *map = ComputeMap::Object();
    char user_des[50];
    for (int i=0; i<map->numComputes(); i++)
    {
        memset(user_des, 0, 50);
        switch ( map->type(i) )
        {
        case computeNonbondedSelfType:
            sprintf(user_des, "computeNonBondedSelfType_%d", i);
            break;
        case computeNonbondedPairType:
            sprintf(user_des, "computeNonBondedPairType_%d", i);
            break;
        case computeBondsType:
            sprintf(user_des, "computeBondsType_%d", i);
            break;
        case computeAnglesType:
            sprintf(user_des, "computeAnglesType_%d", i);
            break;
        case computeDihedralsType:
            sprintf(user_des, "computeDihedralsType_%d", i);
            break;
        case computeImpropersType:
            sprintf(user_des, "computeImpropersType_%d", i);
            break;
        case computeCrosstermsType:
            sprintf(user_des, "computeCrosstermsType_%d", i);
            break;
        case computeSelfBondsType:
            sprintf(user_des, "computeSelfBondsType_%d", i);
            break;
        case computeSelfAnglesType:
            sprintf(user_des, "computeSelfAnglesType_%d", i);
            break;
        case computeSelfDihedralsType:
            sprintf(user_des, "computeSelfDihedralsType_%d", i);
            break;
        case computeSelfImpropersType:
            sprintf(user_des, "computeSelfImpropersType_%d", i);
            break;
        case computeSelfCrosstermsType:
            sprintf(user_des, "computeSelfCrosstermsType_%d", i);
            break;
#ifdef DPMTA
        case computeDPMTAType:
            sprintf(user_des, "computeDPMTAType_%d", i);
            break;
#endif
#ifdef DPME
        case computeDPMEType:
            sprintf(user_des, "computeDPMEType_%d", i);
            break;
#endif
        case computePmeType:
            sprintf(user_des, "computePMEType_%d", i);
            break;
        case computeEwaldType:
            sprintf(user_des, "computeEwaldType_%d", i);
            break;
        case computeFullDirectType:
            sprintf(user_des, "computeFullDirectType_%d", i);
            break;
        case computeGlobalType:
            sprintf(user_des, "computeGlobalType_%d", i);
            break;
        case computeStirType:
            sprintf(user_des, "computeStirType_%d", i);
            break;
        case computeExtType:
            sprintf(user_des, "computeExtType_%d", i);
            break;
        case computeEFieldType:
            sprintf(user_des, "computeEFieldType_%d", i);
            break;
            /* BEGIN gf */
        case computeGridForceType:
            sprintf(user_des, "computeGridForceType_%d", i);
            break;
            /* END gf */
        case computeSphericalBCType:
            sprintf(user_des, "computeSphericalBCType_%d", i);
            break;
        case computeCylindricalBCType:
            sprintf(user_des, "computeCylindricalBCType_%d", i);
            break;
        case computeTclBCType:
            sprintf(user_des, "computeTclBCType_%d", i);
            break;
        case computeRestraintsType:
            sprintf(user_des, "computeRestraintsType_%d", i);
            break;
        case computeConsForceType:
            sprintf(user_des, "computeConsForceType_%d", i);
            break;
        case computeConsTorqueType:
            sprintf(user_des, "computeConsTorqueType_%d", i);
            break;
        default:
            NAMD_bug("Unknown compute type in ComputeMgr::registerUserEventForAllComputeObjs().");
            break;
        }
	int user_des_len = strlen(user_des);
	char *user_des_cst = new char[user_des_len+1];
	memcpy(user_des_cst, user_des, user_des_len);
	user_des_cst[user_des_len] = 0;
	//Since the argument in traceRegisterUserEvent is supposed
	//to be a const string which will not be copied inside the
	//function when a new user event is created, user_des_cst 
	//has to be allocated in heap.
        int reEvenId = traceRegisterUserEvent(user_des_cst, TRACE_COMPOBJ_IDOFFSET+i);
	//printf("Register user event (%s) with id (%d)\n", user_des, reEvenId);
    }
#else
    return;
#endif
}

void
ComputeMgr::createComputes(ComputeMap *map)
{
    Node *node = Node::Object();
    SimParameters *simParams = node->simParameters;
    int myNode = node->myid();

    numNonbondedSelf = 0;
    numNonbondedPair = 0;
    ComputeNonbondedUtil::select();

    if ( simParams->globalForcesOn && !myNode )
    {
        DebugM(4,"Mgr running on Node "<<CkMyPe()<<"\n");
        /* create a master server to allow multiple masters */
        masterServerObject = new GlobalMasterServer(this,
                PatchMap::Object()->numNodesWithPatches());

        /* create the individual global masters */
        // masterServerObject->addClient(new GlobalMasterTest());
        if (simParams->tclForcesOn)
            masterServerObject->addClient(new GlobalMasterTcl());
        if (simParams->IMDon)
            masterServerObject->addClient(new GlobalMasterIMD());

        if (simParams->SMDOn)
            masterServerObject->addClient(
                new GlobalMasterSMD(simParams->SMDk, simParams->SMDk2,
				    simParams->SMDVel,
                                    simParams->SMDDir, simParams->SMDOutputFreq,
                                    simParams->firstTimestep, simParams->SMDFile,
                                    node->molecule->numAtoms)
            );
        if (simParams->TMDOn)
            masterServerObject->addClient(new GlobalMasterTMD());
        if (simParams->miscForcesOn)
            masterServerObject->addClient(new GlobalMasterMisc());
        if ( simParams->freeEnergyOn )
            masterServerObject->addClient(new GlobalMasterFreeEnergy());
    }

    for (int i=0; i < map->nComputes; i++)
    {
        if ( ! ( i % 100 ) )
        {
        }
        if ( map->computeData[i].node != myNode ) continue;
        DebugM(1,"Compute " << i << '\n');
        DebugM(1,"  node = " << map->computeData[i].node << '\n');
        DebugM(1,"  type = " << map->computeData[i].type << '\n');
        DebugM(1,"  patchBased = " << map->computeData[i].patchBased << '\n');
        DebugM(1,"  numPids = " << map->computeData[i].numPids << '\n');
        DebugM(1,"  numPidsAllocated = " << map->computeData[i].numPidsAllocated << '\n');
        for (int j=0; j < map->computeData[i].numPids; j++)
        {
            //      DebugM(1,"  pid " << map->computeData[i].pids[j] << '\n');
            if (!((j+1) % 6))
                DebugM(1,'\n');
        }
        DebugM(1,"\n---------------------------------------");
        DebugM(1,"---------------------------------------\n");

        createCompute(i, map);

    }

}


void ComputeMgr:: sendComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
    (CProxy_ComputeMgr(CkpvAccess(BOCclass_group).computeMgr)).recvComputeGlobalConfig(msg);
}

void ComputeMgr:: recvComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
    if ( computeGlobalObject )
    {
        computeGlobalObject->recvConfig(msg);
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalData(ComputeGlobalDataMsg *msg)
{
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
#if CHARM_VERSION > 050402
    cm[0].recvComputeGlobalData(msg);
#else
    cm.recvComputeGlobalData(msg, 0);
#endif
}

void ComputeMgr:: recvComputeGlobalData(ComputeGlobalDataMsg *msg)
{
    if (masterServerObject)  // make sure it has been initialized
    {
        masterServerObject->recvData(msg);
    }
    else NAMD_die("ComputeMgr::masterServerObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
    (CProxy_ComputeMgr(CkpvAccess(BOCclass_group).computeMgr)).recvComputeGlobalResults(msg);
}

void ComputeMgr:: recvComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
    if ( computeGlobalObject )
    {
        computeGlobalObject->recvResults(msg);
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

/*
 * Begin Ewald messages
 */
void ComputeMgr:: sendComputeEwaldData(ComputeEwaldMsg *msg)
{
    if (computeEwaldObject)
    {
        int node = computeEwaldObject->getMasterNode();
        CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
#if CHARM_VERSION > 050402
        cm[node].recvComputeEwaldData(msg);
#else
        cm.recvComputeEwaldData(msg, node);
#endif
    }
    else if (!PatchMap::Object()->numHomePatches())
    {
        CkPrintf("skipping message on Pe(%d)\n", CkMyPe());
        delete msg;
    }
    else NAMD_die("ComputeMgr::computeEwaldObject is NULL!");
}

void ComputeMgr:: recvComputeEwaldData(ComputeEwaldMsg *msg)
{
    if (computeEwaldObject)
        computeEwaldObject->recvData(msg);
    else NAMD_die("ComputeMgr::computeEwaldObject in recvData is NULL!");
}

void ComputeMgr:: sendComputeEwaldResults(ComputeEwaldMsg *msg)
{
    (CProxy_ComputeMgr(CkpvAccess(BOCclass_group).computeMgr)).recvComputeEwaldResults(msg);
}

void ComputeMgr::recvComputeEwaldResults(ComputeEwaldMsg *msg)
{
    if (computeEwaldObject)
        computeEwaldObject->recvResults(msg);
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeEwaldObject in recvResults is NULL!");
}

void ComputeMgr:: sendComputeDPMEData(ComputeDPMEDataMsg *msg)
{
    if ( computeDPMEObject )
    {
#ifdef DPME
        int node = computeDPMEObject->getMasterNode();
        CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
        cm.recvComputeDPMEData(msg,node);
#endif
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: recvComputeDPMEData(ComputeDPMEDataMsg *msg)
{
    if ( computeDPMEObject )
    {
#ifdef DPME
        computeDPMEObject->recvData(msg);
#endif
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: sendComputeDPMEResults(ComputeDPMEResultsMsg *msg, int node)
{
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
#if CHARM_VERSION > 050402
    cm[node].recvComputeDPMEResults(msg);
#else
    cm.recvComputeDPMEResults(msg, node);
#endif
}

void ComputeMgr:: recvComputeDPMEResults(ComputeDPMEResultsMsg *msg)
{
    if ( computeDPMEObject )
    {
#ifdef DPME
        computeDPMEObject->recvResults(msg);
#endif
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr::recvComputeConsForceMsg(ComputeConsForceMsg *msg)
{
    Molecule *m = Node::Object()->molecule;
    delete [] m->consForceIndexes;
    delete [] m->consForce;
    int n = msg->aid.size();
    if (n > 0)
    {
        m->consForceIndexes = new int32[m->numAtoms];
        m->consForce = new Vector[n];
        int i;
        for (i=0; i<m->numAtoms; i++) m->consForceIndexes[i] = -1;
        for (i=0; i<msg->aid.size(); i++)
        {
            m->consForceIndexes[msg->aid[i]] = i;
            m->consForce[i] = msg->f[i];
        }
    }
    else
    {
        m->consForceIndexes = NULL;
        m->consForce = NULL;
    }
    delete msg;
}

#include "ComputeMgr.def.h"

