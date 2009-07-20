/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEGRIDFORCE_H
#define COMPUTEGRIDFORCE_H

#include "ComputeHomePatches.h"
#include "ReductionMgr.h"
#include "GridForceGrid.h"
#include "ComputeGridForceMgr.decl.h"

#define GRIDSEGSZ 8000000

struct BoxState {
  Box<Patch, CompAtom> *posBox;
  Box<Patch, Results> *forceBox;
  CompAtom *atoms;
  Results *results;
  HomePatch *homePatch;
};

class ComputeGridForce : public ComputeHomePatches
{

public:
	ComputeGridForce(ComputeID c); 	//  Constructor
	virtual ~ComputeGridForce();			//  Destructor

        void doWork();
	void doForce(HomePatch *hp,FullAtom* p);
        void finishWork();
        void finishForce(HomePatch* homePatch,
                         FullAtom* p,
                         Results* r);
//        static int getNumGridIndices();

	SubmitReduction *reduction;
private:
        BoxState *box_state;
        int num_boxes;
};

struct GridVal {
  int index;
//  float val[GRIDBLOCKX*GRIDBLOCKY*GRIDBLOCKZ];
  float val;
};

class GridDepositMsg : public CMessage_GridDepositMsg {
public:
  int gridnum;
  GridforceGrid* grid;
  int num_grids;
};

class SubReqMsg : public CMessage_SubReqMsg {
public:
  int rank;
  ComputeGridForce* gf;
};

class GridSegmentMsg : public CMessage_GridSegmentMsg {
public:
  int gridnum;
  int count;
  int start_index;
  float grid[GRIDSEGSZ];
};

class GridRequestMsg : public CMessage_GridRequestMsg {
public:
  int from_node;
  int num_grid_indices;
  int *gridStartIndex;
  int *gridIndexList;
};

class GridValuesMsg : public CMessage_GridValuesMsg {
public:
  int num_grid_indices;
  int *gridStartIndex;
  GridVal *gridVals;
};

class ComputeGridForceMgr : public BOCclass {
public:
  ComputeGridForceMgr();
  void request(ComputeGridForce *gf);
  void finishWork();
private:
  ComputeGridForce *gf;
};

class ComputeGridForceNodeMgr : public NodeGroup {
public:
  ComputeGridForceNodeMgr();
  void depositInitialGrid(GridDepositMsg *msg);
  void requestInitialGridData();
  void receiveInitialGridData(GridSegmentMsg *msg);
  void submitRequest(SubReqMsg *msg);
  void fetchGridValues(GridRequestMsg *msg);
  void recvGridValues(GridValuesMsg *msg);
  void resume();
  
  private:
    int myNode;
    int nodeSize;
    int num_grids;
    int grids_deposited;
    GridforceGrid **grids;
    int requestGrids_count;
    ComputeGridForce **gf;
    int proc_count;
    int msgSentCount;
    int cur_grid;
    infostream *myout;
};


#endif







