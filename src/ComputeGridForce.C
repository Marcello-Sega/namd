/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeGridForce.h"
#include "GridForceGrid.h"
#include "Node.h"
#include "SimParameters.h"
#include "HomePatch.h"
#include "Molecule.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

#include "MGridforceParams.h"

#define GF_FORCE_OUTPUT
#define GF_FORCE_OUTPUT_FREQ 100

ComputeGridForceNodeMgr::ComputeGridForceNodeMgr() 
{
  __sdag_init();

  CkpvAccess(BOCclass_group).computeGridForceNodeMgr = thisgroup;

  grids_deposited = 0;
  grids = NULL;
  requestGrids_count = 0;
  myNode = CkMyNode();
  numNodes = CkNumNodes();
  int nodeSize = CkNodeSize();
  activeProcCount = -1;
  activeProcs = new int[nodeSize];
  gf = new ComputeGridForce*[nodeSize];
  proc_count = 0;
  cur_grid = 0;
//  myout = new infostream[nodeSize];
//  iout << iINFO << "[" << CkMyPe() << "] GridForceNodeMgr created"
//       << " procs=" << CkNumPes()
//       << " nodes=" << CkNumNodes()
//       << " nodesize=" << CkNodeSize(CkMyNode())
//       << "\n" << endi;
}
  
void ComputeGridForceNodeMgr::depositInitialGrid(GridDepositMsg *msg)
{
  const int gridnum = msg->gridnum;
  GridforceGrid *grid = msg->grid;
  num_grids = msg->num_grids;
  delete msg;

//  iout << iINFO << "[" << CkMyPe() << "] Grid " << gridnum 
//       << " deposited\n" << endi;
  
  if (grids_deposited == 0) {
//    iout << iINFO << "[" << CkMyPe() << "] Allocating for  " 
//         << num_grids << " grids"
//         << "\n" << endi;
    grids = new GridforceGrid*[num_grids];
  }
  grids[gridnum] = grid;
  grids_deposited++;
  if (grids_deposited == num_grids) {
//    iout << iINFO << "[" << CkMyPe() << "]" 
//         << " all grids deposited\n" << endi;
    int i;
    for(i=0;i<num_grids; i++) {
      grids[i]->init2();
    }
//    iout << iINFO << "[" << CkMyPe() << "]" 
//             << " requesting from node 0"
//             << "\n" << endi;
    // Once all grids are deposited here, request grid data from node 0
    CProxy_ComputeGridForceNodeMgr myproxy(thisgroup);
    myproxy[0].requestInitialGridData();
  }
}

void ComputeGridForceNodeMgr::requestInitialGridData()
{
//  iout << iINFO << "[" << CkMyPe() << "] Grids requested " 
//       << requestGrids_count 
//       << "\n" << endi;
  requestGrids_count++;
  
  // When all of the grids on a node have been deposited initially
  // a message will get sent here. After that, each call to 
  // receiveInitialGridData on every node will generate another call
  // here until the data has all been sent
  if (requestGrids_count == CkNumNodes()) {
    CProxy_ComputeGridForceNodeMgr myproxy(thisgroup);
    myproxy[0].broadcastInitialGridData();
  }
#if 0
    requestGrids_count = 0;
    bool sent = false;
    while (!sent && cur_grid < num_grids) {
      GridSegmentMsg *msg = new GridSegmentMsg;
      
      msg->gridnum = cur_grid;
      grids[cur_grid]->init3(msg->grid,&(msg->start_index),
                             &(msg->count));
      
      if (msg->count > 0) {
//        CkPrintf("Msg is %p\n",msg);
//        iout << iINFO << "[" << CkMyPe() << "] Sending grid " 
//          << msg->gridnum 
//          << " of " << num_grids
//          << " start " << msg->start_index 
//          << " count " << msg->count
//          << "\n" << endi;
        CProxy_ComputeGridForceNodeMgr myproxy(thisgroup);
//        iout << iINFO << "[" << CkMyPe() << "] got proxy"
//             << "\n" << endi;
        myproxy.receiveInitialGridData(msg);
//        iout << iINFO << "["  << CkMyPe() << "] sent data from grid " 
//             << cur_grid
//             << "\n" << endi;
        sent = true;
      } else {
//        iout << iINFO << "[" << CkMyPe() << "] no data from grid " 
//             << msg->gridnum 
//             << "\n" << endi;
        cur_grid++;
      }
    }
//    if (!sent) {
//       iout << iINFO << "[" << CkMyPe() << "] Done sending grid data"
//         << "\n" << endi;
//    }
  }
//  iout << iINFO << "Done w/ request\n" << endi;
#endif
}

void ComputeGridForceNodeMgr::receiveInitialGridData(GridSegmentMsg *msg)
{
//  iout << iINFO << "[" << CkMyPe() << "] Grid received " 
//       << "\n" << endi;
//  iout << iINFO << "[" << CkMyPe() << "] Grid num=" << msg->gridnum
//       << " count=" << msg->count
//       << " start=" << msg->start_index
//       << "\n" << endi;
  grids[msg->gridnum]->init4(msg->grid,msg->start_index,msg->count);
  // delete msg;
  
  // CProxy_ComputeGridForceNodeMgr myproxy(thisgroup);
  // myproxy[0].requestInitialGridData();
}



void ComputeGridForceNodeMgr::submitRequest(SubReqMsg *msg)
{
//  CkPrintf("[%d] submitRequest received message %d from %d\n",
//           CkMyPe(),proc_count,msg->rank); 
  const int msgrank = msg->rank;

  // If we haven't figured out how many client procs this node has, do it now
  if (activeProcCount == -1) {
    int nodeSz = CkNodeSize();
    int i, pe;
    activeProcCount = 0;
    PatchMap* map = PatchMap::Object();
    for(i=0,pe=CkNodeFirst(myNode) ; i < nodeSz; pe++, i++) {
      if (map->numPatchesOnNode(pe)) {
	activeProcs[activeProcCount++] = pe;
      }
    }
  }
    
  // Save the data from this request
  gf[msgrank] = msg->gf;
  delete msg;
  proc_count++;

  if (proc_count < activeProcCount) {
//    iout << iINFO << "[" << CkMyPe() << "] " 
//         << " submitRequest rank=" << msgrank
//         << " still waiting"
//         << "\n" << endi;
    return;
  } else proc_count = 0;

  // All requests received. Figure out what we actually need
  int i,rank;
//  CkPrintf("[%d] num_grids=%d nodeSize=%d\n",
//           CkMyPe(),num_grids,nodeSize);

  for(i=0; i < num_grids; i++) {
    grids[i]->consolidateGridRequests();
  }
  
  // Package up grid element requests and send them off
  int num_indices = 0;
  for(i=0; i < num_grids; i++) {
    num_indices += grids[i]->getNumGridRequests();
  }
  
  if (num_indices == 0) {
//    iout << iINFO << "[" << CkMyPe() << "] no vals needed " 
//         << "\n" << endi;
    resume();
  } else {
//    CkPrintf("[%d] requesting %d grid vals\n",CkMyPe(),num_indices);

    msgSentCount=0;
    int *gridStart = new int[num_grids+1];
    int *gridIndices = new int[num_indices];
    int *gridNodeList = new int[num_indices];
    int *gridNodeCount = new int[CkNumNodes()];
    GridforceGrid::
      getGridIndices(gridStart,gridIndices,gridNodeList,gridNodeCount);

    GridRequestMsg **msgs = new GridRequestMsg*[CkNumNodes()];
    int *cur_ind = new int[CkNumNodes()];
    int i;
    for(i=0; i < CkNumNodes(); i++) {
      msgs[i] = new(num_grids+1,gridNodeCount[i]) GridRequestMsg;
      msgs[i]->from_node = myNode;
      msgs[i]->num_grid_indices = gridNodeCount[i];
      cur_ind[i] = 0;
    }

    int grid;
    CProxy_ComputeGridForceNodeMgr me(thisgroup);
    
    for(grid=0; grid < num_grids; grid++) {
      for(i=0; i < CkNumNodes(); i++)
	msgs[i]->gridStartIndex[grid] = cur_ind[i];

      for(i=gridStart[grid]; i < gridStart[grid+1]; i++) {
	int node = gridNodeList[i];
	msgs[node]->gridIndexList[cur_ind[node]] = gridIndices[i];
	cur_ind[node]++;
      }
    }
    for(i=0; i < CkNumNodes(); i++)
      msgs[i]->gridStartIndex[num_grids] = cur_ind[i];

    for(i=0;i<CkNumNodes();i++) {
      if (gridNodeCount[i] > 0) {
//       iout << iINFO << "[" << CkMyPe() << "] fetching from " << i
//             << "\n" << endi;
        me[i].fetchGridValues(msgs[i]);
        msgSentCount++;
      } else {
        delete msgs[i];
      }
    }
    delete [] cur_ind;
    delete [] msgs;
    delete [] gridNodeCount;
    delete [] gridNodeList;
    delete [] gridIndices;
    delete [] gridStart;
  }
//  iout << iINFO << "[" << CkMyPe() << "] " 
//       << "submitRequest DONE rank=" << msgrank
//       << "\n" << endi;
}

void ComputeGridForceNodeMgr::fetchGridValues(GridRequestMsg *msg) 
{
  // Be careful that this isn't running at the same time resume() is.
  // As long as this entry is exclusive, and resume() is only
  // called from an exclusive entry, this should be fine.
//  iout << iINFO << "[" << CkMyPe() << "] Getting grid values from " 
//       << msg->from_node << " for " << msg->num_grid_indices 
//       << " indices "
//       << "\n" << endi;
  int *startIndex = new int[num_grids+1];
  startIndex[0] = 0;
  int gridnum;
  int i;
  for(gridnum=0;gridnum<num_grids;gridnum++) {
    int nvals = 0;
    for(i=msg->gridStartIndex[gridnum];
        i < msg->gridStartIndex[gridnum+1]; i++) {
      const int block = msg->gridIndexList[i];
      nvals += grids[gridnum]->getBlockSize(block);
//      iout << iINFO << "[" << CkMyPe() << "] grid block "
//           << block << " size is " 
//         << nvals << "\n" << endi;
    }
    startIndex[gridnum+1] = startIndex[gridnum] + nvals;
  }
  GridValuesMsg *outmsg 
    = new(num_grids+1,startIndex[num_grids]) GridValuesMsg;
  outmsg->num_grid_indices = startIndex[num_grids];

  for(gridnum=0;gridnum<=num_grids;gridnum++) {
//    iout << iINFO << "[" << CkMyPe() << "] Grid start " 
//         << startIndex[gridnum]
//         << "\n" << endi;
    outmsg->gridStartIndex[gridnum] = startIndex[gridnum];
  }
  delete [] startIndex;
  
//  Molecule *mol = Node::Object()->molecule;
  for (gridnum = 0; gridnum < num_grids; gridnum++) {
    GridforceGrid* grid = grids[gridnum];
//    iout << iINFO << "[" << CkMyPe() << "] filling grid from "
//         << msg->gridStartIndex[gridnum]
//         << " to " 
//         << msg->gridStartIndex[gridnum+1]
//         << "\n" << endi;
    int outindex = outmsg->gridStartIndex[gridnum];
    for(i=msg->gridStartIndex[gridnum];
        i < msg->gridStartIndex[gridnum+1]; i++) {
      const int block = msg->gridIndexList[i];
      int j;
      const int jmax = grid->getBlockSize(block);
//      iout << iINFO << "[" << CkMyPe() << "] grid block[" << i << "] "
//           << block << " size is "
//         << jmax << " offset " << outindex
//         << "\n" << endi;

      for(j=0; j < jmax; j++, outindex++) {
        const int idx = grid->getIndexForBlock(block,j);
//        CkPrintf("PVal %d Idx %d block %d boffset %d\n",
//                  j1,idx,block,j);
        outmsg->gridVals[outindex].index = idx;
        outmsg->gridVals[outindex].val = grid->get_grid(idx);
//        iout << iINFO << "[" << CkMyPe() << "] Got grid value("
//             << i << ")(" << outmsg->gridVals[i].index << ") "
//             << outmsg->gridVals[i].val
//             << "\n" << endi;
      }
    }
  }

  CProxy_ComputeGridForceNodeMgr me(thisgroup);
  me[msg->from_node].recvGridValues(outmsg);
  delete msg;
}

void ComputeGridForceNodeMgr::recvGridValues(GridValuesMsg *msg)
{
//  iout << iINFO << "[" << CkMyPe() << "] Receiving grid values" 
//       << "\n" << endi;

  int gridnum;
  const int myrank = CkMyRank();
  for (gridnum = 0; gridnum < num_grids; gridnum++) {
    GridforceGrid *grid = grids[gridnum];
    int i;
//    CkPrintf("Unpacking from %d to %d\n",
//             msg->gridStartIndex[gridnum],msg->gridStartIndex[gridnum+1]);
    for(i=msg->gridStartIndex[gridnum]; 
        i < msg->gridStartIndex[gridnum+1]; i++) {
        const int idx = msg->gridVals[i].index;
        int bi,boffset;
        grid->getGridBlockIdxOffset(idx,bi,boffset);
//        CkPrintf("UVal %d Idx %d block %d boffset %d\n",
//                  i,idx,bi,boffset);
        grid->addGridValue(idx,msg->gridVals[i].val);
/*         iout << iINFO << "[" << CkMyPe() << "] Unpacking grid value(" */
/*            << i << ")(" << msg->gridVals[i].index << ") " */
/*            << grid->get_real_grid(msg->gridVals[i].index)  */
/*            << " --- " << msg->gridVals[i].val */
/*            << "\n" << endi; */
/*          */
    }
  }
  msgSentCount--;
//  iout << iINFO << "[" << CkMyPe() << "] Done receiving grid values" 
//       << " msgSentCount=" << msgSentCount 
//       << "\n" << endi;
  
  if (msgSentCount == 0) {
    resume();
  }
  delete msg;
}

ComputeGridForceMgr::ComputeGridForceMgr() 
{
  CkpvAccess(BOCclass_group).computeGridForceMgr = thisgroup;
}

void ComputeGridForceMgr::request(ComputeGridForce *mygf) {
//  iout << iINFO << "[" << CkMyPe() << "] Sending off request " 
//       << "\n" << endi;
  gf = mygf;
  SubReqMsg *msg = new SubReqMsg;
  msg->rank = CkMyRank();
  msg->gf = mygf;
  CProxy_ComputeGridForceNodeMgr nodemgr(CkpvAccess(BOCclass_group).
                                         computeGridForceNodeMgr);
  nodemgr[CkMyNode()].submitRequest(msg);
}

void ComputeGridForceNodeMgr::resume() 
{
  // Grid calculations should be cheap, so we won't bother parallelizing
  // them. As long as only one PE accesses the cache at a time, we don't
  // have to synchronize. Another option is to execute each finishWork()
  // on a separate PE, but then we have to make sure the GridCache is
  // thread-safe, or sequentialize access to it
  int rank;
  CProxy_ComputeGridForceMgr pemgr(CkpvAccess(BOCclass_group).
                                   computeGridForceMgr);

  for(rank=0; rank < activeProcCount; rank++) {
    pemgr[activeProcs[rank]].finishWork();
  }
}

void ComputeGridForceMgr::finishWork() {
  gf->finishWork();
}

ComputeGridForce::ComputeGridForce(ComputeID c)
    : ComputeHomePatches(c)
{

    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}
/*			END OF FUNCTION ComputeGridForce	*/


ComputeGridForce::~ComputeGridForce()
{
    delete reduction;
}
/*			END OF FUNCTION ~ComputeGridForce	*/

void ComputeGridForce::doWork() {
  Molecule *mol = Node::Object()->molecule;
  ResizeArrayIter<PatchElem> ap(patchList);
  DebugM(3,patchID << ": doWork() called.\n");

  
  int i = 0;
  num_boxes = patchList.size();
  box_state = new BoxState[num_boxes];
  
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    // Open up positionBox, forceBox, and atomBox
    box_state[i].posBox = (*ap).positionBox;
    box_state[i].atoms = box_state[i].posBox->open();
    box_state[i].forceBox = (*ap).forceBox;
    box_state[i].results = box_state[i].forceBox->open();
    box_state[i].homePatch = (*ap).p;
    FullAtom* a = box_state[i].homePatch->getAtomList().begin();
    doForce(box_state[i].homePatch,a);
  //  CkPrintf("[%d] Opening box %d for patch %d\n",
//             CkMyPe(),i,box_state[i].homePatch);
    // Leave boxes open, close them in finishWork()
    i++;
  }
  ComputeGridForceMgr* mgr = CProxy_ComputeGridForceMgr::
    ckLocalBranch(CkpvAccess(BOCclass_group).computeGridForceMgr);
  mgr->request(this);

//  iout << iINFO << "Printing grid access\n" << endi;
//  for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
//    const GridforceGrid* grid = mol->get_gridfrc_grid(gridnum);
//    grid->print_elems_used();
//  }
//  CkPrintf("[%d] doWork for %d patches\n",CkMyPe(),num_boxes);

  DebugM(2,patchID << ": doWork() completed.\n");
}

void ComputeGridForce::finishWork() 
{
  Molecule *mol = Node::Object()->molecule;
  ResizeArrayIter<PatchElem> ap(patchList);
  DebugM(3,patchID << ": doWork() called.\n");
  
//  CkPrintf("[%d] finishWork for %d patches\n",CkMyPe(),num_boxes);

  int i;
  for (i=0; i < num_boxes; i++) {
    // Open up positionBox, forceBox, and atomBox
    // Have to re-open them to get the pointers
    CompAtom *x = box_state[i].atoms;
    Results *r = box_state[i].results;
    HomePatch *p = box_state[i].homePatch;
    FullAtom* a = p->getAtomList().begin();
    finishForce(p,a,r);
    box_state[i].posBox->close(&box_state[i].atoms);
    box_state[i].forceBox->close(&box_state[i].results);
//    CkPrintf("[%d] Closing box %d for patch %d\n",
//             CkMyPe(),i,box_state[i].homePatch);
  }
//  iout << iINFO << "Printing grid access\n" << endi;
//  for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
//    const GridforceGrid* grid = mol->get_gridfrc_grid(gridnum);
//    grid->print_elems_used();
//  }
  DebugM(2,patchID << ": finishWork() completed.\n");
  reduction->submit();
  delete [] box_state;
  box_state = 0;
}

void ComputeGridForce::doForce(HomePatch* homePatch,FullAtom* p)
{
    SimParameters *simParams = Node::Object()->simParameters;
    Molecule *mol = Node::Object()->molecule;
    //Lattice lattice = homePatch->lattice;
    
//    Force *forces = r->f[Results::normal];
//    BigReal energy = 0;
//    Force extForce = 0.;
//    Tensor extVirial;
    
//    Real scale;			// Scaling factor
//    Charge charge;		// Charge
//    GridforceGrid::Box box;	// Structure with potential info
//    Force f;
//    float v;
    int numAtoms = homePatch->getNumAtoms();
    
    DebugM(4, "numGridforceGrids = " << mol->numGridforceGrids << "\n" << endi);
    
    for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
	GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
	DebugM(4, "scale = " << grid->get_scale() << "\n" << endi);
    
    Position center = grid->get_center();
    Tensor inv = grid->get_inv();
    Vector gfScale = grid->get_scale();
    
    // Check that grid doesn't cross periodic bounds
    if (CkMyPe() == 0) {
	Position origin = grid->get_origin();
	Tensor e = grid->get_e();
	
	for (int i0 = 0; i0 < 2; i0++) {
	    float x = i0 * (grid->get_k0()-1);
	    for (int i1 = 0; i1 < 2; i1++) {
		float y = i1 * (grid->get_k1()-1);
		for (int i2 = 0; i2 < 2; i2++) {
		    float z = i2 * (grid->get_k2()-1);
		    
		    Position pos_orig = origin + e * Position(x, y, z);
                    Position pos = pos_orig + homePatch->lattice.wrap_delta(pos_orig);
		    pos += homePatch->lattice.delta(pos, center) - (pos - center);
		    
		    if ((pos - pos_orig).length() > 1.0) {	// FIX
			char err_msg[128];

			sprintf(err_msg, "Gridforce grid too large for periodic cell: (%f %f %f) != (%f %f %f)",
				pos.x, pos.y, pos.z, pos_orig.x, pos_orig.y, pos_orig.z);
			NAMD_die(err_msg);
                    }

		}
	    }
	}
    }
    
    //  Loop through and check each atom
    for (int i = 0; i < numAtoms; i++) {
//        iout << iINFO << "Gridforced " << p[i].id << "--" << gridnum << "\n" << endi;

	if (mol->is_atom_gridforced(p[i].id, gridnum))
	{

//	    mol->get_gridfrc_params(scale, charge, p[i].id, gridnum);
		
	    // Wrap coordinates using grid center
	    Position pos = p[i].position;
	    pos += homePatch->lattice.wrap_delta(p[i].position);
	    pos += homePatch->lattice.delta(pos, center) - (pos - center);
		
	    grid->request_box(pos);
	}
    }
  }
}

//int ComputeGridForce::getNumGridIndices() const {
//  Molecule *mol = Node::Object()->molecule;
//  int gridnum;
//  int grid_index=0;
    
  // Count how much total data we need, and store the pointers
  // for the start of each grid's requested data
//  for (gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
//    const GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
//    int req_sz = grid->getGridRequests()->size();
//    grid_index += req_sz;
//  }
//  return grid_index;
//}


void ComputeGridForce::finishForce(HomePatch* homePatch,FullAtom* p,
                                   Results* r)
{
    SimParameters *simParams = Node::Object()->simParameters;
    Molecule *mol = Node::Object()->molecule;
    //Lattice lattice = homePatch->lattice;
    
    Force *forces = r->f[Results::normal];
    BigReal energy = 0;
    Force extForce = 0.;
    Tensor extVirial;
    
    Real scale;			// Scaling factor
    Charge charge;		// Charge
    GridforceGrid::Box box;	// Structure with potential info
    Force f;
    float v;
    int numAtoms = homePatch->getNumAtoms();
    
    DebugM(4, "numGridforceGrids = " << mol->numGridforceGrids << "\n" << endi);
    
    for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
	const GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
	DebugM(4, "scale = " << grid->get_scale() << "\n" << endi);
    
    Position center = grid->get_center();
    Tensor inv = grid->get_inv();
    Vector gfScale = grid->get_scale();
    
    DebugM(4, "center = " << center << "\n" << endi);
    
    //  Loop through and check each atom
    for (int i = 0; i < numAtoms; i++) {
//        iout << iINFO << "Gridforced " << p[i].id << "--" << gridnum << "\n" << endi;

	if (mol->is_atom_gridforced(p[i].id, gridnum))
	{

	    mol->get_gridfrc_params(scale, charge, p[i].id, gridnum);
		
	    // Wrap coordinates using grid center
	    Position pos = p[i].position;
	    pos += homePatch->lattice.wrap_delta(p[i].position);
	    pos += homePatch->lattice.delta(pos, center) - (pos - center);
		
	    int err = grid->get_box(&box, pos);
	    if (err) {
		DebugM(4, "ts = " << patch->flags.step << " gridnum = " << gridnum << " force4 = 0 0 0\n" << endi);
		DebugM(4, "v4 = 0\n" << endi);
		    		    
		continue;  // This means the current atom is outside the potential
	    }
	    
	    float a[64];
// 	    for (int j = 0; j < 64; j++) {
// 		a[j] = 0;
// 		for (int k = 0; k < 64; k++) {
// 		    a[j] += A[j][k] * box.b[k];
// 		}
// 	    }
	    
	    // Multiply 'box.b' vector by matrix ... ugly but efficient (?)
	    a[0] = box.b[0];
	    a[1] = box.b[8];
	    a[2] = -3*box.b[0] + 3*box.b[1] - 2*box.b[8] - box.b[9];
	    a[3] = 2*box.b[0] - 2*box.b[1] + box.b[8] + box.b[9];
	    a[4] = box.b[16];
	    a[5] = box.b[32];
	    a[6] = -3*box.b[16] + 3*box.b[17] - 2*box.b[32] - box.b[33];
	    a[7] = 2*box.b[16] - 2*box.b[17] + box.b[32] + box.b[33];
	    a[8] = -3*box.b[0] + 3*box.b[2] - 2*box.b[16] - box.b[18];
	    a[9] = -3*box.b[8] + 3*box.b[10] - 2*box.b[32] - box.b[34];
	    a[10] = 9*box.b[0] - 9*box.b[1] - 9*box.b[2] + 9*box.b[3] + 6*box.b[8] + 3*box.b[9] - 6*box.b[10] - 3*box.b[11]
		+ 6*box.b[16] - 6*box.b[17] + 3*box.b[18] - 3*box.b[19] + 4*box.b[32] + 2*box.b[33] + 2*box.b[34] + box.b[35];
	    a[11] = -6*box.b[0] + 6*box.b[1] + 6*box.b[2] - 6*box.b[3] - 3*box.b[8] - 3*box.b[9] + 3*box.b[10] + 3*box.b[11]
		- 4*box.b[16] + 4*box.b[17] - 2*box.b[18] + 2*box.b[19] - 2*box.b[32] - 2*box.b[33] - box.b[34] - box.b[35];
	    a[12] = 2*box.b[0] - 2*box.b[2] + box.b[16] + box.b[18];
	    a[13] = 2*box.b[8] - 2*box.b[10] + box.b[32] + box.b[34];
	    a[14] = -6*box.b[0] + 6*box.b[1] + 6*box.b[2] - 6*box.b[3] - 4*box.b[8] - 2*box.b[9] + 4*box.b[10] + 2*box.b[11]
		- 3*box.b[16] + 3*box.b[17] - 3*box.b[18] + 3*box.b[19] - 2*box.b[32] - box.b[33] - 2*box.b[34] - box.b[35];
	    a[15] = 4*box.b[0] - 4*box.b[1] - 4*box.b[2] + 4*box.b[3] + 2*box.b[8] + 2*box.b[9] - 2*box.b[10] - 2*box.b[11]
		+ 2*box.b[16] - 2*box.b[17] + 2*box.b[18] - 2*box.b[19] + box.b[32] + box.b[33] + box.b[34] + box.b[35];
	    a[16] = box.b[24];
	    a[17] = box.b[40];
	    a[18] = -3*box.b[24] + 3*box.b[25] - 2*box.b[40] - box.b[41];
	    a[19] = 2*box.b[24] - 2*box.b[25] + box.b[40] + box.b[41];
	    a[20] = box.b[48];
	    a[21] = box.b[56];
	    a[22] = -3*box.b[48] + 3*box.b[49] - 2*box.b[56] - box.b[57];
	    a[23] = 2*box.b[48] - 2*box.b[49] + box.b[56] + box.b[57];
	    a[24] = -3*box.b[24] + 3*box.b[26] - 2*box.b[48] - box.b[50];
	    a[25] = -3*box.b[40] + 3*box.b[42] - 2*box.b[56] - box.b[58];
	    a[26] = 9*box.b[24] - 9*box.b[25] - 9*box.b[26] + 9*box.b[27] + 6*box.b[40] + 3*box.b[41] - 6*box.b[42] - 3*box.b[43]
		+ 6*box.b[48] - 6*box.b[49] + 3*box.b[50] - 3*box.b[51] + 4*box.b[56] + 2*box.b[57] + 2*box.b[58] + box.b[59];
	    a[27] = -6*box.b[24] + 6*box.b[25] + 6*box.b[26] - 6*box.b[27] - 3*box.b[40] - 3*box.b[41] + 3*box.b[42] + 3*box.b[43]
		- 4*box.b[48] + 4*box.b[49] - 2*box.b[50] + 2*box.b[51] - 2*box.b[56] - 2*box.b[57] - box.b[58] - box.b[59];
	    a[28] = 2*box.b[24] - 2*box.b[26] + box.b[48] + box.b[50];
	    a[29] = 2*box.b[40] - 2*box.b[42] + box.b[56] + box.b[58];
	    a[30] = -6*box.b[24] + 6*box.b[25] + 6*box.b[26] - 6*box.b[27] - 4*box.b[40] - 2*box.b[41] + 4*box.b[42] + 2*box.b[43]
		- 3*box.b[48] + 3*box.b[49] - 3*box.b[50] + 3*box.b[51] - 2*box.b[56] - box.b[57] - 2*box.b[58] - box.b[59];
	    a[31] = 4*box.b[24] - 4*box.b[25] - 4*box.b[26] + 4*box.b[27] + 2*box.b[40] + 2*box.b[41] - 2*box.b[42] - 2*box.b[43]
		+ 2*box.b[48] - 2*box.b[49] + 2*box.b[50] - 2*box.b[51] + box.b[56] + box.b[57] + box.b[58] + box.b[59];
	    a[32] = -3*box.b[0] + 3*box.b[4] - 2*box.b[24] - box.b[28];
	    a[33] = -3*box.b[8] + 3*box.b[12] - 2*box.b[40] - box.b[44];
	    a[34] = 9*box.b[0] - 9*box.b[1] - 9*box.b[4] + 9*box.b[5] + 6*box.b[8] + 3*box.b[9] - 6*box.b[12] - 3*box.b[13]
		+ 6*box.b[24] - 6*box.b[25] + 3*box.b[28] - 3*box.b[29] + 4*box.b[40] + 2*box.b[41] + 2*box.b[44] + box.b[45];
	    a[35] = -6*box.b[0] + 6*box.b[1] + 6*box.b[4] - 6*box.b[5] - 3*box.b[8] - 3*box.b[9] + 3*box.b[12] + 3*box.b[13]
		- 4*box.b[24] + 4*box.b[25] - 2*box.b[28] + 2*box.b[29] - 2*box.b[40] - 2*box.b[41] - box.b[44] - box.b[45];
	    a[36] = -3*box.b[16] + 3*box.b[20] - 2*box.b[48] - box.b[52];
	    a[37] = -3*box.b[32] + 3*box.b[36] - 2*box.b[56] - box.b[60];
	    a[38] = 9*box.b[16] - 9*box.b[17] - 9*box.b[20] + 9*box.b[21] + 6*box.b[32] + 3*box.b[33] - 6*box.b[36] - 3*box.b[37]
		+ 6*box.b[48] - 6*box.b[49] + 3*box.b[52] - 3*box.b[53] + 4*box.b[56] + 2*box.b[57] + 2*box.b[60] + box.b[61];
	    a[39] = -6*box.b[16] + 6*box.b[17] + 6*box.b[20] - 6*box.b[21] - 3*box.b[32] - 3*box.b[33] + 3*box.b[36] + 3*box.b[37]
		- 4*box.b[48] + 4*box.b[49] - 2*box.b[52] + 2*box.b[53] - 2*box.b[56] - 2*box.b[57] - box.b[60] - box.b[61];
	    a[40] = 9*box.b[0] - 9*box.b[2] - 9*box.b[4] + 9*box.b[6] + 6*box.b[16] + 3*box.b[18] - 6*box.b[20] - 3*box.b[22]
		+ 6*box.b[24] - 6*box.b[26] + 3*box.b[28] - 3*box.b[30] + 4*box.b[48] + 2*box.b[50] + 2*box.b[52] + box.b[54];
	    a[41] = 9*box.b[8] - 9*box.b[10] - 9*box.b[12] + 9*box.b[14] + 6*box.b[32] + 3*box.b[34] - 6*box.b[36] - 3*box.b[38]
		+ 6*box.b[40] - 6*box.b[42] + 3*box.b[44] - 3*box.b[46] + 4*box.b[56] + 2*box.b[58] + 2*box.b[60] + box.b[62];
	    a[42] = -27*box.b[0] + 27*box.b[1] + 27*box.b[2] - 27*box.b[3] + 27*box.b[4] - 27*box.b[5] - 27*box.b[6] + 27*box.b[7]
		- 18*box.b[8] - 9*box.b[9] + 18*box.b[10] + 9*box.b[11] + 18*box.b[12] + 9*box.b[13] - 18*box.b[14] - 9*box.b[15]
		- 18*box.b[16] + 18*box.b[17] - 9*box.b[18] + 9*box.b[19] + 18*box.b[20] - 18*box.b[21] + 9*box.b[22] - 9*box.b[23]
		- 18*box.b[24] + 18*box.b[25] + 18*box.b[26] - 18*box.b[27] - 9*box.b[28] + 9*box.b[29] + 9*box.b[30] - 9*box.b[31]
		- 12*box.b[32] - 6*box.b[33] - 6*box.b[34] - 3*box.b[35] + 12*box.b[36] + 6*box.b[37] + 6*box.b[38] + 3*box.b[39]
		- 12*box.b[40] - 6*box.b[41] + 12*box.b[42] + 6*box.b[43] - 6*box.b[44] - 3*box.b[45] + 6*box.b[46] + 3*box.b[47]
		- 12*box.b[48] + 12*box.b[49] - 6*box.b[50] + 6*box.b[51] - 6*box.b[52] + 6*box.b[53] - 3*box.b[54] + 3*box.b[55]
		- 8*box.b[56] - 4*box.b[57] - 4*box.b[58] - 2*box.b[59] - 4*box.b[60] - 2*box.b[61] - 2*box.b[62] - box.b[63];
	    a[43] = 18*box.b[0] - 18*box.b[1] - 18*box.b[2] + 18*box.b[3] - 18*box.b[4] + 18*box.b[5] + 18*box.b[6] - 18*box.b[7]
		+ 9*box.b[8] + 9*box.b[9] - 9*box.b[10] - 9*box.b[11] - 9*box.b[12] - 9*box.b[13] + 9*box.b[14] + 9*box.b[15]
		+ 12*box.b[16] - 12*box.b[17] + 6*box.b[18] - 6*box.b[19] - 12*box.b[20] + 12*box.b[21] - 6*box.b[22] + 6*box.b[23]
		+ 12*box.b[24] - 12*box.b[25] - 12*box.b[26] + 12*box.b[27] + 6*box.b[28] - 6*box.b[29] - 6*box.b[30] + 6*box.b[31]
		+ 6*box.b[32] + 6*box.b[33] + 3*box.b[34] + 3*box.b[35] - 6*box.b[36] - 6*box.b[37] - 3*box.b[38] - 3*box.b[39]
		+ 6*box.b[40] + 6*box.b[41] - 6*box.b[42] - 6*box.b[43] + 3*box.b[44] + 3*box.b[45] - 3*box.b[46] - 3*box.b[47]
		+ 8*box.b[48] - 8*box.b[49] + 4*box.b[50] - 4*box.b[51] + 4*box.b[52] - 4*box.b[53] + 2*box.b[54] - 2*box.b[55]
		+ 4*box.b[56] + 4*box.b[57] + 2*box.b[58] + 2*box.b[59] + 2*box.b[60] + 2*box.b[61] + box.b[62] + box.b[63];
	    a[44] = -6*box.b[0] + 6*box.b[2] + 6*box.b[4] - 6*box.b[6] - 3*box.b[16] - 3*box.b[18] + 3*box.b[20] + 3*box.b[22]
		- 4*box.b[24] + 4*box.b[26] - 2*box.b[28] + 2*box.b[30] - 2*box.b[48] - 2*box.b[50] - box.b[52] - box.b[54];
	    a[45] = -6*box.b[8] + 6*box.b[10] + 6*box.b[12] - 6*box.b[14] - 3*box.b[32] - 3*box.b[34] + 3*box.b[36] + 3*box.b[38]
		- 4*box.b[40] + 4*box.b[42] - 2*box.b[44] + 2*box.b[46] - 2*box.b[56] - 2*box.b[58] - box.b[60] - box.b[62];
	    a[46] = 18*box.b[0] - 18*box.b[1] - 18*box.b[2] + 18*box.b[3] - 18*box.b[4] + 18*box.b[5] + 18*box.b[6] - 18*box.b[7]
		+ 12*box.b[8] + 6*box.b[9] - 12*box.b[10] - 6*box.b[11] - 12*box.b[12] - 6*box.b[13] + 12*box.b[14] + 6*box.b[15]
		+ 9*box.b[16] - 9*box.b[17] + 9*box.b[18] - 9*box.b[19] - 9*box.b[20] + 9*box.b[21] - 9*box.b[22] + 9*box.b[23]
		+ 12*box.b[24] - 12*box.b[25] - 12*box.b[26] + 12*box.b[27] + 6*box.b[28] - 6*box.b[29] - 6*box.b[30] + 6*box.b[31]
		+ 6*box.b[32] + 3*box.b[33] + 6*box.b[34] + 3*box.b[35] - 6*box.b[36] - 3*box.b[37] - 6*box.b[38] - 3*box.b[39]
		+ 8*box.b[40] + 4*box.b[41] - 8*box.b[42] - 4*box.b[43] + 4*box.b[44] + 2*box.b[45] - 4*box.b[46] - 2*box.b[47]
		+ 6*box.b[48] - 6*box.b[49] + 6*box.b[50] - 6*box.b[51] + 3*box.b[52] - 3*box.b[53] + 3*box.b[54] - 3*box.b[55]
		+ 4*box.b[56] + 2*box.b[57] + 4*box.b[58] + 2*box.b[59] + 2*box.b[60] + box.b[61] + 2*box.b[62] + box.b[63];
	    a[47] = -12*box.b[0] + 12*box.b[1] + 12*box.b[2] - 12*box.b[3] + 12*box.b[4] - 12*box.b[5] - 12*box.b[6] + 12*box.b[7]
		- 6*box.b[8] - 6*box.b[9] + 6*box.b[10] + 6*box.b[11] + 6*box.b[12] + 6*box.b[13] - 6*box.b[14] - 6*box.b[15]
		- 6*box.b[16] + 6*box.b[17] - 6*box.b[18] + 6*box.b[19] + 6*box.b[20] - 6*box.b[21] + 6*box.b[22] - 6*box.b[23]
		- 8*box.b[24] + 8*box.b[25] + 8*box.b[26] - 8*box.b[27] - 4*box.b[28] + 4*box.b[29] + 4*box.b[30] - 4*box.b[31]
		- 3*box.b[32] - 3*box.b[33] - 3*box.b[34] - 3*box.b[35] + 3*box.b[36] + 3*box.b[37] + 3*box.b[38] + 3*box.b[39]
		- 4*box.b[40] - 4*box.b[41] + 4*box.b[42] + 4*box.b[43] - 2*box.b[44] - 2*box.b[45] + 2*box.b[46] + 2*box.b[47]
		- 4*box.b[48] + 4*box.b[49] - 4*box.b[50] + 4*box.b[51] - 2*box.b[52] + 2*box.b[53] - 2*box.b[54] + 2*box.b[55]
		- 2*box.b[56] - 2*box.b[57] - 2*box.b[58] - 2*box.b[59] - box.b[60] - box.b[61] - box.b[62] - box.b[63];
	    a[48] = 2*box.b[0] - 2*box.b[4] + box.b[24] + box.b[28];
	    a[49] = 2*box.b[8] - 2*box.b[12] + box.b[40] + box.b[44];
	    a[50] = -6*box.b[0] + 6*box.b[1] + 6*box.b[4] - 6*box.b[5] - 4*box.b[8] - 2*box.b[9] + 4*box.b[12] + 2*box.b[13]
		- 3*box.b[24] + 3*box.b[25] - 3*box.b[28] + 3*box.b[29] - 2*box.b[40] - box.b[41] - 2*box.b[44] - box.b[45];
	    a[51] = 4*box.b[0] - 4*box.b[1] - 4*box.b[4] + 4*box.b[5] + 2*box.b[8] + 2*box.b[9] - 2*box.b[12] - 2*box.b[13]
		+ 2*box.b[24] - 2*box.b[25] + 2*box.b[28] - 2*box.b[29] + box.b[40] + box.b[41] + box.b[44] + box.b[45];
	    a[52] = 2*box.b[16] - 2*box.b[20] + box.b[48] + box.b[52];
	    a[53] = 2*box.b[32] - 2*box.b[36] + box.b[56] + box.b[60];
	    a[54] = -6*box.b[16] + 6*box.b[17] + 6*box.b[20] - 6*box.b[21] - 4*box.b[32] - 2*box.b[33] + 4*box.b[36] + 2*box.b[37]
		- 3*box.b[48] + 3*box.b[49] - 3*box.b[52] + 3*box.b[53] - 2*box.b[56] - box.b[57] - 2*box.b[60] - box.b[61];
	    a[55] = 4*box.b[16] - 4*box.b[17] - 4*box.b[20] + 4*box.b[21] + 2*box.b[32] + 2*box.b[33] - 2*box.b[36] - 2*box.b[37]
		+ 2*box.b[48] - 2*box.b[49] + 2*box.b[52] - 2*box.b[53] + box.b[56] + box.b[57] + box.b[60] + box.b[61];
	    a[56] = -6*box.b[0] + 6*box.b[2] + 6*box.b[4] - 6*box.b[6] - 4*box.b[16] - 2*box.b[18] + 4*box.b[20] + 2*box.b[22]
		- 3*box.b[24] + 3*box.b[26] - 3*box.b[28] + 3*box.b[30] - 2*box.b[48] - box.b[50] - 2*box.b[52] - box.b[54];
	    a[57] = -6*box.b[8] + 6*box.b[10] + 6*box.b[12] - 6*box.b[14] - 4*box.b[32] - 2*box.b[34] + 4*box.b[36] + 2*box.b[38]
		- 3*box.b[40] + 3*box.b[42] - 3*box.b[44] + 3*box.b[46] - 2*box.b[56] - box.b[58] - 2*box.b[60] - box.b[62];
	    a[58] = 18*box.b[0] - 18*box.b[1] - 18*box.b[2] + 18*box.b[3] - 18*box.b[4] + 18*box.b[5] + 18*box.b[6] - 18*box.b[7]
		+ 12*box.b[8] + 6*box.b[9] - 12*box.b[10] - 6*box.b[11] - 12*box.b[12] - 6*box.b[13] + 12*box.b[14] + 6*box.b[15]
		+ 12*box.b[16] - 12*box.b[17] + 6*box.b[18] - 6*box.b[19] - 12*box.b[20] + 12*box.b[21] - 6*box.b[22] + 6*box.b[23]
		+ 9*box.b[24] - 9*box.b[25] - 9*box.b[26] + 9*box.b[27] + 9*box.b[28] - 9*box.b[29] - 9*box.b[30] + 9*box.b[31]
		+ 8*box.b[32] + 4*box.b[33] + 4*box.b[34] + 2*box.b[35] - 8*box.b[36] - 4*box.b[37] - 4*box.b[38] - 2*box.b[39]
		+ 6*box.b[40] + 3*box.b[41] - 6*box.b[42] - 3*box.b[43] + 6*box.b[44] + 3*box.b[45] - 6*box.b[46] - 3*box.b[47]
		+ 6*box.b[48] - 6*box.b[49] + 3*box.b[50] - 3*box.b[51] + 6*box.b[52] - 6*box.b[53] + 3*box.b[54] - 3*box.b[55]
		+ 4*box.b[56] + 2*box.b[57] + 2*box.b[58] + box.b[59] + 4*box.b[60] + 2*box.b[61] + 2*box.b[62] + box.b[63];
	    a[59] = -12*box.b[0] + 12*box.b[1] + 12*box.b[2] - 12*box.b[3] + 12*box.b[4] - 12*box.b[5] - 12*box.b[6] + 12*box.b[7]
		- 6*box.b[8] - 6*box.b[9] + 6*box.b[10] + 6*box.b[11] + 6*box.b[12] + 6*box.b[13] - 6*box.b[14] - 6*box.b[15]
		- 8*box.b[16] + 8*box.b[17] - 4*box.b[18] + 4*box.b[19] + 8*box.b[20] - 8*box.b[21] + 4*box.b[22] - 4*box.b[23]
		- 6*box.b[24] + 6*box.b[25] + 6*box.b[26] - 6*box.b[27] - 6*box.b[28] + 6*box.b[29] + 6*box.b[30] - 6*box.b[31]
		- 4*box.b[32] - 4*box.b[33] - 2*box.b[34] - 2*box.b[35] + 4*box.b[36] + 4*box.b[37] + 2*box.b[38] + 2*box.b[39]
		- 3*box.b[40] - 3*box.b[41] + 3*box.b[42] + 3*box.b[43] - 3*box.b[44] - 3*box.b[45] + 3*box.b[46] + 3*box.b[47]
		- 4*box.b[48] + 4*box.b[49] - 2*box.b[50] + 2*box.b[51] - 4*box.b[52] + 4*box.b[53] - 2*box.b[54] + 2*box.b[55]
		- 2*box.b[56] - 2*box.b[57] - box.b[58] - box.b[59] - 2*box.b[60] - 2*box.b[61] - box.b[62] - box.b[63];
	    a[60] = 4*box.b[0] - 4*box.b[2] - 4*box.b[4] + 4*box.b[6] + 2*box.b[16] + 2*box.b[18] - 2*box.b[20] - 2*box.b[22]
		+ 2*box.b[24] - 2*box.b[26] + 2*box.b[28] - 2*box.b[30] + box.b[48] + box.b[50] + box.b[52] + box.b[54];
	    a[61] = 4*box.b[8] - 4*box.b[10] - 4*box.b[12] + 4*box.b[14] + 2*box.b[32] + 2*box.b[34] - 2*box.b[36] - 2*box.b[38]
		+ 2*box.b[40] - 2*box.b[42] + 2*box.b[44] - 2*box.b[46] + box.b[56] + box.b[58] + box.b[60] + box.b[62];
	    a[62] = -12*box.b[0] + 12*box.b[1] + 12*box.b[2] - 12*box.b[3] + 12*box.b[4] - 12*box.b[5] - 12*box.b[6] + 12*box.b[7]
		- 8*box.b[8] - 4*box.b[9] + 8*box.b[10] + 4*box.b[11] + 8*box.b[12] + 4*box.b[13] - 8*box.b[14] - 4*box.b[15]
		- 6*box.b[16] + 6*box.b[17] - 6*box.b[18] + 6*box.b[19] + 6*box.b[20] - 6*box.b[21] + 6*box.b[22] - 6*box.b[23]
		- 6*box.b[24] + 6*box.b[25] + 6*box.b[26] - 6*box.b[27] - 6*box.b[28] + 6*box.b[29] + 6*box.b[30] - 6*box.b[31]
		- 4*box.b[32] - 2*box.b[33] - 4*box.b[34] - 2*box.b[35] + 4*box.b[36] + 2*box.b[37] + 4*box.b[38] + 2*box.b[39]
		- 4*box.b[40] - 2*box.b[41] + 4*box.b[42] + 2*box.b[43] - 4*box.b[44] - 2*box.b[45] + 4*box.b[46] + 2*box.b[47]
		- 3*box.b[48] + 3*box.b[49] - 3*box.b[50] + 3*box.b[51] - 3*box.b[52] + 3*box.b[53] - 3*box.b[54] + 3*box.b[55]
		- 2*box.b[56] - box.b[57] - 2*box.b[58] - box.b[59] - 2*box.b[60] - box.b[61] - 2*box.b[62] - box.b[63];
	    a[63] = 8*box.b[0] - 8*box.b[1] - 8*box.b[2] + 8*box.b[3] - 8*box.b[4] + 8*box.b[5] + 8*box.b[6] - 8*box.b[7]
		+ 4*box.b[8] + 4*box.b[9] - 4*box.b[10] - 4*box.b[11] - 4*box.b[12] - 4*box.b[13] + 4*box.b[14] + 4*box.b[15]
		+ 4*box.b[16] - 4*box.b[17] + 4*box.b[18] - 4*box.b[19] - 4*box.b[20] + 4*box.b[21] - 4*box.b[22] + 4*box.b[23]
		+ 4*box.b[24] - 4*box.b[25] - 4*box.b[26] + 4*box.b[27] + 4*box.b[28] - 4*box.b[29] - 4*box.b[30] + 4*box.b[31]
		+ 2*box.b[32] + 2*box.b[33] + 2*box.b[34] + 2*box.b[35] - 2*box.b[36] - 2*box.b[37] - 2*box.b[38] - 2*box.b[39]
		+ 2*box.b[40] + 2*box.b[41] - 2*box.b[42] - 2*box.b[43] + 2*box.b[44] + 2*box.b[45] - 2*box.b[46] - 2*box.b[47]
		+ 2*box.b[48] - 2*box.b[49] + 2*box.b[50] - 2*box.b[51] + 2*box.b[52] - 2*box.b[53] + 2*box.b[54] - 2*box.b[55]
		+ box.b[56] + box.b[57] + box.b[58] + box.b[59] + box.b[60] + box.b[61] + box.b[62] + box.b[63];
	    
	    for (int j = 0; j < 64; j++) DebugM(2, "a[" << j << "] = " << a[j] << "\n" << endi);
	    for (int j = 0; j < 64; j++) DebugM(2, "b[" << j << "] = " << box.b[j] << "\n" << endi);
	    
	    // Calculate powers of x, y, z for later use
	    // e.g. x[2] = x^2
	    float x[4], y[4], z[4];
	    x[0] = 1; y[0] = 1; z[0] = 1;
	    for (int j = 1; j < 4; j++) {
		x[j] = x[j-1] * box.loc.x;
		y[j] = y[j-1] * box.loc.y;
		z[j] = z[j-1] * box.loc.z;
	    }
	    
	    int ind = 0;
	    f = 0;
	    v = 0;
	    for (int l = 0; l < 4; l++) {
		for (int k = 0; k < 4; k++) {
		    for (int j = 0; j < 4; j++) {
			v += a[ind] * x[j] * y[k] * z[l];
			if (j > 0) f.x -= a[ind] * j * x[j-1] * y[k]   * z[l];
			if (k > 0) f.y -= a[ind] * k * x[j]   * y[k-1] * z[l];
			if (l > 0) f.z -= a[ind] * l * x[j]   * y[k]   * z[l-1];
			ind++;
		    }
		}
	    }
	    
// 	    Force force = scale * Tensor::diagonal(simParams->gridforceScale) * charge * (inv * f);
//	    Force force = scale * gfScale * charge * (f * inv); // Must multiply ON THE RIGHT by inv tensor
	    Force force = scale * Tensor::diagonal(gfScale) * charge * ((box.scale * f) * inv); // Must multiply ON THE RIGHT by inv tensor
	    
	    DebugM(4, "f4 = " << f << "\n" << endi);
	    DebugM(4, "ts = " << patch->flags.step << " gridnum = " << gridnum << " force4 = " << force << "\n" << endi);
	    DebugM(4, "v4 = " << v << "\n" << endi);
	    
	    //#ifdef GF_FORCE_OUTPUT
	    //	    int t = patch->flags.step;
	    //	    if (t % GF_FORCE_OUTPUT_FREQ == 0) {
	    //		iout << "GRIDFORCE (TS ATOM FORCE)  " << t << ' ' << i << ' ' << force*PNPERKCALMOL << '\n' << endi;
	    //	    }
	    //#endif
	    
	    forces[i] += force;
	    extForce += force;
	    Position vpos = homePatch->lattice.reverse_transform(p[i].position, p[i].transform);
	    
	    DebugM(4, "transform = " << (int)p[i].transform.i << " "
		   << (int)p[i].transform.j << " " << (int)p[i].transform.k << "\n" << endi);
	    
	    //energy -= force * (vpos - homePatch->lattice.origin());
	    if (gfScale.x == gfScale.y && gfScale.x == gfScale.z)
	    {
		// only makes sense when scaling is isotropic
		energy += v * charge * scale * gfScale.x;
	    }
	    extVirial += outer(force,vpos);
	}
    }

    }
//    for(int i=0; i < numAtoms; i++) {
//      iout << "ZZZAtom " << p[i].id << " Force " << forces[i].x << "," 
//           << forces[i].y << "," << forces[i].z
//           << endl;
//    }
   reduction->item(REDUCTION_MISC_ENERGY) += energy;
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
    
}
/*			END OF FUNCTION force				*/


#include "ComputeGridForceMgr.def.h"
