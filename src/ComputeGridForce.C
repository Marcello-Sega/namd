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

#define MIN_DEBUG_LEVEL 2
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
  int nodeSize = CkNodeSize(myNode);
  activeProcCount = -1;
  activeProcs = new int[nodeSize];
  gf = new ComputeGridForce*[nodeSize];
  proc_count = 0;
  cur_grid = 0;
}
  
void ComputeGridForceNodeMgr::depositInitialGrid(GridDepositMsg *msg)
{
  //CkPrintf("depositInitialGrid\n");
  const int gridnum = msg->gridnum;
  GridforceGrid *grid = msg->grid;
  num_grids = msg->num_grids;
  DebugM(4, "depositInitialGrid: read message: gridnum = " << gridnum << ", num_grids = " << num_grids << "\n");
  delete msg;

  if (grids_deposited == 0) {
    DebugM(4, "creating grids array\n");
    grids = new GridforceGrid*[num_grids];
  }
  grids[gridnum] = grid;
  grids_deposited++;
  if (grids_deposited == num_grids) {
    DebugM(4, "all grids deposited, calling init2\n");
    int i;
    for(i=0;i<num_grids; i++) {
      grids[i]->init2();
    }
    // Once all grids are deposited here, request grid data from node 0
    CProxy_ComputeGridForceNodeMgr myproxy(thisgroup);
    myproxy[0].requestInitialGridData();
  }
}

void ComputeGridForceNodeMgr::requestInitialGridData()
{
  requestGrids_count++;
  
  // When all of the grids on a node have been deposited initially
  // a message will get sent here. After that, each call to 
  // receiveInitialGridData on every node will generate another call
  // here until the data has all been sent
  if (requestGrids_count == CkNumNodes()) {
    CProxy_ComputeGridForceNodeMgr myproxy(thisgroup);
    myproxy[0].broadcastInitialGridData();
  }
}

void ComputeGridForceNodeMgr::receiveInitialGridData(GridSegmentMsg *msg)
{
  grids[msg->gridnum]->init4(msg->grid,msg->start_index,msg->count);
}


void ComputeGridForceNodeMgr::submitRequest(SubReqMsg *msg)
{
  const int msgrank = msg->rank;

  // If we haven't figured out how many client procs this node has, do it now
  if (activeProcCount == -1) {
    int nodeSz = CkNodeSize(myNode);
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

  for(i=0; i < num_grids; i++) {
    //CkPrintf("consolidating grid requests\n");
    grids[i]->consolidateGridRequests();
  }
  
  // Package up grid element requests and send them off
  int num_indices = 0;
  for(i=0; i < num_grids; i++) {
    num_indices += grids[i]->getNumGridRequests();
    //CkPrintf("[%d] grid %d num_indices = %d\n", CkMyNode(), i, num_indices);
  }

  if (num_indices == 0) {
    resume();
  } else {
    msgSentCount=0;
    int *gridStart = new int[num_grids+1];
    int *gridIndices = new int[num_indices];
    int *gridNodeList = new int[num_indices];
    int *gridNodeCount = new int[CkNumNodes()];
    GridforceGrid::
      getGridIndices(gridStart,gridIndices,gridNodeList,gridNodeCount);
    
    //CkPrintf("[%d] creating messages\n", CkMyNode());
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
    
    //CkPrintf("[%d] filling messages\n", CkMyNode());
    for(grid=0; grid < num_grids; grid++) {
      for(i=0; i < CkNumNodes(); i++)
	msgs[i]->gridStartIndex[grid] = cur_ind[i];
      
      for(i=gridStart[grid]; i < gridStart[grid+1]; i++) {
	int node = gridNodeList[i];
	msgs[node]->gridIndexList[cur_ind[node]] = gridIndices[i];
	cur_ind[node]++;
      }
    }
    //CkPrintf("[%d] filling messages 2\n", CkMyNode());
    for(i=0; i < CkNumNodes(); i++)
      msgs[i]->gridStartIndex[num_grids] = cur_ind[i];

    //CkPrintf("[%d] fetching grid values\n", CkMyNode());
    for(i=0;i<CkNumNodes();i++) {
      if (gridNodeCount[i] > 0) {
        me[i].fetchGridValues(msgs[i]);
        msgSentCount++;
      } else {
        delete msgs[i];
      }
    }
    //CkPrintf("[%d] done\n", CkMyNode());
    delete [] cur_ind;
    delete [] msgs;
    delete [] gridNodeCount;
    delete [] gridNodeList;
    delete [] gridIndices;
    delete [] gridStart;
  }
}

void ComputeGridForceNodeMgr::fetchGridValues(GridRequestMsg *msg) 
{
  // Be careful that this isn't running at the same time resume() is.
  // As long as this entry is exclusive, and resume() is only
  // called from an exclusive entry, this should be fine.
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
    }
    startIndex[gridnum+1] = startIndex[gridnum] + nvals;
  }
  GridValuesMsg *outmsg 
    = new(num_grids+1,startIndex[num_grids]) GridValuesMsg;
  outmsg->num_grid_indices = startIndex[num_grids];

  for(gridnum=0;gridnum<=num_grids;gridnum++) {
    outmsg->gridStartIndex[gridnum] = startIndex[gridnum];
  }
  delete [] startIndex;
  
  for (gridnum = 0; gridnum < num_grids; gridnum++) {
    GridforceGrid* grid = grids[gridnum];
    int outindex = outmsg->gridStartIndex[gridnum];
    for(i=msg->gridStartIndex[gridnum];
        i < msg->gridStartIndex[gridnum+1]; i++) {
      const int block = msg->gridIndexList[i];
      int j;
      const int jmax = grid->getBlockSize(block);

      for(j=0; j < jmax; j++, outindex++) {
        const int idx = grid->getIndexForBlock(block,j);
        outmsg->gridVals[outindex].index = idx;
        outmsg->gridVals[outindex].val = grid->get_grid(idx);
      }
    }
  }

  CProxy_ComputeGridForceNodeMgr me(thisgroup);
  me[msg->from_node].recvGridValues(outmsg);
  delete msg;
}

void ComputeGridForceNodeMgr::recvGridValues(GridValuesMsg *msg)
{
  int gridnum;
  const int myrank = CkMyRank();
  for (gridnum = 0; gridnum < num_grids; gridnum++) {
    GridforceGrid *grid = grids[gridnum];
    int i;
    for(i=msg->gridStartIndex[gridnum]; 
        i < msg->gridStartIndex[gridnum+1]; i++) {
        const int idx = msg->gridVals[i].index;
        int bi,boffset;
        grid->getGridBlockIdxOffset(idx,bi,boffset);
        grid->addGridValue(idx,msg->gridVals[i].val);
    }
  }
  msgSentCount--;
  
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
    // Leave boxes open, close them in finishWork()
    i++;
  }
  ComputeGridForceMgr* mgr = CProxy_ComputeGridForceMgr::
    ckLocalBranch(CkpvAccess(BOCclass_group).computeGridForceMgr);
  mgr->request(this);

  //DebugM(2,patchID << ": doWork() completed.\n");
}

void ComputeGridForce::finishWork() 
{
  Molecule *mol = Node::Object()->molecule;
  ResizeArrayIter<PatchElem> ap(patchList);
  //DebugM(3,patchID << ": doWork() called.\n");
  
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
  }
  //DebugM(2,patchID << ": finishWork() completed.\n");
  reduction->submit();
  delete [] box_state;
  box_state = 0;
}

void ComputeGridForce::doForce(HomePatch* homePatch,FullAtom* p)
{
    SimParameters *simParams = Node::Object()->simParameters;
    Molecule *mol = Node::Object()->molecule;
    int numAtoms = homePatch->getNumAtoms();
    
    for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
	GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
	
	Position center = grid->get_center();
	
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
//	    iout << iINFO << "Gridforced " << p[i].id << "--" << gridnum << "\n" << endi;
	    
	    if (mol->is_atom_gridforced(p[i].id, gridnum))
	    {
//		mol->get_gridfrc_params(scale, charge, p[i].id, gridnum);
		
		// Wrap coordinates using grid center
		Position pos = p[i].position;
		pos += homePatch->lattice.wrap_delta(p[i].position);
		pos += homePatch->lattice.delta(pos, center) - (pos - center);
		
		DebugM(2, "calling request_box() for maingrid " << gridnum << "\n" << endi);
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
    //GridforceGrid::Box box;	// Structure with potential info
    Vector dV;
    float V;
    int numAtoms = homePatch->getNumAtoms();
    
    for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
	const GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
	
	Position center = grid->get_center();
//	Tensor inv = grid->get_inv();
	Vector gfScale = grid->get_scale();
	
	DebugM(1, "center = " << center << "\n" << endi);
	
	//  Loop through and check each atom
	for (int i = 0; i < numAtoms; i++) {
//	    iout << iINFO << "Gridforced " << p[i].id << "--" << gridnum << "\n" << endi;
	    
	    if (mol->is_atom_gridforced(p[i].id, gridnum))
	    {
		DebugM(1, "Atom " << p[i].id << " is gridforced\n" << endi);
		
		mol->get_gridfrc_params(scale, charge, p[i].id, gridnum);
		
		// Wrap coordinates using grid center
		Position pos = p[i].position;
		pos += homePatch->lattice.wrap_delta(p[i].position);
		pos += homePatch->lattice.delta(pos, center) - (pos - center);
		
		DebugM(1, "pos = " << pos << "\n" << endi);
		
		// Here's where the action happens
		int err = grid->compute_VdV(pos, V, dV);
		if (err) {
		    DebugM(2, "V = 0\n" << endi);
		    DebugM(2, "dV = 0 0 0\n" << endi);
		    continue;  // This means the current atom is outside the potential
		}
		
		Force force = scale * Tensor::diagonal(gfScale) * (-charge * dV);
		
		DebugM(2, "V = " << V << "\n" << endi);
		DebugM(2, "dV = " << dV << "\n" << endi);
		DebugM(2, "force = " << force << " pos = " << pos << " V = " << V << " dV = " << dV << " step = " << homePatch->flags.step << " index = " << p[i].id << "\n" << endi);
		
		if (V != V) {
		    iout << iWARN << "V is NaN!\natomid = " << p[i].id << " loc = " << p[i].position << "\n" << endi;
		}
		
		//#ifdef GF_FORCE_OUTPUT
		//	    int t = patch->flags.step;
		//	    if (t % GF_FORCE_OUTPUT_FREQ == 0) {
		//		iout << "GRIDFORCE (TS ATOM FORCE)  " << t << ' ' << i << ' ' << force*PNPERKCALMOL << '\n' << endi;
		//	    }
		//#endif
	    
		forces[i] += force;
		extForce += force;
		Position vpos = homePatch->lattice.reverse_transform(p[i].position, p[i].transform);
		
		DebugM(1, "transform = " << (int)p[i].transform.i << " "
		       << (int)p[i].transform.j << " " << (int)p[i].transform.k << "\n" << endi);
		
		//energy -= force * (vpos - homePatch->lattice.origin());
		if (gfScale.x == gfScale.y && gfScale.x == gfScale.z)
		{
		    // only makes sense when scaling is isotropic
		    energy += scale * gfScale.x * (charge * V);
		    
		    // add something when we're off the grid? I'm thinking no
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
