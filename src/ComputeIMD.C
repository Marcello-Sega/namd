
#include "vmdsock.h"
#include "Node.h"
#include "imd.h"
#include "SimParameters.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "ComputeIMD.h"
#include "ComputeMgr.h"
#include <errno.h>

ComputeIMD::ComputeIMD(ComputeGlobal *h)
: ComputeGlobalMaster(h) {
  configMsg = 0;
  resultsMsg = 0;
  num_vmd_atoms = 0;
  vmd_atoms = NULL;
  vmd_forces = NULL;

  SimParameters *simparams = Node::Object()->simParameters;
  int port = simparams->IMDport;

  sock = vmdsock_create();
  int rc = vmdsock_bind(sock, port);
  if (rc < 0) {
    vmdsock_destroy(sock);
    NAMD_die("Unable to connect to given IMDport\n");
  }
  else {
    rc = vmdsock_listen(sock); 
    // Wait for VMD to connect
    iout << iINFO << "Waiting for VMD to connect to port "<<port<<'\n'<<endi;
    while (!vmdsock_selread(sock));
    rc = vmdsock_accept(sock);
    iout << iDEBUG << "Accept returned " << rc << '\n' << endi;
  }
  Node::Object()->IMDinit(sock);

}

ComputeIMD::~ComputeIMD() {
  if (sock) 
    vmdsock_destroy(sock);
}

void ComputeIMD::initialize() {
  configMsg = new ComputeGlobalConfigMsg;

  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0;
}

void ComputeIMD::calculate() {
  // Assume for now that the only thing we get from VMD is a set of forces.
  // Later we'll want to look for and implement more sophisticated control
  // parameters. 
 
  resultsMsg = new ComputeGlobalResultsMsg;
  resultsMsg->gforce.resize(gmass.size());

  // Check/get new forces from VMD
  if (get_vmd_forces()) {

    // Copy the forces from this class' member data into a ComputeGlobalMsg
    // and send it out.  These resize arrays are lame; if they were just 
    // pointers I could do a memcpy and call it a day.

    resultsMsg->aid.resize(num_vmd_atoms);
    resultsMsg->f.resize(num_vmd_atoms);
    int i;
    AtomIDList::iterator aid_i = resultsMsg->aid.begin();
    ForceList::iterator f_i = resultsMsg->f.begin();
    for ( i = 0; i < num_vmd_atoms; ++i ) {
      aid_i[i] = vmd_atoms[i];
      f_i[i].x = vmd_forces[3*i];
      f_i[i].y = vmd_forces[3*i+1];
      f_i[i].z = vmd_forces[3*i+2];
    }
  }
  // Send results to clients
  host->comm->sendComputeGlobalResults(resultsMsg);
  resultsMsg = 0;
}

int ComputeIMD::get_vmd_forces() {
  char *buf;
  IMDHeaderType htype;
  int hlength;
  int hsize; 
  int retval = 0;
  int rc;
  while (vmdsock_selread(sock))  {     // Drain the socket
    rc = imd_readheader(sock, &htype, &hlength, &hsize); 
    if (rc != sizeof(IMDHeader)) {
      iout << iDEBUG << "header: " << strerror(errno) << '\n' << endi;
      NAMD_die("Error reading header\n"); 
    }
    // interpret header
    switch (htype) {
      case IMD_MDCOMM:
        // Expect the msglength to give number of indicies, and the data
        // message to consist of first the indicies, then the coordinates
        // in xyz1 xyz2... format.
        if (num_vmd_atoms < hlength) { // need to resize
          delete [] vmd_atoms;
          delete [] vmd_forces;
          vmd_atoms = new int[hlength];
          vmd_forces = new float[hlength*3];
          num_vmd_atoms = hlength;
        }
        rc  = imd_readn(sock,(char *)vmd_atoms,num_vmd_atoms * sizeof(int));
        if (rc != num_vmd_atoms*sizeof(int)) {
          iout << iDEBUG << "read indices returned " << rc 
               << ", should have returned " << num_vmd_atoms*sizeof(int)
               << '\n' << endi; 
          NAMD_die("Error reading indices\n");
        }
        rc = imd_readn(sock,(char *)vmd_forces,num_vmd_atoms*3*sizeof(float));
        if (rc != num_vmd_atoms*3*sizeof(float)) {
          iout << iDEBUG << "read forces returned " << rc 
               << ", should have returned " << num_vmd_atoms*3*sizeof(float)
               << '\n' << endi; 
          NAMD_die("Error reading forces\n");
  	}
        retval = 1;
        break;
      default:
        iout << iWARN << "ComputeIMD: Unsupported header received \n "<<endi;
        buf = new char[hsize];
        rc = imd_readn(sock,buf,hsize);
        delete [] buf; 
        if (rc != hsize) {
          NAMD_die("Unable to read full message from VMD\n"); 
        }
    }
  }
  return retval;
}

