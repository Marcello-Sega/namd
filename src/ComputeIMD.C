
#include "vmdsock.h"
#include "Node.h"
#include "imd.h"
#include "SimParameters.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "ComputeIMD.h"
#include "ComputeMgr.h"
#include "UniqueSortedArray.h"
#include <errno.h>

struct vmdforce {
  int index;
  Vector force;
  int operator <(const vmdforce& v) {return index < v.index;}
  int operator ==(const vmdforce& v) {return index == v.index;}
  vmdforce& operator=(const vmdforce& v) {
    index=v.index;
    force=v.force; 
    return *this;
  }  
};

static UniqueSortedArray<vmdforce> vmdforces;

ComputeIMD::ComputeIMD(ComputeGlobal *h)
: ComputeGlobalMaster(h) {
  SimParameters *simparams = Node::Object()->simParameters;
  int port = simparams->IMDport;

  sock = vmdsock_create();
  int rc = vmdsock_bind(sock, port);
  if (rc < 0) {
    vmdsock_destroy(sock);
    NAMD_die("Unable to connect to given IMDport\n");
  }
  rc = vmdsock_listen(sock); 
  // Wait for VMD to connect
  iout << iINFO << "Waiting for VMD to connect to port "<<port<<'\n'<<endi;
  while ( !(rc = vmdsock_selread(sock)));
  if (rc < 0) {
    iout << iDEBUG << "select: " << strerror(errno) << '\n' << endi;
    vmdsock_destroy(sock);
    NAMD_die("Connection error\n");
  }
  rc = vmdsock_accept(sock);
  if (rc < 0) {
    iout << iDEBUG << "Accept: " << strerror(errno) << '\n' << endi;
    NAMD_die("Connection error\n");
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
  get_vmd_forces();

  // Right now I don't check to see if any new forces were obtained.
  // An optimization would be cache the results message.  However, there
  // would still be copying since it looks like the messages get deleted
  // by the receiver.

    int num = vmdforces.size();
    resultsMsg->aid.resize(num);
    resultsMsg->f.resize(num);
    int i;
    AtomIDList::iterator aid_i = resultsMsg->aid.begin();
    ForceList::iterator f_i = resultsMsg->f.begin();
    UniqueSortedArray<vmdforce>::iterator v_i = vmdforces.begin();
    for ( i = 0; i < num; ++i, ++v_i) {
      aid_i[i] = v_i->index; 
      f_i[i] = v_i->force;
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
    vmdforce *vtest, vnew;
    int i;
    int vmd_atoms[100];
    float vmd_forces[300];
    switch (htype) {
      case IMD_MDCOMM:
        // Expect the msglength to give number of indicies, and the data
        // message to consist of first the indicies, then the coordinates
        // in xyz1 xyz2... format.
        rc  = imd_readn(sock,(char *)vmd_atoms,hlength* sizeof(int));
        if (rc != hlength*sizeof(int)) {
          iout << iDEBUG << "read indices returned " << rc 
               << ", should have returned " << hlength*sizeof(int)
               << '\n' << endi; 
          NAMD_die("Error reading indices\n");
        }
        rc = imd_readn(sock,(char *)vmd_forces,hlength*3*sizeof(float));
        if (rc != hlength*3*sizeof(float)) {
          iout << iDEBUG << "read forces returned " << rc 
               << ", should have returned " << hlength*3*sizeof(float)
               << '\n' << endi; 
          NAMD_die("Error reading forces\n");
  	}
        for (i=0; i<hlength; i++) {
          vnew.index = vmd_atoms[i];
          if ( (vtest=vmdforces.find(vnew)) != NULL) {
            // find was successful, so overwrite the old force values
            vtest->force.x = vmd_forces[3*i];
            vtest->force.y = vmd_forces[3*i+1];
            vtest->force.z = vmd_forces[3*i+2];
          }
          else {
            // Create a new entry in the table
            vnew.force.x = vmd_forces[3*i];
            vnew.force.y = vmd_forces[3*i+1];
            vnew.force.z = vmd_forces[3*i+2];
            vmdforces.add(vnew);
          }
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

