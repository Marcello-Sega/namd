#ifndef COMMUNICATE_OBJ_H
#define COMMUNICATE_OBJ_H

#include "Message.h"

// dummy communicate class 



class Communicate 
{

public:
  // send method
  enum SendMethod { NOW, WAIT };

  // error codes
  enum CommError { NOERROR, ERROR, NONODES, NOSEND, NORECEIVE };


public:
  // constructor and destructor
  Communicate(void) {}
  virtual ~Communicate(void) {}


  // return info about connections in general
  int nodes(void) { return 1; }
  int this_node(void) { return 0; }

  int send(Message *msg, int node, int tag, int delmsg = TRUE) { return 0; }

  int send_all(void){ return 0; };

  int send_now(Message *msg, int node, int tag, int delmsg = TRUE) 
  {
    int retval;
    retval = send(msg, node, tag, delmsg);
    return retval;
  }
  
  Message *receive(int& node, int& tag){return 0;}

  int broadcast_all(Message *msg , int i){ return 0; }
  int broadcast_others(Message *msg, int i, int delmsg=TRUE){ return 0; }
};

#endif
