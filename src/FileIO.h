#include "FileIO.decl.h"


class FileIO : public CBase_FileIO {

 private:
  int numAtomsToWriteAvg;

 public:
  FileIO(CkMigrateMessage* msg);
  FileIO();

  void writeDCD(int base);
  void doneDCD(CkReductionMsg *msg);
  void recvCoordinates();

};
