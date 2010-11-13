/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/**
 * AtomMap maps the global atom ID (int) to the atom's assigned
 * patch ID (int) and local index (int) within that patch.
 * An array of "LocalID" of length (number of atoms) is allocated.
 * The total space required is 2*sizeof(int)*(number of atoms).
 */

#ifndef ATOMMAP_H
#define ATOMMAP_H

#include "NamdTypes.h"

#ifdef MEM_OPT_VERSION
#include "ckhashtable.h"

class MyHashtable: public CkHashtableT<int, LocalID>{
public:	
	
	MyHashtable(int initLen, float loadFactor, CkHashFunction hf, CkHashCompare hcmp) : 
		CkHashtableT<int, LocalID>(initLen, loadFactor, hf, hcmp) {}
	
	~MyHashtable() { empty(); }
	
	inline LocalID *getEntry(int key){
		return (LocalID *)CkHashtable::get((const void *)&key);
	}
	
	inline void putEntry(int key, LocalID &ent){
		LocalID *pos = (LocalID *)CkHashtable::put((const void *)&key);
		*pos = ent;		
	}
	
	//return -1 if not found because atom id will never be -1.
	int getCollidedEntry(int key, LocalID **retEnt);	
};
#endif

enum { notUsed = -1 };

class AtomMap
{
public:
  static AtomMap *Instance();
  inline static AtomMap *Object() { return CkpvAccess(AtomMap_instance); }
  ~AtomMap(void);
  void checkMap();

  void allocateMap(int nAtomIDs);

  int registerIDs(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end);
  int registerIDsFullAtom(PatchID pid, const FullAtom *begin, const FullAtom *end);
  int unregisterIDs(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end);
  int unregisterIDs(PatchID pid, const FullAtom *begin, const FullAtom *end);

  LocalID localID(AtomID id);

  void clearMap(void);
  void print(void);


protected:
  AtomMap(void);

private:
#ifdef MEM_OPT_VERSION
	int *keys; //one-to-one mapping to the localIDTable
	LocalID *localIDTable; //the first-level mapped values;
	int tableSz; //it should be less than "1<<MAXBITS"
	bool onlyUseTbl;
	MyHashtable *collidedAtoms; //the second-level mapped values;
	bool *isCollided;
#else
  LocalID *localIDTable;
  int tableSz;
  bool cleared;
#endif

};

#ifndef MEM_OPT_VERSION
//----------------------------------------------------------------------
// LocalID contains patch pid and local patch atom index
// for a given global atom number
inline LocalID AtomMap::localID(AtomID id)
{
  return localIDTable[id];
}
#endif

#endif /* ATOMMAP_H */

