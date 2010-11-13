/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Tracks location of Atoms on node.  Singleton.
*/

#include "ProcessorPrivate.h"
#include "AtomMap.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

#ifdef MEM_OPT_VERSION
#define MAXBITS 20
#define MAXNUMATOMS (1<<MAXBITS)

static inline int firstHashCode(const int key){
	return key & (MAXNUMATOMS-1);
}

CkHashCode CollidedHashFunc(const void *key, size_t keyLen){
	int aid = *((int *)key);
	CkHashCode ret = 0;
	//the oldkey is the reason that two atoms compete for the same
	//slot of the localIDTable
	int oldkey = firstHashCode(aid);
	ret = oldkey*11 + ((MAXNUMATOMS*3)>>3) + 3103;
	return ret;
}

int CollidedCmpFunc(const void *key1, const void *key2, size_t keyLen){
	int k1 = *((int *)key1);
	int k2 = *((int *)key2);
	
	return (k1==k2);
}

//Returns the key of the first entry in the hash table that has the same
//hash code of the "key". -1 is returned if the entry is not found. -1
//is OK here because the atom id will never be -1. -Chao Mei
int MyHashtable::getCollidedEntry(int key, LocalID **retEnt){
	int oldHash = firstHashCode(key);
	int i = hash(&key, sizeof(int))%len;
	int startSpot = i;
	do{
		char *cur = entry(i);
		if(layout.isEmpty(cur)) {
			*retEnt = NULL;
			return -1;
		}
		char *curKey = layout.getKey(cur);
		int retKey = *(int *)curKey;
		int newHash = firstHashCode(retKey);
		if(newHash == oldHash){
			*retEnt = (LocalID *)layout.getObject(cur);
			return *((int *)curKey);
		}		
	}while(inc(i) != startSpot);
	*retEnt = NULL;
	return -1;	
}

#endif

// Singleton method
AtomMap *AtomMap::Instance() {
  if (CkpvAccess(AtomMap_instance) == 0) {
    CkpvAccess(AtomMap_instance) = new AtomMap;	// this is never deleted!
  }
  return CkpvAccess(AtomMap_instance);
}

//----------------------------------------------------------------------
AtomMap::AtomMap(void)
{
#ifndef MEM_OPT_VERSION
  localIDTable = NULL;
  cleared = false;
#else
	keys = NULL;
	localIDTable = NULL;
	onlyUseTbl = false;
	collidedAtoms = NULL;
	isCollided = NULL;
#endif
}

void
AtomMap::checkMap(void)
{ }
  

//----------------------------------------------------------------------
AtomMap::~AtomMap(void)
{
#ifndef MEM_OPT_VERSION
  delete [] localIDTable;  // Delete on a NULL pointer should be ok
#else
	delete [] localIDTable;
	delete [] keys;		
	delete collidedAtoms;
	delete [] isCollided;	
#endif
}

//----------------------------------------------------------------------
// Creates fixed size table
void AtomMap::allocateMap(int nAtomIds)
{
#ifdef MEM_OPT_VERSION
	if(nAtomIds <= MAXNUMATOMS){
		localIDTable = new LocalID[nAtomIds];
		for(int i=0; i < nAtomIds; i++)
			localIDTable[i].pid = localIDTable[i].index = notUsed;
		tableSz = nAtomIds;
		onlyUseTbl = true;	
	}else{
		localIDTable = new LocalID[MAXNUMATOMS];
		keys = new int[MAXNUMATOMS];
		isCollided = new bool[MAXNUMATOMS];
		tableSz = MAXNUMATOMS;
		for(int i=0; i < MAXNUMATOMS; i++){			
			localIDTable[i].pid = localIDTable[i].index = notUsed;
			keys[i] = notUsed;
			isCollided[i] = false;
		}		
		collidedAtoms = new MyHashtable(1<<(32-MAXBITS),0.75, 
												CollidedHashFunc, CollidedCmpFunc);
	}
#else
  localIDTable = new LocalID[nAtomIds];
  tableSz = nAtomIds;
  for(int i=0; i < nAtomIds; i++)
    localIDTable[i].pid = localIDTable[i].index = notUsed;
  cleared = true;
#endif
}

//
int AtomMap::unregisterIDs(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end)
{
#ifndef MEM_OPT_VERSION
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const CompAtomExt *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	if (localIDTable[ali].pid == pid) {
	    localIDTable[ali].pid = notUsed;
	    localIDTable[ali].index = notUsed;
	}
    }
    return 0;
  }
#else
	if(onlyUseTbl){
		for(const CompAtomExt *a = begin; a != end; ++a){
			unsigned int ali = a->id;
			if(localIDTable[ali].pid == pid){
				localIDTable[ali].pid = notUsed;
				localIDTable[ali].index = notUsed;
			}
		}		
	}else{
		for(const CompAtomExt *a = begin; a != end; ++a){
			unsigned int ali = a->id;			
			int firstTry = firstHashCode(ali);
			if(keys[firstTry] == ali){
				//possibly find the atom in the localIDTable
				if(localIDTable[firstTry].pid == pid){										
					//move the atom from hashtable to the localIDTable if there's one
					LocalID *ent;
					int otherKey = -1;
					if(isCollided[firstTry]){
						otherKey = collidedAtoms->getCollidedEntry(ali, &ent);
					}
					if(otherKey!=-1){						
						keys[firstTry] = otherKey;
						localIDTable[firstTry].pid = ent->pid;
						localIDTable[firstTry].index = ent->index;
						collidedAtoms->remove(otherKey);
					}else{
						keys[firstTry] = notUsed;
						localIDTable[firstTry].pid = notUsed;
						localIDTable[firstTry].index = notUsed;
						isCollided[firstTry] = false;						
					}										
				}				
			}else if(keys[firstTry] != notUsed){
				//find the atom in the hash table
				LocalID *ent = collidedAtoms->getEntry(ali);
				CmiAssert(ent != NULL);
				if(ent->pid == pid){					
					collidedAtoms->remove(ali);
				}
			}

			//possible to unregister an atom before it is in AtomMap
		}	
	}
	
	return 0; 
#endif
}

//----------------------------------------------------------------------
int AtomMap::unregisterIDs(PatchID pid, const FullAtom *begin, const FullAtom *end)
{
#ifndef MEM_OPT_VERSION	
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const FullAtom *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	if (localIDTable[ali].pid == pid) {
	    localIDTable[ali].pid = notUsed;
	    localIDTable[ali].index = notUsed;
	}
    }
    return 0;
  }
#else
	if(onlyUseTbl){
		for(const FullAtom *a = begin; a != end; ++a){
			unsigned int ali = a->id;
			if(localIDTable[ali].pid == pid){
				localIDTable[ali].pid = notUsed;
				localIDTable[ali].index = notUsed;
			}
		}		
	}else{
		for(const FullAtom *a = begin; a != end; ++a){
			unsigned int ali = a->id;					
			int firstTry = firstHashCode(ali);
			if(keys[firstTry] == ali){
				//possibly find the atom in the localIDTable
				if(localIDTable[firstTry].pid == pid){										
					//move the atom from hashtable to the localIDTable if there's one
					LocalID *ent;
					int otherKey = -1;
					if(isCollided[firstTry]) {
						otherKey = collidedAtoms->getCollidedEntry(ali, &ent);
					}
					if(otherKey!=-1){						
						keys[firstTry] = otherKey;
						localIDTable[firstTry].pid = ent->pid;
						localIDTable[firstTry].index = ent->index;
						collidedAtoms->remove(otherKey);
					}else{
						keys[firstTry] = notUsed;
						localIDTable[firstTry].pid = notUsed;
						localIDTable[firstTry].index = notUsed;
						isCollided[firstTry] = false;
					}
				}				
			}else if(keys[firstTry] != notUsed){
				//find the atom in the hash table
				LocalID *ent = collidedAtoms->getEntry(ali);
				CmiAssert(ent!=NULL);
				if(ent->pid == pid){
					collidedAtoms->remove(ali);
				}
			}
			//it's possible to unregister an atom id that is not in the atom
			//map, such as the case of calling doAtomMigration from 
			//HomePatch::positionReady.				
		}	
	}
	
	return 0; 
#endif  
}

//It's possible to register the same atom for a new patch before it is moved 
//from the old patch on the same processor!

//----------------------------------------------------------------------
int AtomMap::registerIDs(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end)
{
#ifndef MEM_OPT_VERSION	
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const CompAtomExt *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	localIDTable[ali].pid = pid;
	localIDTable[ali].index = a - begin;
    }
    cleared = false;
    return 0;
  }
#else
	if(localIDTable == NULL) return -1;
	if(onlyUseTbl){
		for(const CompAtomExt *a = begin; a != end; ++a)
		{
			unsigned int ali = a->id;
			localIDTable[ali].pid = pid;
			localIDTable[ali].index = a - begin;
		}    		
	}else{
		LocalID one;
		one.pid = pid;
		for(const CompAtomExt *a = begin; a != end; ++a)
		{
			unsigned int ali = a->id;
			int firstTry = firstHashCode(ali);
			if(keys[firstTry] == notUsed){
				keys[firstTry] = ali;
				localIDTable[firstTry].pid = pid;
				localIDTable[firstTry].index = a-begin;
			}else if(keys[firstTry] == ali){
				localIDTable[firstTry].pid = pid;
				localIDTable[firstTry].index = a-begin;
			}else{
				//collision happens, insert this atom into hashtable
				one.index = a-begin;
				collidedAtoms->putEntry(ali, one);
				isCollided[firstTry] = true;
			}
		}    
	}
	return 0;  
#endif  
}

//----------------------------------------------------------------------
int AtomMap::registerIDsFullAtom(PatchID pid, const FullAtom *begin, const FullAtom *end)
{
#ifndef MEM_OPT_VERSION	
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const FullAtom *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	localIDTable[ali].pid = pid;
	localIDTable[ali].index = a - begin;
    }
    cleared = false;
    return 0;
  }
#else
	if(localIDTable == NULL) return -1;
	if(onlyUseTbl){
		for(const FullAtom *a = begin; a != end; ++a)
		{
			unsigned int ali = a->id;
			localIDTable[ali].pid = pid;
			localIDTable[ali].index = a - begin;			
		}    		
	}else{
		LocalID one;
		one.pid = pid;
		for(const FullAtom *a = begin; a != end; ++a)
		{
			unsigned int ali = a->id;						
			int firstTry = firstHashCode(ali);
			if(keys[firstTry] == notUsed){
				keys[firstTry] = ali;
				localIDTable[firstTry].pid = pid;
				localIDTable[firstTry].index = a-begin;				
			}else if(keys[firstTry] == ali){
				localIDTable[firstTry].pid = pid;
				localIDTable[firstTry].index = a-begin;				
			}else{
				//collision happens, insert this atom into hashtable
				one.index = a-begin;
				collidedAtoms->putEntry(ali, one);	
				isCollided[firstTry] = true;			
			}
		}    
	}
	return 0;  
#endif  
}

//----------------------------------------------------------------------
// resets map to notUsed condition
void AtomMap::clearMap(void)
{
#ifndef MEM_OPT_VERSION	
  if (!cleared && localIDTable != NULL)
  {
    for(int i=0; i < tableSz; i++)
      localIDTable[i].pid = localIDTable[i].index = notUsed;
    cleared = true;
  }
#else
	for(int i=0; i<tableSz;i++){
		keys[i] = notUsed;
		localIDTable[i].pid = localIDTable[i].index = notUsed;
	}
	if(collidedAtoms) collidedAtoms->empty();
	if(isCollided) memset(isCollided, 0, sizeof(bool)*tableSz);
#endif  
}

void AtomMap::print()
{
#if !USE_STL_MAP	
  for (int i=0; i<tableSz; i++) {
    CkPrintf("AtomMap on node %d\n", CkMyPe());
    CkPrintf("AtomID %d -> PatchID %d:Index %d\n", i, localIDTable[i].pid,
      localIDTable[i].index);
  }
#else
	for(int i=0; i<tableSz; i++){
		if(keys[i]==notUsed) continue;
		CkPrintf("AtomMap on node %d ", CkMyPe());
		CkPrintf("AtomID %d -> PatchID %d:Index %d\n", i, localIDTable[i].pid,
			localIDTable[i].index);	
	}
	if(collidedAtoms == NULL) return;
	CkHashtableIterator *iter = collidedAtoms->iterator();
	LocalID *ent = NULL;
	int *aid = NULL;
	while((ent=(LocalID *)iter->next((void **)&aid))!=NULL){
		CkPrintf("AtomMap on node %d in hashtable", CkMyPe());
		if(*aid == 5852){
			printf("print from 5852\n");
		}
		CkPrintf("AtomID %d -> PatchID %d:Index %d\n", *aid, ent->pid, ent->index);
	}
	delete iter;
#endif  
}

#ifdef MEM_OPT_VERSION
LocalID AtomMap::localID(AtomID id)
{
	if(onlyUseTbl){
		return localIDTable[id];
	}else{
		int firstTry = firstHashCode(id);
		if(keys[firstTry] == id){
			return localIDTable[firstTry];
		}else if(keys[firstTry] != notUsed){
			//possibly in the hashtable
			LocalID *ent = collidedAtoms->getEntry(id);
			if(ent){
				return *ent;
			}else{
				LocalID one;
				one.pid = notUsed;
				one.index = notUsed;
				return one;				
			}
		}else {
			//not found
			LocalID one;
			one.pid = notUsed;
			one.index = notUsed;
			return one;				
		}
	}
}
#endif

