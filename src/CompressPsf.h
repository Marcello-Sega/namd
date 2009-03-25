#ifndef COMPRESSPSF_H
#define COMPRESSPSF_H

#include "structures.h"
#include <string>
#include <deque>
using namespace std;

#define COMPRESSED_PSF_VER 1.60

//used to detemine big-endian or little-endian for 
//the per-atom binary file
#define COMPRESSED_PSF_MAGICNUM 1234

#define BINARY_PERATOM_OUTPUT 1

class Molecule;
class Parameters;
class SimParameters;
class ConfigList;

void compress_psf_file(Molecule *mol, char *psfFileName, Parameters *param, SimParameters *simParam, ConfigList* cfgList);

void compress_molecule_info(Molecule *mol, char *psfFileName, Parameters *param, SimParameters *simParam, ConfigList* cfgList);

void flipNum(char *elem, int elemSize, int numElems);

template <typename T>
int lookupCstPool(const vector<T>& pool, const T& val)
{
    for(int i=0; i<pool.size(); i++)
    {
        if(pool[i]==val)
            return i;
    }
    return -1;
}
#endif
