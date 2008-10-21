#ifndef COMPRESSPSF_H
#define COMPRESSPSF_H

#include "structures.h"
#include <string>
#include <deque>
using namespace std;

#define COMPRESSED_PSF_VER 1.3

class Molecule;
class Parameters;
class SimParameters;
class ConfigList;

void compress_psf_file(Molecule *mol, char *psfFileName, Parameters *param, SimParameters *simParam, ConfigList* cfgList);

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
