#include <algorithm>
#include "CompressPsf.h"
#include "strlib.h"
#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "InfoStream.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"


/**
 * Very tricky part: 
 * 
 * Generating the tuple type affects the Parameter object.
 * Particularly, the "multiplicity" field in the TUPLE_array 
 * will be changed in the function assign_TUPLE_index. 
 * TUPLE=dihedral, improper.
 * 
 * And when we read the compressed psf files, assign_TUPLE_index
 * functions will not be called so that the value of "multiplicity"
 * will not be updated. This results in the incorrect energy from bonded
 * interaction.
 * 
 * Therefore, we need to store its value in the compressed psf
 * file. When it comes to read them, we need to update the Parameter object
 * with these values.
 */

//global variables recording compressed psf file information
Molecule *g_mol = NULL;
Parameters *g_param = NULL;
SimParameters *g_simParam = NULL;
ConfigList *g_cfgList = NULL;

struct BasicAtomInfo
{
    short segNameIdx;
    short resNameIdx;
    short atomNameIdx;
    short atomTypeIdx;
    short chargeIdx;
    short massIdx;
    int resID;

    int atomSigIdx;
    int exclSigIdx;
};

struct AtomSigInfo
{
    vector<short> bondSigIndices;
    vector<short> angleSigIndices;
    vector<short> dihedralSigIndices;
    vector<short> improperSigIndices;
    vector<short> crosstermSigIndices;

    AtomSigInfo()
    {}
    AtomSigInfo(const AtomSigInfo& sig)
    {
        bondSigIndices.clear();
        for(int i=0; i<sig.bondSigIndices.size(); i++)
            bondSigIndices.push_back(sig.bondSigIndices[i]);

        angleSigIndices.clear();
        for(int i=0; i<sig.angleSigIndices.size(); i++)
            angleSigIndices.push_back(sig.angleSigIndices[i]);

        dihedralSigIndices.clear();
        for(int i=0; i<sig.dihedralSigIndices.size(); i++)
            dihedralSigIndices.push_back(sig.dihedralSigIndices[i]);

        improperSigIndices.clear();
        for(int i=0; i<sig.improperSigIndices.size(); i++)
            improperSigIndices.push_back(sig.improperSigIndices[i]);

        crosstermSigIndices.clear();
        for(int i=0; i<sig.crosstermSigIndices.size(); i++)
            crosstermSigIndices.push_back(sig.crosstermSigIndices[i]);
    }

    ~AtomSigInfo()
    {
        bondSigIndices.clear();
        angleSigIndices.clear();
        dihedralSigIndices.clear();
        improperSigIndices.clear();
        crosstermSigIndices.clear();
    }

    void sortTupleSigIndices()
    {
        sort(bondSigIndices.begin(), bondSigIndices.end());
        sort(angleSigIndices.begin(), angleSigIndices.end());
        sort(dihedralSigIndices.begin(), dihedralSigIndices.end());
        sort(improperSigIndices.begin(), improperSigIndices.end());
        sort(crosstermSigIndices.begin(), crosstermSigIndices.end());
    }
};

int operator==(const AtomSigInfo &s1, const AtomSigInfo& s2)
{
    if(s1.bondSigIndices.size() != s2.bondSigIndices.size())
        return 0;
    if(s1.angleSigIndices.size() != s2.angleSigIndices.size())
        return 0;
    if(s1.dihedralSigIndices.size() != s2.dihedralSigIndices.size())
        return 0;
    if(s1.improperSigIndices.size() != s2.improperSigIndices.size())
        return 0;
    if(s1.crosstermSigIndices.size() != s2.crosstermSigIndices.size())
        return 0;

    int equalCnt;
    equalCnt=0;
    int bondSigCnt = s1.bondSigIndices.size();
    for(int i=0; i<bondSigCnt; i++)
        equalCnt += (s1.bondSigIndices[i]==s2.bondSigIndices[i]);
    if(equalCnt!=bondSigCnt)
        return 0;

    equalCnt=0;
    int angleSigCnt = s1.angleSigIndices.size();
    for(int i=0; i<angleSigCnt; i++)
        equalCnt += (s1.angleSigIndices[i]==s2.angleSigIndices[i]);
    if(equalCnt!=angleSigCnt)
        return 0;

    equalCnt=0;
    int dihedralSigCnt = s1.dihedralSigIndices.size();
    for(int i=0; i<dihedralSigCnt; i++)
        equalCnt += (s1.dihedralSigIndices[i]==s2.dihedralSigIndices[i]);
    if(equalCnt!=dihedralSigCnt)
        return 0;

    equalCnt=0;
    int improperSigCnt = s1.improperSigIndices.size();
    for(int i=0; i<improperSigCnt; i++)
        equalCnt += (s1.improperSigIndices[i]==s2.improperSigIndices[i]);
    if(equalCnt!=improperSigCnt)
        return 0;

    equalCnt=0;
    int crosstermSigCnt = s1.crosstermSigIndices.size();
    for(int i=0; i<crosstermSigCnt; i++)
        equalCnt += (s1.crosstermSigIndices[i]==s2.crosstermSigIndices[i]);
    if(equalCnt!=crosstermSigCnt)
        return 0;

    return 1;
}

struct ExclSigInfo
{
    vector<int> fullExclOffset;
    vector<int> modExclOffset;

    ExclSigInfo()
    {}
    ExclSigInfo(const ExclSigInfo& sig)
    {
        fullExclOffset.clear();
        for(int i=0; i<sig.fullExclOffset.size(); i++)
            fullExclOffset.push_back(sig.fullExclOffset[i]);

        modExclOffset.clear();
        for(int i=0; i<sig.modExclOffset.size(); i++)
            modExclOffset.push_back(sig.modExclOffset[i]);
    }

    ~ExclSigInfo()
    {
        fullExclOffset.clear();
        modExclOffset.clear();
    }

    void sortExclOffset()
    {
        sort(fullExclOffset.begin(), fullExclOffset.end());
        sort(modExclOffset.begin(), modExclOffset.end());
    }
};
int operator==(const ExclSigInfo &s1, const ExclSigInfo &s2)
{
    if(s1.fullExclOffset.size()!=s2.fullExclOffset.size())
        return 0;
    if(s1.modExclOffset.size()!=s2.modExclOffset.size())
        return 0;

    for(int i=0; i<s1.fullExclOffset.size(); i++)
    {
        if(s1.fullExclOffset[i] != s2.fullExclOffset[i])
            return 0;
    }

    for(int i=0; i<s1.modExclOffset.size(); i++)
    {
        if(s1.modExclOffset[i] != s2.modExclOffset[i])
            return 0;
    }
    return 1;
}


vector<string> segNamePool;
vector<string> resNamePool;
vector<string> atomNamePool;
vector<string> atomTypePool;
vector<Real> chargePool;
vector<Real> massPool;
vector<AtomSigInfo> atomSigPool;
BasicAtomInfo *atomData;

//Recording cluster information after reading all bonds info
int *eachAtomClusterID = NULL;
vector<int> eachClusterSize;
int g_numClusters = 0;
//indicate whether the atom ids in a cluster are contiguous or not
int g_isClusterContiguous = 0;

vector<TupleSignature> sigsOfBonds;
vector<TupleSignature> sigsOfAngles;
vector<TupleSignature> sigsOfDihedrals;
vector<TupleSignature> sigsOfImpropers;
vector<TupleSignature> sigsOfCrossterms;
AtomSigInfo *eachAtomSigs;

vector<ExclSigInfo> sigsOfExclusions;
ExclSigInfo *eachAtomExclSigs;

//Structures for extraBond options
vector<Bond> extraBonds;
vector<Angle> extraAngles;
vector<Dihedral> extraDihedrals;
vector<Improper> extraImpropers;
vector<BondValue> extraBondParams;
vector<AngleValue> extraAngleParams;
vector<DihedralValue> extraDihedralParams;
vector<ImproperValue> extraImproperParams;

int operator==(const BondValue &b1, const BondValue &b2)
{
    return (b1.k==b2.k) && (b1.x0==b2.x0);
}

int operator==(const AngleValue &a1, const AngleValue &a2)
{
    return (a1.k==a2.k) && (a1.k_ub==a2.k_ub) && (a1.r_ub==a2.r_ub) && (a1.theta0==a2.theta0);
}

int operator!=(const FourBodyConsts& f1, const FourBodyConsts& f2)
{
    return (f1.delta!=f2.delta) || (f1.k!=f2.k) || (f1.n!=f2.n);
}

int operator==(const DihedralValue &d1, const DihedralValue &d2)
{
    if(d1.multiplicity != d2.multiplicity)
        return 0;
    for(int i=0; i<MAX_MULTIPLICITY; i++)
    {
        if(d1.values[i] != d2.values[i])
            return 0;
    }
    return 1;
}

int operator==(const ImproperValue &d1, const ImproperValue &d2)
{
    if(d1.multiplicity != d2.multiplicity)
        return 0;
    for(int i=0; i<MAX_MULTIPLICITY; i++)
    {
        if(d1.values[i] != d2.values[i])
            return 0;
    }
    return 1;
}

void readPsfFile(char *psfFileName);
void integrateAllAtomSigs();
void outputPsfFile(FILE *ofp);

//reading extraBond's information
void getExtraBonds(StringList *file);

//reading atom's basic information
void getAtomData(FILE *ifp);
void getBondData(FILE *ifp);
void getAngleData(FILE *ifp);
void getDihedralData(FILE *ifp);
void getImproperData(FILE *ifp);
void getDonorData(FILE *ifp);
void getAcceptorData(FILE *ifp);
void getExclusionData(FILE *ifp);
void getCrosstermData(FILE *ifp);

//Functions related with building exclusions
void buildExclusions();
void build12Excls(UniqueSet<Exclusion>&, vector<int> *);
void build13Excls(UniqueSet<Exclusion>&, vector<int> *);
void build14Excls(UniqueSet<Exclusion>&, vector<int> *, int);

void clearGlobalVectors()
{
    segNamePool.clear();
    resNamePool.clear();
    atomNamePool.clear();
    atomTypePool.clear();
    chargePool.clear();
    massPool.clear();
    sigsOfBonds.clear();
    sigsOfAngles.clear();
    sigsOfDihedrals.clear();
    sigsOfImpropers.clear();
    sigsOfExclusions.clear();

    eachClusterSize.clear();
}

void compress_psf_file(Molecule *mol, char *psfFileName, Parameters *param, SimParameters *simParam, ConfigList *cfgList)
{

    g_mol = mol;
    g_param = param;
    g_simParam = simParam; //used for building exclusions
    g_cfgList = cfgList; //used for integrating extra bonds

    //read psf files
    readPsfFile(psfFileName);

    integrateAllAtomSigs();

    buildExclusions();


    char *outFileName = new char[strlen(psfFileName)+10];
    sprintf(outFileName, "%s.inter", psfFileName);
    FILE *ofp = fopen(outFileName, "w");
    delete [] outFileName;
    //output compressed psf file
    outputPsfFile(ofp);
    fclose(ofp);
}

//Almost same with Molecule::read_psf_file
void readPsfFile(char *fname)
{
    char err_msg[512];  //  Error message for NAMD_die
    char buffer[512];  //  Buffer for file reading
    int i;      //  Loop counter
    int NumTitle;    //  Number of Title lines in .psf file
    int ret_code;    //  ret_code from NAMD_read_line calls
    FILE *psf_file;
    Parameters *params = g_param;

    /* Try and open the .psf file           */
    if ( (psf_file = fopen(fname, "r")) == NULL)
    {
        sprintf(err_msg, "UNABLE TO OPEN .psf FILE %s", fname);
        NAMD_die(err_msg);
    }

    /*  Read till we have the first non-blank line of file    */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to see if we dropped out of the loop because of a     */
    /*  read error.  This shouldn't happen unless the file is empty */
    if (ret_code!=0)
    {
        sprintf(err_msg, "EMPTY .psf FILE %s", fname);
        NAMD_die(err_msg);
    }

    /*  The first non-blank line should contain the word "psf".    */
    /*  If we can't find it, die.               */
    if (!NAMD_find_word(buffer, "psf"))
    {
        sprintf(err_msg, "UNABLE TO FIND \"PSF\" STRING IN PSF FILE %s",
                fname);
        NAMD_die(err_msg);
    }

    /*  Read until we find the next non-blank line      */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to see if we dropped out of the loop because of a     */
    /*  read error.  This shouldn't happen unless there is nothing  */
    /*  but the PSF line in the file        */
    if (ret_code!=0)
    {
        sprintf(err_msg, "MISSING EVERYTHING BUT PSF FROM %s", fname);
        NAMD_die(err_msg);
    }

    /*  This line should have the word "NTITLE" in it specifying    */
    /*  how many title lines there are        */
    if (!NAMD_find_word(buffer, "NTITLE"))
    {
        sprintf(err_msg,"CAN NOT FIND \"NTITLE\" STRING IN PSF FILE %s",
                fname);
        NAMD_die(err_msg);
    }

    sscanf(buffer, "%d", &NumTitle);

    /*  Now skip the next NTITLE non-blank lines and then read in the*/
    /*  line which should contain NATOM        */
    i=0;

    while ( ((ret_code=NAMD_read_line(psf_file, buffer)) == 0) &&
            (i<NumTitle) )
    {
        if (!NAMD_blank_string(buffer))
            i++;
    }

    /*  Make sure we didn't exit because of a read error    */
    if (ret_code!=0)
    {
        sprintf(err_msg, "FOUND EOF INSTEAD OF NATOM IN PSF FILE %s",
                fname);
        NAMD_die(err_msg);
    }

    while (NAMD_blank_string(buffer))
    {
        NAMD_read_line(psf_file, buffer);
    }

    /*  Check to make sure we have the line we want      */
    if (!NAMD_find_word(buffer, "NATOM"))
    {
        sprintf(err_msg, "DIDN'T FIND \"NATOM\" IN PSF FILE %s",
                fname);
        NAMD_die(err_msg);
    }

    /*  Read in the number of atoms, and then the atoms themselves  */
    sscanf(buffer, "%d", &g_mol->numAtoms);

    getAtomData(psf_file);

    //read extra bonds/angles/dihedrals/impropers information first
    //and then integrate them into the following reading procedures
    if(g_simParam->extraBondsOn)
    {
        getExtraBonds(g_cfgList->find("extraBondsFile"));
    }

    //initialize eachAtomSigs
    eachAtomSigs = new AtomSigInfo[g_mol->numAtoms];

    /*  Read until we find the next non-blank line      */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0)
    {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NBONDS IN PSF");
    }

    /*  Look for the string "NBOND"          */
    if (!NAMD_find_word(buffer, "NBOND"))
    {
        NAMD_die("DID NOT FIND NBOND AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of bonds and then the bonds themselves  */
    sscanf(buffer, "%d", &g_mol->numBonds);

    if (g_mol->numBonds)
        getBondData(psf_file);

    /*  Read until we find the next non-blank line      */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0)
    {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NTHETA IN PSF");
    }

    /*  Look for the string "NTHETA"        */
    if (!NAMD_find_word(buffer, "NTHETA"))
    {
        NAMD_die("DID NOT FIND NTHETA AFTER BOND LIST IN PSF");
    }

    /*  Read in the number of angles and then the angles themselves */
    sscanf(buffer, "%d", &g_mol->numAngles);

    if (g_mol->numAngles)
        getAngleData(psf_file);

    /*  Read until we find the next non-blank line      */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0)
    {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NPHI IN PSF");
    }

    /*  Look for the string "NPHI"          */
    if (!NAMD_find_word(buffer, "NPHI"))
    {
        NAMD_die("DID NOT FIND NPHI AFTER ANGLE LIST IN PSF");
    }

    /*  Read in the number of dihedrals and then the dihedrals      */
    sscanf(buffer, "%d", &g_mol->numDihedrals);

    if (g_mol->numDihedrals)
        getDihedralData(psf_file);

    /*  Read until we find the next non-blank line      */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0)
    {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NIMPHI IN PSF");
    }

    /*  Look for the string "NIMPHI"        */
    if (!NAMD_find_word(buffer, "NIMPHI"))
    {
        NAMD_die("DID NOT FIND NIMPHI AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of Impropers and then the impropers  */
    sscanf(buffer, "%d", &g_mol->numImpropers);

    if (g_mol->numImpropers)
        getImproperData(psf_file);

    /*  Read until we find the next non-blank line      */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0)
    {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NDON IN PSF");
    }

    /*  Look for the string "NDON"        */
    if (!NAMD_find_word(buffer, "NDON"))
    {
        NAMD_die("DID NOT FIND NDON AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of hydrogen bond donors and then the donors */
    sscanf(buffer, "%d", &g_mol->numDonors);

    if (g_mol->numDonors)
        getDonorData(psf_file);

    /*  Read until we find the next non-blank line      */
    ret_code = NAMD_read_line(psf_file, buffer);

    while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
    {
        ret_code = NAMD_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0)
    {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NACC IN PSF");
    }

    /*  Look for the string "NACC"        */
    if (!NAMD_find_word(buffer, "NACC"))
    {
        NAMD_die("DID NOT FIND NACC AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of hydrogen bond donors and then the donors */
    sscanf(buffer, "%d", &g_mol->numAcceptors);

    if (g_mol->numAcceptors)
        getAcceptorData(psf_file);

    /*  look for the explicit non-bonded exclusion section.     */
    while (!NAMD_find_word(buffer, "NNB"))
    {
        ret_code = NAMD_read_line(psf_file, buffer);

        if (ret_code != 0)
        {
            NAMD_die("EOF ENCOUNTERED LOOKING FOR NNB IN PSF FILE");
        }
    }

    /*  Read in the number of exclusions and then the exclusions    */
    sscanf(buffer, "%d", &g_mol->numExclusions);

    if (g_mol->numExclusions)
        getExclusionData(psf_file);

    /*  look for the cross-term section.     */
    int crossterms_present = 1;
    while (!NAMD_find_word(buffer, "NCRTERM"))
    {
        ret_code = NAMD_read_line(psf_file, buffer);

        if (ret_code != 0)
        {
            // hit EOF before finding cross-term section
            crossterms_present = 0;
            break;
        }
    }

    if ( crossterms_present)
    {

        /*  Read in the number of cross-terms and then the cross-terms*/
        sscanf(buffer, "%d", &g_mol->numCrossterms);

        if (g_mol->numCrossterms)
            getCrosstermData(psf_file);

    }

    //  analyze the data and find the status of each atom
    //build_atom_status();

    fclose(psf_file);
}

void integrateAllAtomSigs()
{
    printf("Bond sigs:  %d\n", (int)sigsOfBonds.size());
    printf("Angle sigs:  %d\n", (int)sigsOfAngles.size());
    printf("Dihedral sigs:  %d\n", (int)sigsOfDihedrals.size());
    printf("Improper sigs:  %d\n", (int)sigsOfImpropers.size());
    printf("Crossterm sigs:  %d\n", (int)sigsOfCrossterms.size());


    for(int i=0; i<g_mol->numAtoms; i++)
    {
        eachAtomSigs[i].sortTupleSigIndices();
        int poolIndex = lookupCstPool(atomSigPool, eachAtomSigs[i]);
        if(poolIndex==-1)
        {
            atomSigPool.push_back(eachAtomSigs[i]);
            poolIndex = atomSigPool.size()-1;
        }
        atomData[i].atomSigIdx = poolIndex;
    }

    printf("Atom's sigs: %d\n", (int)atomSigPool.size());

    delete[] eachAtomSigs;
}

void outputPsfFile(FILE *ofp)
{
    fprintf(ofp, "FORMAT VERSION: %f\n", COMPRESSED_PSF_VER);

    fprintf(ofp, "%d !NSEGMENTNAMES\n", segNamePool.size());
    for(int i=0; i<segNamePool.size(); i++)
    {
        fprintf(ofp, "%s\n", segNamePool[i].c_str());
    }

    fprintf(ofp, "%d !NRESIDUENAMES\n", resNamePool.size());
    for(int i=0; i<resNamePool.size(); i++)
    {
        fprintf(ofp, "%s\n", resNamePool[i].c_str());
    }

    fprintf(ofp, "%d !NATOMNAMES\n", atomNamePool.size());
    for(int i=0; i<atomNamePool.size(); i++)
    {
        fprintf(ofp, "%s\n", atomNamePool[i].c_str());
    }

    fprintf(ofp, "%d !NATOMTYPES\n", atomTypePool.size());
    for(int i=0; i<atomTypePool.size(); i++)
    {
        fprintf(ofp, "%s\n", atomTypePool[i].c_str());
    }

    fprintf(ofp, "%d !NCHARGES\n", chargePool.size());
    for(int i=0; i<chargePool.size(); i++)
    {
        fprintf(ofp, "%f\n", chargePool[i]);
    }

    fprintf(ofp, "%d !NMASSES\n", massPool.size());
    for(int i=0; i<massPool.size(); i++)
    {
        fprintf(ofp, "%f\n", massPool[i]);
    }


    fprintf(ofp, "%d !NATOMSIGS\n", atomSigPool.size());
    for(int i=0; i<atomSigPool.size(); i++)
    {
        AtomSigInfo& oneAtomSig = atomSigPool[i];
        int oneTypeCnt = oneAtomSig.bondSigIndices.size();
        fprintf(ofp, "%d !%sSIGS\n", oneTypeCnt, "NBOND");
        for(int j=0; j<oneTypeCnt; j++)
        {
            short idx = oneAtomSig.bondSigIndices[j];
            TupleSignature& tSig = sigsOfBonds[idx];
            tSig.output(ofp);
        }

        oneTypeCnt = oneAtomSig.angleSigIndices.size();
        fprintf(ofp, "%d !%sSIGS\n", oneTypeCnt, "NTHETA");
        for(int j=0; j<oneTypeCnt; j++)
        {
            short idx = oneAtomSig.angleSigIndices[j];
            TupleSignature& tSig = sigsOfAngles[idx];
            tSig.output(ofp);
        }

        oneTypeCnt = oneAtomSig.dihedralSigIndices.size();
        fprintf(ofp, "%d !%sSIGS\n", oneTypeCnt, "NPHI");
        for(int j=0; j<oneTypeCnt; j++)
        {
            short idx = oneAtomSig.dihedralSigIndices[j];
            TupleSignature& tSig = sigsOfDihedrals[idx];
            tSig.output(ofp);
        }

        oneTypeCnt = oneAtomSig.improperSigIndices.size();
        fprintf(ofp, "%d !%sSIGS\n", oneTypeCnt, "NIMPHI");
        for(int j=0; j<oneTypeCnt; j++)
        {
            short idx = oneAtomSig.improperSigIndices[j];
            TupleSignature& tSig = sigsOfImpropers[idx];
            tSig.output(ofp);
        }

        oneTypeCnt = oneAtomSig.crosstermSigIndices.size();
        fprintf(ofp, "%d !%sSIGS\n", oneTypeCnt, "NCRTERM");
        for(int j=0; j<oneTypeCnt; j++)
        {
            short idx = oneAtomSig.crosstermSigIndices[j];
            TupleSignature& tSig = sigsOfCrossterms[idx];
            tSig.output(ofp);
        }
    }

    //2. Output exclusion signatures
    int exclSigCnt = sigsOfExclusions.size();
    fprintf(ofp, "%d !NEXCLSIGS\n", exclSigCnt);
    for(int i=0; i<exclSigCnt; i++)
    {
        ExclSigInfo *sig = &sigsOfExclusions[i];
        //first line is for full exclusions (1-2, 1-3) in the format of count offset1 offset2 offset3 ...
        fprintf(ofp, "%d", sig->fullExclOffset.size());
        for(int j=0; j<sig->fullExclOffset.size(); j++)
            fprintf(ofp, " %d", sig->fullExclOffset[j]);
        fprintf(ofp, "\n");

        //second line is for modified exclusions (1-4)
        fprintf(ofp, "%d", sig->modExclOffset.size());
        for(int j=0; j<sig->modExclOffset.size(); j++)
            fprintf(ofp, " %d", sig->modExclOffset[j]);
        fprintf(ofp, "\n");
    }

    //3. Output the cluster information
    fprintf(ofp, "%d !NCLUSTERS\n", g_numClusters);
    fprintf(ofp, "%d !CLUSTERCONTIGUITY\n", g_isClusterContiguous);

    //4. Output atom info
    fprintf(ofp, "%d !NATOM\n", g_mol->numAtoms);
    for(int i=0; i<g_mol->numAtoms; i++)
    {
        BasicAtomInfo &one = atomData[i];
        fprintf(ofp, "%d %d %d %d %d %d %d %d %d %d\n", one.segNameIdx, one.resID, one.resNameIdx, one.atomNameIdx,
                one.atomTypeIdx, one.chargeIdx, one.massIdx, one.atomSigIdx, one.exclSigIdx, eachAtomClusterID[i]);
    }
    //fprintf(ofp, "\n");

    delete[] atomData;
    delete[] eachAtomClusterID;

    //4.Output the parameter new values if extraBonds are present.
    if(g_simParam->extraBondsOn)
    {
        // 1) output extra params for bonds
        int numExtraParams = extraBondParams.size();
        fprintf(ofp, "%d!NEXTRABONDPARAMS\n", numExtraParams);
        if(numExtraParams>0)
        {
            vector<BondValue>::iterator bpIter;
            for(bpIter=extraBondParams.begin(); bpIter!=extraBondParams.end(); bpIter++)
            {
                fprintf(ofp, "%f %f\n", bpIter->k, bpIter->x0);
            }
        }

        // 2) output extra params for angles
        numExtraParams = extraAngleParams.size();
        fprintf(ofp, "%d!NEXTRAANGLEPARAMS\n", numExtraParams);
        if(numExtraParams>0)
        {
            NAMD_die("Output extra angle params not implemented!");
        }

        // 3) output extra params for dihedrals
        numExtraParams = extraDihedralParams.size();
        fprintf(ofp, "%d!NEXTRADIHEDRALPARAMS\n", numExtraParams);
        if(numExtraParams>0)
        {
            NAMD_die("Output extra dihedral params not implemented!");
        }

        // 4) output extra params for impropers
        numExtraParams = extraImproperParams.size();
        fprintf(ofp, "%d!NEXTRAIMPROPERPARAMS\n", numExtraParams);
        if(numExtraParams>0)
        {
            NAMD_die("Output extra improper params not implemented!");
        }
    }


    //5. Output the "multiplicity" field TUPLE_array of the Parameter object
    fprintf(ofp, "!DIHEDRALPARAMARRAY\n");
    for(int i=0; i<g_param->NumDihedralParams; i++)
    {
        fprintf(ofp, "%d ", g_param->dihedral_array[i].multiplicity);
    }
    fprintf(ofp, "\n");
    fprintf(ofp, "!IMPROPERPARAMARRAY\n");
    for(int i=0; i<g_param->NumImproperParams; i++)
    {
        fprintf(ofp, "%d ", g_param->improper_array[i].multiplicity);
    }
    fprintf(ofp, "\n");
}

void getAtomData(FILE *ifp)
{
    char buffer[512]; //Buffer for reading from file
    int atomID=0;
    int lastAtomID=0; //to ensure atoms are in order
    char segmentName[11];
    int residueID;
    char residueName[11];
    char atomName[11];
    char atomType[11];
    Real charge;
    Real mass;
    int fieldsCnt; //Number of fields read by sscanf

    int numAtoms = g_mol->numAtoms;

    //1. parse atom data to build constant pool (atom name, mass, charge etc.)
    atomData = new BasicAtomInfo[numAtoms];
    while(atomID < numAtoms)
    {
        NAMD_read_line(ifp, buffer);
        //If it is blank or a comment, skip it
        if(NAMD_blank_string(buffer) || buffer[0]=='!')
            continue;

        //Fields are arranged as follows:
        //atomID, segName, resID, resName, atomName, atomType, charge, mass
        fieldsCnt = sscanf(buffer, "%d %s %d %s %s %s %f %f",
                           &atomID, segmentName, &residueID, residueName,
                           atomName, atomType, &charge, &mass);

        //check to make sure we found what we were expecting
        if(fieldsCnt!=8)
        {
            char errMsg[128];
            sprintf(errMsg, "BAD ATOM LINE FORMAT IN PSF FILE FOR ATOM %d (%s)", lastAtomID+1, buffer);
            NAMD_die(errMsg);
        }

        // check if this is in XPLOR format
        int atomTypeNum;
        if(sscanf(atomType, "%d", &atomTypeNum)>0)
            NAMD_die("Structure (psf) file is either in CHARMM format (with numbers for atoms types, the X-PLOR format using names is required) or the segment name field is empty.");

        //make sure the atoms were in sequence
        if(atomID!=lastAtomID+1)
        {
            char errMsg[128];
            sprintf(errMsg, "ATOM NUMBERS OUT OF ORDER AT ATOM #%d IN FSF FILE",
                    lastAtomID+1);
            NAMD_die(errMsg);
        }

        lastAtomID++;

        //building constant pool
        int poolIndex;
        string fieldName;
        fieldName.assign(segmentName);
        poolIndex = lookupCstPool(segNamePool, fieldName);
        if(poolIndex==-1)
        {
            segNamePool.push_back(fieldName);
            poolIndex = segNamePool.size()-1;
        }
        atomData[atomID-1].segNameIdx = poolIndex;

        //poolIndex = lookupCstPool(resIDPool, residueID);
        //if(poolIndex==-1)
        //    resIDPool.push_back(residueID);
        atomData[atomID-1].resID = residueID;

        fieldName.assign(residueName);
        poolIndex = lookupCstPool(resNamePool, fieldName);
        if(poolIndex==-1)
        {
            resNamePool.push_back(fieldName);
            poolIndex = resNamePool.size()-1;
        }
        atomData[atomID-1].resNameIdx = poolIndex;

        fieldName.assign(atomName);
        poolIndex = lookupCstPool(atomNamePool, fieldName);
        if(poolIndex==-1)
        {
            atomNamePool.push_back(fieldName);
            poolIndex = atomNamePool.size()-1;
        }
        atomData[atomID-1].atomNameIdx = poolIndex;

        fieldName.assign(atomType);
        poolIndex = lookupCstPool(atomTypePool, fieldName);
        if(poolIndex==-1)
        {
            atomTypePool.push_back(fieldName);
            poolIndex = atomTypePool.size()-1;
        }
        atomData[atomID-1].atomTypeIdx = poolIndex;

        poolIndex = lookupCstPool(chargePool, charge);
        if(poolIndex==-1)
        {
            chargePool.push_back(charge);
            poolIndex = chargePool.size()-1;
        }
        atomData[atomID-1].chargeIdx = poolIndex;

        poolIndex = lookupCstPool(massPool, mass);
        if(poolIndex==-1)
        {
            massPool.push_back(mass);
            poolIndex = massPool.size()-1;
        }
        atomData[atomID-1].massIdx = poolIndex;
    }
}

void getExtraBonds(StringList *file)
{
    char err_msg[512];
    int a1, a2, a3, a4;
    float k, ref;

    int numAtoms = g_mol->numAtoms;

    if(!file)
    {
        NAMD_die("NO EXTRA BONDS FILES SPECIFIED");
    }

    for ( ; file; file = file->next )
    {  // loop over files
        FILE *f = fopen(file->data,"r");
        if ( ! f )
        {
            sprintf(err_msg, "UNABLE TO OPEN EXTRA BONDS FILE %s", file->data);
            NAMD_err(err_msg);
        }
        else
        {
            iout << iINFO << "READING EXTRA BONDS FILE " << file->data <<"\n"<<endi;
        }

        while ( 1 )
        {
            char buffer[512];
            int ret_code;
            int poolIndex;
            do
            {
                ret_code = NAMD_read_line(f, buffer);
            }
            while ( (ret_code==0) && (NAMD_blank_string(buffer)) );
            if (ret_code!=0)
                break;

            char type[512];
            sscanf(buffer,"%s",type);

#define CHECKATOMID(ATOMID) if ( ATOMID < 0 || ATOMID >= numAtoms ) badatom = 1;

            int badline = 0;
            int badatom = 0;
            if ( ! strncasecmp(type,"bond",4) )
            {
                if ( sscanf(buffer, "%s %d %d %f %f %s",
                            type, &a1, &a2, &k, &ref, err_msg) != 5 )
                    badline = 1;
                else
                {
                    CHECKATOMID(a1)
                    CHECKATOMID(a2)
                }
                Bond tmp;                
                tmp.atom1 = a1;
                tmp.atom2 = a2;

                BondValue tmpv;
                tmpv.k = k;
                tmpv.x0 = ref;

                poolIndex = lookupCstPool(extraBondParams, tmpv);
                if(poolIndex==-1){               
                    extraBondParams.push_back(tmpv);
                    poolIndex = extraBondParams.size()-1;
                }                 
                tmp.bond_type = g_param->NumBondParams + poolIndex;
                extraBonds.push_back(tmp);
            }
            else if ( ! strncasecmp(type,"angle",4) )
            {
                if ( sscanf(buffer, "%s %d %d %d %f %f %s",
                            type, &a1, &a2, &a3, &k, &ref, err_msg) != 6 )
                    badline = 1;
                else
                {
                    CHECKATOMID(a1)
                    CHECKATOMID(a2)
                    CHECKATOMID(a3)
                }
                Angle tmp;
                tmp.atom1 = a1;
                tmp.atom2 = a2;
                tmp.atom3 = a3;
                
                AngleValue tmpv;
                tmpv.k = k;
                tmpv.theta0 = ref / 180. * PI;
                tmpv.k_ub = 0;
                tmpv.r_ub = 0;

                poolIndex = lookupCstPool(extraAngleParams, tmpv);
                if(poolIndex==-1){               
                    extraAngleParams.push_back(tmpv);
                    poolIndex = extraAngleParams.size()-1;
                }                 
                tmp.angle_type = g_param->NumAngleParams + poolIndex;
                extraAngles.push_back(tmp);
            }
            else if ( ! strncasecmp(type,"dihedral",4) )
            {
                if ( sscanf(buffer, "%s %d %d %d %d %f %f %s",
                            type, &a1, &a2, &a3, &a4, &k, &ref, err_msg) != 7 )
                    badline = 1;
                else
                {
                    CHECKATOMID(a1)
                    CHECKATOMID(a2)
                    CHECKATOMID(a3)
                    CHECKATOMID(a4)
                }
                Dihedral tmp;
                tmp.atom1 = a1;
                tmp.atom2 = a2;
                tmp.atom3 = a3;
                tmp.atom4 = a4;
                
                DihedralValue tmpv;
                tmpv.multiplicity = 1;
                tmpv.values[0].n = 0;
                tmpv.values[0].k = k;
                tmpv.values[0].delta = ref / 180. * PI;

                poolIndex = lookupCstPool(extraDihedralParams, tmpv);
                if(poolIndex==-1){               
                    extraDihedralParams.push_back(tmpv);
                    poolIndex = extraDihedralParams.size()-1;
                }                 
                tmp.dihedral_type = g_param->NumDihedralParams + poolIndex;
                extraDihedrals.push_back(tmp);
            }
            else if ( ! strncasecmp(type,"improper",4) )
            {
                if ( sscanf(buffer, "%s %d %d %d %d %f %f %s",
                            type, &a1, &a2, &a3, &a4, &k, &ref, err_msg) != 7 )
                    badline = 1;
                else
                {
                    CHECKATOMID(a1)
                    CHECKATOMID(a2)
                    CHECKATOMID(a3)
                    CHECKATOMID(a4)
                }
                Improper tmp;
                tmp.atom1 = a1;
                tmp.atom2 = a2;
                tmp.atom3 = a3;
                tmp.atom4 = a4;
                
                ImproperValue tmpv;
                tmpv.multiplicity = 1;
                tmpv.values[0].n = 0;
                tmpv.values[0].k = k;
                tmpv.values[0].delta = ref / 180. * PI;

                poolIndex = lookupCstPool(extraImproperParams, tmpv);
                if(poolIndex==-1){               
                    extraImproperParams.push_back(tmpv);
                    poolIndex = extraImproperParams.size()-1;
                }                 
                tmp.improper_type = g_param->NumImproperParams + poolIndex;
                extraImpropers.push_back(tmp);
            }
            else if ( ! strncasecmp(type,"#",1) )
            {
                continue;  // comment
            }
            else
            {
                badline = 1;
            }
#undef CHECKATOMID
            if ( badline )
            {
                sprintf(err_msg, "BAD LINE IN EXTRA BONDS FILE %s: %s",
                        file->data, buffer);
                NAMD_die(err_msg);
            }
            if ( badatom )
            {
                sprintf(err_msg, "BAD ATOM ID IN EXTRA BONDS FILE %s: %s",
                        file->data, buffer);
                NAMD_die(err_msg);
            }
        }
        fclose(f);
    }  // loop over files
}


//In this function, safety check will also be performed:
//1. bond itself 2.duplicate bonds such as (a, b) & (b, a)
void getBondData(FILE *fd)
{

    //first read bonds
    int atom_nums[2];  // Atom indexes for the bonded atoms
    char atom1name[11];  // Atom type for atom #1
    char atom2name[11];  // Atom type for atom #2
    register int j;      // Loop counter
    int num_read=0;    // Number of bonds read so far

    int numBonds = g_mol->numBonds;
    int origNumBonds = numBonds;   // number of bonds in file header

    /*  Allocate the array to hold the bonds      */
    Bond *bonds=new Bond[numBonds + extraBonds.size()];

    if (bonds == NULL)
    {
        NAMD_die("memory allocations failed in Molecule::read_bonds");
    }

    /*  Loop through and read in all the bonds      */
    while (num_read < numBonds)
    {
        /*  Loop and read in the two atom indexes    */
        for (j=0; j<2; j++)
        {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=NAMD_read_int(fd, "BONDS")-1;

            /*  Check to make sure the index isn't too big  */
            if (atom_nums[j] >= g_mol->numAtoms)
            {
                char err_msg[128];

                sprintf(err_msg, "BOND INDEX %d GREATER THAN NATOM %d IN BOND # %d IN PSF FILE", atom_nums[j]+1, g_mol->numAtoms, num_read+1);
                NAMD_die(err_msg);
            }
        }

        if(atom_nums[0] == atom_nums[1])
        { //bond to itself
            char err_msg[128];
            sprintf(err_msg, "ATOM %d is bonded to itself!", atom_nums[0]+1);
            NAMD_die(err_msg);
        }

        /*  Get the atom type for the two atoms.  When we query */
        /*  the parameter object, we need to send the atom type */
        /*  that is alphabetically first as atom 1.    */
        const char *atom1Type = atomTypePool[atomData[atom_nums[0]].atomTypeIdx].c_str();
        const char *atom2Type = atomTypePool[atomData[atom_nums[1]].atomTypeIdx].c_str();
        if (strcasecmp(atom1Type,atom2Type) < 0)
        {
            strcpy(atom1name, atom1Type);
            strcpy(atom2name, atom2Type);
        }
        else
        {
            strcpy(atom2name, atom1Type);
            strcpy(atom1name, atom2Type);
        }

        /*  Query the parameter object for the constants for    */
        /*  this bond            */
        Bond *b = &(bonds[num_read]);
        g_param->assign_bond_index(atom1name, atom2name, b);

        /*  Assign the atom indexes to the array element  */
        b->atom1=atom_nums[0];
        b->atom2=atom_nums[1];

        /*  Make sure this isn't a fake bond meant for shake in x-plor.  */
        Real k, x0;
        g_param->get_bond_params(&k,&x0,b->bond_type);
        if ( k == 0. )
            --numBonds;  // fake bond
        else
            ++num_read;  // real bond
    }

    /*  Tell user about our subterfuge  */
    if ( numBonds != origNumBonds )
    {
        iout << iWARN << "Ignored " << origNumBonds - numBonds <<
        " bonds with zero force constants.\n" << endi;
        iout << iWARN <<
        "Will get H-H distance in rigid H2O from H-O-H angle.\n" << endi;
    }

    //copy extra bonds to bonds structure
    int numRealBonds = numBonds;
    for(int i=0; i<extraBonds.size(); i++)
        bonds[numBonds+i] = extraBonds[i];
    numBonds += extraBonds.size();
    extraBonds.clear();

    g_mol->numBonds = numBonds;

    //then creating bond's tupleSignature
    for(int i=0; i<numBonds; i++)
    {
        Bond *b = bonds+i;
        TupleSignature oneSig(1,BOND,b->bond_type);
        oneSig.offset[0] = b->atom2 - b->atom1;
        oneSig.isReal = (i<numRealBonds );

        int poolIndex = lookupCstPool(sigsOfBonds, oneSig);
        int newSig=0;
        if(poolIndex == -1)
        {
            sigsOfBonds.push_back(oneSig);
            poolIndex = (short)sigsOfBonds.size()-1;
            newSig=1;
        }

        if(!newSig)
        {//check duplicate bonds in the form of (a, b) && (a, b);
            int dupIdx = lookupCstPool(eachAtomSigs[b->atom1].bondSigIndices, (short)poolIndex);
            if(dupIdx!=-1)
            {
                char err_msg[128];
                sprintf(err_msg, "Duplicate bond %d-%d!", b->atom1+1, b->atom2+1);
                NAMD_die(err_msg);
            }
        }
        eachAtomSigs[b->atom1].bondSigIndices.push_back(poolIndex);
    }

    //check duplicate bonds in the form of (a, b) && (b, a)
    for(int i=0; i<numBonds; i++)
    {
        Bond *b=bonds+i;
        int atom2 = b->atom2;
        int thisOffset = atom2 - b->atom1;
        for(int j=0; j<eachAtomSigs[atom2].bondSigIndices.size(); j++)
        {
            short atom2BondId = eachAtomSigs[atom2].bondSigIndices[j];
            TupleSignature *secSig = &(sigsOfBonds[atom2BondId]);
            if(thisOffset== -(secSig->offset[0]))
            {
                char err_msg[128];
                sprintf(err_msg, "Duplicate bond %d-%d because two atoms are just reversed!", b->atom1+1, atom2+1);
                NAMD_die(err_msg);
            }
        }
    }      

    //building clusters for this simulation system in two steps
    //1. create a list for each atom where each atom in the list is bonded with that atom
    vector<int> *atomListOfBonded = new vector<int>[g_mol->numAtoms];

    for(int i=0; i<numRealBonds; i++)
    {
        Bond *b=bonds+i;
        int atom1 = b->atom1;
        int atom2 = b->atom2;
        atomListOfBonded[atom1].push_back(atom2);
        atomListOfBonded[atom2].push_back(atom1);
    }

    delete [] bonds;

    //2. using breadth-first-search to build the clusters. Here, we avoid recursive call
    // because the depth of calls may be of thousands which will blow up the stack, and
    //recursive call is slower than the stack-based BFS.
    //Considering such structure
    //1->1245; 7->1243; 1243->1245
    eachAtomClusterID = new int[g_mol->numAtoms];
    for(int i=0; i<g_mol->numAtoms; i++)
        eachAtomClusterID[i] = -1;

    for(int i=0; i<g_mol->numAtoms; i++)
    {
        int curClusterID=eachAtomClusterID[i];
        if(curClusterID==-1)
        {
            curClusterID=i;
        }

        deque<int> toVisitAtoms;
        toVisitAtoms.push_back(i);
        while(!toVisitAtoms.empty())
        {
            int visAtomID = toVisitAtoms.front();
            toVisitAtoms.pop_front();
            eachAtomClusterID[visAtomID] = curClusterID;
            for(int j=0; j<atomListOfBonded[visAtomID].size(); j++)
            {
                int otherAtom = atomListOfBonded[visAtomID][j];
                if(eachAtomClusterID[otherAtom]!=curClusterID)
                    toVisitAtoms.push_back(otherAtom);
            }
        }
    }

    //Now the clusterID of each atom should be usually in the non-decreasing order.
    //In other words, the atom ids of a cluster are generally contiguous.
    //If this is the case, the temporary memory usage of output IO during
    //the simulation can be dramatically reduced. So g_isClusterContiguous
    //is used to differentiate the two cases

    int curClusterID;
    int prevClusterID=eachAtomClusterID[0];
    int curClusterSize=1;
    g_isClusterContiguous = 1;
    for(int i=1; i<g_mol->numAtoms; i++)
    {
        curClusterID = eachAtomClusterID[i];
        if(curClusterID > prevClusterID)
        {
            eachClusterSize.push_back(curClusterSize);
            curClusterSize=1;
        }
        else if(curClusterID == prevClusterID)
        {
            curClusterSize++;
        }
        else
        { 
            g_isClusterContiguous = 0;
            break;
        }
        prevClusterID = curClusterID;
    }


    //Now iterate over the cluster size again to filter out the repeating cluster size.
    //There are two cases:
    //1. if the cluster id of atoms is monotonically increasing, the size of the cluster can be used
    //as this cluster's signature.
    //After this, eachAtomClusterID will store the cluster signature. Only the atom whose id is "clustersize"
    //will store the size. Others will be -1
    //2. if the cluster id of atoms is not monotonically increasing, in other words,
    //the atom ids of a cluster are not contiguous, then eachAtomClusterID remains
    //unchanged.

    if(g_isClusterContiguous) {
        eachClusterSize.push_back(curClusterSize); //record the last sequence of cluster
        int aid=0;
        for(int clusterIdx=0; clusterIdx<eachClusterSize.size(); clusterIdx++)
        {
            int curSize = eachClusterSize[clusterIdx];
            eachAtomClusterID[aid] = curSize;
            for(int i=aid+1; i<aid+curSize; i++)
                eachAtomClusterID[i] = -1;
            aid += curSize;
        }
        g_numClusters = eachClusterSize.size();
        eachClusterSize.clear();
    }else{
        eachClusterSize.clear();
        //reiterate over the atoms to figure out the number of clusters
        char *clusters = new char[g_mol->numAtoms];
        memset(clusters, 0, sizeof(char)*g_mol->numAtoms);
        for(int i=0; i<g_mol->numAtoms; i++) {
            clusters[eachAtomClusterID[i]] = 1;
        }
        
        for(int zeroPos=0; zeroPos<g_mol->numAtoms; zeroPos++) {
            if(clusters[zeroPos]==0) {
                //since the cluster ids are contiguous integers starting from 0,
                //if it becomes zero, we know the number of clusters is zeroPos
                g_numClusters = zeroPos;
                break;
            }
        }
        delete [] clusters;
    }

    //check whether cluster is built correctly
    /*printf("num clusters: %d\n", eachClusterSize.size());
    FILE *checkFile = fopen("cluster.opt", "w");
    for(int i=0; i<g_mol->numAtoms; i++)  fprintf(checkFile, "%d\n", eachAtomClusterID[i]);
    fclose(checkFile);*/

    
    for(int i=0; i<g_mol->numAtoms; i++)
        atomListOfBonded[i].clear();
    delete [] atomListOfBonded;
}

void getAngleData(FILE *fd)
{

    int atom_nums[3];  //  Atom numbers for the three atoms
    char atom1name[11];  //  Atom type for atom 1
    char atom2name[11];  //  Atom type for atom 2
    char atom3name[11];  //  Atom type for atom 3
    register int j;      //  Loop counter
    int num_read=0;    //  Number of angles read so far

    int numAngles = g_mol->numAngles;
    int origNumAngles = numAngles;  // Number of angles in file
    /*  Alloc the array of angles          */
    Angle *angles=new Angle[numAngles + extraAngles.size()];

    if (angles == NULL)
    {
        NAMD_die("memory allocation failed in Molecule::read_angles");
    }

    /*  Loop through and read all the angles      */
    while (num_read < numAngles)
    {
        /*  Loop through the 3 atom indexes in the current angle*/
        for (j=0; j<3; j++)
        {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=NAMD_read_int(fd, "ANGLES")-1;

            /*  Check to make sure the atom index doesn't   */
            /*  exceed the Number of Atoms      */
            if (atom_nums[j] >= g_mol->numAtoms)
            {
                char err_msg[128];

                sprintf(err_msg, "ANGLES INDEX %d GREATER THAN NATOM %d IN ANGLES # %d IN PSF FILE", atom_nums[j]+1, g_mol->numAtoms, num_read+1);
                NAMD_die(err_msg);
            }
        }

        /*  Place the bond name that is alphabetically first  */
        /*  in the atom1name.  This is OK since the order of    */
        /*  atom1 and atom3 are interchangable.  And to search  */
        /*  the tree of angle parameters, we need the order     */
        /*  to be predictable.          */

        const char *atom1Type = atomTypePool[atomData[atom_nums[0]].atomTypeIdx].c_str();
        const char *atom2Type = atomTypePool[atomData[atom_nums[1]].atomTypeIdx].c_str();
        const char *atom3Type = atomTypePool[atomData[atom_nums[2]].atomTypeIdx].c_str();

        if (strcasecmp(atom1Type,atom2Type) < 0)
        {
            strcpy(atom1name, atom1Type);
            strcpy(atom2name, atom2Type);
            strcpy(atom3name, atom3Type);
        }
        else
        {
            strcpy(atom1name, atom1Type);
            strcpy(atom2name, atom2Type);
            strcpy(atom3name, atom3Type);
        }

        /*  Get the constant values for this bond from the  */
        /*  parameter object          */
        g_param->assign_angle_index(atom1name, atom2name,
                                    atom3name, &(angles[num_read]));

        /*  Assign the three atom indices      */
        angles[num_read].atom1=atom_nums[0];
        angles[num_read].atom2=atom_nums[1];
        angles[num_read].atom3=atom_nums[2];

        /*  Make sure this isn't a fake angle meant for shake in x-plor.  */
        Real k, t0, k_ub, r_ub;
        g_param->get_angle_params(&k,&t0,&k_ub,&r_ub,angles[num_read].angle_type);
        if ( k == 0. && k_ub == 0. )
            --numAngles;  // fake angle
        else
            ++num_read;  // real angle
    }

    /*  Tell user about our subterfuge  */
    if ( numAngles != origNumAngles )
    {
        iout << iWARN << "Ignored " << origNumAngles - numAngles <<
        " angles with zero force constants.\n" << endi;
    }

    //copy extra angles to angles structure
    for(int i=0; i<extraAngles.size(); i++)
        angles[numAngles+i] = extraAngles[i];
    numAngles += extraAngles.size();
    extraAngles.clear();

    g_mol->numAngles = numAngles;

    //create angles' tupleSignature
    for(int i=0; i<numAngles; i++)
    {
        Angle *tuple = angles+i;
        TupleSignature oneSig(2,ANGLE,tuple->angle_type);
        int offset[2];
        offset[0] = tuple->atom2 - tuple->atom1;
        offset[1] = tuple->atom3 - tuple->atom1;
        oneSig.setOffsets(offset);

        int poolIndex = lookupCstPool(sigsOfAngles, oneSig);
        if(poolIndex == -1)
        {
            sigsOfAngles.push_back(oneSig);
            poolIndex = (short)sigsOfAngles.size()-1;
        }
        eachAtomSigs[tuple->atom1].angleSigIndices.push_back(poolIndex);
    }

    delete [] angles;

}

void getDihedralData(FILE *fd)
{

    int atom_nums[4];  // The 4 atom indexes
    int last_atom_nums[4];  // Atom numbers from previous bond
    char atom1name[11];  // Atom type for atom 1
    char atom2name[11];  // Atom type for atom 2
    char atom3name[11];  // Atom type for atom 3
    char atom4name[11];  // Atom type for atom 4
    register int j;      // loop counter
    int num_read=0;    // number of dihedrals read so far
    int multiplicity=1;  // multiplicity of the current bond
    Bool duplicate_bond;  // Is this a duplicate of the last bond
    int num_unique=0;   // Number of unique dihedral bonds

    //  Initialize the array used to check for duplicate dihedrals
    for (j=0; j<4; j++)
        last_atom_nums[j] = -1;

    int numDihedrals = g_mol->numDihedrals;
    int numAtoms = g_mol->numAtoms;

    /*  Allocate an array to hold the Dihedrals      */
    Dihedral *dihedrals = new Dihedral[numDihedrals + extraDihedrals.size()];

    if (dihedrals == NULL)
    {
        NAMD_die("memory allocation failed in Molecule::read_dihedrals");
    }

    /*  Loop through and read all the dihedrals      */
    while (num_read < numDihedrals)
    {
        duplicate_bond = TRUE;

        /*  Loop through and read the 4 indexes for this bond   */
        for (j=0; j<4; j++)
        {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=NAMD_read_int(fd, "DIHEDRALS")-1;

            /*  Check for an atom index that is too large  */
            if (atom_nums[j] >= numAtoms)
            {
                char err_msg[128];

                sprintf(err_msg, "DIHEDRALS INDEX %d GREATER THAN NATOM %d IN DIHEDRALS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
                NAMD_die(err_msg);
            }

            //  Check to see if this atom matches the last bond
            if (atom_nums[j] != last_atom_nums[j])
            {
                duplicate_bond = FALSE;
            }

            last_atom_nums[j] = atom_nums[j];
        }

        /*  Get the atom types for the 4 atoms so we can look  */
        /*  up the constants in the parameter object    */
        const char *atom1Type = atomTypePool[atomData[atom_nums[0]].atomTypeIdx].c_str();
        const char *atom2Type = atomTypePool[atomData[atom_nums[1]].atomTypeIdx].c_str();
        const char *atom3Type = atomTypePool[atomData[atom_nums[2]].atomTypeIdx].c_str();
        const char *atom4Type = atomTypePool[atomData[atom_nums[3]].atomTypeIdx].c_str();
        strcpy(atom1name, atom1Type);
        strcpy(atom2name, atom2Type);
        strcpy(atom3name, atom3Type);
        strcpy(atom4name, atom4Type);

        //  Check to see if this is really a new bond or just
        //  a repeat of the last one
        if (duplicate_bond)
        {
            //  This is a duplicate, so increase the multiplicity
            multiplicity++;

            if (multiplicity == 2)
            {
                g_mol->numMultipleDihedrals++;
            }
        }
        else
        {
            multiplicity=1;
            num_unique++;
        }

        /*  Get the constants for this dihedral bond    */
        g_param->assign_dihedral_index(atom1name, atom2name,
                                       atom3name, atom4name, &(dihedrals[num_unique-1]),
                                       multiplicity);

        /*  Assign the atom indexes        */
        dihedrals[num_unique-1].atom1=atom_nums[0];
        dihedrals[num_unique-1].atom2=atom_nums[1];
        dihedrals[num_unique-1].atom3=atom_nums[2];
        dihedrals[num_unique-1].atom4=atom_nums[3];

        num_read++;
    }

    //copy extra dihedrals to dihedrals structure
    if(extraDihedrals.size()>0){
        iout << iWARN << "It's your responsibility to ensure there is no duplicates among these extra dihedrals\n" << endi;
    }
    for(int i=0; i<extraDihedrals.size(); i++)
        dihedrals[num_unique+i] = extraDihedrals[i];
    num_unique += extraDihedrals.size();
    extraDihedrals.clear();

    g_mol->numDihedrals = num_unique;

    //create dihedrals' tupleSignature
    for(int i=0; i<num_unique; i++)
    {
        Dihedral *tuple = dihedrals+i;
        TupleSignature oneSig(3,DIHEDRAL,tuple->dihedral_type);
        int offset[3];
        offset[0] = tuple->atom2 - tuple->atom1;
        offset[1] = tuple->atom3 - tuple->atom1;
        offset[2] = tuple->atom4 - tuple->atom1;
        oneSig.setOffsets(offset);

        int poolIndex = lookupCstPool(sigsOfDihedrals, oneSig);
        if(poolIndex == -1)
        {
            sigsOfDihedrals.push_back(oneSig);
            poolIndex = (short)sigsOfDihedrals.size()-1;
        }
        eachAtomSigs[tuple->atom1].dihedralSigIndices.push_back(poolIndex);
    }

    delete[] dihedrals;

}

void getImproperData(FILE *fd)
{
    int atom_nums[4];  //  Atom indexes for the 4 atoms
    int last_atom_nums[4];  //  Atom indexes from previous bond
    char atom1name[11];  //  Atom type for atom 1
    char atom2name[11];  //  Atom type for atom 2
    char atom3name[11];  //  Atom type for atom 3
    char atom4name[11];  //  Atom type for atom 4
    register int j;      //  Loop counter
    int num_read=0;    //  Number of impropers read so far
    int multiplicity=1;  // multiplicity of the current bond
    Bool duplicate_bond;  // Is this a duplicate of the last bond
    int num_unique=0;   // Number of unique dihedral bonds

    //  Initialize the array used to look for duplicate improper
    //  entries.  Set them all to -1 so we know nothing will match
    for (j=0; j<4; j++)
        last_atom_nums[j] = -1;

    int numImpropers = g_mol->numImpropers;
    int numAtoms = g_mol->numAtoms;

    /*  Allocate the array to hold the impropers      */
    Improper *impropers=new Improper[numImpropers];

    if (impropers == NULL)
    {
        NAMD_die("memory allocation failed in Molecule::read_impropers");
    }

    /*  Loop through and read all the impropers      */
    while (num_read < numImpropers)
    {
        duplicate_bond = TRUE;

        /*  Loop through the 4 indexes for this improper  */
        for (j=0; j<4; j++)
        {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=NAMD_read_int(fd, "IMPROPERS")-1;

            /*  Check to make sure the index isn't too big  */
            if (atom_nums[j] >= numAtoms)
            {
                char err_msg[128];

                sprintf(err_msg, "IMPROPERS INDEX %d GREATER THAN NATOM %d IN IMPROPERS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
                NAMD_die(err_msg);
            }

            if (atom_nums[j] != last_atom_nums[j])
            {
                duplicate_bond = FALSE;
            }

            last_atom_nums[j] = atom_nums[j];
        }

        /*  Get the atom types so we can look up the parameters */
        const char *atom1Type = atomTypePool[atomData[atom_nums[0]].atomTypeIdx].c_str();
        const char *atom2Type = atomTypePool[atomData[atom_nums[1]].atomTypeIdx].c_str();
        const char *atom3Type = atomTypePool[atomData[atom_nums[2]].atomTypeIdx].c_str();
        const char *atom4Type = atomTypePool[atomData[atom_nums[3]].atomTypeIdx].c_str();
        strcpy(atom1name, atom1Type);
        strcpy(atom2name, atom2Type);
        strcpy(atom3name, atom3Type);
        strcpy(atom4name, atom4Type);

        //  Check to see if this is a duplicate improper
        if (duplicate_bond)
        {
            //  This is a duplicate improper.  So we don't
            //  really count this entry, we just update
            //  the parameters object
            multiplicity++;

            if (multiplicity == 2)
            {
                //  Count the number of multiples.
                g_mol->numMultipleImpropers++;
            }
        }
        else
        {
            //  Not a duplicate
            multiplicity = 1;
            num_unique++;
        }

        /*  Look up the constants for this bond      */
        g_param->assign_improper_index(atom1name, atom2name,
                                       atom3name, atom4name, &(impropers[num_unique-1]),
                                       multiplicity);

        /*  Assign the atom indexes        */
        impropers[num_unique-1].atom1=atom_nums[0];
        impropers[num_unique-1].atom2=atom_nums[1];
        impropers[num_unique-1].atom3=atom_nums[2];
        impropers[num_unique-1].atom4=atom_nums[3];

        num_read++;
    }

    //copy extra impropers to impropers structure
    if(extraImpropers.size()>0){
        iout << iWARN << "It's your responsibility to ensure there is no duplicates among these extra impropers\n" << endi;
    }
    for(int i=0; i<extraImpropers.size(); i++)
        impropers[num_unique+i] = extraImpropers[i];
    num_unique += extraImpropers.size();
    extraImpropers.clear();

    //  Now reset the numImpropers value to the number of UNIQUE
    //  impropers.  Sure, we waste a few entries in the improper_array
    //  on the master node, but it is very little space . . .
    numImpropers = num_unique;

    //create improper's tupleSignature
    for(int i=0; i<num_unique; i++)
    {
        Improper *tuple = impropers+i;
        TupleSignature oneSig(3,IMPROPER,tuple->improper_type);
        int offset[3];
        offset[0] = tuple->atom2 - tuple->atom1;
        offset[1] = tuple->atom3 - tuple->atom1;
        offset[2] = tuple->atom4 - tuple->atom1;
        oneSig.setOffsets(offset);

        int poolIndex = lookupCstPool(sigsOfImpropers, oneSig);
        if(poolIndex == -1)
        {
            sigsOfImpropers.push_back(oneSig);
            poolIndex = (short)sigsOfImpropers.size()-1;
        }
        eachAtomSigs[tuple->atom1].improperSigIndices.push_back(poolIndex);
    }

    delete[] impropers;

}

void getDonorData(FILE *fd)
{
    int d[2];               // temporary storage of donor atom index
    register int j;      // Loop counter
    int num_read=0;    // Number of bonds read so far
    int num_no_hydr=0;      // Number of bonds with no hydrogen given

    /*  Allocate the array to hold the bonds      */
    int numDonors = g_mol->numDonors;
    int numAtoms = g_mol->numAtoms;
    Bond *donors=new Bond[numDonors];

    if (donors == NULL)
    {
        NAMD_die("memory allocations failed in Molecule::read_donors");
    }

    /*  Loop through and read in all the donors      */
    while (num_read < numDonors)
    {
        /*  Loop and read in the two atom indexes    */
        for (j=0; j<2; j++)
        {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            d[j]=NAMD_read_int(fd, "DONORS")-1;

            /*  Check to make sure the index isn't too big  */
            if (d[j] >= numAtoms)
            {
                char err_msg[128];

                sprintf(err_msg,
                        "DONOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE",
                        d[j]+1, numAtoms, num_read+1);
                NAMD_die(err_msg);
            }

            /*  Check if there is a hydrogen given */
            if (d[j] < 0)
                num_no_hydr++;
        }

        /*  Assign the atom indexes to the array element  */
        Bond *b = &(donors[num_read]);
        b->atom1=d[0];
        b->atom2=d[1];

        num_read++;
    }

    delete [] donors;
}

void getAcceptorData(FILE *fd)
{
    int d[2];               // temporary storage of atom index
    register int j;      // Loop counter
    int num_read=0;    // Number of bonds read so far
    int num_no_ante=0;      // number of pairs with no antecedent

    int numAcceptors = g_mol->numAcceptors;

    /*  Allocate the array to hold the bonds      */
    Bond *acceptors=new Bond[numAcceptors];

    if (acceptors == NULL)
    {
        NAMD_die("memory allocations failed in Molecule::read_acceptors");
    }

    /*  Loop through and read in all the acceptors      */
    while (num_read < numAcceptors)
    {
        /*  Loop and read in the two atom indexes    */
        for (j=0; j<2; j++)
        {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            d[j]=NAMD_read_int(fd, "ACCEPTORS")-1;

            /*  Check to make sure the index isn't too big  */
            if (d[j] >= g_mol->numAtoms)
            {
                char err_msg[128];

                sprintf(err_msg, "ACCEPTOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE", d[j]+1, g_mol->numAtoms, num_read+1);
                NAMD_die(err_msg);
            }

            /*  Check if there is an antecedent given */
            if (d[j] < 0)
                num_no_ante++;
        }

        /*  Assign the atom indexes to the array element  */
        Bond *b = &(acceptors[num_read]);
        b->atom1=d[0];
        b->atom2=d[1];

        num_read++;
    }

    delete [] acceptors;
}

void getExclusionData(FILE *fd)
{
    //reading explicit exclusions from PSF file
    //TODO: Implement it
    //currently just abort saying it is not supported
    printf("ERROR: The current compression doesn't support explicit exclusions!\n");
    NAMD_die("Compressing .psf file is not finished!\n");  
}

void getCrosstermData(FILE *fd)
{
  int atom_nums[8];  //  Atom indexes for the 4 atoms
  int last_atom_nums[8];  //  Atom indexes from previous bond
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  char atom4name[11];  //  Atom type for atom 4
  char atom5name[11];  //  Atom type for atom 5
  char atom6name[11];  //  Atom type for atom 6
  char atom7name[11];  //  Atom type for atom 7
  char atom8name[11];  //  Atom type for atom 8
  int j;      //  Loop counter
  int num_read=0;    //  Number of items read so far
  Bool duplicate_bond;  // Is this a duplicate of the last bond

  //  Initialize the array used to look for duplicate crossterm
  //  entries.  Set them all to -1 so we know nothing will match
  for (j=0; j<8; j++)
    last_atom_nums[j] = -1;

  int numCrossterms = g_mol->numCrossterms;
  int numAtoms = g_mol->numAtoms;

  /*  Allocate the array to hold the cross-terms */
  Crossterm* crossterms=new Crossterm[numCrossterms];

  if (crossterms == NULL)
  {
    NAMD_die("memory allocation failed in getCrosstermData when compressing psf file");
  }

  /*  Loop through and read all the cross-terms      */
  while (num_read < numCrossterms)
  {
    duplicate_bond = TRUE;

    /*  Loop through the 8 indexes for this cross-term */
    for (j=0; j<8; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "CROSS-TERMS")-1;

      /*  Check to make sure the index isn't too big  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "CROSS-TERM INDEX %d GREATER THAN NATOM %d IN CROSS-TERMS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }
      
      if (atom_nums[j] != last_atom_nums[j]){
        duplicate_bond = FALSE;
      }
      
      last_atom_nums[j] = atom_nums[j];
    }
    
    /*  Get the atom types so we can look up the parameters */
    const char *atom1Type = atomTypePool[atomData[atom_nums[0]].atomTypeIdx].c_str();
    const char *atom2Type = atomTypePool[atomData[atom_nums[1]].atomTypeIdx].c_str();
    const char *atom3Type = atomTypePool[atomData[atom_nums[2]].atomTypeIdx].c_str();
    const char *atom4Type = atomTypePool[atomData[atom_nums[3]].atomTypeIdx].c_str();
    const char *atom5Type = atomTypePool[atomData[atom_nums[4]].atomTypeIdx].c_str();
    const char *atom6Type = atomTypePool[atomData[atom_nums[5]].atomTypeIdx].c_str();
    const char *atom7Type = atomTypePool[atomData[atom_nums[6]].atomTypeIdx].c_str();
    const char *atom8Type = atomTypePool[atomData[atom_nums[7]].atomTypeIdx].c_str();
    strcpy(atom1name, atom1Type);
    strcpy(atom2name, atom2Type);
    strcpy(atom3name, atom3Type);
    strcpy(atom4name, atom4Type);
    strcpy(atom5name, atom5Type);
    strcpy(atom6name, atom6Type);
    strcpy(atom7name, atom7Type);
    strcpy(atom8name, atom8Type);

    //  Check to see if this is a duplicate term
    if (duplicate_bond)
    {
      iout << iWARN << "Duplicate cross-term detected.\n" << endi;
    }

    /*  Look up the constants for this bond      */
    g_param->assign_crossterm_index(atom1name, atom2name, 
       atom3name, atom4name, atom5name, atom6name,
       atom7name, atom8name, &(crossterms[num_read]));

    /*  Assign the atom indexes        */
    crossterms[num_read].atom1=atom_nums[0];
    crossterms[num_read].atom2=atom_nums[1];
    crossterms[num_read].atom3=atom_nums[2];
    crossterms[num_read].atom4=atom_nums[3];
    crossterms[num_read].atom5=atom_nums[4];
    crossterms[num_read].atom6=atom_nums[5];
    crossterms[num_read].atom7=atom_nums[6];
    crossterms[num_read].atom8=atom_nums[7];

    num_read++;
  }

  numCrossterms = num_read;

  //create crossterm's tupleSignature
  for(int i=0; i<numCrossterms; i++)
  {
    Crossterm *tuple = crossterms+i;
    TupleSignature oneSig(7, CROSSTERM, tuple->crossterm_type);
    int offset[7];
    offset[0] = tuple->atom2 - tuple->atom1;
    offset[1] = tuple->atom3 - tuple->atom1;
    offset[2] = tuple->atom4 - tuple->atom1;
    offset[3] = tuple->atom5 - tuple->atom1;
    offset[4] = tuple->atom6 - tuple->atom1;
    offset[5] = tuple->atom7 - tuple->atom1;
    offset[6] = tuple->atom8 - tuple->atom1;
    oneSig.setOffsets(offset);
   
    int poolIndex = lookupCstPool(sigsOfCrossterms, oneSig);
    if(poolIndex == -1)
    {
      sigsOfCrossterms.push_back(oneSig);
      poolIndex = (short)sigsOfCrossterms.size()-1;
    }
    eachAtomSigs[tuple->atom1].crosstermSigIndices.push_back(poolIndex);
  }
  
  delete[] crossterms;
}

void buildExclusions()
{
    //1. Build exclusions: mainly accomplish the function of
    //Molecule::build_exclusions (based on the bonds)
    UniqueSet<Exclusion> allExclusions;

    int exclude_flag; //Exclusion policy
    exclude_flag = g_simParam->exclude;
    //int stripHGroupExclFlag = (simParams->splitPatch == SPLIT_PATCH_HYDROGEN);

    //Commented now since no explicit exclusions are read
    //  Go through the explicit exclusions and add them to the arrays
    //for(i=0; i<numExclusions; i++){
    //	exclusionSet.add(exclusions[i]);
    //}

    // If this is AMBER force field, and readExclusions is TRUE,
    // then all the exclusions were read from parm file, and we
    // shouldn't generate any of them.
    // Comment on stripHGroupExcl:
    // 1. Inside this function, hydrogenGroup is initialized in
    // build_atom_status, therefore, not available when reading psf files
    // 2. this function's main purpose is to reduce memory usage. Since exclusion
    // signatures are used, this function could be overlooked  --Chao Mei

    vector<int> *eachAtomNeighbors = new vector<int>[g_mol->numAtoms];   
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        AtomSigInfo *aSig = &atomSigPool[atomData[atom1].atomSigIdx];
        for(int j=0; j<aSig->bondSigIndices.size(); j++)
        {
            TupleSignature *tSig = &sigsOfBonds[aSig->bondSigIndices[j]];
            if(!tSig->isReal) continue;
            int atom2 = atom1+tSig->offset[0];
            eachAtomNeighbors[atom1].push_back(atom2);
            eachAtomNeighbors[atom2].push_back(atom1);
        }
    }

    if (!g_simParam->amberOn || !g_simParam->readExclusions)
    { //  Now calculate the bonded exlcusions based on the exclusion policy
        switch (exclude_flag)
        {
        case NONE:
            break;
        case ONETWO:
            build12Excls(allExclusions, eachAtomNeighbors);
            break;
        case ONETHREE:
            build12Excls(allExclusions, eachAtomNeighbors);
            build13Excls(allExclusions, eachAtomNeighbors);
            //if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        case ONEFOUR:
            build12Excls(allExclusions, eachAtomNeighbors);
            build13Excls(allExclusions, eachAtomNeighbors);
            build14Excls(allExclusions, eachAtomNeighbors, 0);
            //if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        case SCALED14:
            build12Excls(allExclusions, eachAtomNeighbors);
            build13Excls(allExclusions, eachAtomNeighbors);
            build14Excls(allExclusions, eachAtomNeighbors, 1);
            //if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        }
    }
    //else if (stripHGroupExclFlag && exclude_flag!=NONE && exclude_flag!=ONETWO)
    //  stripHGroupExcl();

    //Commented since atomFepFlags information is not available when reading psf file
    //stripFepExcl();

    for(int i=0; i<g_mol->numAtoms; i++)
        eachAtomNeighbors[i].clear();
    delete [] eachAtomNeighbors;

    //2. Build each atom's list of exclusions
    UniqueSetIter<Exclusion> exclIter(allExclusions);
    eachAtomExclSigs = new ExclSigInfo[g_mol->numAtoms];
    for(exclIter=exclIter.begin(); exclIter!=exclIter.end(); exclIter++)
    {
        int atom1 = exclIter->atom1;
        int atom2 = exclIter->atom2;
        int offset21 = atom2-atom1;
        if(exclIter->modified)
        {
            eachAtomExclSigs[atom1].modExclOffset.push_back(offset21);
            eachAtomExclSigs[atom2].modExclOffset.push_back(-offset21);
        }
        else
        {
            eachAtomExclSigs[atom1].fullExclOffset.push_back(offset21);
            eachAtomExclSigs[atom2].fullExclOffset.push_back(-offset21);
        }
    }
    allExclusions.clear();

    //3. Build up exclusion signatures and determine each atom's
    //exclusion signature index
    for(int i=0; i<g_mol->numAtoms; i++)
    {
        eachAtomExclSigs[i].sortExclOffset();
        int poolIndex = lookupCstPool(sigsOfExclusions, eachAtomExclSigs[i]);
        if(poolIndex==-1)
        {
            poolIndex = sigsOfExclusions.size();
            sigsOfExclusions.push_back(eachAtomExclSigs[i]);
        }
        atomData[i].exclSigIdx = poolIndex;
    }
    delete [] eachAtomExclSigs;
    eachAtomExclSigs = NULL;
    printf("Exclusion signatures: %d\n", (int)sigsOfExclusions.size());
}

void build12Excls(UniqueSet<Exclusion>& allExcls, vector<int> *eachAtomNeighbors)
{
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        vector<int> *atom1List = &eachAtomNeighbors[atom1];
        for(int j=0; j<atom1List->size(); j++)
        {
            int atom2 = atom1List->at(j);
            if(atom1<atom2)
                allExcls.add(Exclusion(atom1, atom2));
            else
                allExcls.add(Exclusion(atom2, atom1));
        }
    }
}

void build13Excls(UniqueSet<Exclusion>& allExcls, vector<int> *eachAtomNeighbors)
{
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        vector<int> *atom1List = &eachAtomNeighbors[atom1];
        for(int j=0; j<atom1List->size(); j++)
        {
            int atom2 = atom1List->at(j);
            vector<int> *atom2List = &eachAtomNeighbors[atom2];
            for(int k=0; k<atom2List->size(); k++)
            {
                int atom3 = atom2List->at(k);
                //atom1-atom2, so atom2List contains atom1 which should not be considered
                if(atom3 == atom1)
                    continue;
                if(atom1<atom3)
                    allExcls.add(Exclusion(atom1, atom3));
                else
                    allExcls.add(Exclusion(atom3, atom1));
            }
        }
    }
}

void build14Excls(UniqueSet<Exclusion>& allExcls, vector<int> *eachAtomNeighbors, int modified)
{
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        vector<int> *atom1List = &eachAtomNeighbors[atom1];
        for(int j=0; j<atom1List->size(); j++)
        {
            int atom2 = atom1List->at(j);
            vector<int> *atom2List = &eachAtomNeighbors[atom2];
            for(int k=0; k<atom2List->size(); k++)
            {
                int atom3 = atom2List->at(k);
                //atom1-atom2, so atom2List contains atom1 which should not be considered
                if(atom3 == atom1)
                    continue;
                vector<int> *atom3List = &eachAtomNeighbors[atom3];
                for(int l=0; l<atom3List->size(); l++)
                {
                    int atom4 = atom3List->at(l);
                    //atom1-atom2, so atom2List contains atom1 which should not be considered
                    if(atom4 == atom2)
                        continue;
                    if(atom1<atom4)
                        allExcls.add(Exclusion(atom1, atom4, modified));
                    else
                        allExcls.add(Exclusion(atom4, atom1, modified));
                }
            }
        }
    }
}
