
/*
 * COPYRIGHT 1992, REGENTS OF THE UNIVERSITY OF CALIFORNIA
 *
 *  prm.c - read information from an amber PARM topology file: 
 *	atom/residue/bond/charge info, plus force field data. 
 *	This file and the accompanying prm.h may be distributed 
 *	provided this notice is retained unmodified and provided 
 *	that any modifications to the rest of the file are noted 
 *	in comments.
 *
 *	Bill Ross, UCSF 1994
 */

// A few changes were made to the original code in order to make it
// easy to use in NAMD. Functions were integrated into the object.

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#include "common.h"
#include "InfoStream.h"
#include "parm.h"

#if 0
extern int 	errno;
#endif

static int	debug = 0;	/* set it if you want */

/*	fortran formats 
 *	 9118 FORMAT(12I6)
 *	 9128 FORMAT(5E16.8)
 */
const char	*f9118 = "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n";


// This function was borrowed from VMD code in "ReadPARM.C". Here it
// is used to replace the "skipeoln()" in the original code, which
// seems to be undefined.
static int readtoeoln(FILE *f) {
  char c;

  /* skip to eoln */
  while((c = getc(f)) != '\n') {
    if (c == EOF) 
      return -1;
  }

  return 0;
}  

/***********************************************************************
 							GENOPEN()
************************************************************************/

/*
 *  genopen() - fopen regular or popen compressed file for reading
 */
//  replaced with NAMD implementation

FILE *Ambertoppar::genopen(const char *name)
{
	return(Fopen(name,"r"));
}

/***********************************************************************
 							GENCLOSE()
************************************************************************/

/*
 *  genclose() - close fopened or popened file
 */
//  replaced with NAMD implementation

void Ambertoppar::genclose(FILE *fileptr)
{
	Fclose(fileptr);
}


/***********************************************************************
 							GET()
************************************************************************/

char *Ambertoppar::get(int size)
{
	char	*ptr;

#ifdef DEBUG
	printf("malloc %d\n", size);
	fflush(stdout);
#endif
	if (size ==0)
		return((char *) NULL);

	if ((ptr = (char *) malloc((unsigned)size)) == NULL) {
		printf("malloc %d", size);
		fflush(stdout);
		NAMD_die("Memory allocation error in Ambertoppar::get()");
		exit(1);
	}
	return(ptr);
}

/***********************************************************************
 							PREADLN()
************************************************************************/

void Ambertoppar::preadln(FILE *file, const char *name, char *string)
{
	int 	i, j;

	for (i=0; i<81; i++) {
		if ((j = getc(file)) == EOF) {
			printf("Error: unexpected EOF in %s\n", name);
			exit(1);
		}
		string[i] = (char) j;
		if (string[i] == '\n') {
			break;
		}
	}
	if (i == 80  &&  string[i] != '\n') {
		printf("Error: line too long in %s:\n%.80s", name, string);
		exit(1);
	}
}

/***************************************************************************
								READPARM()
****************************************************************************/

/*
 * readparm() - instantiate a given Ambertoppar
 */

int Ambertoppar::readparm(char *name)
{
	_REAL 		*H;
	int		i, idum, res, ifpert;
	FILE 		*file;

	if (data_read)
	{ printf("Duplicate parm data in one object!\n");
	  return(0);}

//	printf("Reading parm file (%s)\n", name);
	iout << "Reading parm file (" << name << ") ...\n" << endi;

//	if ((file = genopen(name, "parm")) == NULL) 
	if ((file = genopen(name)) == NULL)
		return(0);

	/* READ TITLE */

	preadln(file, name, ititl);
// "ititle" doesn't guarantee to have '\0' (as the end of a string),
// so the following is disabled in order to avoid strange output
//	printf("%s title:\n%s", name, ititl);

	/* READ CONTROL INTEGERS */

	fscanf(file, f9118, 
		&Natom,  &Ntypes, &Nbonh, &Mbona, 
		&Ntheth, &Mtheta, &Nphih, &Mphia, 
		&Nhparm, &Nparm,  &Nnb,   &Nres);

	fscanf(file, f9118, 
		&Nbona,  &Ntheta, &Nphia, &Numbnd, 
		&Numang, &Nptra,  &Natyp, &Nphb, 
		&ifpert,      &idum,        &idum,       &idum);

	if (ifpert) {
		printf("not equipped to read perturbation prmtop\n");
		return(0);
	}
	fscanf(file, " %d %d %d %d %d %d", 
		&idum, &idum,&idum,&IfBox,&Nmxrs,&IfCap);

//	skipeoln(file);
	readtoeoln(file);	

	/* ALLOCATE MEMORY */

	Nat3 = 3 * Natom;
	Ntype2d = Ntypes * Ntypes;
	Nttyp = Ntypes*(Ntypes+1)/2;

	/*
	 * get most of the indirect stuff; some extra allowed for char arrays
	 */

	AtomNames = (char *) get(4*Natom+81);
	Charges = (_REAL *) get(sizeof(_REAL)*Natom);
	Masses = (_REAL *) get(sizeof(_REAL)*Natom);
	Iac = (int *) get(sizeof(int)*Natom);
	Iblo = (int *) get(sizeof(int)*Natom);
	Cno = (int *) get(sizeof(int)* Ntype2d);
	ResNames = (char *) get(4* Nres+81);
	Ipres = (int *) get(sizeof(int)*( Nres+1));
        Rk = (_REAL *) get(sizeof(_REAL)* Numbnd);
        Req = (_REAL *) get(sizeof(_REAL)* Numbnd);
        Tk = (_REAL *) get(sizeof(_REAL)* Numang);
        Teq = (_REAL *) get(sizeof(_REAL)* Numang);
        Pk = (_REAL *) get(sizeof(_REAL)* Nptra);
        Pn = (_REAL *) get(sizeof(_REAL)* Nptra);
        Phase = (_REAL *) get(sizeof(_REAL)* Nptra);
        Solty = (_REAL *) get(sizeof(_REAL)* Natyp);
        Cn1 = (_REAL *) get(sizeof(_REAL)* Nttyp);
        Cn2 = (_REAL *) get(sizeof(_REAL)* Nttyp);
	BondHAt1 = (int *) get(sizeof(int)* Nbonh);
	BondHAt2 = (int *) get(sizeof(int)* Nbonh);
	BondHNum = (int *) get(sizeof(int)* Nbonh);
	BondAt1 = (int *) get(sizeof(int)* Nbona);
	BondAt2 = (int *) get(sizeof(int)* Nbona);
	BondNum = (int *) get(sizeof(int)* Nbona);
	AngleHAt1 = (int *) get(sizeof(int)* Ntheth);
	AngleHAt2 = (int *) get(sizeof(int)* Ntheth);
	AngleHAt3 = (int *) get(sizeof(int)* Ntheth);
	AngleHNum = (int *) get(sizeof(int)* Ntheth);
	AngleAt1 = (int *) get(sizeof(int)* Ntheta);
	AngleAt2 = (int *) get(sizeof(int)*Ntheta);
	AngleAt3 = (int *) get(sizeof(int)*Ntheta);
	AngleNum = (int *) get(sizeof(int)*Ntheta);
	DihHAt1 = (int *) get(sizeof(int)*Nphih);
	DihHAt2 = (int *) get(sizeof(int)*Nphih);
	DihHAt3 = (int *) get(sizeof(int)*Nphih);
	DihHAt4 = (int *) get(sizeof(int)*Nphih);
	DihHNum = (int *) get(sizeof(int)*Nphih);
	DihAt1 = (int *) get(sizeof(int)*Nphia);
	DihAt2 = (int *) get(sizeof(int)*Nphia);
	DihAt3 = (int *) get(sizeof(int)*Nphia);
	DihAt4 = (int *) get(sizeof(int)*Nphia);
	DihNum = (int *) get(sizeof(int)*Nphia);
	ExclAt = (int *) get(sizeof(int)*Nnb);
	HB12 = (_REAL *) get(sizeof(_REAL)*Nphb);
	HB6 = (_REAL *) get(sizeof(_REAL)*Nphb);
	AtomSym = (char *) get(4*Natom+81);
	AtomTree = (char *) get(4*Natom+81);
	TreeJoin = (int *) get(sizeof(int)*Natom);
	AtomRes = (int *) get(sizeof(int)*Natom);

	/* 
	 * READ ATOM NAMES -IH(M04)
	 */

	for (i=0; i<(Natom/20 + (Natom%20 ? 1 : 0)); i++)
		preadln(file, "", &AtomNames[i*80]);

	/* 
	 * READ ATOM CHARGES -X(L15)
	 *	(pre-multiplied by an energy factor of 18.2223 == sqrt(332)
	 *	 for faster force field calculations)
	 */

	for (i=0; i<Natom; i++)
#ifdef DOUBLE
		fscanf(file, " %lf", &Charges[i]);
#else
		fscanf(file, " %f", &Charges[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ ATOM MASSES -X(L20)
	 */

	for (i=0; i<Natom; i++)
#ifdef DOUBLE
		fscanf(file, " %le", &Masses[i]);
#else
		fscanf(file, " %e", &Masses[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ ATOM L-J TYPES -IX(I04)
	 */

	for (i=0; i<Natom; i++)
		fscanf(file, " %d", &Iac[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ ATOM INDEX TO 1st IN EXCLUDED ATOM LIST "NATEX" -IX(I08)
	 */

	for (i=0; i<Natom; i++)
		fscanf(file, " %d", &Iblo[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ TYPE INDEX TO N-B TYPE -IX(I06)
	 */

	for (i=0; i<Ntype2d; i++)
		fscanf(file, " %d", &Cno[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ RES NAMES (4 chars each, 4th blank) -IH(M02)
	 */

	for (i=0; i<(Nres/20 + (Nres%20 ? 1 : 0)); i++)
		preadln(file, "", &ResNames[i*80]);

	/* 
	 * READ RES POINTERS TO 1st ATOM 		-IX(I02)
	 */

	for (i=0; i<Nres; i++) 
		fscanf(file, " %d", &Ipres[i]);
	Ipres[Nres] = Natom + 1;
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ BOND FORCE CONSTANTS 			-RK()
	 */

	for (i=0; i< Numbnd; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Rk[i]);
#else
		fscanf(file, " %f", &Rk[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ BOND LENGTH OF MINIMUM ENERGY  		-REQ()
	 */

	for (i=0; i< Numbnd; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Req[i]);
#else
		fscanf(file, " %f", &Req[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ BOND ANGLE FORCE CONSTANTS (following Rk nomen) -TK()
	 */

	for (i=0; i< Numang; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Tk[i]);
#else
		fscanf(file, " %f", &Tk[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ BOND ANGLE OF MINIMUM ENERGY (following Req nomen) -TEQ()
	 */

	for (i=0; i< Numang; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Teq[i]);
#else
		fscanf(file, " %f", &Teq[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ DIHEDRAL PEAK MAGNITUDE 		-PK()
	 */

	for (i=0; i< Nptra; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Pk[i]);
#else
		fscanf(file, " %f", &Pk[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ DIHEDRAL PERIODICITY 			-PN()
	 */

	for (i=0; i< Nptra; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Pn[i]);
#else
		fscanf(file, " %f", &Pn[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ DIHEDRAL PHASE  			-PHASE()
	 */

	for (i=0; i< Nptra; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Phase[i]);
#else
		fscanf(file, " %f", &Phase[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * ?? "RESERVED" 				-SOLTY()
	 */

	for (i=0; i< Natyp; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Solty[i]);
#else
		fscanf(file, " %f", &Solty[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ L-J R**12 FOR ALL PAIRS OF ATOM TYPES  	-CN1()
	 *	(SHOULD BE 0 WHERE H-BONDS)
	 */

	for (i=0; i< Nttyp; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Cn1[i]);
#else
		fscanf(file, " %f", &Cn1[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ L-J R**6 FOR ALL PAIRS OF ATOM TYPES 	-CN2()
	 *	(SHOULD BE 0 WHERE H-BONDS)
	 */

	for (i=0; i< Nttyp; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &Cn2[i]);
#else
		fscanf(file, " %f", &Cn2[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ COVALENT BOND W/ HYDROGEN (3*(atnum-1)): 
	 *	IBH = ATOM1 		-IX(I12)
	 *	JBH = ATOM2 		-IX(I14)
	 *	ICBH = BOND ARRAY PTR	-IX(I16)
	 */

	for (i=0; i<Nbonh; i++) 
		fscanf(file, " %d %d %d", 
		    &BondHAt1[i], &BondHAt2[i], &BondHNum[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ COVALENT BOND W/OUT HYDROGEN (3*(atnum-1)):
	 *	IB = ATOM1		-IX(I18)
	 *	JB = ATOM2		-IX(I20)
	 *	ICB = BOND ARRAY PTR	-IX(I22)
	 */

	for (i=0; i<Nbona; i++)
		fscanf(file, " %d %d %d", 
			&BondAt1[i], &BondAt2[i], &BondNum[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ ANGLE W/ HYDROGEN: 
	 *	ITH = ATOM1			-IX(I24)
	 *	JTH = ATOM2			-IX(I26)
	 *	KTH = ATOM3			-IX(I28)
	 *	ICTH = ANGLE ARRAY PTR		-IX(I30)
	 */

	for (i=0; i<Ntheth; i++)
		fscanf(file, " %d %d %d %d", 
		    		&AngleHAt1[i], &AngleHAt2[i], 
				&AngleHAt3[i], &AngleHNum[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ ANGLE W/OUT HYDROGEN: 
	 *	IT = ATOM1			-IX(I32)
	 *	JT = ATOM2			-IX(I34)
	 *	KT = ATOM3			-IX(I36)
	 *	ICT = ANGLE ARRAY PTR		-IX(I38)
	 */

	for (i=0; i<Ntheta; i++)
		fscanf(file, " %d %d %d %d", 
		    		&AngleAt1[i], &AngleAt2[i], 
				&AngleAt3[i], &AngleNum[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ DIHEDRAL W/ HYDROGEN: 
	 *	ITH = ATOM1			-IX(40)
	 *	JTH = ATOM2			-IX(42)
	 *	KTH = ATOM3			-IX(44)
	 *	LTH = ATOM4			-IX(46)
	 *	ICTH = DIHEDRAL ARRAY PTR	-IX(48)
	 */

	for (i=0; i<Nphih; i++)
		fscanf(file, " %d %d %d %d %d", 
		    	&DihHAt1[i], &DihHAt2[i], &DihHAt3[i], 
			&DihHAt4[i], &DihHNum[i]);
//	skipeoln(file);
	readtoeoln(file);

	/* 
	 * READ DIHEDRAL W/OUT HYDROGEN: 
	 *	IT = ATOM1
	 *	JT = ATOM2
	 *	KT = ATOM3
	 *	LT = ATOM4
	 *	ICT = DIHEDRAL ARRAY PTR
	 */

	for (i=0; i<Nphia; i++) {
		fscanf(file, " %d %d %d %d %d", 
		    	&DihAt1[i], &DihAt2[i], &DihAt3[i], 
			&DihAt4[i], &DihNum[i]);
	}
//	skipeoln(file);
	readtoeoln(file);

	/*
	 * READ EXCLUDED ATOM LIST	-IX(I10)
	 */
	for (i=0; i<Nnb; i++)
		fscanf(file, " %d", &ExclAt[i]);
//	skipeoln(file);
	readtoeoln(file);

	/*
	 * READ H-BOND R**12 TERM FOR ALL N-B TYPES	-ASOL()
	 */

	for (i=0; i<Nphb; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &HB12[i]);
#else
		fscanf(file, " %f", &HB12[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);

	/*
	 * READ H-BOND R**6 TERM FOR ALL N-B TYPES	-BSOL()
	 */

	for (i=0; i<Nphb; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &HB6[i]);
#else
		fscanf(file, " %f", &HB6[i]);
#endif
//	skipeoln(file);
	readtoeoln(file);
      
	/*
	 * READ H-BOND CUTOFF (NOT USED) ??		-HBCUT()
	 */

	H = (_REAL *) get(Nphb * sizeof(_REAL));
	for (i=0; i<Nphb; i++) 
#ifdef DOUBLE
		fscanf(file, " %lf", &H[i]);
#else
		fscanf(file, " %f", &H[i]);
#endif
	free((char *)H);
	
//	skipeoln(file);
	readtoeoln(file);
      
	/*
	 * READ ATOM SYMBOLS (FOR ANALYSIS PROGS)	-IH(M06)
	 */

	for (i=0; i<(Natom/20 + (Natom%20 ? 1 : 0)); i++)
		preadln(file, "", &AtomSym[i*80]);

	/*
	 * READ TREE SYMBOLS (FOR ANALYSIS PROGS)	-IH(M08)
	 */

	for (i=0; i<(Natom/20 + (Natom%20 ? 1 : 0)); i++)
		preadln(file, "", &AtomTree[i*80]);
      
	/*
	 * READ TREE JOIN INFO (FOR ANALYSIS PROGS)	-IX(I64)
	 */

	for (i=0; i<Natom; i++)
		fscanf(file, " %d", &TreeJoin[i]);
//	skipeoln(file);
	readtoeoln(file);
      
	/*
	 * READ PER-ATOM RES NUMBER			-IX(I66)
	 *	NOTE: this appears to be something entirely different
	 *	NOTE: overwriting this with correct PER-ATOM RES NUMBERs
	 */

	for (i=0; i<Natom; i++)
		fscanf(file, " %d", &AtomRes[i]);
	res = 0;
	for (i=0; i<Natom; i++) {
		if (i+1 == Ipres[res+1])	/* atom is 1st of next res */
			res++;
		AtomRes[i] = res;
	}
      
	/*
	 * BOUNDARY CONDITION STUFF
	 */

	if (!IfBox) {
		Nspm = 1;
		Boundary = (int *) get(sizeof(int)*Nspm);
		Boundary[0] = Natom;
	} else {
//		skipeoln(file);
		readtoeoln(file);
		fscanf(file, " %d %d %d", &Iptres, &Nspm, 
								&Nspsol);
//		skipeoln(file);
		readtoeoln(file);
		Boundary = (int *) get(sizeof(int)*Nspm);
		for (i=0; i<Nspm; i++)
			fscanf(file, " %d", &Boundary[i]);
//		skipeoln(file);
		readtoeoln(file);
#ifdef DOUBLE
		fscanf(file, " %lf %lf %lf", 
#else
		fscanf(file, " %f %f %f", 
#endif
				&Box[0], &Box[1], &Box[2]);
//		skipeoln(file);
		readtoeoln(file);
		if (Iptres)
			Ipatm = Ipres[Iptres] - 1; 
      		/* IF(IPTRES.GT.0) IPTATM = IX(I02+IPTRES-1+1)-1 */
	}

	/*
	 * ----- LOAD THE CAP INFORMATION IF NEEDED -----
	 */

	if (IfCap) {
		/* if (IfBox) 
			skipeoln(file); */
#ifdef DOUBLE
		fscanf(file, " %d %lf %lf %lf %lf", 
#else
		fscanf(file, " %d %f %f %f %f", 
#endif
				&Natcap, &Cutcap, 
				&Xcap, &Ycap, &Zcap);
	}
	genclose(file);
        if (debug) {
		printf("rdprm done\n");
		fflush(stdout);
	}
	data_read = 1;
	return(1);
}

int Ambertoppar::firstwat()
{
	char	*restr = ResNames; 
	char	*lastres = ResNames + Nres * 4 + 1;
	int	res = 0;

	/*
	 *  find 1st water residue
	 */

	for (; restr<lastres; restr+=4) {
		if (!strncmp(restr, "WAT ", 4)) {
		  printf("first water: res = %d, atom = %d (%.4s)\n", 
					res+1, Ipres[res],
					&AtomNames[Ipres[res]]);
		  fflush(stdout);
		  return(Ipres[res]-1);
		}
		res++;
	}
	return(0);
}

// Constructer: simply set all the pointers Null
Ambertoppar::parm()
{
  data_read = 0;	// No data are read yet
  AtomNames = ResNames = AtomSym = AtomTree = NULL;
  Charges = Masses = Rk = Req = Tk = Teq = Pk = Pn = Phase = NULL;
  Solty = Cn1 = Cn2 = HB12 = HB6 = NULL;
  Iac = Iblo = Cno = Ipres = ExclAt = TreeJoin = AtomRes = NULL;
  BondHAt1 = BondHAt2 = BondHNum = BondAt1 = BondAt2 = NULL;
  BondNum = AngleHAt1 = AngleHAt2 = AngleHAt3 = AngleHNum = NULL;
  AngleAt1 = AngleAt2 = AngleAt3 = AngleNum = DihHAt1 = NULL;
  DihHAt2 = DihHAt3 = DihHAt4 = DihHNum = DihAt1 = DihAt2 = NULL;
  DihAt3 = DihAt4 = DihNum = Boundary = NULL;}

// Destructer: free all the allocated memory for arrays
Ambertoppar::~parm()
{ free(AtomNames);
  free(Charges);
  free(Masses);
  free(Iac);
  free(Iblo);
  free(Cno);
  free(ResNames);
  free(Ipres);
  free(Rk);
  free(Req);
  free(Tk);
  free(Teq);
  free(Pk);
  free(Pn);
  free(Phase);
  free(Solty);
  free(Cn1);
  free(Cn2);
  free(BondHAt1);
  free(BondHAt2);
  free(BondHNum);
  free(BondAt1);
  free(BondAt2);
  free(BondNum);
  free(AngleHAt1);
  free(AngleHAt2);
  free(AngleHAt3);
  free(AngleHNum);
  free(AngleAt1);
  free(AngleAt2);
  free(AngleAt3);
  free(AngleNum);
  free(DihHAt1);
  free(DihHAt2);
  free(DihHAt3);
  free(DihHAt4);
  free(DihHNum);
  free(DihAt1);
  free(DihAt2);
  free(DihAt3);
  free(DihAt4);
  free(DihNum);
  free(ExclAt);
  free(HB12);
  free(HB6);
  free(AtomSym);
  free(AtomTree);
  free(TreeJoin);
  free(AtomRes);}

