/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   dcdlib contains C routines used for reading and writing binary
   dcd files.  The format of these files is from binary FORTRAN output,
   so its pretty ugly.  If you're squimish, don't look!!
*/

#ifndef DCDLIB_H
#define DCDLIB_H

/*  DEFINE ERROR CODES THAT MAY BE RETURNED BY DCD ROUTINES		*/
#define DCD_DNE		-2	/*  DCD file does not exist		*/
#define DCD_OPENFAILED	-3	/*  Open of DCD file failed		*/
#define DCD_BADREAD 	-4	/*  read call on DCD file failed	*/
#define DCD_BADEOF	-5	/*  premature EOF found in DCD file	*/
#define DCD_BADFORMAT	-6	/*  format of DCD file is wrong		*/
#define DCD_FILEEXISTS  -7	/*  output file already exists		*/
#define DCD_BADMALLOC   -8	/*  malloc failed			*/


/*			FUNCTION ALLUSIONS				*/
int open_dcd_read(char *);      /*  Open a DCD file for reading 	*/
int read_dcdheader(int, int*, int*, int*, int*, double*, int*, int**);	
				/*  Read the DCD header			*/
int read_dcdstep(int, int, float*, float*, float*, int, int, int*);	
				/*  Read a timestep's values		*/
int open_dcd_write(char *);     /*  Open a DCD file for writing		*/
int write_dcdstep(int, int, float *, float *, float *);
				/*  Write out a timesteps values	*/
int write_dcdheader(int, char*, int, int, int, int, int, double);	
				/*  Write a dcd header			*/
void close_dcd_read(int, int, int *);
				/*  Close a dcd file open for reading   */
void close_dcd_write(int);	/*  Close a dcd file open for writing   */

#endif /* ! DCDLIB_H */

