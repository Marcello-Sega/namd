/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
/**********************************************************
* Abdulnour (NOUR) Toukmaji, Duke University, ECE Dept 1995 
*
* This file has all the define calls for PVM calls in DPME
***********************************************************/
#define MAXPROC         16   /* the max number of slave processors (exclude master) */

#define PvmDataType     PvmDataRaw       /* define the pvm pack method 
					     {PvmDataInPlace, PvmDataRaw, PvmDataDefault} */
/********************************************************/
/* PVM DEFINES */
#define GROUP "direct_solvers2"
#define MSG_100 1000 /* used in merge_i */ 
#define MSG_200 2000 /* used in swap */
#define MSG_300 3000 /* used in merge_d : dir_ene*/
#define MSG_320 3200 /* used in merge_d : adjust_rcp_ene*/
#define MSG_350 3500 /* used in merge_d : adjust_dir_ene*/
#define MSG_380 3800 /* used in merge_d : self_ene */
#define MSG_400 4000 /* used in setup_recip_sum broadcast */
#define MSG_500 5000 /* used in merge_v */
#define MSG_510 5100 /* used in merge_v */
#define MSG_520 5200 /* used in merge_v */
/***************************************************
* define the working directory where executables reside
* and the data directory where input data file is
*/
#if 0
#define  WORKINGDIR  "/home/scicomp0/ayt/pvm3/bin/SUN4SOL2/PPME/LJS/"
#define  DATADIR  "/home/scicomp0/ayt/pvm3/bin/SUN4SOL2/PPME/LJS/"
#else
#define  WORKINGDIR  "./"
#define  DATADIR  "./"
#endif
