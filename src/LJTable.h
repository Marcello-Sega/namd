/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: LJTable.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/08 04:49:01 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: LJTable.h,v $
 * Revision 1.1  1996/10/08 04:49:01  brunner
 * Initial revision
 *
 ***************************************************************************/
#include "common.h"

class LJTable
{
public:
  struct TableEntry
  {
    BigReal A;
    BigReal B;
  };

  LJTable();
  ~LJTable();
  inline TableEntry *table_entry(int index);
  inline int table_index(int i,int j, int scaled14);
  inline void get_LJ_params(int index, BigReal *A, BigReal *B);
private:
  void compute_vdw_params(int i, int j, 
			  TableEntry *cur, TableEntry *cur_scaled);

  TableEntry *table;
  int half_table_sz;
  int table_dim;

};

//======================================================================
inline LJTable::TableEntry *LJTable::table_entry(int index)
{
  return &(table[index]);
}

//======================================================================
inline int LJTable::table_index(int i, int j, int scaled14)
{
  if ((i<0) || (i>=table_dim) || (j<0) || (j>table_dim))
  {
    NAMD_die("Unexpected LJ table value requested in LJTable::table_entry()");
  }
  int index = i * table_dim + j;
  if (scaled14)
    index += half_table_sz;

  return index;
}

//----------------------------------------------------------------------
inline void LJTable::get_LJ_params(int index, BigReal *A, BigReal *B)
{
  TableEntry *cur = table + index;

  *A = cur->A;
  *B = cur->B;
  
return;
}


