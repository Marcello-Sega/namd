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
 *	$Revision: 1.3 $	$Date: 1996/10/12 21:42:18 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: LJTable.h,v $
 * Revision 1.3  1996/10/12 21:42:18  brunner
 * Fixed NAMD_DEBUG test, and added table_val
 *
 * Revision 1.2  1996/10/12 00:11:27  brunner
 * Optimized table_index a bit, and took out bounds checking
 *
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
  inline TableEntry *table_entry(const int index) const;
  inline int table_index(const int i,const int j, const int scaled14) const;
  inline void get_LJ_params(const int index, BigReal *A, BigReal *B) const;
  inline TableEntry *table_val(const int gi, const int gj, 
			       const int scaled14) const;

private:
  void compute_vdw_params(int i, int j, 
			  TableEntry *cur, TableEntry *cur_scaled);

  TableEntry *table;
  int half_table_sz;
  int table_dim;

};

//======================================================================
inline LJTable::TableEntry *LJTable::table_entry(const int index) const
{
  return &(table[index]);
}

//----------------------------------------------------------------------
inline int 
LJTable::table_index(const int i, const int j, const int scaled14) const
{
#ifdef NAMD_DEBUG
  if ((i<0) || (i>=table_dim) || (j<0) || (j>table_dim))
  {
    NAMD_die("Unexpected LJ table value requested in LJTable::table_index()");
  }
#endif

  return i * table_dim + j + (scaled14 ? half_table_sz : 0);
}

//----------------------------------------------------------------------
inline LJTable::TableEntry *
LJTable::table_val(const int i, const int j, const int scaled14) const
{
#if NAMD_DEBUG
  if ((i<0) || (i>=table_dim) || (j<0) || (j>table_dim))
  {
    NAMD_die("Unexpected LJ table value requested in LJTable::table_val()");
  }
#endif
  return table + i * table_dim + j + (scaled14 ? half_table_sz : 0);
}

//----------------------------------------------------------------------
inline void 
LJTable::get_LJ_params(const int index, BigReal *A, BigReal *B) const
{
  TableEntry *cur = table + index;

  *A = cur->A;
  *B = cur->B;
  
return;
}


