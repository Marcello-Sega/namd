//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/
/* DESCRIPTION:                                                            */
/*                                                                         */
/***************************************************************************/
#ifndef LJTABLE_H
#define LJTABLE_H

#include "common.h"
#include "ProcessorPrivate.h"

class LJTable
{
public:
  struct TableEntry
  {
    BigReal exclcut2;
    BigReal A;
    BigReal B;
  };

  static LJTable *Instance(void);
  inline static LJTable *Object(void) {return CpvAccess(LJTable_instance);}
  ~LJTable(void);

  inline TableEntry *table_entry(const int index) const;
  inline int table_index(const int i,const int j, const int scaled14) const;
  inline void get_LJ_params(const int index, BigReal *A, BigReal *B) const;
  inline TableEntry *table_val(const int gi, const int gj, 
			       const int scaled14) const;
  inline TableEntry *table_val(const int gi, const int gj) const;

protected:
  LJTable(void);

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
inline LJTable::TableEntry *
LJTable::table_val(const int i, const int j) const
{
#if NAMD_DEBUG
  if ((i<0) || (i>=table_dim) || (j<0) || (j>table_dim))
  {
    NAMD_die("Unexpected LJ table value requested in LJTable::table_val()");
  }
#endif
  return table + i * table_dim + j;
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

#endif

