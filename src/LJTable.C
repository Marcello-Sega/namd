/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
/* DESCRIPTION:
/*
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/LJTable.C,v 1.777 1997/01/17 19:36:18 ari Exp $";
#include "LJTable.h"
#include "Node.h"
#include "Parameters.h"
#include "InfoStream.h"

LJTable *LJTable::_instance = 0;

LJTable *LJTable::Instance() {
  if (_instance == 0) {
    _instance = new LJTable;
  }
  return _instance;
}

//----------------------------------------------------------------------  
LJTable::LJTable()
{
  table_dim = Node::Object()->parameters->get_num_vdw_params();
  half_table_sz = table_dim * table_dim;

  iout << iINFO << "Allocating LJ Table: size = " << table_dim << "\n" << endi;
  
  table = new TableEntry[half_table_sz * 2];

  for (int i=0; i < table_dim; i++)
    for (int j=i; j < table_dim; j++)
    {
      TableEntry *curij = &(table[i*table_dim+j]);
      TableEntry *curji = &(table[j*table_dim+i]);
      compute_vdw_params(i,j,curij,curij+half_table_sz);

      // Copy to transpose entry
      *curji = *curij;
      *(curji + half_table_sz) = *(curij + half_table_sz);
    }
  iout << iINFO << "LJ Table done\n" << endi;
}

//----------------------------------------------------------------------  
LJTable::~LJTable()
{
  delete [] table;
}

//----------------------------------------------------------------------
void LJTable::compute_vdw_params(int i, int j,
				 LJTable::TableEntry *cur, 
				 LJTable::TableEntry *cur_scaled)
{
  Parameters *params = Node::Object()->parameters;

  Real A, B, A14, B14;
  BigReal sigma_max;
  //  We need the A and B parameters for the Van der Waals.  These can
  //  be explicitly be specified for this pair or calculated from the
  //  sigma and epsilon values for the two atom types
  if (params->get_vdw_pair_params(i,j, &A, &B, &A14, &B14))
  {
    cur->A = A;
    cur->B = B;
    cur_scaled->A = A14;
    cur_scaled->B = B14;

    BigReal sigma_ij = pow(A/B,1./6.);
    BigReal sigma_ij14 = pow(A14/B14,1./6.);

    sigma_max = ( sigma_ij > sigma_ij14 ? sigma_ij : sigma_ij14 );
  }
  else
  {
    //  We didn't find explicit parameters for this pair. So instead,
    //  get the parameters for each atom type separately and use them
    //  to calculate the values we need
    Real sigma_i, sigma_i14, epsilon_i, epsilon_i14;
    Real sigma_j, sigma_j14, epsilon_j, epsilon_j14;

    params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
				       &epsilon_i14,i);
    params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14, 
				       &epsilon_j14,j);
  	
    BigReal sigma_ij = 0.5 * (sigma_i+sigma_j);
    BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);
    BigReal sigma_ij14 = 0.5 * (sigma_i14+sigma_j14);
    BigReal epsilon_ij14 = sqrt(epsilon_i14*epsilon_j14);

    sigma_max = ( sigma_ij > sigma_ij14 ? sigma_ij : sigma_ij14 );

    //  Calculate sigma^6
    sigma_ij *= sigma_ij*sigma_ij;
    sigma_ij *= sigma_ij;
    sigma_ij14 *= sigma_ij14*sigma_ij14;
    sigma_ij14 *= sigma_ij14;
    
    //  Calculate LJ constants A & B
    cur->B = 4.0 * sigma_ij * epsilon_ij;
    cur->A = cur->B * sigma_ij;
    cur_scaled->B = 4.0 * sigma_ij14 * epsilon_ij14;
    cur_scaled->A = cur_scaled->B * sigma_ij14;
  }
  //  Calculate exclcut2
  cur_scaled->exclcut2 = cur->exclcut2 = 0.64 * sigma_max * sigma_max;

}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: LJTable.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:36:18 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: LJTable.C,v $
 * Revision 1.777  1997/01/17 19:36:18  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.4  1996/11/11 19:54:09  nealk
 * Modified to use InfoStream instead of Inform.
 *
 * Revision 1.3  1996/11/05 04:59:56  jim
 * Added exclcut2 to table.
 *
 * Revision 1.2  1996/10/31 22:19:51  jim
 * first incarnation in NAMD 2.0, added singleton pattern
 *
 * Revision 1.1  1996/10/08 04:49:01  brunner
 * Initial revision
 *
 ***************************************************************************/
