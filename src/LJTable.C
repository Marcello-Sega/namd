/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: LJTable.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/08 04:49:01 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: LJTable.C,v $
 * Revision 1.1  1996/10/08 04:49:01  brunner
 * Initial revision
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/LJTable.C,v 1.1 1996/10/08 04:49:01 brunner Exp $";
#include "LJTable.h"
#include "Node.h"
#include "Parameters.h"
#include "Inform.h"

//----------------------------------------------------------------------  
LJTable::LJTable()
{
  table_dim = namdMyNode->params->get_num_vdw_params();
  half_table_sz = table_dim * table_dim;

  namdInfo << "Allocating LJ Table: size = " << table_dim << "\n" << sendmsg;
  
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
  namdInfo << "LJ Table done\n" << sendmsg;
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

  Real A, B, A14, B14;
  //  We need the A and B parameters for the Van der Waals.  These can
  //  be explicitly be specified for this pair or calculated from the
  //  sigma and epsilon values for the two atom types
  if (namdMyNode->params->get_vdw_pair_params(i,j, &A, &B, &A14, &B14))
  {
    cur->A = A;
    cur->B = B;
    cur_scaled->A = A14;
    cur_scaled->B = B14;
  }
  else
  {
    //  We didn't find explicit parameters for this pair. So instead,
    //  get the parameters for each atom type separately and use them
    //  to calculate the values we need
    Real sigma_i, sigma_i14, epsilon_i, epsilon_i14;
    Real sigma_j, sigma_j14, epsilon_j, epsilon_j14;

    namdMyNode->params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
				       &epsilon_i14,i);
    namdMyNode->params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14, 
				       &epsilon_j14,j);
  	
	    
    BigReal sigma_ij = 0.5 * (sigma_i+sigma_j);
    BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);
    BigReal sigma_ij14 = 0.5 * (sigma_i14+sigma_j14);
    BigReal epsilon_ij14 = sqrt(epsilon_i14*epsilon_j14);

    
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
}


