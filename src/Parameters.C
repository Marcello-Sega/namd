/***************************************************************************/
/*      (C) Copyright 1995,1996,1997 The Board of Trustees of the          */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *  The class Parameters is used to hold all of the parameters read
 * in from the parameter files.  The class provides a routine to read in
 * parameter files (as many parameter files as desired can be read in) and
 * a series of routines that allow the parameters that have been read in
 * to be queried.
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Parameters.C,v 1.1005 1997/10/01 16:47:00 milind Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "Parameters.h"
#include "InfoStream.h"
#include "Communicate.h"
#include "ConfigList.h"

typedef Real *RealPtr;
typedef int  *intPtr;

//  struct bond_params is used to form a binary tree of bond parameters.
//  The two atom names are used to determine the order of the nodes in the
//  tree.  atom1name should ALWAYS be lexically before atom2name

struct bond_params
{
  char atom1name[11];
  char atom2name[11];
  Real forceconstant;
  Real distance;
  Index index;
  struct bond_params *left;
  struct bond_params *right;
};

//  struct angle_params is used to form a binary tree of bond parameters.
//  The three atom names are used to determine the order of the nodes in
//  the tree.  atom1name should ALWAYS be lexically before atom3name

struct angle_params
{
  char atom1name[11];
  char atom2name[11];
  char atom3name[11];
  Real forceconstant;
  Real angle;
  Real k_ub;
  Real r_ub;
  Index index;
  struct angle_params *left;
  struct angle_params *right;
};

//  struct dihedral_params is used to form a linked list of the dihedral
//  parameters.  The linked list is arranged in such a way that any
//  bonds with wildcards are at the end of the list so that a linear
//  search can be done but we will still find exact matches before
//  wildcard matches

struct dihedral_params
{
  char atom1name[11];
  char atom2name[11];
  char atom3name[11];
  char atom4name[11];
  int multiplicity;
  FourBodyConsts values[4];
  Index index;
  struct dihedral_params *next;
};

//  struct improper_params is used to form a linked list of the improper
//  parameters.  The linked list is arranged in such a way that any
//  bonds with wildcards are at the end of the list so that a linear
//  search can be done but we will still find exact matches before
//  wildcard matches

struct improper_params
{
  char atom1name[11];
  char atom2name[11];
  char atom3name[11];
  char atom4name[11];
  int multiplicity;
  FourBodyConsts values[4];
  Index index;
  struct improper_params *next;
};

//  struct vdw_params is used to form a binary serach tree of the
//  vdw paramters for a single atom.

struct vdw_params
{
  char atomname[11];
  Real sigma;
  Real epsilon;
  Real sigma14;
  Real epsilon14;
  Index index;
  struct vdw_params *left;
  struct vdw_params *right;
};

//  struct vdw_pair_params is used to form a linked list of the
//  vdw parameters for a pair of atoms

struct vdw_pair_params
{
  char atom1name[11];
  char atom2name[11];
  Real A;
  Real B;
  Real A14;
  Real B14;
  struct vdw_pair_params *next;
};


/************************************************************************/
/*                  */
/*      FUNCTION Parameters        */
/*                  */
/*  This is the constructor for the class.  It simlpy sets the      */
/*  pointers to the list and trees to NULL and the count of all the     */
/*  parameters to 0.              */
/*                  */
/************************************************************************/

Parameters::Parameters(StringList *f)
{
  /*  Set all the pointers to NULL        */
        atomTypeNames=NULL;
  bondp=NULL;
  anglep=NULL;
  improperp=NULL;
  dihedralp=NULL;
  vdwp=NULL;
  vdw_pairp=NULL;
  bond_array=NULL;
  angle_array=NULL;
  dihedral_array=NULL;
  improper_array=NULL;
  vdw_pair_tree=NULL;
  maxDihedralMults=NULL;
  maxImproperMults=NULL;

  /*  Set all the counts to 0          */
  NumBondParams=0;
  NumAngleParams=0;
  NumDihedralParams=0;
  NumImproperParams=0;
  NumVdwParams=0;
  NumVdwPairParams=0;

  /* Set up AllFilesRead flag to FALSE.  Once all of the files    */
  /* have been read in, then this will be set to true and the     */
  /* arrays of parameters will be set up        */
  AllFilesRead = FALSE;

  if (NULL != f) {
       do
      {
    read_parameter_file(f->data);
    f = f->next;
      } while ( f != NULL );
      done_reading_files();
  }
}
/*      END OF FUNCTION Parameters      */

/************************************************************************/
/*                  */
/*      FUNCTION ~Parameters        */
/*                  */
/*  This is the destructor for this class.  It basically just       */
/*  frees all of the memory allocated for the parameters.    */
/*                  */
/************************************************************************/

Parameters::~Parameters()

{
        if (atomTypeNames)
          delete [] atomTypeNames;

  if (bondp != NULL)
    free_bond_tree(bondp);

  if (anglep != NULL)
    free_angle_tree(anglep);

  if (dihedralp != NULL)
    free_dihedral_list(dihedralp);

  if (improperp != NULL)
    free_improper_list(improperp);

  if (vdwp != NULL)
    free_vdw_tree(vdwp);

  if (vdw_pairp != NULL)
    free_vdw_pair_list();

  if (bond_array != NULL)
    delete [] bond_array;

  if (angle_array != NULL)
    delete [] angle_array;

  if (dihedral_array != NULL)
    delete [] dihedral_array;

  if (improper_array != NULL)
    delete [] improper_array;

  if (vdw_array != NULL)
    delete [] vdw_array;
  
  if (vdw_pair_tree != NULL)
    free_vdw_pair_tree(vdw_pair_tree);

  if (maxDihedralMults != NULL)
    delete [] maxDihedralMults;

  if (maxImproperMults != NULL)
    delete [] maxImproperMults;
}
/*      END OF FUNCTION ~Parameters      */

/************************************************************************/
/*                  */
/*      FUNCTION read_paramter_file      */
/*                  */
/*   INPUTS:                */
/*  fname - name of the parameter file to read      */
/*                  */
/*  This function reads in a parameter file and adds the parameters */
/*   from this file to the current group of parameters.  The basic      */
/*   logic of the routine is to read in a line from the file, looks at  */
/*   the first word of the line to determine what kind of parameter we  */
/*   have, and then call the appropriate routine to add the parameter   */
/*   to the parameters that we have.          */
/*                  */
/************************************************************************/

void Parameters::read_parameter_file(char *fname)

{
  char buffer[512];  //  Buffer to store each line of the file
  char first_word[512];  //  First word of the current line
  FILE *pfile;    //  File descriptor for the parameter file

  /*  Check to make sure that we haven't previously been told     */
  /*  that all the files were read        */
  if (AllFilesRead)
  {
    NAMD_die("Tried to read another parameter file after being told that all files were read!");
  }

  /*  Try and open the file          */
  if ( (pfile = Fopen(fname, "r")) == NULL)
  {
    char err_msg[256];

    sprintf(err_msg, "UNABLE TO OPEN PARAMETER FILE %s\n", fname);
    NAMD_die(err_msg);
  }

  /*  Keep reading in lines until we hit the EOF      */
  while (NAMD_read_line(pfile, buffer) != -1)
  {
    /*  Get the first word of the line      */
    NAMD_find_first_word(buffer, first_word);

    /*  First, screen out things that we ignore, such as    */
    /*  blank lines, lines that start with '!', lines that  */
    /*  start with "REMARK", lines that start with set",    */
    /*  and most of the HBOND parameters which include      */
    /*  AEXP, REXP, HAEX, AAEX, but not the HBOND statement */
    /*  which is parsed.                                    */
    if ((buffer[0] != '!') && 
        !NAMD_blank_string(buffer) &&
        (strncasecmp(first_word, "REMARK", 6) != 0) &&
        (strcasecmp(first_word, "set")!=0) &&
        (strncasecmp(first_word, "AEXP", 4) != 0) &&
        (strncasecmp(first_word, "REXP", 4) != 0) &&
        (strncasecmp(first_word, "HAEX", 4) != 0) &&
        (strncasecmp(first_word, "AAEX", 4) != 0) &&
        (strncasecmp(first_word, "NBOND", 5) != 0) &&
        (strncasecmp(first_word, "CUTNB", 5) != 0) &&
        (strncasecmp(first_word, "END", 3) != 0) &&
        (strncasecmp(first_word, "CTONN", 5) != 0) &&
        (strncasecmp(first_word, "EPS", 3) != 0) &&
        (strncasecmp(first_word, "VSWI", 4) != 0) &&
        (strncasecmp(first_word, "NBXM", 4) != 0) &&
        (strncasecmp(first_word, "INHI", 4) != 0) )
    {
      /*  Now, call the appropriate function based    */
      /*  on the type of parameter we have    */
      if (strncasecmp(first_word, "bond", 4)==0)
      {
        add_bond_param(buffer);
        NumBondParams++;
      }
      else if (strncasecmp(first_word, "angl", 4)==0)
      {
        add_angle_param(buffer);
        NumAngleParams++;
      }
      else if (strncasecmp(first_word, "dihe", 4)==0)
      {
        add_dihedral_param(buffer, pfile);
        NumDihedralParams++;
      }
      else if (strncasecmp(first_word, "impr", 4)==0)
      {
        add_improper_param(buffer, pfile);
        NumImproperParams++;
      }
      else if (strncasecmp(first_word, "nonb", 4)==0)
      {
        add_vdw_param(buffer);
        NumVdwParams++;
      }
      else if (strncasecmp(first_word, "nbfi", 4)==0)
      {
        add_vdw_pair_param(buffer);
        NumVdwPairParams++;
      }
      else if (strncasecmp(first_word, "hbon", 4)==0)
      {
        add_hb_pair_param(buffer);
      }
      else
      {
        /*  This is an unknown paramter.        */
        /*  This is BAD        */
        char err_msg[512];

        sprintf(err_msg, "UNKNOWN PARAMETER IN PARAMETER FILE %s\nLINE=*%s*",
           fname, buffer);
        NAMD_die(err_msg);
      }
    }
  }

  /*  Close the file            */
  Fclose(pfile);

  return;
}
/*      END OF FUNCTION read_paramter_file    */

/************************************************************************/
/*                  */
/*      FUNCTION add_bond_param        */
/*                  */
/*   INPUTS:                */
/*  buf - Line from parameter file containing bond parameters  */
/*                  */
/*  This function adds a new bond paramter to the binary tree of    */
/*   angle paramters that we have.  If a duplicate is found, a warning  */
/*   message is printed and the new parameters are used.    */
/*                  */
/************************************************************************/

void Parameters::add_bond_param(char *buf)

{
  char atom1name[11];    //  Atom type for atom 1
  char atom2name[11];    //  Atom type for atom 2
  Real forceconstant;    //  Force constant for bond
  Real distance;      //  Rest distance for bond
  int read_count;      //  Count from sscanf
  struct bond_params *new_node;  //  New node in tree

  /*  Use sscanf to parse up the input line      */
  read_count=sscanf(buf, "%*s %s %s %f %f\n", atom1name, atom2name, 
     &forceconstant, &distance);

  /*  Check to make sure we found everything we expeceted    */
  if (read_count != 4)
  {
    char err_msg[512];

    sprintf(err_msg, "BAD BOND FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buf);
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new bond_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_bond_param\n");
  }

  /*  Order the atoms so that the atom that comes alphabetically  */
  /*  first is atom 1.  Since the bond is symmetric, it doesn't   */
  /*  matter physically which atom is first.  And this allows the */
  /*  search of the binary tree to be done in a logical manner    */
  if (strcasecmp(atom1name, atom2name) < 0)
  {
    strcpy(new_node->atom1name, atom1name);
    strcpy(new_node->atom2name, atom2name);
  }
  else
  {
    strcpy(new_node->atom2name, atom1name);
    strcpy(new_node->atom1name, atom2name);
  }

  /*  Assign force constant and distance        */
  new_node->forceconstant = forceconstant;
  new_node->distance = distance;

  /*  Set pointers to null          */
  new_node->left = NULL;
  new_node->right = NULL;

  /*  Make call to recursive call to actually add the node to the */
  /*  tree              */
  bondp=add_to_bond_tree(new_node, bondp);

  return;
}
/*      END OF FUNCTION add_bond_param      */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_bond_tree      */
/*                  */
/*   INPUTS:                */
/*  new_node - Node to add to the tree        */
/*  tree - tree to add the node to          */
/*                  */
/*   OUTPUTS:                */
/*  ths function returns a pointer to the new tree with the node    */
/*   added to it.  Most of the time it will be the same pointer as was  */
/*   passed in, but not if the current tree is empty.      */
/*                  */
/*  this is a receursive function that adds a node to the binary    */
/*   tree used to store bond parameters.        */
/*                  */
/************************************************************************/

struct bond_params *Parameters::add_to_bond_tree(struct bond_params *new_node,
             struct bond_params *tree)

{
  int compare_code;  //  Results from strcasecmp

  /*  If the tree is currently empty, then the new tree consists  */
  /*  only of the new node          */
  if (tree == NULL)
    return(new_node);

  /*  Compare the atom1 name from the new node and the head of    */
  /*  the tree              */
  compare_code = strcasecmp(new_node->atom1name, tree->atom1name);

  /*  Check to see if they are the same        */
  if (compare_code == 0)
  {
    /*  The atom 1 names are the same, compare atom 2  */
    compare_code = strcasecmp(new_node->atom2name, tree->atom2name);

    /*  If atom 1 AND atom 2 are the same, we have a duplicate */
    if (compare_code == 0)
    {
      /*  We have a duplicate.  So print out a warning*/
      /*  message.  Then assign the new values to the */
      /*  tree and free the new_node      */
      iout << iWARN << "DUPLICATE BOND ENTRY FOR "
        << new_node->atom1name << "-"
        << new_node->atom2name
        << "\nPREVIOUS VALUES  k=" << tree->forceconstant
        << "  x0=" << tree->distance
        << "\n   USING VALUES  k=" << new_node->forceconstant
        << "  x0=" << new_node->distance
        << "\n" << endi;

      tree->forceconstant=new_node->forceconstant;
      tree->distance=new_node->distance;

      delete new_node;

      return(tree);
    }
  }

  /*  We don't have a duplicate, so if the new value is less      */
  /*  than the head of the tree, add it to the left child,   */
  /*  otherwise add it to the right child        */
  if (compare_code < 0)
  {
    tree->left = add_to_bond_tree(new_node, tree->left);
  }
  else
  {
    tree->right = add_to_bond_tree(new_node, tree->right);
  }

  return(tree);
}
/*    END OF FUNCTION add_to_bond_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION add_angle_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line from paramter file with angle parameters    */
/*                  */
/*  this function adds an angle parameter.  It parses up the input  */
/*   line and then adds it to the binary tree used to store the angle   */
/*   parameters.              */
/*                  */
/************************************************************************/

void Parameters::add_angle_param(char *buf)

{
  char atom1name[11];    // Type for atom 1
  char atom2name[11];    // Type for atom 2
  char atom3name[11];    // Type for atom 3
  Real forceconstant;    // Force constant
  Real angle;      // Theta 0
  Real k_ub;      // Urey-Bradley force constant
  Real r_ub;      // Urey-Bradley distance
  int read_count;      // count from sscanf
  struct angle_params *new_node;  // new node in tree

  /*  parse up the input line with sscanf        */
  read_count=sscanf(buf, "%*s %s %s %s %f %f UB %f %f\n", 
     atom1name, atom2name, atom3name, &forceconstant, &angle,
     &k_ub, &r_ub);

  /*  Check to make sure we got what we expected      */
  if ( (read_count != 5) && (read_count != 7) )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD ANGLE FORMAT IN PARAMETER FILE\nLINE=*%s*\n",
       buf);
    NAMD_die(err_msg);
  }

  /*  Allocate the new node          */
  new_node = new angle_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_angle_param");
  }

  /*  As with the bond, we want the atom type is comes first  */
  /*  alphbetically first between atom 1 and atom 3 to be in      */
  /*  atom 1 so that we can search the tree reliably.    */
  if (strcasecmp(atom1name, atom3name) < 0)
  {
    strcpy(new_node->atom1name, atom1name);
    strcpy(new_node->atom2name, atom2name);
    strcpy(new_node->atom3name, atom3name);
  }
  else
  {
    strcpy(new_node->atom3name, atom1name);
    strcpy(new_node->atom2name, atom2name);
    strcpy(new_node->atom1name, atom3name);
  }

  /*  Assign the constants and pointer values      */
  new_node->forceconstant = forceconstant;
  new_node->angle = angle;

  if (read_count == 7)
  {
    //  Urey-Bradley constants
    new_node->k_ub = k_ub;
    new_node->r_ub = r_ub;
  }
  else
  {
    new_node->k_ub = 0.0;
    new_node->r_ub = 0.0;
  }

  new_node->left = NULL;
  new_node->right = NULL;

  /*  Insert it into the tree          */
  anglep = add_to_angle_tree(new_node, anglep);

  return;
}
/*      END OF FUNCTION add_angle_param      */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_angle_tree      */
/*                  */
/*   INPUTS:                */
/*  new_node - new node to add to the angle tree      */
/*  tree - tree to add the node to          */
/*                  */
/*   OUTPUTS:                */
/*  the function returns a pointer to the new tree with the node    */
/*   added.  Most of the time, this will be the same as the value passed*/
/*   in, but not in the case where the tree is empty.      */
/*                  */
/*  this is a recursive function that adds an angle parameter  */
/*   to the binary tree storing the angle parameters.  If a duplicate   */
/*   is found, a warning message is printed, the current values in the  */
/*   tree are replaced with the new values, and the new node is free'd  */
/*                  */
/************************************************************************/

struct angle_params *Parameters::add_to_angle_tree(struct angle_params *new_node,
             struct angle_params *tree)

{
  int compare_code;  //  Return code from strcasecmp

  /*  If the tree is empty, then the new_node is the tree    */
  if (tree == NULL)
    return(new_node);

  /*  Compare atom 1 from the new node and the head of the tree   */
  compare_code = strcasecmp(new_node->atom1name, tree->atom1name);

  if (compare_code == 0)
  {
    /*  Atom 1 is the same, compare atom 2      */
    compare_code = strcasecmp(new_node->atom2name, tree->atom2name);

    if (compare_code == 0)
    {
      /*  Atoms 1 & 2 are the same, compare atom 3  */
      compare_code = strcasecmp(new_node->atom3name, 
            tree->atom3name);

      if (compare_code == 0)
      {
        /*  All three atoms were the same, this */
        /*  is a duplicate.  Print a warning    */
        /*  message, replace the current values,*/
        /*  and free the new node    */
        iout << iWARN << "DUPLICATE ANGLE ENTRY FOR "
          << new_node->atom1name << "-"
          << new_node->atom2name << "-"
          << new_node->atom3name
          << "\nPREVIOUS VALUES  k="
          << tree->forceconstant << "  theta0="
          << tree->angle << " k_ub="
          << tree->k_ub << " r_ub="
          << tree->r_ub
          << "\n   USING VALUES  k="
          << new_node->forceconstant << "  theta0="
          << new_node->angle
          << new_node->angle << " k_ub="
          << new_node->k_ub << " r_ub="
          << "\n" << endi;

        tree->forceconstant=new_node->forceconstant;
        tree->angle=new_node->angle;
        tree->k_ub=new_node->k_ub;
        tree->r_ub=new_node->r_ub;

        delete new_node;

        return(tree);
      }
    }
  }

  /*  Didn't find a duplicate, so if the new_node is smaller  */
  /*  than the current head, add it to the left child.  Otherwise */
  /*  add it to the right child.          */
  if (compare_code < 0)
  {
    tree->left = add_to_angle_tree(new_node, tree->left);
  }
  else
  {
    tree->right = add_to_angle_tree(new_node, tree->right);
  }

  return(tree);
}
/*      END OF FUNCTION add_to_angle_tree    */

/************************************************************************/
/*                  */
/*      FUNCTION add_dihedral_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line from paramter file with dihedral parameters    */
/*                  */
/*  this function adds an dihedral parameter.  It parses up the     */
/*   input line and then adds it to the binary tree used to store the   */
/*   dihedral parameters.            */
/*                  */
/************************************************************************/

void Parameters::add_dihedral_param(char *buf, FILE *fd)

{
  char atom1name[11];       //  Type of atom 1
  char atom2name[11];       //  Type of atom 2
  char atom3name[11];       //  Type of atom 3
  char atom4name[11];       //  Type of atom 4
  Real forceconstant;       //  Force constant
  int periodicity;       //  Periodicity
  Real phase_shift;       //  Phase shift
  int read_count;         //  Count from sscanf
  struct dihedral_params *new_node;  //  New node
  int multiplicity;       //  Multiplicity for bonds
  int i;           //  Loop counter
  char buffer[513];       //  Buffer for new line
  int ret_code;         //  Return code

  /*  Parse up the input line using sscanf      */
  read_count=sscanf(buf, "%*s %s %s %s %s MULTIPLE= %d %f %d %f\n", 
     atom1name, atom2name, atom3name, atom4name, &multiplicity,
     &forceconstant, &periodicity, &phase_shift);

  if ( (read_count != 4) && (read_count != 8) )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD DIHEDRAL FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buf);
    NAMD_die(err_msg);
  }

  if (read_count == 4)
  {
    read_count=sscanf(buf, "%*s %*s %*s %*s %*s %f %d %f\n", 
          &forceconstant, &periodicity, &phase_shift);

    /*  Check to make sure we got what we expected    */
    if (read_count != 3)
    {
      char err_msg[512];

      sprintf(err_msg, "BAD DIHEDRAL FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buf);
      NAMD_die(err_msg);
    }

    multiplicity = 1;
  }

  if (multiplicity > MAX_MULTIPLICITY)
  {
    char err_msg[181];

    sprintf(err_msg, "Multiple dihedral with multiplicity of %d greater than max of %d",
       multiplicity, MAX_MULTIPLICITY);
    NAMD_die(err_msg);
  }

  /*  Allocate new node            */
  new_node = new dihedral_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_dihedral_param\n");
  }

  /*  Assign all of the values for this node.  Notice that since  */
  /*  the dihedrals and impropers are implemented with a linked   */
  /*  list rather than a binary tree, we don't really care about  */
  /*  the order of the atoms any more        */
  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);
  strcpy(new_node->atom3name, atom3name);
  strcpy(new_node->atom4name, atom4name);
  new_node->multiplicity = multiplicity;
  new_node->values[0].k = forceconstant;
  new_node->values[0].n = periodicity;
  new_node->values[0].delta = phase_shift;

  new_node->next = NULL;

  //  If the multiplicity is greater than 1, then read in other parameters
  if (multiplicity > 1)
  {
    for (i=1; i<multiplicity; i++)
    {
      ret_code = NAMD_read_line(fd, buffer);

      //  Get rid of comments at the end of a line
      if (ret_code == 0)
      {
        NAMD_remove_comment(buffer);
      }

      //  Keep reading lines until we get one that isn't blank
      while ( (ret_code == 0) && (NAMD_blank_string(buffer)) )
      {
        ret_code = NAMD_read_line(fd, buffer);
      }

      if (ret_code != 0)
      {
        NAMD_die("EOF encoutner in middle of multiple dihedral");
      }

      read_count=sscanf(buffer, "%f %d %f\n", 
            &forceconstant, &periodicity, &phase_shift);

      if (read_count != 3)
      {
        char err_msg[512];

        sprintf(err_msg, "BAD MULTIPLE FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buffer);
        NAMD_die(err_msg);
      }

      new_node->values[i].k = forceconstant;
      new_node->values[i].n = periodicity;
      new_node->values[i].delta = phase_shift;
    }
  }

  /*  Add this node to the list          */
  add_to_dihedral_list(new_node);

  return;
}
/*      END OF FUNCTION add_dihedral_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_dihedral_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node that is to be added to dihedral_list    */
/*                  */
/*  this function adds a new dihedral parameter to the linked list  */
/*   of dihedral parameters.  First, it checks for duplicates.  If a    */
/*   duplicate is found, a warning message is printed, the old values   */
/*   are replaced with the new values, and the new node is freed.  If   */
/*   Otherwise, the node is added to the list.  This list is arranged   */
/*   so that bods with wildcards are placed at the tail of the list.    */
/*   This will guarantee that if we just do a linear search, we will    */
/*   always find an exact match before a wildcard match.    */
/*                  */
/************************************************************************/

void Parameters::add_to_dihedral_list(
        struct dihedral_params *new_node)

{
  static struct dihedral_params *ptr;   //  position within list
  static struct dihedral_params *tail;  //  Pointer to the end of 
                //  the list so we can add
                //  entries to the end of the
                //  list in constant time
  int i;              //  Loop counter

  /*  If the list is currently empty, then the new node is the list*/
  if (dihedralp == NULL)
  {
    dihedralp=new_node;
    tail=new_node;

    return;
  }

  /*  The list isn't empty, so check for a duplicate    */
  ptr=dihedralp;

  while (ptr != NULL)
  {
    if ( ( (strcasecmp(new_node->atom1name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom4name, ptr->atom4name) == 0) ) ||
         ( (strcasecmp(new_node->atom4name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom1name, ptr->atom4name) == 0) ) )
    {
      /*  Found a duplicate        */
      iout << iWARN << "DUPLICATE DIHEDRAL ENTRY FOR "
        << ptr->atom1name << "-"
        << ptr->atom2name << "-"
        << ptr->atom3name << "-"
        << ptr->atom4name
        << "\nPREVIOUS VALUES MULTIPLICITY " << ptr->multiplicity << "\n";
        
      for (i=0; i<ptr->multiplicity; i++)
      {
        iout << iWARN << " k=" << ptr->values[i].k
                 << "  n=" << ptr->values[i].n
                 << "  delta=" << ptr->values[i].delta;
      }

      iout << iWARN << "\nUSING VALUES MULTIPLICITY " << new_node->multiplicity << "\n";

      for (i=0; i<new_node->multiplicity; i++)
      {
        iout << iWARN << " k=" << new_node->values[i].k
                 << "  n=" << new_node->values[i].n
                 << "  delta=" << new_node->values[i].delta;
      }

      iout << iWARN << endi;

      ptr->multiplicity = new_node->multiplicity;

      for (i=0; i<new_node->multiplicity; i++)
      {
        ptr->values[i].k = new_node->values[i].k;
        ptr->values[i].n = new_node->values[i].n;
        ptr->values[i].delta = new_node->values[i].delta;
      }

      delete new_node;

      return;
    }

    ptr=ptr->next;
  }

  /*  Check to see if we have any wildcards.  Since specific  */
  /*  entries are to take precedence, we'll put anything without  */
  /*  wildcards at the begining of the list and anything with     */
  /*  wildcards at the end of the list.  Then, we can just do a   */
  /*  linear search for a bond and be guaranteed to have specific */
  /*  entries take precendence over over wildcards          */
  if ( (strcasecmp(new_node->atom1name, "X") == 0) ||
       (strcasecmp(new_node->atom2name, "X") == 0) ||
       (strcasecmp(new_node->atom3name, "X") == 0) ||
       (strcasecmp(new_node->atom4name, "X") == 0) )
  {
    /*  add to the end of the list        */
    tail->next=new_node;
    tail=new_node;

    return;
  }
  else
  {
    /*  add to the head of the list        */
    new_node->next=dihedralp;
    dihedralp=new_node;

    return;
  }

}
/*    END OF FUNCTION add_to_dihedral_list      */

/************************************************************************/
/*                  */
/*      FUNCTION add_improper_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line from paramter file with improper parameters    */
/*                  */
/*  this function adds an improper parameter.  It parses up the     */
/*   input line and then adds it to the binary tree used to store the   */
/*   improper parameters.            */
/*                  */
/************************************************************************/

void Parameters::add_improper_param(char *buf, FILE *fd)

{
  char atom1name[11];       //  Atom 1 type
  char atom2name[11];       //  Atom 2 type
  char atom3name[11];       //  Atom 3 type
  char atom4name[11];       //  Atom 4 type
  Real forceconstant;       //  Force constant 
  int periodicity;       //  Periodicity
  Real phase_shift;       //  Phase shift
  int read_count;         //  Count from sscanf
  struct improper_params *new_node;  //  New node
  int multiplicity;       //  Multiplicity for bonds
  int i;           //  Loop counter
  char buffer[513];       //  Buffer for new line
  int ret_code;         //  Return code

  /*  Parse up the line with sscanf        */
  read_count=sscanf(buf, "%*s %s %s %s %s MULTIPLE= %d %f %d %f\n", 
     atom1name, atom2name, atom3name, atom4name, &multiplicity, 
     &forceconstant, &periodicity, &phase_shift);

  if ( (read_count != 4) && (read_count != 8) )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD IMPROPER FORMAT IN PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }

  if (read_count == 4)
  {
    read_count=sscanf(buf, "%*s %*s %*s %*s %*s %f %d %f\n", 
          &forceconstant, &periodicity, &phase_shift);

    /*  Check to make sure we got what we expected    */
    if (read_count != 3)
    {
      char err_msg[512];

      sprintf(err_msg, "BAD IMPROPER FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buf);
      NAMD_die(err_msg);
    }

    multiplicity = 1;
  }

  if (multiplicity > MAX_MULTIPLICITY)
  {
    char err_msg[181];

    sprintf(err_msg, "Multiple improper with multiplicity of %d greater than max of %d",
       multiplicity, MAX_MULTIPLICITY);
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new improper_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_improper_param");
  }

  /*  Assign the values for this bond.  As with the dihedrals,    */
  /*  the atom order doesn't matter        */
  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);
  strcpy(new_node->atom3name, atom3name);
  strcpy(new_node->atom4name, atom4name);
  new_node->multiplicity = multiplicity;
  new_node->values[0].k = forceconstant;
  new_node->values[0].n = periodicity;
  new_node->values[0].delta = phase_shift;

  new_node->next = NULL;

  //  Check to see if this improper has multiple values
  if (multiplicity > 1)
  {
    //  Loop through and read the other values
    for (i=1; i<multiplicity; i++)
    {
      ret_code = NAMD_read_line(fd, buffer);

      //  Strip off comments at the end of the line
      if (ret_code == 0)
      {
        NAMD_remove_comment(buffer);
      }

      //  Skip blank lines
      while ( (ret_code == 0) && (NAMD_blank_string(buffer)) )
      {
        ret_code = NAMD_read_line(fd, buffer);
      }

      if (ret_code != 0)
      {
        NAMD_die("EOF encoutner in middle of multiple improper");
      }

      //  Get the values from the line
      read_count=sscanf(buffer, "%f %d %f\n", 
            &forceconstant, &periodicity, &phase_shift);

      if (read_count != 3)
      {
        char err_msg[512];

        sprintf(err_msg, "BAD MULTIPLE FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buffer);
        NAMD_die(err_msg);
      }

      new_node->values[i].k = forceconstant;
      new_node->values[i].n = periodicity;
      new_node->values[i].delta = phase_shift;
    }
  }

  /*  Add the paramter to the list        */
  add_to_improper_list(new_node);

  return;
}
/*      END OF FUNCTION add_improper_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_improper_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node that is to be added to imporper_list    */
/*                  */
/*  this function adds a new dihedral parameter to the linked list  */
/*   of improper parameters.  First, it checks for duplicates.  If a    */
/*   duplicate is found, a warning message is printed, the old values   */
/*   are replaced with the new values, and the new node is freed.  If   */
/*   Otherwise, the node is added to the list.  This list is arranged   */
/*   so that bods with wildcards are placed at the tail of the list.    */
/*   This will guarantee that if we just do a linear search, we will    */
/*   always find an exact match before a wildcard match.    */
/*                  */
/************************************************************************/

void Parameters::add_to_improper_list(struct improper_params *new_node)

{
  int i;              //  Loop counter
  static struct improper_params *ptr;   //  position within list
  static struct improper_params *tail;  //  Pointer to the end of 
                //  the list so we can add
                //  entries to the end of the
                //  list in constant time

  /*  If the list is currently empty, then the new node is the list*/
  if (improperp == NULL)
  {
    improperp=new_node;
    tail=new_node;

    return;
  }

  /*  The list isn't empty, so check for a duplicate    */
  ptr=improperp;

  while (ptr != NULL)
  {
    if ( ( (strcasecmp(new_node->atom1name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom4name, ptr->atom4name) == 0) ) ||
         ( (strcasecmp(new_node->atom4name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom1name, ptr->atom4name) == 0) ) )
    {
      /*  Found a duplicate        */
      iout << iWARN << "DUPLICATE IMPROPER DIHEDRAL ENTRY FOR "
        << ptr->atom1name << "-"
        << ptr->atom2name << "-"
        << ptr->atom3name << "-"
        << ptr->atom4name
        << "\nPREVIOUS VALUES MULTIPLICITY " << ptr->multiplicity << "\n";
        
      for (i=0; i<ptr->multiplicity; i++)
      {
        iout << iWARN << " k=" << ptr->values[i].k
                 << "  n=" << ptr->values[i].n
                 << "  delta=" << ptr->values[i].delta;
      }

      iout << iWARN << "\nUSING VALUES MULTIPLICITY " << new_node->multiplicity << "\n";

      for (i=0; i<new_node->multiplicity; i++)
      {
        iout << iWARN << " k=" << new_node->values[i].k
                 << "  n=" << new_node->values[i].n
                 << "  delta=" << new_node->values[i].delta;
      }

      iout << iWARN << endi;

      ptr->multiplicity = new_node->multiplicity;

      for (i=0; i<new_node->multiplicity; i++)
      {
        ptr->values[i].k = new_node->values[i].k;
        ptr->values[i].n = new_node->values[i].n;
        ptr->values[i].delta = new_node->values[i].delta;
      }


      delete new_node;

      return;
    }

    ptr=ptr->next;
  }

  /*  Check to see if we have any wildcards.  Since specific  */
  /*  entries are to take precedence, we'll put anything without  */
  /*  wildcards at the begining of the list and anything with     */
  /*  wildcards at the end of the list.  Then, we can just do a   */
  /*  linear search for a bond and be guaranteed to have specific */
  /*  entries take precendence over over wildcards          */
  if ( (strcasecmp(new_node->atom1name, "X") == 0) ||
       (strcasecmp(new_node->atom2name, "X") == 0) ||
       (strcasecmp(new_node->atom3name, "X") == 0) ||
       (strcasecmp(new_node->atom4name, "X") == 0) )
  {
    /*  add to the end of the list        */
    tail->next=new_node;
    tail=new_node;

    return;
  }
  else
  {
    /*  add to the head of the list        */
    new_node->next=improperp;
    improperp=new_node;

    return;
  }
}
/*    END OF FUNCTION add_to_improper_list      */

/************************************************************************/
/*                  */
/*      FUNCTION add_vdw_param        */
/*                  */
/*  INPUTS:                */
/*  buf - line containing the vdw information      */
/*                  */
/*  add_vdw_param adds a vdw parameter for an atom to the current   */
/*  binary tree of values.            */
/*                  */
/************************************************************************/

void Parameters::add_vdw_param(char *buf)

{
  char atomname[11];    //  atom type of paramter
  Real sigma;      //  sigma value for this atom
  Real epsilon;      //  epsilon value for this atom
  Real sigma14;      //  sigma value for 1-4 interactions
  Real epsilon14;      //  epsilon value for 1-4 interactions
  int read_count;      //  count returned by sscanf
  struct vdw_params *new_node;  //  new node for tree

  /*  Parse up the line with sscanf        */
  read_count=sscanf(buf, "%*s %s %f %f %f %f\n", atomname, 
     &epsilon, &sigma, &epsilon14, &sigma14);

  /*  Check to make sure we got what we expected      */
  if (read_count != 5)
  {
    char err_msg[512];

    sprintf(err_msg, "BAD vdW FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buf);
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new vdw_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_vdw_param");
  }

  /*  Assign the values to the new node        */
  strcpy(new_node->atomname, atomname);
  new_node->sigma = sigma;
  new_node->sigma14 = sigma14;
  new_node->epsilon = epsilon;
  new_node->epsilon14 = epsilon14;

  new_node->left = NULL;
  new_node->right = NULL;

  /*  Add the new node into the tree        */
  vdwp=add_to_vdw_tree(new_node, vdwp);

  return;
}
/*      END OF FUNCTION add_vdw_param      */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_vdw_tree      */
/*                  */
/*   INPUTS:                */
/*  new_node - node to add to tree          */
/*  tree - tree to add the node to          */
/*                  */
/*   OUTPUTS:                */
/*  the function returns a pointer to the tree with the node added  */
/*                  */
/*  this function adds a vdw to the binary tree containing the      */
/*   parameters.              */
/*                  */
/************************************************************************/

struct vdw_params *Parameters::add_to_vdw_tree(struct vdw_params *new_node,
             struct vdw_params *tree)

{
  int compare_code;  //  Return code from strcasecmp

  /*  If the tree is currently empty, the new node is the tree    */
  if (tree == NULL)
    return(new_node);

  compare_code = strcasecmp(new_node->atomname, tree->atomname);

  /*  Check to see if we have a duplicate        */
  if (compare_code==0)
  {
    /*  We have a duplicate.  So print out a warning   */
    /*  message, copy the new values into the current node  */
    /*  of the tree, and then free the new_node    */
    iout << iWARN << "DUPLICATE vdW ENTRY FOR " << tree->atomname
      << "\nPREVIOUS VALUES  sigma=" << tree->sigma
      << " epsilon=" << tree->epsilon
      << " sigma14=" << tree->sigma14
      << " epsilon14" << tree->epsilon14
      << "\n   USING VALUES  sigma=" << new_node->sigma
      << " epsilon=" << new_node->epsilon
      << " sigma14=" << new_node->sigma14
      << " epsilon14" << new_node->epsilon14
      << "\n" << endi;

    tree->sigma=new_node->sigma;
    tree->epsilon=new_node->epsilon;
    tree->sigma14=new_node->sigma14;
    tree->epsilon14=new_node->epsilon14;

    delete new_node;

    return(tree);
  }

  /*  Otherwise, if the new node is less than the head of    */
  /*  the tree, add it to the left child, and if it is greater  */
  /*  add it to the right child          */
  if (compare_code < 0)
  {
    tree->left = add_to_vdw_tree(new_node, tree->left);
  }
  else
  {
    tree->right = add_to_vdw_tree(new_node, tree->right);
  }

  return(tree);
}
/*      END OF FUNCTION add_to_vdw_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION add_vdw_pair_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line containing the vdw_pair information      */
/*                  */
/*  this function adds a vdw_pair parameter to the current          */
/*   parameters.              */
/*                  */
/************************************************************************/

void Parameters::add_vdw_pair_param(char *buf)

{
  char atom1name[11];      //  Atom 1 name
  char atom2name[11];      //  Atom 2 name
  Real A;          //  A value for pair
  Real B;          //  B value for pair
  Real A14;        //  A value for 1-4 ints
  Real B14;        //  B value for 1-4 ints
  int read_count;        //  count from sscanf
  struct vdw_pair_params *new_node;  //  new node

  /*  Parse up the input line using sscanf      */
  read_count=sscanf(buf, "%*s %s %s %f %f %f %f\n", atom1name, 
     atom2name, &A, &B, &A14, &B14);

  /*  Check to make sure we got what we expected      */
  if (read_count != 6)
  {
    char err_msg[512];

    sprintf(err_msg, "BAD vdW PAIR FORMAT IN PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new vdw_pair_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_vdw_pair_param\n");
  }

  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);

  /*  Assign values to this node          */
  new_node->A = A;
  new_node->A14 = A14;
  new_node->B = B;
  new_node->B14 = B14;

  new_node->next = NULL;

  /*  Add this node to the tree          */
  add_to_vdw_pair_list(new_node);

  return;
}
/*      END OF FUNCTION add_vdw_par_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_hb_pair_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line containing the hydrogen bond information    */
/*                  */
/*  this function adds data for a hydrogen bond interaction pair    */
/*   to the hbondParams object.                                         */
/*                  */
/************************************************************************/

void Parameters::add_hb_pair_param(char *buf)

{
  char a1n[11];      //  Atom 1 name
  char a2n[11];      //  Atom 2 name
  Real A, B;      //  A, B value for pair

  /*  Parse up the input line using sscanf      */
  if (sscanf(buf, "%*s %s %s %f %f\n", a1n, a2n, &A, &B) != 4)
  {
    char err_msg[512];
    sprintf(err_msg, "BAD HBOND PAIR FORMAT IN PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }

  /*  add data */
  if (hbondParams.add_hbond_pair(a1n, a2n, A, B) == FALSE) {
    iout << iWARN << "Duplicate HBOND parameters for types " << a1n
    << " and " << a2n << " found; using latest values." << "\n" << endi;
  }
}
/*      END OF FUNCTION add_hb_par_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_vdw_pair_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node to be added to list        */
/*                  */
/*  This function adds a link to the end of the vdw_pair_list list  */
/*                  */
/************************************************************************/

void Parameters::add_to_vdw_pair_list(struct vdw_pair_params *new_node)

{
     static struct vdw_pair_params *tail=NULL;
  struct vdw_pair_params *ptr;
  int compare_code;
  

  //  If the list was empty, then just make the new node the list
  if (vdw_pairp == NULL)
  {
     vdw_pairp = new_node;
     tail = new_node;
     return;
  }
  
  ptr = vdw_pairp;

  //  Now check the list to see if we have a duplicate entry
  while (ptr!=NULL)
  {
      /*  Compare atom 1            */
      compare_code = strcasecmp(new_node->atom1name, ptr->atom1name);
      
      if (compare_code == 0)
      {
    /*  Atom 1 is the same, compare atom 2      */
    compare_code = strcasecmp(new_node->atom2name, ptr->atom2name);

    if (compare_code==0)
    {
      /*  Found a duplicate.  Print out a warning   */
      /*  message, assign the values to the current   */
      /*  node in the tree, and then free the new_node*/
      iout << iWARN << "DUPLICATE vdW PAIR ENTRY FOR "
        << new_node->atom1name << "-"
        << new_node->atom2name
        << "\nPREVIOUS VALUES  A=" << ptr->A
        << " B=" << ptr->B
        << " A14=" << ptr->A14
        << " B14" << ptr->B14
        << "\n   USING VALUES  A=" << new_node->A
        << " B=" << new_node->B
        << " A14=" << new_node->A14
        << " B14" << new_node->B14
        << "\n" << endi;

      ptr->A=new_node->A;
      ptr->B=new_node->B;
      ptr->A14=new_node->A14;
      ptr->B14=new_node->B14;

      delete new_node;

      return;
    }
      }
      
      ptr = ptr->next;
  }

  //  We didn't find a duplicate, so add this node to the end
  //  of the list
  tail->next = new_node;
  tail = new_node;
}
/*      END OF FUNCTION add_to_vdw_pair_list    */

/************************************************************************/
/*                  */
/*      FUNCTION done_reading_files      */
/*                  */
/*  This function is used to signal the Parameters object that all  */
/*  of the parameter files have been read.  Once the object knows this, */
/*  it can set un indexes for all the parameters and transfer the values*/
/*  to linear arrays.  This will allow constant time access from this   */
/*  point on.                */
/*                  */
/************************************************************************/

void Parameters::done_reading_files()

{
  AllFilesRead = TRUE;

  //  Allocate space for all of the arrays
  if (NumBondParams)
  {
    bond_array = new BondValue[NumBondParams];

    if (bond_array == NULL)
    {
      NAMD_die("memory allocation of bond_array failed!");
    }
  }

  if (NumAngleParams)
  {
    angle_array = new AngleValue[NumAngleParams];

    if (angle_array == NULL)
    {
      NAMD_die("memory allocation of angle_array failed!");
    }
  }

  if (NumDihedralParams)
  {
    dihedral_array = new DihedralValue[NumDihedralParams];

    if (dihedral_array == NULL)
    {
      NAMD_die("memory allocation of dihedral_array failed!");
    }
  }

  if (NumImproperParams)
  {
    improper_array = new ImproperValue[NumImproperParams];

    if (improper_array == NULL)
    {
      NAMD_die("memory allocation of improper_array failed!");
    }
  }

  if (NumVdwParams)
  {
          atomTypeNames = new char[NumVdwParams*(MAX_ATOMTYPE_CHARS+1)];
    vdw_array = new VdwValue[NumVdwParams];
    
    if (vdw_array == NULL)
    {
      NAMD_die("memory allocation of vdw_array failed!");
    }
  }

  //  Assign indexes to each of the parameters and populate the
  //  arrays using the binary trees and linked lists that we have
  //  already read in
  index_bonds(bondp, 0);
  index_angles(anglep, 0);
  NumVdwParamsAssigned = index_vdw(vdwp, 0);
  index_dihedrals();
  index_impropers();
  
  //  Convert the vdw pairs
  convert_vdw_pairs();
}
/*      END OF FUNCTION done_reading_files    */

/************************************************************************/
/*                  */
/*      FUNCTION index_bonds        */
/*                  */
/*   INPUTS:                */
/*  tree - The tree that is to be indexed        */
/*  index - index to start with          */
/*                  */
/*  This is a recursive routine that will traverse the binary tree  */
/*   of bond parameters, assigning an index to each one, and copying    */
/*   the data from the binary tree to the array that will be used from  */
/*   here on.                */
/*                  */
/************************************************************************/

Index Parameters::index_bonds(struct bond_params *tree, Index index)

{
  //  Tree is empty, do nothing
  if (tree==NULL)
    return(index);

  //  If I have a left subtree, index it first
  if (tree->left != NULL)
  {
    index=index_bonds(tree->left, index);
  }

  //  Now assign an index to top node and populate array
  tree->index = index;
  bond_array[index].k = tree->forceconstant;
  bond_array[index].x0 = tree->distance;
  index++;

  //  If I have a right subtree, index it
  if (tree->right != NULL)
  {
    index=index_bonds(tree->right, index);
  }

  return(index);
}
/*      END OF FUNCTION index_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION index_angles        */
/*                  */
/*   INPUTS:                */
/*  tree - The tree that is to be indexed        */
/*  index - index to start with          */
/*                  */
/*  This is a recursive routine that will traverse the binary tree  */
/*   of angle parameters, assigning an index to each one, and copying   */
/*   the data from the binary tree to the array that will be used from  */
/*   here on.                */
/*                  */
/************************************************************************/

Index Parameters::index_angles(struct angle_params *tree, Index index)

{
  //  Tree is empty, do nothing
  if (tree==NULL)
    return(index);

  //  If I have a left subtree, index it first
  if (tree->left != NULL)
  {
    index=index_angles(tree->left, index);
  }

  //  Now assign an index to top node and populate array
  tree->index = index;

  angle_array[index].k = tree->forceconstant;
  angle_array[index].k_ub = tree->k_ub;
  angle_array[index].r_ub = tree->r_ub;

  //  Convert the angle to radians before storing it
  angle_array[index].theta0 = (tree->angle*PI)/180.0;
  index++;

  //  If I have a right subtree, index it
  if (tree->right != NULL)
  {
    index=index_angles(tree->right, index);
  }

  return(index);
}
/*      END OF FUNCTION index_angles      */

/************************************************************************/
/*                  */
/*      FUNCTION index_dihedrals      */
/*                  */
/*  This function walks down the linked list of dihedral parameters */
/*  and assigns an index to each one.  It also copies the data from this*/
/*  linked list to the arrays that will be used from here on out  */
/*                  */
/************************************************************************/

void Parameters::index_dihedrals()

{
  struct dihedral_params *ptr;  //  Current location in list
  Index index=0;      //  Current index value
  int i;        //  Loop counter

  //  Allocate an array to hold the multiplicity present in the
  //  parameter file for each bond.  This will be used to check
  //  the multiplicities that are detected in the psf file

  //  This is kind of ugly, but necessary because of the way that
  //  X-PLOR psf files deal with Charmm22 parameters.  The way
  //  that multiple periodicities are specified is by having
  //  the bonds appear multiple times in the psf file.  This even
  //  if a bond type has multiple parameters defined, they
  //  will be used if the bond appears multiple times in the
  //  psf file.  So we need to store the number of parameters
  //  we have to make sure the psf file doesn't ask for more
  //  parameters than we really have, and we also need to track
  //  how many times the bond appears in the psf file so that
  //  we can decide how many parameters to actually use.
  maxDihedralMults = new int[NumDihedralParams];

  if (maxDihedralMults == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::index_dihedrals()");
  }
  
  //  Start at the head
  ptr = dihedralp;

  while (ptr != NULL)
  {
    //  Copy data to array and assign index

    //  Save the multiplicity in another array
    maxDihedralMults[index] = ptr->multiplicity;

    //  Assign the multiplicity in the actual structure a bogus value
    //  that we will update in assign_dihedral_index
    dihedral_array[index].multiplicity = -1;

    for (i=0; i<ptr->multiplicity; i++)
    {
      dihedral_array[index].values[i].k = ptr->values[i].k;
      dihedral_array[index].values[i].n = ptr->values[i].n;

      //  Convert the angle to radians before storing it
      dihedral_array[index].values[i].delta = ptr->values[i].delta*PI/180.0;
    }

    ptr->index = index;

    index++;
    ptr=ptr->next;
  }
}
/*      END OF FUNCTION index_dihedrals      */

/************************************************************************/
/*                  */
/*      FUNCTION index_impropers      */
/*                  */
/*  This function walks down the linked list of improper parameters */
/*  and assigns an index to each one.  It also copies the data from this*/
/*  linked list to the arrays that will be used from here on out  */
/*                  */
/************************************************************************/

void Parameters::index_impropers()

{
  struct improper_params *ptr;  //  Current place in list
  Index index=0;      //  Current index value
  int i;        //  Loop counter

  //  Allocate an array to hold the multiplicity present in the
  //  parameter file for each bond.  This will be used to check
  //  the multiplicities that are detected in the psf file

  //  This is kind of ugly, but necessary because of the way that
  //  X-PLOR psf files deal with Charmm22 parameters.  The way
  //  that multiple periodicities are specified is by having
  //  the bonds appear multiple times in the psf file.  This even
  //  if a bond type has multiple parameters defined, they
  //  will be used if the bond appears multiple times in the
  //  psf file.  So we need to store the number of parameters
  //  we have to make sure the psf file doesn't ask for more
  //  parameters than we really have, and we also need to track
  //  how many times the bond appears in the psf file so that
  //  we can decide how many parameters to actually use.
  maxImproperMults = new int[NumImproperParams];

  if (maxImproperMults == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::index_impropers()");
  }
  
  //  Start at the head
  ptr = improperp;

  while (ptr != NULL)
  {
    //  Copy data to array and assign index

    //  Save the multiplicity in another array
    maxImproperMults[index] = ptr->multiplicity;

    //  Assign the multiplicity in the actual structure a bogus value
    //  that we will update in assign_dihedral_index
    improper_array[index].multiplicity = -1;

    for (i=0; i<ptr->multiplicity; i++)
    {
      improper_array[index].values[i].k = ptr->values[i].k;
      improper_array[index].values[i].n = ptr->values[i].n;

      //  Convert the angle to radians before storing it
      improper_array[index].values[i].delta = ptr->values[i].delta*PI/180.0;
    }

    ptr->index=index;

    index++;
    ptr=ptr->next;
  }
}
/*      END OF FUNCTION index_impropers      */

/************************************************************************/
/*                  */
/*      FUNCTION index_vdw        */
/*                  */
/*   INPUTS:                */
/*  tree - The tree that is to be indexed        */
/*  index - index to start with          */
/*                  */
/*  This is a recursive routine that will traverse the binary tree  */
/*   of vdw parameters, assigning an index to each one, and copying     */
/*   the data from the binary tree to the array that will be used from  */
/*   here on.                */
/*                  */
/************************************************************************/

Index Parameters::index_vdw(struct vdw_params *tree, Index index)

{
  //  If the tree is empty, do nothing
  if (tree==NULL)
    return(index);

  //  If I have a left subtree, populate it first
  if (tree->left != NULL)
  {
    index=index_vdw(tree->left, index);
  }

  //  Assign the index and copy the data to the array
  tree->index = index;

  vdw_array[index].sigma = tree->sigma;
  vdw_array[index].epsilon = tree->epsilon;
  vdw_array[index].sigma14 = tree->sigma14;
  vdw_array[index].epsilon14 = tree->epsilon14;

  char *nameloc = atom_type_name(index);
  strncpy(nameloc, tree->atomname, MAX_ATOMTYPE_CHARS);
  nameloc[MAX_ATOMTYPE_CHARS] = '\0';

//  iout << iWARN << "Parameters: Stored name for type " << index << ": '";
//      iout << iWARN << nameloc << "'" << "\n" << endi;

  index++;

  //  If I have a right subtree, index it
  if (tree->right != NULL)
  {
    index=index_vdw(tree->right, index);
  }

  return(index);
}
/*      END OF FUNCTION index_vdw      */

/************************************************************************/
/*                  */
/*      FUNCTION assign_vdw_index      */
/*                  */
/*   INPUTS:                */
/*  atomtype - atom type to find          */
/*  atom_ptr - pointer to the atom structure to find vdw paramters  */
/*       for              */
/*                  */
/*   OUTPUTS:                */
/*  the vdw_index field of the atom structure is populated    */
/*                  */
/*  This function searches the binary tree of vdw parameters so     */
/*   that an index can be assigned to this atom.  If the parameter is   */
/*   is found, then the index is assigned.  If the parameter is not     */
/*   found, then NAMD terminates.          */
/*                  */
/************************************************************************/

void Parameters::assign_vdw_index(char *atomtype, Atom *atom_ptr)

{
  struct vdw_params *ptr;    //  Current position in trees
  int found=0;      //  Flag 1->found match
  int comp_code;      //  return code from strcasecmp

  /*  Check to make sure the files have all been read    */
  if (!AllFilesRead)
  {
    NAMD_die("Tried to assign vdw index before all parameter files were read");
  }

  /*  Start at the top            */
  ptr=vdwp;

  /*  While we haven't found a match, and we haven't reached      */
  /*  the bottom of the tree, compare the atom passed in with     */
  /*  the current value and decide if we have a match, or if not, */
  /*  which way to go            */
  while (!found && (ptr!=NULL))
  {
    comp_code = strcasecmp(atomtype, ptr->atomname);

    if (comp_code == 0)
    {
      /*  Found a match!        */
      atom_ptr->vdw_type=ptr->index;
      found=1;
    }
    else if (comp_code < 0)
    {
      /*  Go to the left        */
      ptr=ptr->left;
    }
    else
    {
      /*  Go to the right        */
      ptr=ptr->right;
    }
  }

  /*  Make sure we found it          */
  if (!found)
  {
    char err_msg[100];

    sprintf(err_msg, "DIDN'T FIND vdW PARAMETER FOR ATOM TYPE %s",
       atomtype);
    NAMD_die(err_msg);
  }

  return;
}
/*      END OF FUNCTION assign_vdw_index    */

/************************************************************************/
/*                  */
/*      FUNCTION get_vdw_pair_params      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  A - A value to populate            */
/*  B - B value to populate            */
/*  A14 - A 1-4 value to populate          */
/*  B14 - B 1-4 value to populate          */
/*                  */
/*   OUTPUTS:                */
/*  If a match is found, A, B, A14, and B14 are all populated and a */
/*   1 is returned.  Otherwise, a 0 is returned.      */
/*                    */
/*  This function finds a set of vdw_pair paramters.  It is given   */
/*   the two types of atoms involved.  This is the only paramter for    */
/*   which a match is NOT guaranteed.  There will only be a match if    */
/*   there are specific van der waals parameters for the two atom types */
/*   involved.                */
/*                  */
/************************************************************************/

int Parameters::get_vdw_pair_params(Index ind1, Index ind2, Real *A, 
        Real *B, Real *A14, Real *B14)

{
  IndexedVdwPair *ptr;    //  Current location in tree
  Index temp;      //  Temporary value for swithcing
          // values
  int found=FALSE;    //  Flag 1-> found a match

  ptr=vdw_pair_tree;

  //  We need the smaller type in ind1, so if it isn't already that 
  //  way, switch them        */
  if (ind1 > ind2)
  {
    temp = ind1;
    ind1 = ind2;
    ind2 = temp;
  }

  /*  While we haven't found a match and we're not at the end  */
  /*  of the tree, compare the bond passed in with the tree  */
  while (!found && (ptr!=NULL))
  {
    if ( (ind1 == ptr->ind1) && (ind2 == ptr->ind2) )
    {
       found = TRUE;
    }
    else if ( (ind1 < ptr->ind1) || 
        ( (ind1==ptr->ind1) && (ind2 < ptr->ind2) ) )
    {
      /*  Go left          */
      ptr=ptr->left;
    }
    else
    {
      /*  Go right          */
      ptr=ptr->right;
    }
  }

  /*  If we found a match, assign the values      */
  if (found)
  {
    *A = ptr->A;
    *B = ptr->B;
    *A14 = ptr->A14;
    *B14 = ptr->B14;

    return(TRUE);
  }
  else
  {
    return(FALSE);
  }
}
/*      END OF FUNCTION get_vdw_pair_params    */


/************************************************************************/
/*                  */
/*        FUNCTION assign_bond_index    */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  bond_ptr - pointer to bond structure to populate    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by bond_ptr is populated    */
/*                  */
/*  This function finds a bond in the binary tree of bond values    */
/*   and assigns its index.  If the bond is found, than the bond_type   */
/*   field of the bond structure is populated.  If the parameter is     */
/*   not found, NAMD will terminate.          */
/*                  */
/************************************************************************/

void Parameters::assign_bond_index(char *atom1, char *atom2, Bond *bond_ptr)

{
  struct bond_params *ptr;  //  Current location in tree
  int found=0;      //  Flag 1-> found a match
  int cmp_code;      //  return code from strcasecmp
  char tmp_name[15];    //  Temporary atom name

  /*  Check to make sure the files have all been read    */
  if (!AllFilesRead)
  {
    NAMD_die("Tried to assign bond index before all parameter files were read");
  }

  /*  We need atom1 < atom2, so if that's not the way they        */
  /*  were passed, flip them          */
  if (strcasecmp(atom1, atom2) > 0)
  {
    strcpy(tmp_name, atom1);
    strcpy(atom1, atom2);
    strcpy(atom2, tmp_name);
  }

  /*  Start at the top            */
  ptr=bondp;

  /*  While we haven't found a match and we're not at the end  */
  /*  of the tree, compare the bond passed in with the tree  */
  while (!found && (ptr!=NULL))
  {
    cmp_code=strcasecmp(atom1, ptr->atom1name);

    if (cmp_code == 0)
    {
      cmp_code=strcasecmp(atom2, ptr->atom2name);
    }

    if (cmp_code == 0)
    {
      /*  Found a match        */
      found=1;
      bond_ptr->bond_type = ptr->index;
    }
    else if (cmp_code < 0)
    {
      /*  Go left          */
      ptr=ptr->left;
    }
    else
    {
      /*  Go right          */
      ptr=ptr->right;
    }
  }

  /*  Check to see if we found anything        */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "CAN'T FIND BOND PARAMETERS FOR BOND %s - %s IN PARAMETER FILES", atom1, atom2);
    NAMD_die(err_msg);
  }

  return;
}
/*      END OF FUNCTION assign_bond_index    */

/************************************************************************/
/*                  */
/*      FUNCTION assign_angle_index      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  atom3 - atom type for atom 3          */
/*  angle_ptr - pointer to angle structure to populate    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by angle_ptr is populated    */
/*                  */
/*  This function assigns an angle index to a specific angle.  */
/*   It searches the binary tree of angle parameters for the appropriate*/
/*   values.  If they are found, the index is assigned.  If they are    */
/*   not found, then NAMD will terminate.        */
/*                  */
/************************************************************************/

void Parameters::assign_angle_index(char *atom1, char *atom2, char*atom3,
          Angle *angle_ptr)

{
  struct angle_params *ptr;  //  Current position in tree
  int comp_val;      //  value from strcasecmp
  int found=0;      //  flag 1->found a match
  char tmp_name[15];    //  Temporary atom name

  /*  Check to make sure the files have all been read    */
  if (!AllFilesRead)
  {
    NAMD_die("Tried to assign angle index before all parameter files were read");
  }

  /*  We need atom1 < atom3.  If that was not what we were   */
  /*  passed, switch them            */
  if (strcasecmp(atom1, atom3) > 0)
  {
    strcpy(tmp_name, atom1);
    strcpy(atom1, atom3);
    strcpy(atom3, tmp_name);
  }

  /*  Start at the top            */
  ptr=anglep;

  /*  While we don't have a match and we haven't reached the  */
  /*  bottom of the tree, compare values        */
  while (!found && (ptr != NULL))
  {
    comp_val = strcasecmp(atom1, ptr->atom1name);

    if (comp_val == 0)
    {
      /*  Atom 1 matches, so compare atom 2    */
      comp_val = strcasecmp(atom2, ptr->atom2name);
      
      if (comp_val == 0)
      {
        /*  Atoms 1&2 match, try atom 3    */
        comp_val = strcasecmp(atom3, ptr->atom3name);
      }
    }

    if (comp_val == 0)
    {
      /*  Found a match        */
      found = 1;
      angle_ptr->angle_type = ptr->index;
    }
    else if (comp_val < 0)
    {
      /*  Go left          */
      ptr=ptr->left;
    }
    else
    {
      /*  Go right          */
      ptr=ptr->right;
    }
  }

  /*  Make sure we found a match          */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "UNABLE TO FIND ANGLE PARAMETERS FOR %s %s %s",
       atom1, atom2, atom3);
    NAMD_die(err_msg);
  }

  return;
}
/*      END OF FUNCTION assign_angle_index    */

/************************************************************************/
/*                  */
/*      FUNCTION assign_dihedral_index      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  atom3 - atom type for atom 3          */
/*  atom4 - atom type for atom 4          */
/*  dihedral_ptr - pointer to dihedral structure to populate  */
/*  multiplicity - Multiplicity to assign to this bond    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by dihedral_ptr is populated    */
/*                  */
/*  This function searchs the linked list of dihedral parameters for*/
/*   a given bond.  If a match is found, the dihedral type is assigned. */
/*   If no match is found, NAMD terminates        */
/*                  */
/************************************************************************/

void Parameters::assign_dihedral_index(char *atom1, char *atom2, char *atom3,
        char *atom4, Dihedral *dihedral_ptr,
        int multiplicity)

{
  struct dihedral_params *ptr;  //  Current position in list
  int found=0;      //  Flag 1->found a match

  /*  Start at the begining of the list        */
  ptr=dihedralp;

  /*  While we haven't found a match and we haven't reached       */
  /*  the end of the list, keep looking        */
  while (!found && (ptr!=NULL))
  {
    /*  Do a linear search through the linked list of   */
    /*  dihedral parameters.  Since the list is arranged    */
    /*  with wildcard paramters at the end of the list, we  */
    /*  can simply do a linear search and be guaranteed that*/
    /*  we will find exact matches before wildcard matches. */
    /*  Also, we must check for an exact match, and a match */
    /*  in reverse, since they are really the same          */
    /*  physically.            */
    if ( ( (strcasecmp(ptr->atom1name, atom1)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom2)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom3)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom4name, atom4)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) )
    {
      /*  Found an exact match      */
      found=1;
    }
    else if ( ( (strcasecmp(ptr->atom4name, atom1)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom2)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom3)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom1name, atom4)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) )
    {
      /*  Found a reverse match      */
      found=1;
    }
    else
    {
      /*  Didn't find a match, go to the next node  */
      ptr=ptr->next;
    }
  }

  /*  Make sure we found a match          */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "CAN'T FIND DIHEDRAL PARAMETERS FOR %s  %s  %s  %s",
       atom1, atom2, atom3, atom4);
    
    NAMD_die(err_msg);
  }

  //  Check to make sure the number of multiples specified in the psf
  //  file doesn't exceed the number of parameters in the parameter
  //  files
  if (multiplicity > maxDihedralMults[ptr->index])
  {
    char err_msg[257];

    sprintf(err_msg, "Multiplicity of Paramters for diehedral bond %s %s %s %s of %d exceeded", atom1, atom2, atom3, atom4, maxDihedralMults[ptr->index]);
    NAMD_die(err_msg);
  }

  //  If the multiplicity from the current bond is larger than that
  //  seen in the past, increase the multiplicity for this bond
  if (multiplicity > dihedral_array[ptr->index].multiplicity)
  {
    dihedral_array[ptr->index].multiplicity = multiplicity;
  }

  dihedral_ptr->dihedral_type = ptr->index;

  return;
}
/*      END OF FUNCTION assign_dihedral_index    */

/************************************************************************/
/*                  */
/*      FUNCTION assign_improper_index      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  atom3 - atom type for atom 3          */
/*  atom4 - atom type for atom 4          */
/*  improper_ptr - pointer to improper structure to populate  */
/*   multiplicity - Multiplicity to assign to this bond    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by improper_ptr is populated    */
/*                  */
/*  This function searchs the linked list of improper parameters for*/
/*   a given bond.  If a match is found, the improper_type is assigned. */
/*   If no match is found, NAMD will terminate.        */
/*                  */
/************************************************************************/

void Parameters::assign_improper_index(char *atom1, char *atom2, char *atom3,
        char *atom4, Improper *improper_ptr,
        int multiplicity)

{
  struct improper_params *ptr;  //  Current position in list
  int found=0;      //  Flag 1->found a match

  /*  Start at the head of the list        */
  ptr=improperp;

  /*  While we haven't fuond a match and haven't reached the end  */
  /*  of the list, keep looking          */
  while (!found && (ptr!=NULL))
  {
    /*  Do a linear search through the linked list of   */
    /*  improper parameters.  Since the list is arranged    */
    /*  with wildcard paramters at the end of the list, we  */
    /*  can simply do a linear search and be guaranteed that*/
    /*  we will find exact matches before wildcard matches. */
    /*  Also, we must check for an exact match, and a match */
    /*  in reverse, since they are really the same          */
    /*  physically.            */
    if ( ( (strcasecmp(ptr->atom1name, atom1)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom2)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom3)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom4name, atom4)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) )
    {
      /*  Found an exact match      */
      found=1;
    }
    else if ( ( (strcasecmp(ptr->atom4name, atom1)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom2)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom3)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom1name, atom4)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) )
    {
      /*  Found a reverse match      */
      found=1;
    }
    else
    {
      /*  Didn't find a match, go to the next node  */
      ptr=ptr->next;
    }
  }

  /*  Make sure we found a match          */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "CAN'T FIND IMPROPER PARAMETERS FOR %s  %s  %s  %s",
       atom1, atom2, atom3, atom4);
    
    NAMD_die(err_msg);
  }

  //  Check to make sure the number of multiples specified in the psf
  //  file doesn't exceed the number of parameters in the parameter
  //  files
  if (multiplicity > maxImproperMults[ptr->index])
  {
    char err_msg[257];

    sprintf(err_msg, "Multiplicity of Paramters for improper bond %s %s %s %s of %d exceeded", atom1, atom2, atom3, atom4, maxImproperMults[ptr->index]);
    NAMD_die(err_msg);
  }

  //  If the multiplicity from the current bond is larger than that
  //  seen in the past, increase the multiplicity for this bond
  if (multiplicity > improper_array[ptr->index].multiplicity)
  {
    improper_array[ptr->index].multiplicity = multiplicity;
  }

  /*  Assign the constants          */
  improper_ptr->improper_type = ptr->index;

  return;
}
/*      END OF FUNCTION assign_improper_index    */

/************************************************************************/
/*                  */
/*      FUNCTION free_bond_tree        */
/*                  */
/*   INPUTS:                */
/*  bond_ptr - pointer to bond tree to free        */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a bond paramter tree.  It makes recursive calls to   */
/*   free the left an right subtress, and then frees the head.  It is   */
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_bond_tree(struct bond_params *bond_ptr)

{
  if (bond_ptr->left != NULL)
  {
    free_bond_tree(bond_ptr->left);
  }

  if (bond_ptr->right != NULL)
  {
    free_bond_tree(bond_ptr->right);
  }

  delete bond_ptr;

  return;
}
/*      END OF FUNCTION free_bond_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION free_angle_tree      */
/*                  */
/*   INPUTS:                */
/*  angle_ptr - pointer to angle tree to free      */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a angle paramter tree.  It makes recursive calls to  */
/*   free the left an right subtress, and then frees the head.  It is   */
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_angle_tree(struct angle_params *angle_ptr)

{
  if (angle_ptr->left != NULL)
  {
    free_angle_tree(angle_ptr->left);
  }

  if (angle_ptr->right != NULL)
  {
    free_angle_tree(angle_ptr->right);
  }

  delete angle_ptr;

  return;
}
/*      END OF FUNCTION free_angle_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION free_dihedral_list      */
/*                  */
/*   INPUTS:                */
/*  dih_ptr - pointer to the list to free        */
/*                  */
/*  this function frees a linked list of dihedral parameters.  It   */
/*   is only called by the destructor.          */
/*                  */
/************************************************************************/

void Parameters::free_dihedral_list(struct dihedral_params *dih_ptr)

{
  struct dihedral_params *ptr;  //  Current position in list
  struct dihedral_params *next; //  Next position in list

  ptr=dih_ptr;

  while (ptr != NULL)
  {
    next=ptr->next;
    delete ptr;
    ptr=next;
  }

  return;
}
/*      END OF FUNCTION free_dihedral_list    */

/************************************************************************/
/*                  */
/*      FUNCTION free_improper_list      */
/*                  */
/*   INPUTS:                */
/*  imp_ptr - pointer to the list to free        */
/*                  */
/*  this function frees a linked list of improper parameters.  It   */
/*   is only called by the destructor.          */
/*                  */
/************************************************************************/

void Parameters::free_improper_list(struct improper_params *imp_ptr)

{
  struct improper_params *ptr;  //  Current position in list
  struct improper_params *next; //  Next position in list

  ptr=imp_ptr;

  while (ptr != NULL)
  {
    next=ptr->next;
    delete ptr;
    ptr=next;
  }

  return;
}
/*      END OF FUNCTION free_improper_list    */
    

/************************************************************************/
/*                  */
/*      FUNCTION free_vdw_tree        */
/*                  */
/*   INPUTS:                */
/*  vdw_ptr - pointer to vdw tree to free        */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a vdw paramter tree.  It makes recursive calls to    */
/*   free the left an right subtress, and then frees the head.  It is   */
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_vdw_tree(struct vdw_params *vdw_ptr)

{
  if (vdw_ptr->left != NULL)
  {
    free_vdw_tree(vdw_ptr->left);
  }

  if (vdw_ptr->right != NULL)
  {
    free_vdw_tree(vdw_ptr->right);
  }

  delete vdw_ptr;

  return;
}
/*      END OF FUNCTION free_vdw_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION free_vdw_pair_list      */
/*                  */
/*  This function frees the vdw_pair_list        */
/*                  */
/************************************************************************/

void Parameters::free_vdw_pair_list()
{
   struct vdw_pair_params *ptr, *next;
   
   ptr=vdw_pairp;
   
   while (ptr != NULL)
   {
      next = ptr->next;
      
      delete ptr;
      
      ptr = next;
   }
   
   vdw_pairp = NULL;
}
/*      END OF FUNCTION free_vdw_pair_list    */

/************************************************************************/
/*                  */
/*      FUNCTION free_vdw_pair_tree      */
/*                  */
/*   INPUTS:                */
/*  vdw_pair_ptr - pointer to vdw_pair tree to free      */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a vdw_pair paramter tree.  It makes recursive calls  */
/*   to free the left an right subtress, and then frees the head.  It is*/
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_vdw_pair_tree(IndexedVdwPair *vdw_pair_ptr)

{
  if (vdw_pair_ptr->left != NULL)
  {
    free_vdw_pair_tree(vdw_pair_ptr->left);
  }

  if (vdw_pair_ptr->right != NULL)
  {
    free_vdw_pair_tree(vdw_pair_ptr->right);
  }

  delete vdw_pair_ptr;

  return;
}
/*      END OF FUNCTION free_vdw_pair_tree    */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_bond_params      */
/*                  */
/*   INPUTS:                */
/*  tree - the bond binary tree to traverse        */
/*                  */
/*  This is a recursive call used for debugging purposes that  */
/*   prints out all the bond paramters in the bond parameter binary  */
/*   search tree. It is only called by print_bond_params    */
/*                  */
/************************************************************************/

void Parameters::traverse_bond_params(struct bond_params *tree)

{
  if (tree==NULL)
    return;

  if (tree->left != NULL)
  {
    traverse_bond_params(tree->left);
  }

  DEBUG_MSG("BOND " <<  tree->atom1name << "  " << tree->atom2name \
      << " index=" << tree->index << " k=" << tree->forceconstant \
      << " x0=" << tree->distance);

  if (tree->right != NULL)
  {
    traverse_bond_params(tree->right);
  }
}
/*      END OF FUNCTION traverse_bond_params    */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_angle_params      */
/*                  */
/*   INPUTS:                */
/*  tree - the angle binary tree to traverse      */
/*                  */
/*  This is a recursive call used for debugging purposes that  */
/*   prints out all the angle paramters in the angle parameter binary  */
/*   search tree. It is only called by print_angle_params    */
/*                  */
/************************************************************************/

void Parameters::traverse_angle_params(struct angle_params *tree)

{
  if (tree==NULL)
    return;

  if (tree->left != NULL)
  {
    traverse_angle_params(tree->left);
  }

  DEBUG_MSG("ANGLE " << tree->atom1name << "  " << tree->atom2name \
      << "  " << tree->atom3name << " index=" << tree->index \
      << " k=" << tree->forceconstant << " theta0=" << tree->angle \
      );

  if (tree->right != NULL)
  {
    traverse_angle_params(tree->right);
  }
}
/*      END OF FUNCTION traverse_angle_params    */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_dihedral_params    */
/*                  */
/*   INPUTS:                */
/*  list - the dihedral linked list to traverse      */
/*                  */
/*  This is a call used for debugging purposes that prints out all  */
/*   the bond paramters in the dihedral parameter linked list. It is    */
/*   only called by print_dihedral_params.        */
/*                  */
/************************************************************************/

void Parameters::traverse_dihedral_params(struct dihedral_params *list)

{
  int i;

  while (list != NULL)
  {
    DEBUG_MSG("DIHEDRAL  " << list->atom1name << "  " \
        << list->atom2name << "  " << list->atom3name \
        << "  " << list->atom4name << " index=" \
        << list->index \
        << " multiplicity=" << list->multiplicity << "\n");
        
    for (i=0; i<list->multiplicity; i++)
    {
      DEBUG_MSG("k=" << list->values[i].k \
          << " n=" << list->values[i].n  \
          << " delta=" << list->values[i].delta);
    }

    list=list->next;
  }
}
/*      END OF FUNCTION traverse_dihedral_params  */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_improper_params    */
/*                  */
/*   INPUTS:                */
/*  list - the improper linked list to traverse      */
/*                  */
/*  This is a call used for debugging purposes that prints out all  */
/*   the improper paramters in the improper parameter linked list. It is*/
/*   only called by print_improper_params.        */
/*                  */
/************************************************************************/

void Parameters::traverse_improper_params(struct improper_params *list)

{
  int i;

  while (list != NULL)
  {
    DEBUG_MSG("Improper  " << list->atom1name << "  " \
        << list->atom2name << "  " << list->atom3name \
        << "  " << list->atom4name << " index="  \
        << list->index  \
        << " multiplicity=" << list->multiplicity << "\n");

    for (i=0; i<list->multiplicity; i++)
    {
       DEBUG_MSG("k=" << list->values[i].k \
           << " n=" << list->values[i].n \
           << " delta=" << list->values[i].delta);
    }

    list=list->next;
  }
}
/*      END OF FUNCTION traverse_improper_params  */


/************************************************************************/
/*                  */
/*      FUNCTION traverse_vdw_params      */
/*                  */
/*   INPUTS:                */
/*  tree - the vw binary tree to traverse        */
/*                  */
/*  This is a recursive call used for debugging purposes that  */
/*   prints out all the vdw paramters in the vdw parameter binary  */
/*   search tree. It is only called by print_vdw_params      */
/*                  */
/************************************************************************/

void Parameters::traverse_vdw_params(struct vdw_params *tree)

{
  if (tree==NULL)
    return;

  if (tree->left != NULL)
  {
    traverse_vdw_params(tree->left);
  }

  DEBUG_MSG("vdW " << tree->atomname << " index=" << tree->index \
      << " sigma=" << tree->sigma << " epsilon=" << \
      tree->epsilon << " sigma 1-4=" << tree->sigma14 \
      << " epsilon 1-4=" << tree->epsilon14);

  if (tree->right != NULL)
  {
    traverse_vdw_params(tree->right);
  }
}
/*      END OF FUNCTION traverse_vdw_params    */


/************************************************************************/
/*                  */
/*      FUNCTION traverse_vdw_pair_params    */
/*                  */
/*   INPUTS:                */
/*  list - the vdw_pair list to traverse        */
/*                  */
/*  This call simply prints out the vdw_pair list      */
/*                  */
/************************************************************************/

void Parameters::traverse_vdw_pair_params(struct vdw_pair_params *list)

{
  if (list==NULL)
    return;

  DEBUG_MSG("vdW PAIR  " << list->atom1name << "  "  \
      << list->atom2name << " A=" << list->A \
      << " B=" << list->B << " A 1-4=" \
      << list->A14 << " B 1-4=" << list->B14 \
      );

  traverse_vdw_pair_params(list->next);
}
/*      END OF FUNCTION traverse_vdw_pair_params  */

/************************************************************************/
/*                  */
/*      FUNCTION print_bond_params      */
/*                  */
/*  This is a debugging routine used to print out all the bond  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_bond_params()
{
  DEBUG_MSG(NumBondParams << " BOND PARAMETERS\n" \
      << "*****************************************"  \
      );

  traverse_bond_params(bondp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_angle_params      */
/*                  */
/*  This is a debugging routine used to print out all the angle  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_angle_params()
{
  DEBUG_MSG(NumAngleParams << " ANGLE PARAMETERS\n"
      << "*****************************************" );

  traverse_angle_params(anglep);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_dihedral_params      */
/*                  */
/*  This is a debugging routine used to print out all the dihedral  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_dihedral_params()
{
  DEBUG_MSG(NumDihedralParams << " DIHEDRAL PARAMETERS\n" \
      << "*****************************************" );

  traverse_dihedral_params(dihedralp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_improper_params      */
/*                  */
/*  This is a debugging routine used to print out all the improper  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_improper_params()
{
  DEBUG_MSG(NumImproperParams << " IMPROPER PARAMETERS\n" \
      << "*****************************************" );

  traverse_improper_params(improperp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_vdw_params      */
/*                  */
/*  This is a debugging routine used to print out all the vdw  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_vdw_params()
{
  DEBUG_MSG(NumVdwParams << " vdW PARAMETERS\n" \
      << "*****************************************" );

  traverse_vdw_params(vdwp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_vdw_pair_params      */
/*                  */
/*  This is a debugging routine used to print out all the vdw_pair  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_vdw_pair_params()
{
  DEBUG_MSG(NumVdwPairParams << " vdW PAIR PARAMETERS\n" \
      << "*****************************************" );

  traverse_vdw_pair_params(vdw_pairp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_param_summary      */
/*                  */
/*  This function just prints out a brief summary of the paramters  */
/*  that have been read in.  It is intended for debugging purposes  */
/*                  */
/************************************************************************/

void Parameters::print_param_summary()
{
  iout << iINFO << "SUMMARY OF PARAMETERS:\n" 
     << NumBondParams << " BONDS\n" 
           << NumAngleParams << " ANGLES\n"
           << NumDihedralParams << " DIHEDRAL\n"
           << NumImproperParams << " IMPROPER\n"
           << NumVdwParams << " VDW\n"
           << NumVdwPairParams << " VDW_PAIRS\n"
     << hbondParams.num() << " HBOND_PAIRS\n" << endi;
}


/************************************************************************/
/*                  */
/*      FUNCTION done_reading_structure      */
/*                  */
/*  This function is used to tell the Parameters object that the    */
/*  structure has been read in.  This is so that the Parameters object  */
/*  can now release the binary trees and linked lists that it was using */
/*  to search for parameters based on the atom type.  From this point   */
/*  on, only the arrays of parameter data will be used.  If this object */
/*  resides on any node BUT the master node, it will never even have    */
/*  these trees and lists.  For the master node, this just frees up     */
/*  some memory for better uses.          */
/*                  */
/************************************************************************/

void Parameters::done_reading_structure()

{
  if (bondp != NULL)
    free_bond_tree(bondp);

  if (anglep != NULL)
    free_angle_tree(anglep);

  if (dihedralp != NULL)
    free_dihedral_list(dihedralp);

  if (improperp != NULL)
    free_improper_list(improperp);

  if (vdwp != NULL)
    free_vdw_tree(vdwp);

  //  Free the arrays used to track multiplicity for dihedrals
  //  and impropers
  if (maxDihedralMults != NULL)
    delete [] maxDihedralMults;

  if (maxImproperMults != NULL)
    delete [] maxImproperMults;

  bondp=NULL;
  anglep=NULL;
  dihedralp=NULL;
  improperp=NULL;
  vdwp=NULL;
  maxImproperMults=NULL;
  maxDihedralMults=NULL;
}
/*      END OF FUNCTION done_reading_structure    */

/************************************************************************/
/*                  */
/*      FUNCTION send_Parameters      */
/*                  */
/*  This function is used by the master node to broadcast the       */
/*   structure Parameters to all the other nodes.        */
/*                  */
/************************************************************************/

void Parameters::send_Parameters(Communicate *comm_obj)
{
  Real *a1, *a2, *a3, *a4;        //  Temporary arrays for sending messages
  int *i1, *i2;      //  Temporary int array
  int i, j;      //  Loop counters
  Real **kvals;      //  Force constant values for dihedrals and impropers
  int **nvals;      //  Periodicity values for  dihedrals and impropers
  Real **deltavals;    //  Phase shift values for  dihedrals and impropers
  MOStream *msg=comm_obj->newOutputStream(ALLBUTME, STATICPARAMSTAG, BUFSIZE);
  if ( msg == NULL )
  {
    NAMD_die("memory allocation failed in Parameters::send_Parameters");
  }

  //  Send the bond parameters
  msg->put(NumBondParams);

  if (NumBondParams)
  {
    a1 = new Real[NumBondParams];
    a2 = new Real[NumBondParams];

    if ( (a1 == NULL) || (a2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<NumBondParams; i++)
    {
      a1[i] = bond_array[i].k;
      a2[i] = bond_array[i].x0;
    }

    msg->put(NumBondParams, a1)->put(NumBondParams, a2);

    delete [] a1;
    delete [] a2;
  }

  //  Send the angle parameters
  msg->put(NumAngleParams);

  if (NumAngleParams)
  {
    a1 = new Real[NumAngleParams];
    a2 = new Real[NumAngleParams];
    a3 = new Real[NumAngleParams];
    a4 = new Real[NumAngleParams];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) ||
         (a4 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<NumAngleParams; i++)
    {
      a1[i] = angle_array[i].k;
      a2[i] = angle_array[i].theta0;
      a3[i] = angle_array[i].k_ub;
      a4[i] = angle_array[i].r_ub;
    }

    msg->put(NumAngleParams, a1)->put(NumAngleParams, a2);
    msg->put(NumAngleParams, a3)->put(NumAngleParams, a4);

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }

  //  Send the dihedral parameters
  msg->put(NumDihedralParams);

  if (NumDihedralParams)
  {
    i1 = new int[NumDihedralParams];
    kvals = new RealPtr[MAX_MULTIPLICITY];
    nvals = new intPtr[MAX_MULTIPLICITY];
    deltavals = new RealPtr[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumDihedralParams];
      nvals[i] = new int[NumDihedralParams];
      deltavals[i] = new Real[NumDihedralParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    for (i=0; i<NumDihedralParams; i++)
    {
      i1[i] = dihedral_array[i].multiplicity;

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        kvals[j][i] = dihedral_array[i].values[j].k;
        nvals[j][i] = dihedral_array[i].values[j].n;
        deltavals[j][i] = dihedral_array[i].values[j].delta;
      }
    }

    msg->put(NumDihedralParams, i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->put(NumDihedralParams, kvals[i]);
      msg->put(NumDihedralParams, nvals[i]);
      msg->put(NumDihedralParams, deltavals[i]);

      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Send the improper parameters
  msg->put(NumImproperParams);

  if (NumImproperParams)
  {
    i1 = new int[NumImproperParams];
    kvals = new RealPtr[MAX_MULTIPLICITY];
    nvals = new intPtr[MAX_MULTIPLICITY];
    deltavals = new RealPtr[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumImproperParams];
      nvals[i] = new int[NumImproperParams];
      deltavals[i] = new Real[NumImproperParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    for (i=0; i<NumImproperParams; i++)
    {
      i1[i] = improper_array[i].multiplicity;

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        kvals[j][i] = improper_array[i].values[j].k;
        nvals[j][i] = improper_array[i].values[j].n;
        deltavals[j][i] = improper_array[i].values[j].delta;
      }
    }

    msg->put(NumImproperParams, i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->put(NumImproperParams, kvals[i]);
      msg->put(NumImproperParams, nvals[i]);
      msg->put(NumImproperParams, deltavals[i]);

      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Send the vdw parameters
  msg->put(NumVdwParams);
  msg->put(NumVdwParamsAssigned);

  if (NumVdwParams)
  {
    a1 = new Real[NumVdwParams];
    a2 = new Real[NumVdwParams];
    a3 = new Real[NumVdwParams];
    a4 = new Real[NumVdwParams];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<NumVdwParams; i++)
    {
      a1[i] = vdw_array[i].sigma;
      a2[i] = vdw_array[i].epsilon;
      a3[i] = vdw_array[i].sigma14;
      a4[i] = vdw_array[i].epsilon14;
    }

    msg->put(NumVdwParams * (MAX_ATOMTYPE_CHARS+1), atomTypeNames);
    msg->put(NumVdwParams, a1);
    msg->put(NumVdwParams, a2);
    msg->put(NumVdwParams, a3);
    msg->put(NumVdwParams, a4);

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }
  
  //  Send the vdw pair parameters
  msg->put(NumVdwPairParams);
  
  if (NumVdwPairParams)
  {
    a1 = new Real[NumVdwPairParams];
    a2 = new Real[NumVdwPairParams];
    a3 = new Real[NumVdwPairParams];
    a4 = new Real[NumVdwPairParams];
    i1 = new int[NumVdwPairParams];
    i2 = new int[NumVdwPairParams];    

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) || (a4 == NULL) || 
         (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }
    
    vdw_pair_to_arrays(i1, i2, a1, a2, a3, a4, 0, vdw_pair_tree);
    
    msg->put(NumVdwPairParams, i1)->put(NumVdwPairParams, i2);
    msg->put(NumVdwPairParams, a1);
    msg->put(NumVdwPairParams, a2)->put(NumVdwPairParams, a3);
    msg->put(NumVdwPairParams, a4);
  }

  // send the hydrogen bond parameters
  hbondParams.create_message(msg);
  msg->end();
  delete msg;
}

/************************************************************************/
/*                  */
/*      FUNCTION receive_Parameters      */
/*                  */
/*  This function is used by all the client processes to receive    */
/*   the structure parameters from the master node.      */
/*                  */
/************************************************************************/

void Parameters::receive_Parameters(MIStream *msg)

{
  int i, j;      //  Loop counters
  Real *a1, *a2, *a3, *a4;  //  Temporary arrays to get data from message in
  int *i1, *i2;      //  Temporary int array to get data from message in
  IndexedVdwPair *new_node;  //  New node for vdw pair param tree
  Real **kvals;      //  Force constant values for dihedrals and impropers
  int **nvals;      //  Periodicity values for dihedrals and impropers
  Real **deltavals;    //  Phase shift values for dihedrals and impropers

  //  Get the bonded parameters
  msg->get(NumBondParams);

  if (NumBondParams)
  {
    bond_array = new BondValue[NumBondParams];
    a1 = new Real[NumBondParams];
    a2 = new Real[NumBondParams];

    if ( (bond_array == NULL) || (a1 == NULL) || (a2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(NumBondParams, a1);
    msg->get(NumBondParams, a2);

    for (i=0; i<NumBondParams; i++)
    {
      bond_array[i].k = a1[i];
      bond_array[i].x0 = a2[i];
    }

    delete [] a1;
    delete [] a2;
  }

  //  Get the angle parameters
  msg->get(NumAngleParams);

  if (NumAngleParams)
  {
    angle_array = new AngleValue[NumAngleParams];
    a1 = new Real[NumAngleParams];
    a2 = new Real[NumAngleParams];
    a3 = new Real[NumAngleParams];
    a4 = new Real[NumAngleParams];

    if ( (angle_array == NULL) || (a1 == NULL) || (a2 == NULL) ||
         (a3 == NULL) || (a4 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(NumAngleParams, a1);
    msg->get(NumAngleParams, a2);
    msg->get(NumAngleParams, a3);
    msg->get(NumAngleParams, a4);

    for (i=0; i<NumAngleParams; i++)
    {
      angle_array[i].k = a1[i];
      angle_array[i].theta0 = a2[i];
      angle_array[i].k_ub = a3[i];
      angle_array[i].r_ub = a4[i];
    }

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }

  //  Get the dihedral parameters
  msg->get(NumDihedralParams);

  if (NumDihedralParams)
  {
    dihedral_array = new DihedralValue[NumDihedralParams];

    i1 = new int[NumDihedralParams];
    kvals = new RealPtr[MAX_MULTIPLICITY];
    nvals = new intPtr[MAX_MULTIPLICITY];
    deltavals = new RealPtr[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) || (dihedral_array == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumDihedralParams];
      nvals[i] = new int[NumDihedralParams];
      deltavals[i] = new Real[NumDihedralParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    msg->get(NumDihedralParams, i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->get(NumDihedralParams, kvals[i]);
      msg->get(NumDihedralParams, nvals[i]);
      msg->get(NumDihedralParams, deltavals[i]);
    }

    for (i=0; i<NumDihedralParams; i++)
    {
      dihedral_array[i].multiplicity = i1[i];

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        dihedral_array[i].values[j].k = kvals[j][i];
        dihedral_array[i].values[j].n = nvals[j][i];
        dihedral_array[i].values[j].delta = deltavals[j][i];
      }
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Get the improper parameters
  msg->get(NumImproperParams);

  if (NumImproperParams)
  {
    improper_array = new ImproperValue[NumImproperParams];
    i1 = new int[NumImproperParams];
    kvals = new RealPtr[MAX_MULTIPLICITY];
    nvals = new intPtr[MAX_MULTIPLICITY];
    deltavals = new RealPtr[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) || (improper_array==NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumImproperParams];
      nvals[i] = new int[NumImproperParams];
      deltavals[i] = new Real[NumImproperParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    msg->get(NumImproperParams,i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->get(NumImproperParams,kvals[i]);
      msg->get(NumImproperParams,nvals[i]);
      msg->get(NumImproperParams,deltavals[i]);
    }

    for (i=0; i<NumImproperParams; i++)
    {
      improper_array[i].multiplicity = i1[i];

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        improper_array[i].values[j].k = kvals[j][i];
        improper_array[i].values[j].n = nvals[j][i];
        improper_array[i].values[j].delta = deltavals[j][i];
      }
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Get the vdw parameters
  msg->get(NumVdwParams);
  msg->get(NumVdwParamsAssigned);

  if (NumVdwParams)
  {
          atomTypeNames = new char[NumVdwParams*(MAX_ATOMTYPE_CHARS+1)];
    vdw_array = new VdwValue[NumVdwParams];
    a1 = new Real[NumVdwParams];
    a2 = new Real[NumVdwParams];
    a3 = new Real[NumVdwParams];
    a4 = new Real[NumVdwParams];

    if ( (vdw_array==NULL) || (a1==NULL) || (a2==NULL) || (a3==NULL)
             || (a4==NULL) || (atomTypeNames==NULL))
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(NumVdwParams * (MAX_ATOMTYPE_CHARS+1), atomTypeNames);
    msg->get(NumVdwParams, a1);
    msg->get(NumVdwParams, a2);
    msg->get(NumVdwParams, a3);
    msg->get(NumVdwParams, a4);

    for (i=0; i<NumVdwParams; i++)
    {
      vdw_array[i].sigma = a1[i];
      vdw_array[i].epsilon = a2[i];
      vdw_array[i].sigma14 = a3[i];
      vdw_array[i].epsilon14 = a4[i];
    }

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }
  
  //  Get the vdw_pair_parameters
  msg->get(NumVdwPairParams);
  
  if (NumVdwPairParams)
  {
    a1 = new Real[NumVdwPairParams];
    a2 = new Real[NumVdwPairParams];
    a3 = new Real[NumVdwPairParams];
    a4 = new Real[NumVdwPairParams];
    i1 = new int[NumVdwPairParams];
    i2 = new int[NumVdwPairParams];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) || (a4 == NULL) || 
         (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }
    
    msg->get(NumVdwPairParams, i1);
    msg->get(NumVdwPairParams, i2);
    msg->get(NumVdwPairParams, a1);
    msg->get(NumVdwPairParams, a2);
    msg->get(NumVdwPairParams, a3);
    msg->get(NumVdwPairParams, a4);
    
    for (i=0; i<NumVdwPairParams; i++)
    {
      new_node = (IndexedVdwPair *) malloc(sizeof(IndexedVdwPair));
      
      if (new_node == NULL)
      {
         NAMD_die("memory allocation failed in Parameters::receive_Parameters");
      }
      
      new_node->ind1 = i1[i];
      new_node->ind2 = i2[i];
      new_node->A = a1[i];
      new_node->A14 = a2[i];
      new_node->B = a3[i];
      new_node->B14 = a4[i];
      new_node->left = NULL;
      new_node->right = NULL;
      
      vdw_pair_tree = add_to_indexed_vdw_pairs(new_node, vdw_pair_tree);
    }
    
    delete [] i1;
    delete [] i2;
    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }
  
  // receive the hydrogen bond parameters
  hbondParams.receive_message(msg);

  AllFilesRead = TRUE;

  delete msg;
}
/*      END OF FUNCTION receive_Parameters    */

/************************************************************************/
/*                  */
/*      FUNCTION convert_vdw_pairs      */
/*                  */
/*  This function converts the linked list of vdw_pairs indexed by  */
/*  atom name into a binary search tree of parameters stored by vdw     */
/*  type index.  This tree is what will be used for real when searching */
/*  for parameters during the simulation.        */
/*                  */
/************************************************************************/

void Parameters::convert_vdw_pairs()
   
{
   Atom atom_struct;    //  Dummy structure for getting indexes
   Index index1, index2;  //  Indexes for the two atoms
   IndexedVdwPair *new_node;  //  New node for tree
   struct vdw_pair_params *ptr, *next;  //  Pointers for traversing list
   
   ptr = vdw_pairp;
   
   //  Go down then entire list and insert each node into the 
   //  binary search tree
   while (ptr != NULL)
   {
      new_node = (IndexedVdwPair *) malloc(sizeof(IndexedVdwPair));
      
      if (new_node == NULL)
      {
   NAMD_die("memory allocation failed in Parameters::convert_vdw_pairs");
      }
      
      //  Get the vdw indexes for the two atoms.  This is kind of a hack
      //  using the goofy Atom structure, but hey, it works
      assign_vdw_index(ptr->atom1name, &atom_struct);
      index1 = atom_struct.vdw_type;
      assign_vdw_index(ptr->atom2name, &atom_struct);
      index2 = atom_struct.vdw_type;
      
      if (index1 > index2)
      {
   new_node->ind1 = index2;
   new_node->ind2 = index1;
      }
      else
      {
   new_node->ind1 = index1;
   new_node->ind2 = index2;
      }
           
      new_node->A = ptr->A;
      new_node->B = ptr->B;
      new_node->A14 = ptr->A14;
      new_node->B14 = ptr->B14;
      
      new_node->left = NULL;
      new_node->right = NULL;
      
      //  Add it to the tree
      vdw_pair_tree = add_to_indexed_vdw_pairs(new_node, vdw_pair_tree);
      
      //  Free the current node and move to the next
      next = ptr->next;
      
      delete ptr;
      
      ptr = next;
   }
   
   vdw_pairp = NULL;
}
/*      END OF FUNCTION convert_vdw_pairs    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_indexed_vdw_pairs    */
/*                  */
/*   INPUTS:                */
/*  new_node - new node to be added to the tree      */
/*  tree - tree to add the node to          */
/*                  */
/*  This is a recursive function that adds a node to the    */
/*   binary search tree of vdw_pair parameters        */
/*                  */
/************************************************************************/

IndexedVdwPair *Parameters::add_to_indexed_vdw_pairs(IndexedVdwPair *new_node,
                 IndexedVdwPair *tree)
   
{
   if (tree == NULL)
      return(new_node);
   
   if ( (new_node->ind1 < tree->ind1) || 
        ((new_node->ind1 == tree->ind1) && (new_node->ind2 < tree->ind2)) )
   {
      tree->left = add_to_indexed_vdw_pairs(new_node, tree->left);
   }
   else
   {
      tree->right = add_to_indexed_vdw_pairs(new_node, tree->right);
   }
   
   return(tree);
}
/*      END OF FUNCTION add_to_indexed_vdw_pairs  */

/************************************************************************/
/*                  */
/*      FUNCTION vdw_pair_to_arrays      */
/*                  */
/*   INPUTS:                */
/*  ind1_array - Array of index 1 values        */
/*  ind2_array - Array of index 2 values        */
/*  A - Array of A values            */
/*  A14 - Array of A 1-4 values          */
/*  B - Array of B values            */
/*  B14 - Array of B 1-4 values          */
/*  arr_index - current position in arrays        */
/*  tree - tree to traverse            */
/*                  */
/*  This is a recursive function that places all the entries of     */
/*   the tree passed in into arrays of values.  This is done so that    */
/*   the parameters can be sent from the master node to the other       */
/*   nodes.                */
/*                  */
/************************************************************************/

int Parameters::vdw_pair_to_arrays(int *ind1_array, int *ind2_array,
          Real *A, Real *A14,
          Real *B, Real *B14,
          int arr_index, IndexedVdwPair *tree)
      
{
   if (tree == NULL)
      return(arr_index);
   
   ind1_array[arr_index] = tree->ind1;
   ind2_array[arr_index] = tree->ind2;
   A[arr_index] = tree->A;
   A14[arr_index] = tree->A14;
   B[arr_index] = tree->B;
   B14[arr_index] = tree->B14;
   
   arr_index++;
   
   arr_index = vdw_pair_to_arrays(ind1_array, ind2_array, A, A14, B, B14,
          arr_index, tree->left);
   arr_index = vdw_pair_to_arrays(ind1_array, ind2_array, A, A14, B, B14,
          arr_index, tree->right);
   
   return(arr_index);
}
/*      END OF FUNCTION vdw_pair_to_arrays    */

/***************************************************************************
 * RCS INFORMATION:
 *
 *  $RCSfile: Parameters.C,v $
 *  $Author: milind $  $Locker:  $    $State: Exp $
 *  $Revision: 1.1005 $  $Date: 1997/10/01 16:47:00 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Parameters.C,v $
 * Revision 1.1005  1997/10/01 16:47:00  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1004  1997/04/07 14:54:33  nealk
 * Changed fclose() to Fclose() (found in common.[Ch]) to use with popen().
 * Also corrected compilation warnings in Set.[Ch].
 *
 * Revision 1.1003  1997/04/03 19:59:10  nealk
 * 1) New Fopen() which handles .Z and .gz files.
 * 2) localWaters and localNonWaters lists on each patch.
 *
 * Revision 1.1002  1997/03/20 15:38:30  nealk
 * Added "\n" before endi in iout statements.
 * Removed endi from DebugM() calls.
 *
 * Revision 1.1001  1997/03/19 11:54:42  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:58:57  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:05  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/24 02:29:53  jim
 * Fixed bug where only first parameter file was read!
 * Added files for hydrogen bond parameter reading.
 *
 * Revision 1.777  1997/01/17 19:36:40  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.5  1996/11/11 19:54:09  nealk
 * Modified to use InfoStream instead of Inform.
 *
 * Revision 1.4  1996/10/31 20:42:40  jim
 * small changes to support LJTable
 *
 * Revision 1.3  1996/08/16 04:55:30  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 04:39:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.15  1996/04/18 18:46:18  billh
 * Updated to read hydrogen bond information.
 *
 * Revision 1.14  1995/10/10 02:57:28  hazen
 * Updated memory allocation to use C++ new/delete
 *
 * Revision 1.13  1995/04/06  13:23:16  nelson
 * Fixed mistype of "struct IndexedVdw ..."
 *
 * Revision 1.12  95/03/08  14:31:11  14:31:11  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.11  95/02/01  14:15:53  14:15:53  nelson (Mark T. Nelson)
 * Replaced references to namdDebug with DEBUG_MSG
 * 
 * Revision 1.10  95/01/26  14:31:14  14:31:14  nelson (Mark T. Nelson)
 * Added Charm22 parameters
 * 
 * Revision 1.9  94/09/28  14:35:19  14:35:19  nelson (Mark T. Nelson)
 * Changed so that delta values for Impropers and Dihedrals are
 * converted to radians
 * 
 * Revision 1.8  94/09/27  12:05:55  12:05:55  nelson (Mark T. Nelson)
 * Added checks for 0 number of any parameters
 * 
 * Revision 1.7  94/09/24  20:12:17  20:12:17  nelson (Mark T. Nelson)
 * added routines to deal with vdw pairs correctly
 * 
 * Revision 1.6  94/09/22  09:55:49  09:55:49  nelson (Mark T. Nelson)
 * Changed it so that angle parameter theta0 was converted to radians
 * 
 * Revision 1.5  94/09/12  17:57:02  17:57:02  gursoy (Attila Gursoy)
 * receive-message is moved to the Node object for charm++ integration
 * 
 * Revision 1.4  94/08/02  16:28:50  16:28:50  nelson (Mark T. Nelson)
 * Implemented send_Parameters and receive_Parameters
 * 
 * Revision 1.3  94/07/07  13:30:48  13:30:48  nelson (Mark T. Nelson)
 * changes to make parameter access index based
 * 
 * Revision 1.2  94/06/24  03:10:27  03:10:27  billh (Bill Humphrey)
 * Replace NAMD_warn calls with usage of namdWarn object.
 * 
 * Revision 1.1  94/06/22  15:05:02  15:05:02  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/
