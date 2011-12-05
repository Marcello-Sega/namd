
This directory contains Tcl scripts that implement replica exchange.
This replaces the old Tcl server and socket connections driving a
separate NAMD process for every replica used in the simulation.

*** NAMD based on a patched MPI build of Charm++ is required! ***

Replica exchanges and energies are recorded in the .history files
written in the output directories.  These can be viewed with, e.g.,
"xmgrace output/*/*.history" and processed via awk or other tools.
There is also a script to load the output into VMD and color each
frame according to replica index.  An example simulation folds
a 66-atom model of a deca-alanine helix in about 10 ns.

replica.namd - master script for replica exchange simulations
  to run: mkdir output
          (cd output; mkdir 0 1 2 3 4 5 6 7)
          mpirun namd2 +replicas 8 fold_alanin.conf --source ../replica.namd +stdout output/%d/job0.%d.log
          mpirun namd2 +replicas 8 restart_1.conf --source ../replica.namd +stdout output/%d/job1.%d.log

The number of MPI ranks must be a multiple of the number of replicas (+replicas).
Increase num_runs below if job completed, or use latest restartXXX.tcl file available.
Be sure to increment jobX for +stdout option on command line.

show_replicas.vmd - script for loading replicas into VMD, first source
    the replica exchange conf file and then this script, repeat for
    restart conf file or for example just do vmd -e load_all.vmd

clone_reps.vmd - provides "clone_reps" commmand to copy graphical
  representation from top molecule to all others

example subdirectory:

alanin_base.namd - basic config options for NAMD
alanin.params - parameters
alanin.psf - structure
unfolded.pdb - initial coordinates
alanin.pdb - folded structure for fitting in show_replicas.vmd

fold_alanin.conf - config file for replica_exchange.tcl script
restart_1.conf - config file to continue alanin folding another 10 ns

load_all.vmd - load all output into VMD and color by replica index

