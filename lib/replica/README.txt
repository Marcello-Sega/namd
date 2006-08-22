
This directory contains Tcl scripts that implement replica exchange
for NAMD, using a Tcl server and socket connections to drive a
separate NAMD process for every replica used in the simulation.
There is also a script to load the output into VMD and color each
frame according to target temperature.  An example simulation folds
a 66-atom model of a deca-alanine helix in about 10 ns.

replica_exchange.tcl - master script for replica exchange simulations
  to run: tclsh ../replica_exchange.tcl fold_alanin.conf
          tclsh ../replica_exchange.tcl restart_1.conf
  uses:
    namd_replica_server.tcl - general script for driving NAMD slaves
    spawn_namd.tcl - variety of methods for launching NAMD slaves

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

load_all.vmd - load all output into VMD and color by target temperature

