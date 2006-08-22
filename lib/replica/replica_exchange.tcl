
if {! [info exists libdir]} { set libdir [file dirname [info script]] }
# If it fails, try the local directory
if { $libdir == "" } { set libdir "." }

foreach file $argv { source $file }

if {! [info exists i_job]} { set i_job 0 }
set job_output_root "$output_root.job$i_job"

if {! [array exists sorted_data]} {
  set i_run 0
  set i_step 0
  for {set i 0} {$i < $num_replicas} {incr i} {
    set sorted_data($i.replica_id) $i
    set temp [expr round($min_temp * \
	 exp( log(1.0*$max_temp/$min_temp)*(1.0*$i/($num_replicas-1)) ) )]
    puts "TEMPERATURE $i: $temp"
    set sorted_data($i.temperature) $temp
    set sorted_data($i.exchanges_attempted) 0
    set sorted_data($i.exchanges_accepted) 0
  }
}

source [file join $libdir namd_replica_server.tcl]
namespace import namd_replica_server::*
upvar 0 namd_replica_server::replica_data replica_data

set replogfile ${job_output_root}.%d.log
if {! [info exists server_port]} { set server_port 3169 }

source [file join $libdir spawn_namd.tcl]
namespace import spawn_namd::*

start_replicas $spawn_namd_command \
	[info hostname] $server_port \
	$num_replicas $replogfile

replica_eval dir "cd [pwd]; pwd"

for {set i 0} {$i < $num_replicas} {incr i} {
  set rid $sorted_data($i.replica_id)
  set replica_data($rid.temperature) $sorted_data($i.temperature)
}
replica_push temperature NEWTEMP

replica_eval null "set base_seed [expr int(100000*rand())]"
replica_eval null "set job_output_root $job_output_root"
replica_eval null {

  proc save_callback {labels values} {
    global saved_labels saved_values
    set saved_labels $labels
    set saved_values $values
  }
  callback save_callback

  proc save_array {} {
    global saved_labels saved_values saved_array
    foreach label $saved_labels value $saved_values {
      set saved_array($label) $value
    }
  }

  seed [expr int(0*srand($base_seed + 100*$replica_id) + 100000*rand())]
  langevin on
  langevinDamping 10.0
  langevinTemp $NEWTEMP
  outputname $job_output_root.$replica_id
}

if {$i_run} { #restart
  replica_eval null "set restart_root $restart_root"
  replica_eval null "firsttimestep $i_step"
  replica_eval null {
    bincoordinates $restart_root.$replica_id.coor
    binvelocities $restart_root.$replica_id.vel
    extendedSystem $restart_root.$replica_id.xsc
    unset restart_root
  }
} else {
  replica_eval null {
    temperature $NEWTEMP
  }
}

replica_eval null "outputEnergies [expr $steps_per_run / 10]"
replica_eval null "dcdFreq [expr $steps_per_run * $runs_per_frame]"
replica_eval null "source $namd_config_file; run 0; save_array"

replica_eval TEMP {expr $saved_array(TEMP)}
replica_eval POTENTIAL {expr $saved_array(TOTAL) - $saved_array(KINETIC)}

set potenergy_file [open "${job_output_root}.potenergy.dat" "w"]
fconfigure $potenergy_file -buffering line
set targtemp_file [open "${job_output_root}.targtemp.dat" "w"]
fconfigure $targtemp_file -buffering line
set realtemp_file [open "${job_output_root}.realtemp.dat" "w"]
fconfigure $realtemp_file -buffering line

proc array_to_list {arrayname field} {
  upvar $arrayname array
  global num_replicas
  for {set i 0} {$i < $num_replicas} {incr i} {
    lappend list $array($i.$field)
  }
  return $list
}

while {$i_run < $num_runs} {
  replica_eval null "run $steps_per_run; save_array"
  incr i_step $steps_per_run
  replica_eval TEMP {expr $saved_array(TEMP)}
  replica_eval POTENTIAL {expr $saved_array(TOTAL) - $saved_array(KINETIC)}
  puts $potenergy_file "$i_step [array_to_list replica_data POTENTIAL]"
  puts $targtemp_file "$i_step [array_to_list replica_data temperature]"
  puts $realtemp_file "$i_step [array_to_list replica_data TEMP]"
  for {set i [expr $i_run % 2]} {$i+1 < $num_replicas} {incr i 2} {
    set i2 [expr $i + 1]
    set rid $sorted_data($i.replica_id)
    set rid2 $sorted_data($i2.replica_id)
    set temp $sorted_data($i.temperature)
    set temp2 $sorted_data($i2.temperature)
    set BOLTZMAN 0.001987191
    set dbeta [expr ((1.0/$temp) - (1.0/$temp2)) / $BOLTZMAN]
    set pot $replica_data($rid.POTENTIAL)
    set pot2 $replica_data($rid2.POTENTIAL)
    set delta [expr $dbeta * ($pot2 - $pot)]
    if {$delta < 0. || exp(-1. * $delta) > rand()} {
      puts "EXCHANGE_ACCEPT $rid ($temp) $rid2 ($temp2)"
      set sorted_data($i.replica_id) $rid2
      set sorted_data($i2.replica_id) $rid
      set replica_data($rid.temperature) $temp2
      set replica_data($rid2.temperature) $temp
      incr sorted_data($i.exchanges_accepted)
    }
    incr sorted_data($i.exchanges_attempted)
  }
  replica_eval null {set OLDTEMP $NEWTEMP}
  replica_push temperature NEWTEMP
  replica_eval null {
    rescalevels [expr sqrt(1.0*$NEWTEMP/$OLDTEMP)]
    langevinTemp $NEWTEMP
  }

  incr i_run

  if { $i_run % ($runs_per_frame * $frames_per_restart) == 0 ||
	$i_run == $num_runs } {  # restart
    set restart_root "$job_output_root.restart$i_run"
    replica_eval null "output $restart_root.\$replica_id"
    set rfile [open "$restart_root.tcl" "w"]
    puts $rfile [list set i_job [expr $i_job + 1]]
    puts $rfile [list set i_run $i_run]
    puts $rfile [list set i_step $i_step]
    puts $rfile [list set restart_root $restart_root]
    puts $rfile [list array set sorted_data [array get sorted_data]]
    close $rfile
    if [info exists old_restart_root] {
      file delete $old_restart_root.tcl
      replica_eval null "
        file delete $old_restart_root.\$replica_id.coor
        file delete $old_restart_root.\$replica_id.vel
        file delete $old_restart_root.\$replica_id.xsc"
    }
    set old_restart_root $restart_root
  }
}

for {set i 0} {$i < $num_replicas} {incr i} {
  set attempts $sorted_data($i.exchanges_attempted)
  if $attempts {
    set i2 [expr $i + 1]
    set temp $sorted_data($i.temperature)
    set temp2 $sorted_data($i2.temperature)
    set accepts $sorted_data($i.exchanges_accepted)
    set ratio [expr 1.0*$accepts/$attempts]
    puts "EXCHANGE_RATIO $temp $temp2 $accepts $attempts $ratio"
  }
}

replica_eval exit {exit}

