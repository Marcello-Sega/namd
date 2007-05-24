
namespace eval spawn_namd {

proc spawn_namd_simple {namd conffiles logfiles} {
  foreach conf $conffiles log $logfiles {
    lappend channels [open \
      "| csh -c \"$namd $conf > $log\" << {} |& cat" "r"]
  }
  return $channels
}

proc spawn_namd_ssh {namd hosts conffiles logfiles} {
  set i 0
  set nhosts [llength $hosts]
  foreach conf $conffiles log $logfiles {
    set host [lindex $hosts $i]
    incr i
    if {$i >= $nhosts} {set i 0}
    puts "SPAWNING $conf on $host"
    lappend channels [open \
      "| ssh $host \"$namd $conf > $log\" << {} |& cat" "r"]
  }
  return $channels
}

proc split_hosts {njobs hosts procs_per_host} {
  # put host that script is running on at end of host list
  # and put an entry in list for every processors on a host
  set myhost [lindex [split [info hostname] "."] 0]
  set newhosts {}
  set selfhosts {}
  foreach host $hosts {
    set shorthost [lindex [split $host "."] 0]
    if [string compare $shorthost $myhost] {
      set list newhosts
    } else {
      set list selfhosts
    }
    for {set i 0} {$i < $procs_per_host} {incr i} {
      lappend $list $host
    }
  }
  foreach host $selfhosts { lappend newhosts $host }
  set hosts $newhosts

  # break new host list evenly among jobs
  set nhosts [llength $hosts]
  set hosts_per_job [expr $nhosts/$njobs]
  for {set i 0} {$i < $njobs} {incr i} {
    set first [expr $i*$hosts_per_job]
    set last [expr $first + $hosts_per_job - 1]
    lappend host_lists [lrange $hosts $first $last]
  }

  return $host_lists
}

proc create_nodelist_file {filename hosts} {
  set file [open $filename "w"]
  puts $file "group main"
  foreach host $hosts {
    puts $file "host $host"
  }
  close $file
}

proc spawn_namd_parallel {namd hosts procs_per_host conffiles logfiles} {
  set njobs [llength $conffiles]
  set host_lists [split_hosts $njobs $hosts $procs_per_host]
  foreach hlist $host_lists conf $conffiles log $logfiles {
    set host [lindex $hlist 0]
    set nfile "$conf.nodelist"
    create_nodelist_file $nfile $hlist
    set nprocs [llength $hlist]
    puts "SPAWNING $conf on $host"
    lappend channels [open \
      "| ssh $host \"$namd ++nodelist $nfile +p$nprocs $conf > $log\" << {} |& cat" "r"]
  }
  return $channels
}

proc pdb_nidlist_to_tcl {pbslist} { ;# convert 15:17..19 to {15 17 18 19}
  foreach range [split $pbslist :] {
    switch [scan $range {%d..%d} lower upper] {
      2 { for { } { $lower <= $upper } { incr lower } { lappend nids $lower } }
      1 { lappend nids $lower }
      default { error "bad nidlist element $range" }
    }
  }
  return $nids
}

proc tcl_nidlist_to_yod {tcllist} { ;# convert {15 17 18 19} to 15,17..19
  set lower [lindex $tcllist 0]
  set upper $lower
  foreach nid [lrange $tcllist 1 end] {
    if { $nid == $upper + 1 } {
      set upper $nid
    } else {
      if { $lower == $upper } {
        lappend ranges $lower
      } else {
        lappend ranges "$lower..$upper"
      }
      set lower $nid
      set upper $lower
    }
  }
  if { $lower == $upper } {
    lappend ranges $lower
  } else {
    lappend ranges "$lower..$upper"
  }
  return [join $ranges ,]
}

proc spawn_namd_crayxt {namd conffiles logfiles} {
  global env
  set njobs [llength $conffiles]
  set procs_per_node [expr $env(PBS_NPROCS) / $env(PBS_NNODES)]
  set nids [pdb_nidlist_to_tcl $env(YOD_NIDLIST)]
  set nnodes [llength $nids]
  if { $nnodes != $env(PBS_NNODES) } {
    error "PBS_NNODES=$env(PBS_NNODES) does not match node count in YOD_NIDLIST=$env(YOD_NIDLIST)"
  }
  set nodes_per_job [expr $nnodes/$njobs]
  if { $nodes_per_job < 1 } {
    error "only $nnodes nodes for $njobs jobs"
  }
  for {set i 0} {$i < $njobs} {incr i} {
    set first [expr $i*$nodes_per_job]
    set last [expr $first + $nodes_per_job - 1]
    lappend node_lists [lrange $nids $first $last]
  }
  # set procs_per_job [expr $procs_per_node * $nodes_per_job]
  foreach nidlist $node_lists conf $conffiles log $logfiles {
    set yodlist [tcl_nidlist_to_yod $nidlist]
    puts "SPAWNING $conf on $yodlist"
    lappend channels [open \
      "| csh -c \"yod -list $yodlist $namd $conf > $log\" << {} |& cat" "r"]
  }
  return $channels
}

namespace export spawn_namd_simple spawn_namd_ssh spawn_namd_parallel \
	 spawn_namd_crayxt

}

# examples below

# set replogfile ${filebase}.%d.log
# set port_number 3169

# start_replicas \
# 	[list spawn_namd_simple "namd2 ++local +netpoll +p1"] \
# 	localhost $port_number \
# 	$num_replicas $replogfile

# start_replicas \
# 	[list spawn_namd_ssh "cd [pwd]; namd2 ++local +netpoll +p1" \
# 		$env(LSB_HOSTS) ] \
# 	[info hostname] $port_number \
# 	$num_replicas $replogfile

# start_replicas \
# 	[list spawn_namd_ssh "cd [pwd]; namd2 ++local +netpoll +p1" \
# 		[read [open $env(HOST_FILE) "r"]] ] \
# 	[info hostname] $port_number \
# 	$num_replicas $replogfile

# start_replicas \
# 	[list spawn_namd_ssh "cd [pwd]; namd2 ++local +netpoll +p1" \
# 		[list budapest istanbul kiev lyon] ] \
# 	[info hostname] $port_number \
# 	$num_replicas $replogfile

# set bindir $env(HOME)/NAMD/
# start_replicas \
# 	[list spawn_namd_parallel
# 		"cd [pwd]; $bindir/charmrun $bindir/namd2 +netpoll" \
# 		$env(LSB_HOSTS) 1 ] \
# 	[info hostname] $port_number \
# 	$num_replicas $replogfile

# run on Cray XT3
# set namd_bin_dir /usr/users/7/jphillip/NAMD_2.6_CRAY-XT3/
# set server_port 0  ;# use file-based transport
# start_replicas \
#	[list spawn_namd_crayxt \
#		"-small_pages [file join $namd_bin_dir namd2]" ] \
#	localhost 0 \
#	$num_replicas $replogfile

