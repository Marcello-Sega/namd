
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
    crease_nodelist_file $nfile $hlist
    set nprocs [llength $hlist]
    puts "SPAWNING $conf on $host"
    lappend channels [open \
      "| ssh $host \"$namd ++nodelist $nfile +p$nprocs $conf > $log\" << {} |& cat" "r"]
  }
  return $channels
}

namespace export spawn_namd_simple spawn_namd_ssh spawn_namd_parallel

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

