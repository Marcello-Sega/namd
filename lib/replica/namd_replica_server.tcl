
namespace eval namd_replica_server {

  variable verbose 1
  variable use_sockets 1

  proc start_replicas {cmd host port num filebase} {
    variable replica_response_field null
    variable num_replicas $num
    if { $port } { start_server $port } else { variable use_sockets 0 }
    for { set i 0 } { $i < $num_replicas } { incr i } {
      foreach {confname logname} [write_config $filebase $i $host $port] {}
      lappend conflist $confname
      lappend loglist $logname
    }
    set cidlist [eval $cmd [list $conflist $loglist]]
    set i 0
    foreach cid $cidlist {
      fileevent $cid readable [namespace code "stderr_handle $i $cid"]
      incr i
    }
    if {$i != $num_replicas} {
      error "$cmd returned fewer channels than requested"
    }
    wait_all
  }

  proc stderr_handle {replica_id cid} {
    if {[gets $cid msg] < 0} {
      puts "disconnect errpipe $replica_id"
      fileevent $cid readable {}
    } else {
      puts "errpipe $replica_id: $msg"
    }
  }


  proc write_config {filebase replica_id server_host server_port} {
    set logname "[format $filebase $replica_id]"
    if [file exists $logname] {error "log file $logname already exists"}
    if { ! $server_port } {
      variable replica_data
      set replica_data($replica_id.ichannel) "${logname}.ichannel"
      poll_file_channel $replica_id "${logname}.ochannel"
    }
    set confname "${logname}.nrc"
    set f [open $confname "w"]
    puts $f "set replica_id $replica_id"
    puts $f "set verbose 1"
    if { ! $server_port } {
    puts $f "set logname [file tail $logname]"
    puts $f {
      set logname [file join [pwd] $logname]
      set f [open "${logname}.ochannel.tmp" "w"]
      puts $f "replica $replica_id"
      close $f
      file rename "${logname}.ochannel.tmp" "${logname}.ochannel"
      print "${logname}.ochannel written"
      while { 1 } {
        while { ! [file exists "${logname}.ichannel"] } { after 100 }
        set f [open "${logname}.ichannel" "r"]
        set cmd [read -nonewline $f]
        close $f
        file delete "${logname}.ichannel"
        if {$verbose > 0} { print "SERVER: $cmd" }
        set response [eval $cmd]
        if {$verbose > 0} { print "CLIENT: $response" }
        set f [open "${logname}.ochannel.tmp" "w"]
        puts $f $response
        close $f
        file rename "${logname}.ochannel.tmp" "${logname}.ochannel"
      }
    }
    } else {
    puts $f "set server_host $server_host"
    puts $f "set server_port $server_port"
    puts $f {
      set server_channel [socket $server_host $server_port]
      fconfigure $server_channel -buffering line
      print "SERVER CONNECTED"
      puts $server_channel "replica $replica_id"
      set server_buffer {}
      while { 1 } {
        while [catch {gets $server_channel cmd} rval] {
          if {![string match "POSIX EINTR*" $errorCode]} {
            error $rval $errorInfo $errorCode
          }
        }
        if {$rval < 0} {
          error "disconnect server"
        }
        append server_buffer $cmd "\n"
        if [info complete $server_buffer] {
          set cmd $server_buffer
          set server_buffer {}
          if {$verbose > 0} { print "SERVER: $cmd" }
          set response [eval $cmd]
          if {$verbose > 0} { print "CLIENT: $response" }
          puts $server_channel $response
        }
      }
    }
    }
    close $f
    return [list $confname $logname]
  }

  variable replica_data

  proc replica_ready {replica_id} {
    variable num_replicas
    variable replica_data
    variable num_replicas_ready
    if $replica_data($replica_id.ready) {
      error "duplicate replica_ready $replica_id"
    }
    set replica_data($replica_id.ready) 1
    incr num_replicas_ready
    if { $num_replicas_ready == $num_replicas } {
      variable all_replicas_barrier_var 0 ; # ends vwait call
    }
  }

  proc wait_all { } {
    variable num_replicas
    variable replica_data
    variable num_replicas_ready 0
    for { set i 0 } { $i < $num_replicas } { incr i } {
      set replica_data($i.ready) 0
    }
    variable all_replicas_barrier_var
    vwait [namespace which -variable all_replicas_barrier_var]
  }


  proc poll_file_channel {replica_id filename} {
    if { [file exists $filename] } {
      variable replica_data
      variable replica_response_field
      set f [open $filename "r"]
      set response [read -nonewline $f]
      close $f
      file delete $filename
      if [catch {replica_ready $replica_id}] {
        error "replica $replica_id unexpected data: $response"
      }
      variable verbose
      if {$verbose > 1} {
        puts "CLIENT: [list $replica_id.$replica_response_field $response]"
      }
      set replica_data($replica_id.$replica_response_field) $response
    }
    after 100 [namespace code "poll_file_channel $replica_id $filename"]
  }

  proc start_server {port} {
    socket -server [namespace code "server_accept"] $port
  }

  proc server_accept {cid addr port} {
    variable num_replicas
    fconfigure $cid -buffering line
    if {[gets $cid request] < 0 ||
        [scan $request "replica %d%c" replica_id garbage] != 1 ||
        $replica_id < 0 || $replica_id >= $num_replicas } {
      puts "bad client initialization: $cid $addr $port $request"
      close $cid
      return
    }
    variable replica_data
    if [info exists replica_data($replica_id.channel)] {
      error "duplicate replica $replica_id $cid $addr $port"
    }
    puts "starting replica $replica_id $cid $addr $port"
    set replica_data($replica_id.channel) $cid
    set replica_data($replica_id.address) $addr
    set replica_data($replica_id.port) $port
    replica_ready $replica_id
    fileevent $cid readable [namespace code "server_handle $replica_id $cid"]
  }

  proc server_handle {replica_id cid} {
    variable replica_data
    set oldcid $replica_data($replica_id.channel)
    if [string compare $cid $oldcid] {
      error "replica $replica_id channel altered: $cid vs $oldcid"
    }
    variable replica_response_field
    if {[gets $cid response] < 0} {
      if {! [string compare $replica_response_field "exit"] } {
        set response {}
        fileevent $cid readable {}
      } else {
        puts "disconnect replica $replica_id"
        fileevent $cid readable {}
        after 1000 [list error "disconnect replica $replica_id"]
        return
      }
    }
    if [catch {replica_ready $replica_id}] {
      error "replica $replica_id unexpected data: $response"
    }
    variable verbose
    if {$verbose > 1} {
      puts "CLIENT: [list $replica_id.$replica_response_field $response]"
    }
    set replica_data($replica_id.$replica_response_field) $response
  }

  proc replica_send {replica_id cmd} {
    variable replica_data
    variable use_sockets
    if { $use_sockets } {
      puts $replica_data($replica_id.channel) $cmd
    } else {
      set filename $replica_data($replica_id.ichannel)
      set f [open "${filename}.tmp" "w"]
      puts $f $cmd
      close $f
      file rename "${filename}.tmp" $filename
    }
  }

  proc replica_eval {field script} {
    variable replica_response_field $field
    variable num_replicas
    variable replica_data
    variable verbose
    if {$verbose > 0} {
      puts "SERVER: $script"
    }
    for { set i 0 } { $i < $num_replicas } { incr i } {
      replica_send $i "eval \{$script\n\}"
    }
    wait_all
  }

  proc replica_push {field var} {
    variable replica_response_field null
    variable num_replicas
    variable replica_data
    variable verbose
    if {$verbose > 0} {
      puts "SERVER: replica_push $field $var"
    }
    for { set i 0 } { $i < $num_replicas } { incr i } {
      replica_send $i "eval \{set $var $replica_data($i.$field)\n\}"
    }
    wait_all
  }

  namespace export start_replicas replica_eval replica_push

}

proc bgerror {error} {
  global errorInfo
  puts "FATAL ERROR: $error"
  foreach {name value} [array get namd_replica_server::replica_data] {
    puts [list $name $value]
  }
  puts $errorInfo
  puts stderr "FATAL ERROR: $error"
  exit -1
}

