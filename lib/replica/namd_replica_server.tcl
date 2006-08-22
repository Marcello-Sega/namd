
namespace eval namd_replica_server {

  variable verbose 1

  proc start_replicas {cmd host port num filebase} {
    variable replica_response_field null
    variable num_replicas $num
    start_server $port
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
    set confname "${logname}.nrc"
    set f [open $confname "w"]
    puts $f "set replica_id $replica_id"
    puts $f "set server_host $server_host"
    puts $f "set server_port $server_port"
    puts $f {
      set verbose 1
      set server_channel [socket $server_host $server_port]
      fconfigure $server_channel -buffering line
      print "SERVER CONNECTED"
      puts $server_channel "replica $replica_id"
      set server_buffer {}
      while { 1 } {
        while [catch {gets $server_channel cmd} rval] {
puts $errorCode
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

  proc replica_eval {field script} {
    variable replica_response_field $field
    variable num_replicas
    variable replica_data
    variable verbose
    if {$verbose > 0} {
      puts "SERVER: $script"
    }
    for { set i 0 } { $i < $num_replicas } { incr i } {
      puts $replica_data($i.channel) "eval \{$script\n\}"
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
      puts $replica_data($i.channel) "eval \{set $var $replica_data($i.$field)\n\}"
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

