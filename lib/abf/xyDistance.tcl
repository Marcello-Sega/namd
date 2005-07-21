#########################################################
# ABF procedures for a distance between two atom groups #
#########################################################

set ABFcoordID "2D-distance between COM of two atom groups"

# Define coordinate-specific optional parameters with default values
array set ABFoptions {
temp		300.0
dxi		0.1
dSmooth		0.2
}

###############################################################
# ABFstartup : declares atoms whose coordinates are requested #
###############################################################

proc ABFstartup {} {

    namespace eval ABFcoord {
	set abf1 $::ABF::abf1
	set abf2 $::ABF::abf2

	# we need this for 'loadtotalforces'
	foreach a $abf1 {addatom $a}
	foreach a $abf2 {addatom $a}

	# this one is convenient for 'loadcoords' and 'addforce'
	set g1 [ addgroup $abf1 ]
	set g2 [ addgroup $abf2 ]
    }
}

	
################################################################
# ABFcoord : reads coord, returns value of reaction coordinate #
################################################################

proc ABFcoord {} {
    namespace eval ABFcoord {
	loadcoords coords

	foreach {vx vy vz} [vecsub $coords($g2) $coords($g1)] {}

	set r [expr {sqrt($vx*$vx + $vy*$vy)} ]

	return $r
    }
}

############################################################
# ABForce : returns force along reaction coordinate        #
############################################################

proc ABForce {} {
    namespace eval ABFcoord {
	
	# Assume that ABFcoord has already been called
	# vx vy vz is set

	set nv	[vecnorm "$vx $vy 0."] ;# unity vector

	loadtotalforces forces

	set f1 0.0
	foreach a $abf1 { set f1 [expr {$f1 + [vecdot $forces($a) $nv]} ]}
	set f2 0.0
	foreach a $abf2 { set f2 [expr {$f2 + [vecdot $forces($a) $nv]} ]}

	return [expr ($f2 - $f1) / 2.0]
    }
}


###############################################################################
# ABFapply : applies the force given as a parameter along reaction coordinate #
###############################################################################

proc ABFapply {type force} {
    set ABFcoord::type $type
    set ABFcoord::force $force

    namespace eval ABFcoord {
	# We wouldn't need type, if we didn't cheat on the Jacobian term

	# Assume that ABFcoord and ABFthermoForce have been called
	# i.e. {vx vy vz} and r are set

	set nv	[vecnorm "$vx $vy 0."] ;# unity vector

	if { $type == "bias" } {
		# compensate for the Jacobian term kT/r
		set force [expr {$force - 0.001986 * $::ABF::temp / $r}]
	}

	set F2 [vecscale $force $nv]

	addforce $g1 [vecinvert $F2]
	addforce $g2 $F2

	return
    }
}
