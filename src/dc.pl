#!/usr/local/bin/perl

# input set of strings to match
@dontprint = @ARGV;
@ARGV=();

$n = @dontprint;

# read from stdin

while (<>) {
# if line ends with : it is the start of a dependency
  chop;
  if ( ($target,$other) = /([a-zA-Z_]*\.o:)(.*)$/ ) {
    print $target;
    $go=1;
    $first=1;
    while ($go) {
	if ($first) {
	  $_ = $other;
	} else {
	  $_ = <> || last;
	  chop;
	}
	$first = 0;

	if ( /\\$/ ) {
	  chop;
	  $go = 1;
	} else {
	  $go = 0;
	}

	@files = split;
	foreach $word (@files) {
	  $bad = 0;
	  foreach $notword (@dontprint) {
	    if ( $word =~ /$notword/ ) {
	      $bad = 1;
	      last;
	    }
	  }
	  if ( ! $bad ) {
	      print " \\\n";
	      print "	",$word;
	  }
	}
    }
  }
  print "\n";
}

# chop off end of line

# split line into filenames

# go thru filenames and print ones that don't match

# continue to next line

# continue to next dependency

# end


