#!/usr/bin/env perl
use strict;
foreach my $module ( @ARGV ) {
  eval "require $module";
  printf( "%-20s: %s\n", $module, $module->VERSION ) unless ( $@ );
}
