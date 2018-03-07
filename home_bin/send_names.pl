#!/usr/bin/env perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use lib '/data/web/isga.cgb/docs/lib/perl5';

use strict;
use warnings;
use ISGA;
use Data::Dumper;
my %names = (Rick => 'Rick', Debbie => 'Debbie', Nicki => 'Nicki', Rob => 'Rob', Kayce => 'Kayce', Aaron => 'Aaron');
my %emails = ( 'Aaron' => 'aaron.buechlein@gmail.com', 
               'Rick' => 'rbuechlein@gmail.com',
               'Debbie' => 'dbuechlein@myvu.vinu.edu',
               'Nicki' => 'tela78@yahoo.com',
               'Rob' => 'robert.m.hoyt@gmail.com',
               'Kayce' => 'kayce.reed@gmail.com'
);
my $support_email = 'abuechle@cgb.indiana.edu';
my @hat = ('Rick', 'Debbie', 'Nicki', 'Rob', 'Kayce', 'Aaron');

while( $names{Rick} eq 'Rick' or $names{Rick} eq 'Debbie' or
       $names{Debbie} eq 'Rick' or $names{Debbie} eq 'Debbie' or
       $names{Nicki} eq 'Rob' or $names{Nicki} eq 'Nicki' or
       $names{Rob} eq 'Rob' or $names{Rob} eq 'Nicki' or
       $names{Kayce} eq 'Aaron' or $names{Kayce} eq 'Kayce' or
       $names{Aaron} eq 'Aaron' or $names{Aaron} eq 'Kayce' ){

 for(my $i = 0; $i < 100; $i++){
   &fisher_yates_shuffle(\@hat);
 }
 $names{Rick} = $hat[0];
 $names{Debbie} = $hat[1];
 $names{Nicki} = $hat[2];
 $names{Rob} = $hat[3];
 $names{Kayce} = $hat[4];
 $names{Aaron} = $hat[5];

print "Rick  $names{Rick}\n" if $names{Rick} eq 'Rick' or $names{Rick} eq 'Debbie';
print "Debbie  $names{Debbie}\n" if $names{Debbie} eq 'Rick' or $names{Debbie} eq 'Debbie';
print "Nicki  $names{Nicki}\n" if $names{Nicki} eq 'Rob' or $names{Nicki} eq 'Nicki';
print "Rob  $names{Rob}\n" if $names{Rob} eq 'Rob' or $names{Rob} eq 'Nicki';
print "Kayce  $names{Kayce}\n" if $names{Kayce} eq 'Aaron' or $names{Kayce} eq 'Kayce';
print "Aaron  $names{Aaron}\n" if $names{Aaron} eq 'Aaron' or $names{Aaron} eq 'Kayce';
print "\n\n";
}

#print Dumper(%names);

foreach my $email (keys %emails){
print "$email  $emails{$email}\n";
  eval {
    my %mail =
               ( To => $emails{$email},
                 From => "Aaron Buechlein <$support_email>",
                 Subject => 'Buechlein Christmas Name Draw',
                 Message => "Hey $email,
The name you drew in the gift draw was $names{$email}." );
    Mail::Sendmail::sendmail(%mail) 
      or X::Mail::Send->throw( text => $Mail::Sendmail::error, message => \%mail ); 
  };
  
  if ( $@ ) {
    print "Failed to send email because: $@";
  }
}
exit;
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

