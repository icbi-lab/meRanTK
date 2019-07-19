#!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.18.2/bin/perl -w
#
#  meRanCompare
#
#  Copyright 2019 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

use strict;
use warnings;

use Pod::Usage;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use IO::File;
use File::Basename;

my $countsFile = "";
my $condNames  = "";

my $help    = undef;
my $man     = undef;
my $version = undef;
my $VERSION = "0.1";

my $DEBUG = 0;
my $debug;

GetOptions(
            'count-file|cf=s'      => \$countsFile,
            'condition-names|cn=s' => \$condNames,
            'debug|d'              => \$DEBUG,
            'help|h'               => \$help,
            'man|m'                => \$man,
            'version'              => \$version,
            );

if ($DEBUG) {
    eval "use Data::Dumper";
    $debug = sub { print STDERR "\n" . join( " ", @_ ) . "\n"; }
}
else {
    $debug = sub { return; };
}

usage() and exit(0) if ($help);
pod2usage( -verbose => 3 ) if $man;

say $VERSION and exit(0) if ($version);

if ( !-e $countsFile ) {
    say STDERR "Need a counts file";
    exit(1);
}

my @conditions = split( ',', $condNames );
my $condcount = scalar @conditions;

my $i = 0;
my @scaledCounts;

my $cdataFH = IO::File->new( $countsFile, O_RDONLY ) || die( $countsFile . ": " . $! );

while ( my $cdataRow = $cdataFH->getline() ) {
    next if $cdataRow =~ /^\#|^_/;
    chomp($cdataRow);
    my @fields = ( split( '\t', $cdataRow ) )[ 1 .. $condcount ];

    next if sum(@fields) == 0;

    my $loggeomean = loggeomean( \@fields );
    next if ( $loggeomean == 0 );

    my @tmp = map { ( ( $_ == 0 ) ? 0 : log($_) ) - $loggeomean } @fields;
    my $j = 0;
    for ( 0 .. $condcount - 1 ) {
        $scaledCounts[$j][$i] = $tmp[$j];
        $j++;
    }

    $i++;
}
$cdataFH->close();

my $j = 0;
say STDOUT "Library size factors estimation for:";
for ( 0 .. $condcount - 1 ) {
    print $conditions[$j] . ": " . sprintf( "%01.4f", exp( median( $scaledCounts[$j] ) ) ) . "\n";
    $j++;
}

sub loggeomean {
    my ($data) = @_;

    my $loggeomean = 1 / ( scalar @$data ) * logsum(@$data);

    return ($loggeomean);
}

sub geomean {
    my ($data) = @_;

    my $product = 1;
    map { ( $_ > 0 ) ? $product *= $_ : $product *= 1 } @$data;
    my $geomean = $product**( 1 / ( scalar @$data ) );

    return ($geomean);
}

sub sum {
    my $sum;
    map { $sum += $_ } @_;
    return ($sum);
}

sub logsum {
    my $sum;
    map { $sum += ( $_ == 0 ) ? 0 : log($_) } @_;
    return ($sum);
}

sub median {
    my ($data) = @_;

    my $median;

    my $n   = scalar @$data;
    my $mid = int $n / 2;

    my @sorted_values = sort { $a <=> $b } @$data;
    if ( @$data % 2 ) {
        $median = $sorted_values[$mid];
    }
    else {
        $median = ( $sorted_values[ $mid - 1 ] + $sorted_values[$mid] ) / 2;
    }

    return ($median);
}

sub usage {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self [options] [-h] [-man] [--version]

Required options all of:
    -count-file|-cf        : read count file. Tabulator separated text file containing
                             the read counts for all conditions and replicates.
                             
                             Format:
                             featureID_1<tab>[counts<tab>...]
                             
                             You may create such a count file by using htseq-count on all
                             your samples and then merging the results:
                             
                             e.g.:
                             
                             paste CondA_rep1_Counts.txt CondA_rep2_Counts.txt CondA_rep3_Counts.txt \
                                   CondB_rep1_Counts.txt CondB_rep2_Counts.txt CondB_rep3_Counts.txt \
                                   | cut -f1-2,4,6,8,10,12 > allSamples_Counts.txt

    -condition-names|-cn   : Condition names as comma separated list.
                             
                             e.g.: -cn CondArep1,CondArep2,CondArep3,CondBrep1,CondBrep2,CondBrep3

                             (default: not set)

Options:

    --version              : Print the program version and exit.
    -h|help                : Print the program help information.
    -man                   : Print a detailed documentation.
    
    -debug|-d              : Print some debugging information.
    
EOF
}

__END__

=head1 NAME

estimateSizeFactors - estimate sequencing library size factors

=head1 SYNOPSIS

=head2  estimate the size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010)

=over 2
 
### Example with 3 replicates and 2 conditions

 estimateSizeFactors \
 -cf allSamples_Counts.txt \
 -cn CondArep1,CondArep2,CondArep3,CondBrep1,CondBrep2,CondBrep3


 The command estimate the size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010)
 
=back

=head1 DEPENDENCIES

 None

=head1 TESTED WITH:

=over

=item *
Perl 5.18.1 (Centos 6.5)

=item *
Perl 5.18.2 (Centos 6.5)

=item *
Perl 5.18.2 (RHEL 6.5)

=back 

=head1 OPTIONS

=head2 general:

=over 2

=item -version

 Print the program version and exit.

=item -h|-help

 Print the program help information with a detailed describtion of options

=item -man

 Print this documentation.

=back

=head1 LICENSE

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
  MA 02110-1301, USA.

  Copyright (C) 2015 D.Rieder

=head1 COPYRIGHT

  Copyright (C) 2015 Dietmar Rieder. All Rights Reserved.

=head1 AUTHOR

Dietmar Rieder

=head1 CONTACT

dietmar . rieder (at) i-med . ac . at

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.


exit;
