#!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.18.2/bin/perl -w
#
#  bedISgff.pl
#
#  Copyright 2015 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
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

use strict;
use warnings;
use Carp;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;
use IO::File;
use File::Basename;

my $bedFile;
my $gffFile  = '';
my $features = 'gene|mRNA|transcript|ncRNA';
my $outFile  = "";

my $VERSION = '0.6';
my $version;
my $help;

GetOptions(
            'bed|b=s'     => \$bedFile,
            'gff|g=s'     => \$gffFile,
            'feature|f=s' => \$features,
            'outfile|o=s' => \$outFile,
            'help|h'      => \$help,
            'version'     => \$version,
            );

my %bedData;
my %gff;

usage() and exit(0) if ($help);
say $VERSION and exit(0) if ($version);

if ( !$bedFile || !$gffFile ) {
    usage() and exit(1);
}

# read the gff file (must be sorted)
readGFF( $features, $gffFile, \%gff );

# read the bed file
readBED( $bedFile, \%bedData );

# search for intersecting features
my $resFH;
if ($outFile) {
    $resFH = IO::File->new( $outFile, O_RDWR | O_CREAT | O_TRUNC ) || die( $outFile . ": " . $! );
}
else {
    open( $resFH, ">-" );
}
foreach my $chr ( keys(%bedData) ) {

    my $i = 0;

    if ( !exists( $gff{$chr} ) ) {
        warn "Chromosome/SeqFeature not found in GFF: " . $chr;
        next;
    }

    foreach my $range ( @{ $bedData{$chr}->{range} } ) {
        my $brs_iterator = binary_range_search( { queries => [$range], database => $gff{$chr}->{gffData} } );

        my $bedStrand = ( split( '\t', $bedData{$chr}->{line}->[$i] ) )[5];
        $bedStrand = ( defined($bedStrand) ) ? $bedStrand : ".";
        my $matchLinesOut = [ [ $bedData{$chr}->{line}->[$i], "genomic region", ".", ".", ".", $bedStrand ] ];

        my $matches = 0;
        while ( my $gff_line = $brs_iterator->() ) {

            my @keyval = split( /;/, $gff_line->{attribute} );
            my %attribs = map { split( /=/, $_ ) } @keyval;

            if ( ( $bedStrand ne "." ) && ( $gff_line->{strand} eq $bedStrand ) ) {
                if ( $gff_line->{feature} eq 'intron' ) { $attribs{Name} = $attribs{Parent}; }
                $matchLinesOut->[$matches] = [
                    $bedData{$chr}->{line}->[$i], ( defined( $attribs{Name} ) ? $attribs{Name} : $attribs{gene} ),
                    $attribs{gene},   $gff_line->{start},
                    $gff_line->{end}, $gff_line->{strand}
                    ];
                $matches++;
            }
            elsif ( $bedStrand eq "." ) {
                if ( $gff_line->{feature} eq 'intron' ) { $attribs{Name} = $attribs{Parent}; }
                $matchLinesOut->[$matches] = [
                    $bedData{$chr}->{line}->[$i], ( defined( $attribs{Name} ) ? $attribs{Name} : $attribs{gene} ),
                    $attribs{gene},   $gff_line->{start},
                    $gff_line->{end}, $gff_line->{strand}
                    ];
                $matches++;
            }
        }
        foreach my $matchLine (@$matchLinesOut) {
            $resFH->print( ( join( "\t", @$matchLine ) . "\n" ) );
        }
        $matchLinesOut = undef;
        $i++;
    }
}
$resFH->close();
undef($resFH);

sub readBED {
    my $bedFile = shift;
    my $bedData = shift;

    my $bedFH = IO::File->new( $bedFile, O_RDONLY | O_EXCL ) || die( $bedFile . ": " . $! );

    while ( defined( my $bedLine = $bedFH->getline() ) ) {
        next if ( $bedLine =~ /^#/ );

        chomp($bedLine);
        my ( $chrom, $chromStart, $chromEnd ) = split( '\t', $bedLine, 4 );
        push( @{ $bedData->{$chrom}->{range} }, [ $chromStart, $chromEnd ] );
        push( @{ $bedData->{$chrom}->{line} }, $bedLine );
    }
    $bedFH->close();
    undef($bedFH);

}

sub readGFF {
    my $reqFeature = shift;
    my $gffFile    = shift;
    my $gff        = shift;

    my $gffFH = IO::File->new( $gffFile, O_RDONLY | O_EXCL ) || die( $gffFile . ": " . $! );

    my $i           = 0;
    my $old_seqname = "";
    my $lastStart   = 0;
    my $old_gffLine = "";

    while ( defined( my $gffLine = $gffFH->getline() ) ) {
        next if ( $gffLine =~ /^\#/ );

        chomp($gffLine);
        my ( $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute ) = split( '\t', $gffLine );

        if ( $old_seqname ne $seqname ) {
            $i         = 0;
            $lastStart = 0;
        }
        $old_seqname = $seqname;

        if ( $feature =~ m/$reqFeature/ ) {
            if ( ( $start < $lastStart ) && ( $lastStart != 0 ) ) {
                say STDERR
"\nGFF file is not sorted, please use sort to sort your gff file, e.g:\n\n\tsort -S 1G -t \$\'\\t\' -k1,1 -k4,4n "
                  . $gffFile
                  . " > mysortedGFF.gff\n\n";
                say STDERR $old_gffLine;
                say STDERR $gffLine;
                exit(1);
            }
            $lastStart = $start;

            $gff->{$seqname}->{gffData}->[$i] = {
                                                  seqname   => $seqname,
                                                  source    => $source,
                                                  feature   => $feature,
                                                  start     => $start,
                                                  end       => $end,
                                                  range     => [ $start, $end ],
                                                  score     => $score,
                                                  strand    => $strand,
                                                  frame     => $frame,
                                                  attribute => $attribute,
                                                  };
            $i++;
            $old_gffLine = $gffLine;
        }

    }
    $gffFH->close();
    undef($gffFH);

}

sub binary_range_search {
    my $options = shift;

    my $queries  = $options->{queries}  || croak 'Need a queries parameter';
    my $database = $options->{database} || croak 'Need a database parameter';

    my ( $low, $high ) = ( 0, $#{$database} );
    my @iterators = ();

  TARGET:
    for my $query (@$queries) {

      RANGE_CHECK:
        while ( $low <= $high ) {

            # middle
            my $try = int( ( $low + $high ) / 2 );

            if ( range_before( $database->[$try]{range}, $query ) ) {
                $low = $try + 1;
                next RANGE_CHECK;
            }
            elsif ( range_before( $query, $database->[$try]{range} ) ) {
                $high = $try - 1;
                next RANGE_CHECK;
            }

            # if we make it here, we've found a source $try which overlaps with the
            # target $query

            my ( $down, $up ) = ($try) x 2;
            my %seen = ();

            # create an iterator which, on every call, returns an overlapping
            # source range, starting with $query and then going down or up the
            # list.
            my $brs_iterator = sub {

                # if $query overlaps with $up + 1, and it's new, return it
                if (     $up + 1 <= $high
                     and !exists $seen{ $up + 1 }
                     and range_overlap( $database->[ $up + 1 ]{range}, $query ) )
                {
                    $seen{ $up + 1 } = undef;
                    return $database->[ ++$up ];
                }

                # if $query overlaps with $down - 1, and it's new, return it
                elsif (     $down - 1 >= $low
                        and !exists $seen{ $down - 1 }
                        and range_overlap( $database->[ $down - 1 ]{range}, $query ) )
                {
                    $seen{ $down - 1 } = undef;
                    return $database->[ --$down ];
                }

                # we already know $try overlaps, so return it too.
                elsif ( !exists $seen{$try} ) {
                    $seen{$try} = undef;
                    return $database->[$try];
                }
                else {
                    return;
                }
            };
            push @iterators, $brs_iterator;
            next TARGET;
        }
    }

    # In scalar context return master iterator that iterates over the list of range iterators.
    # In list context returns a list of range iterators.
    return wantarray
      ? @iterators
      : sub {
        while (@iterators) {
            if ( my $range = $iterators[0]->() ) {
                return $range;
            }
            shift @iterators;
        }
        return;
      };
}

sub range_before {
    my ( $range_a, $range_b ) = @_;

    return $range_a->[1] < $range_b->[0];
}

sub range_after {
    my ( $range_a, $range_b ) = @_;

    return $range_a->[0] > $range_b->[1];
}

sub range_overlap {
    my ( $range_a, $range_b ) = @_;

    return !( range_before( $range_a, $range_b ) || range_after( $range_a, $range_b ) );
}

sub usage {
    my $self = basename($0);

    print STDERR<<EOF
bedISGgff is a tool to search BED regions overlapping with selected GFF
features. It uses a fast binary search algorithm initally developed by
Pedro Silva <psilva (at) pedrosilva (dot) pt>.

USAGE: $self [options] ][-h] [--version]

Required options any of:
    -bed|-b               : BED file to intersect with gff.

    -gff|-g               : Sorted GFF3 file.
                            You can sort your GFF file using the unix "sort"
                            command,. e.g.:

                            sort -S 1G -t \$'\\t' -k1,1 -k4,4n gffFile.gff > mysortedGFF.gff

Options:

    -feautre|-f           : GFF features you want to intersect with your BED file,
                            e.g. if you are interested in mRNA and ncRNA use

                                 -f 'mRNA|ncRNA'

                            (default: gene|mRNA|transcript|ncRNA)

    -outfile|-o           : Result file where to store the intersecting/overlapping
                            features.
                            (default: STDOUT)

    --version             : Print the program version and exit.
    
    -h|help               : Print the program help information.
    
EOF
}
