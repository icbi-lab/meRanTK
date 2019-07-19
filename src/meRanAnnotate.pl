#!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.24.0/bin/perl -w
#
#  meRanAnnotate.pl
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

use strict;
use warnings;
use Carp;

use Module::Load::Conditional qw(check_install);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;
use IO::File;
use File::Basename;

my $MCE;

BEGIN {
    if ( check_install( module => 'MCE::Loop' ) ) {
        $MCE = 1;
    }
    else {
        $MCE = 0;
    }
}

my $DEBUG = 0;
my $tabFile;
my $bedFile;
my $gffFile  = '';
my $outFile  = "";
my $isGTF;
my $isBED;
my $chrPrefixTab;
my $chrPrefixGFF;

my $VERSION = '1.2.1b';
my $version;
my $help;

my $expandResults;
my $maxWorkers = 1;

my $supportedFeaturesGFF = 'gene|mRNA|transcript|ncRNA|rRNA|tRNA|exon|intron|CDS';
my $supportedFeaturesGTF = 'gene|transcript|exon|intron|CDS|UTR';
my $featuresGFF = 'gene|mRNA|transcript|ncRNA';
my $featuresGTF = 'gene|transcript';
my $features;
my $reportDist;
my $relativeDist;

GetOptions(
            'tab|t=s'            => \$tabFile,
            'gff|g=s'            => \$gffFile,
            'feature|f=s'        => \$features,
            'expandResults|er'   => \$expandResults,
            'outfile|o=s'        => \$outFile,
            'parallel|p=i'       => \$maxWorkers,
            'BED|b=s'            => \$bedFile,
            'ensGTF|gtf'         => \$isGTF,
            'chrPrefix|cp=s'     => \$chrPrefixTab,
            'chrPrefixG|cpg=s'   => \$chrPrefixGFF,
            'reportDist|rd'      => \$reportDist,
            'relativeDist|reld'  => \$relativeDist,
            'help|h'             => \$help,
            'version'            => \$version,
            'debug|d'            => \$DEBUG,
            );

my %tabData;
my %gff;
my $header = "";

usage() and exit(0) if ($help);
say $VERSION and exit(0) if ($version);

my $debug;
if ($DEBUG) {
    eval "use Data::Dumper";
    eval "use diagnostics";
    $debug = sub { print STDERR "\n" . join( " ", @_ ) . "\n"; }
}
else {
    $debug = sub { return; };
}

if ( !$tabFile && !$bedFile || !$gffFile ) {
    usage() and exit(1);
}

if ( $maxWorkers > 1 ) {
    if ($MCE) {
        use MCE::Loop;
    }
    else {
        warn "MCE module not available, running in single core mode\n";
    }
}

# Force separate result lines for each overlapping feature if we need to report distances
if ($reportDist) { $expandResults = 1 };

my $supportedFeatures = ($isGTF) ? $supportedFeaturesGTF : $supportedFeaturesGFF;

# lets see if we use the default search features
$features = ($features) ? $features : ( ($isGTF) ? $featuresGTF : $featuresGFF );

my $searchFeatures;
if ( index( $features, '|' ) != -1 ) {
    my $i = 0;
    foreach my $f ( split( /\|/, $features ) ) {
        say STDERR $f . " is not supported" and next if $f !~ /\b$supportedFeatures\b/;
        $searchFeatures->[ $i++ ] = $f;
    }
    $features = join('|', @$searchFeatures);
}
else {
    say STDERR $features . " is not supported" if $features !~ /\b$supportedFeatures\b/;
    $searchFeatures->[0] = $features;
}
if (! defined($searchFeatures->[0])) {
    say STDERR "Not enough features to search for";
    exit(1);
}

# read the gff file (must be sorted)
readGFF( $features, $gffFile, \%gff, $chrPrefixGFF );

# read the tab/bed file
if ( $tabFile && (!$bedFile) ) {
    readTAB( $tabFile, \%tabData, $header, $chrPrefixTab );
}
elsif ( $bedFile && (!$tabFile) ) {
    readBED( $bedFile, \%tabData, $header, $chrPrefixTab );
}
else {
    say STDERR "Specifiy either a meRanCall/meRanCompare or a BED file";
    exit(1);
}

# search for intersecting features
my $resFH;
if ($outFile) {
    $resFH = IO::File->new( $outFile, O_RDWR | O_CREAT | O_TRUNC ) || die( $outFile . ": " . $! );
}
else {
    open( $resFH, ">-" );
}

if ($header) {
    $header .= "\t" . join( ';', ( sort @$searchFeatures ) );
    if ($reportDist) {
        $header .= "\t" . join("\t", ("5\' end", "3\' ends", "5\' to 3\'", "3\' to 5\'", "center to 5\'", "center to 3\'"));
    }
    $resFH->print($header . "\n");
}

my $report = sub { $resFH->print( $_[0] ) };

# check if we use multi processing
if ( $MCE && ( $maxWorkers > 1 ) ) {
    MCE::Loop::init {
                      chunk_size  => 1,
                      max_workers => $maxWorkers,
                      use_threads => 0,
                      tmp_dir     => '/dev/shm',
                      };

    $report = sub { MCE->print( $resFH, $_[0] ) };
    mce_loop { getIntersections($_) } ( keys(%tabData) );
}
else {
    $MCE = undef;
    foreach my $chr ( keys(%tabData) ) {
        getIntersections($chr);
    }
}

# Done lets close the result file handle
$resFH->close();
undef($resFH);


sub getIntersections {
    my $chr = $_[0];

    my $i = 0;

    if ( !exists( $gff{$chr} ) ) {
        say STDERR "Chromosome/SeqFeature not found in GFF3: " . $chr;
        return;
    }

    foreach my $range ( @{ $tabData{$chr}->{range} } ) {
        my $tabStrand = ( defined( $tabData{$chr}->{strand}->[$i] ) ) ? $tabData{$chr}->{strand}->[$i] : ".";

        my $overlappingFeaturesString;
        my %overlappingFeatures = ();
        my $featureDists;

        foreach my $sf ( keys %{ $gff{$chr}->{gffData} } ) {
            my $brs_iterator = binary_range_search(
                              { queries => [$range], database => $gff{$chr}->{gffData}->{$sf}, strand => $tabStrand } );

            while ( my $gff_line = $brs_iterator->() ) {
                my @keyval = split( /;\s?/, $gff_line->{attribute} );
                my %attribs;
                if ($isGTF) {
                    %attribs = map { my ($t, $v) = split( /\ /, $_, 2 ); ($v && ($t =~ /transcript_id|gene_id|gene_name/)) ? ($t => $v) : () } @keyval;
                }
                else {
                    %attribs = map { my ($t, $v) = split( /=/, $_, 2 ); ($v && ($t =~ /Name|gene/)) ? ($t => $v) : () } @keyval;
                }
                if (! %attribs) {
                    say STDERR "No tag/value pair found for: "
                        . $chr
                        . ":"
                        . $gff_line->{start}
                        . "-"
                        . $gff_line->{end}
                        . "... not a "
                        . ( ($isGTF) ? "GTF file? ... try to omit -gtf " : "GFF3 file? ... try with -gtf" );
                        next;
                }
                # &$debug(Dumper(\%attribs));
                if ( $gff_line->{feature} =~ /intron|exon/ ) { $attribs{Name} = $attribs{transcript_id}; }
                if ( ( $tabStrand ne "." ) && ( $gff_line->{strand} eq $tabStrand ) ) {
                    if ($isGTF) {
                        push(
                              @{ $overlappingFeatures{ $gff_line->{feature} } },
                              ( defined( $attribs{transcript_id} ) ? $attribs{transcript_id} : $attribs{gene_name} )
                              );
                    }
                    else {
                        push(
                              @{ $overlappingFeatures{ $gff_line->{feature} } },
                              ( defined( $attribs{Name} ) ? $attribs{Name} : $attribs{gene} )
                              );
                    }
                    if ($reportDist) {
                        my $dists = calcDist( {
                                                qRange  => $range,
                                                mRange  => [$gff_line->{start}, $gff_line->{end}],
                                                strand  => $tabStrand,
                                                relDist => $relativeDist
                                              } );
                        push( @{ $featureDists->{ $gff_line->{feature} } }, $dists );
                    }
                }
                elsif ( $tabStrand eq "." ) {
                    if ($isGTF) {
                        push(
                              @{ $overlappingFeatures{ $gff_line->{feature} } },
                              ( defined( $attribs{transcript_id} ) ? $attribs{transcript_id} : $attribs{gene_name} )
                              );
                    }
                    else {
                        push(
                              @{ $overlappingFeatures{ $gff_line->{feature} } },
                              ( defined( $attribs{Name} ) ? $attribs{Name} : $attribs{gene} )
                              );
                    }
                }
            }
        }    # end search feature

        if ( !%overlappingFeatures ) {
            if ($expandResults) {
                foreach my $f (@$searchFeatures) {
                    $overlappingFeatures{$f} = ['no_feature'];
                }
            }
            else {
                $overlappingFeatures{$features} = ['no_feature'];
            }
        }

        elsif ($expandResults) {
            foreach my $f (@$searchFeatures) {
                $overlappingFeatures{$f} = ['no_feature'] unless $overlappingFeatures{$f}[0];
            }
            
        }

        if ( ! $expandResults ) {
            foreach my $f ( sort ( keys %overlappingFeatures ) ) {
                $overlappingFeaturesString .= $f . ":" . join( ',', ( sort @{ $overlappingFeatures{$f} } ) ) . ';';
            }
            
            &$report( ( join( "\t", ( $tabData{$chr}->{line}->[$i], $overlappingFeaturesString ) ) . "\n" ) );
        }
        elsif ( $expandResults && ( ! $reportDist ) ) {
            foreach my $f ( sort ( keys %overlappingFeatures ) ) {
                foreach my $fm (sort @{ $overlappingFeatures{$f} }) {
                    &$report( ( join( "\t", ( $tabData{$chr}->{line}->[$i], $f . ":" . $fm ) ) . "\n" ) );
                }
            }
        }
        elsif ($reportDist) {
            foreach my $f ( sort ( keys %overlappingFeatures ) ) {
                my $c = 0;
                foreach my $fm (sort @{ $overlappingFeatures{$f} }) {
                    &$report( ( join( "\t", ( $tabData{$chr}->{line}->[$i], $f . ":" . $fm ) ) . "\t"
                             . ($featureDists->{$f}->[$c]->{q5m5} // ".") . "\t"
                             . ($featureDists->{$f}->[$c]->{q3m3} // ".") . "\t"
                             . ($featureDists->{$f}->[$c]->{q5m3} // ".") . "\t"
                             . ($featureDists->{$f}->[$c]->{q3m5} // ".") . "\t"
                             . ($featureDists->{$f}->[$c]->{qCm5} // ".") . "\t"
                             . ($featureDists->{$f}->[$c]->{qCm3} // ".")
                             . "\n" ) );
                    $c++;
                }
            }
        }
        
        $i++;
    }
}

sub readTAB {
    my $tabFile   = $_[0];
    my $tabData   = $_[1];
    my $hl        = $_[2];
    my $chrPrefix = $_[3] // "";
    
    my $haveHeader = undef;

    my $tabFH = IO::File->new( $tabFile, O_RDONLY | O_EXCL ) || die( $tabFile . ": " . $! );

    while ( defined( my $tabLine = $tabFH->getline() ) ) {
        chomp($tabLine);

        if ( $tabLine =~ /^#/ ) {
            $hl = $tabLine unless $haveHeader;
            $haveHeader = 1;
            next;
        }

        my ( $chrom, $chromPos, $strand, $fields ) = split( '\t', $tabLine, 4 );
        if ($chrPrefix) {
            $chrom = $chrPrefix . $chrom;
        }
        if ( $chromPos !~ /\d/ ) {
            say STDERR "This does not look like a meRanTK result file: " . $tabFile;
            exit(1);
        }
        if ( ($strand !~ /[+-\.]/) ) {
            say STDERR "This does not look like a BED file: strand can not be " . $strand . " :" . $tabFile;
            exit(1);
        }

        push( @{ $tabData->{$chrom}->{range} }, [ ( $chromPos - 1 ), $chromPos ] );
        push( @{ $tabData->{$chrom}->{strand} }, $strand );
        push( @{ $tabData->{$chrom}->{line} },   $tabLine );
    }
    $tabFH->close();
    undef($tabFH);
    $_[2] = $hl;
}

sub readBED {
    my $bedFile   = $_[0];
    my $bedData   = $_[1];
    my $hl        = $_[2];
    my $chrPrefix = $_[3] // "";
    
    my $haveHeader = undef;

    my $bedFH = IO::File->new( $bedFile, O_RDONLY | O_EXCL ) || die( $bedFile . ": " . $! );

    while ( defined( my $bedLine = $bedFH->getline() ) ) {
        chomp($bedLine);

        if ( $bedLine =~ /^#/ ) {
            $hl = $bedLine unless $haveHeader;
            $haveHeader = 1;
            next;
        }

        my ( $chrom, $chromStart, $chromEnd, $name, $score, $strand ) = split( '\t', $bedLine );
        if ($chrPrefix) {
            $chrom = $chrPrefix . $chrom;
        }
        if ( ($chromStart !~ /\d/) || ($chromEnd !~ /\d/) || ($chromEnd < $chromStart) ) {
            say STDERR "This does not look like a BED file: " . $bedFile;
            exit(1);
        }
        if ( ($strand !~ /[+-\.]/) ) {
            say STDERR "This does not look like a BED file: strand can not be " . $strand . " :" . $bedFile;
            exit(1);
        }

        push( @{ $bedData->{$chrom}->{range} }, [ ( $chromStart ), $chromEnd ] );
        push( @{ $bedData->{$chrom}->{strand} }, $strand );
        push( @{ $bedData->{$chrom}->{line} },   $bedLine );
    }
    $bedFH->close();
    undef($bedFH);
    $_[2] = $hl;
}

sub readGFF {
    my $reqFeature = shift;
    my $gffFile    = shift;
    my $gff        = shift;
    my $chrPrefix  = shift // "";

    my $gffFH = IO::File->new( $gffFile, O_RDONLY | O_EXCL ) || die( $gffFile . ": " . $! );

    my $i           = 0;
    my $old_seqname = "";
    my $lastEnd     = 0;
    my $old_gffLine = "";

    my $needToSort;

    my %featureCounter = map { $_ => 0 } @$searchFeatures;

    while ( defined( my $gffLine = $gffFH->getline() ) ) {
        next if ( $gffLine =~ /^\#/ );

        chomp($gffLine);
        my ( $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute ) = split( '\t', $gffLine );

        if ( $old_seqname ne $seqname ) {
            $i              = 0;
            $lastEnd        = 0;
            %featureCounter = map { $_ => 0 } @$searchFeatures;
        }
        $old_seqname = $seqname;

        if ( $feature =~ m/\b($reqFeature)\b/ ) {
            if ( ( $end < $lastEnd ) && ( $lastEnd != 0 ) ) {
                #
                # Have now internal sorting function :-)
                #
                #say STDERR
                #    "\nGFF file is not sorted, please use sort to sort your gff file, e.g:\n\n\tsort -S 1G -t \$\'\\t\' -k1,1 -k5,5n "
                #  . $gffFile
                #  . " > mysortedGFF.gff\n\n";
                #say STDERR $old_gffLine;
                #say STDERR $gffLine;
                #exit(1);
                
                # Add to sort jobs
                $needToSort->{$seqname}->{$feature} = 1;
            }
            $lastEnd = $end;

            if ($chrPrefix) {
                $seqname = $chrPrefix . $seqname;
            }

            $gff->{$seqname}->{gffData}->{$feature}->[ $featureCounter{$feature} ] = {
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

            $featureCounter{$feature}++;
            $old_gffLine = $gffLine;
        }

    }
    $gffFH->close();
    undef($gffFH);

    if ($needToSort) {
        foreach my $s (keys %$needToSort) {
            foreach my $f (keys %{$needToSort->{$s}}) {
                &$debug("Sorting:", $s, "-", $f, "-", $0, "Line:", __LINE__);
                sortGFF($gff->{$s}->{gffData}->{$f});
            }
        }
    }

}

sub sortGFF {
    my $toSort = $_[0];

    @$toSort = sort { $a->{end} <=> $b->{end} } @$toSort;
}

sub calcDist {
    my $options = shift;
    
    my $qRange  = $options->{qRange};
    my $mRange  = $options->{mRange};
    my $strand  = $options->{strand};
    my $relDist = $options->{relDist};

    my $dists;
    my ($fivePE, $threePE) = (0, 1);
    
    my $lf = 1;
    if ($relDist) {
        my $mLength = $mRange->[$threePE] - $mRange->[$fivePE] + 1;
        $lf = 1 / $mLength * 100;
    }
        
    if ($strand eq '-') {
        ($fivePE, $threePE) = (1, 0);
    }
    
    $dists->{q5m5} = abs( $qRange->[$fivePE]  - $mRange->[$fivePE]  ) * $lf;
    $dists->{q3m3} = abs( $mRange->[$threePE] - $qRange->[$threePE] ) * $lf;
    
    $dists->{q5m3} = abs( $mRange->[$threePE] - $qRange->[$fivePE] ) * $lf;
    $dists->{q3m5} = abs( $qRange->[$threePE] - $mRange->[$fivePE] ) * $lf;
    
    my $qC = int( ($qRange->[$fivePE] + $qRange->[$threePE]) / 2 );

    $dists->{qCm5} =  abs( $qC - $mRange->[$fivePE] ) * $lf;
    $dists->{qCm3} =  abs( $mRange->[$threePE] - $qC ) * $lf;

 
    while ( my($d, $v) = each %$dists) {
        if($relDist) {
            $v = sprintf("%01.3f" , $v);
            if ( ($v > 100) || ($v < 0) ) {
                say STDERR "Invalid relative distance calculated: " . $d . " -> " . $v;
                &$debug( "Got these data:", Dumper($options), "Line:", __LINE__ );
                $dists->{$d} = undef;
            }
            else {
                $dists->{$d} = $v;
            }
        }
        else {
            if ( $v < 0 ) {
                say STDERR "Invalid distance calculated: " . $d . " -> " . $v;
                &$debug( "Got these data:", Dumper($options), "Line:", __LINE__ );
                $dists->{$d} = undef;
            }
        }
    }

    return($dists);
}

sub binary_range_search {
    my $options = shift;

    my $queries  = $options->{queries}  || croak 'Need a queries parameter';
    my $database = $options->{database} || croak 'Need a database parameter';
    my $strand   = $options->{strand}   || croak 'Need a strand parameter';

    my ( $low, $high ) = ( 0, $#{$database} );
    my @iterators = ();

    my $maxtry = $high;

  TARGET:
    for my $query (@$queries) {
        my $end = 0;

      RANGE_CHECK:
        while ( $low <= $high ) {

            # middle2
            my $try = int( ( $low + $high ) / 2 );

            if ( range_before( $database->[$try]{range}, $query ) ) {
                $low = $try + 1;
                if ( $low == $high - 1 ) {
                    $low = $try;
                }
                if ( $high == $low ) {

                    # go linearely up
                    while ( $try < $maxtry ) {
                        $try++;
                        if ( range_overlap( $query, $database->[$try]{range} )
                             && ( $strand eq $database->[$try]{strand} ) )
                        {
                            goto RANGE_FOUND;
                        }
                    }
                }

                next RANGE_CHECK;
            }
            elsif ( range_before( $query, $database->[$try]{range} ) ) {
                $high = $try - 1;
                if ( $high == $low + 1 ) {
                    $high = $try;
                }
                if ( $high == $low ) {

                    # go linearely down
                    while ( $try < $maxtry ) {
                        $try++;
                        if ( range_overlap( $query, $database->[$try]{range} )
                             && ( $strand eq $database->[$try]{strand} ) )
                        {
                            goto RANGE_FOUND;
                        }
                    }
                }
                next RANGE_CHECK;
            }
            elsif ( $strand ne $database->[$try]{strand} ) {

                # go linearely down
                my $downtry = $try;
                my $uptry   = $try;
                while ( $downtry > 0 && $uptry < $maxtry ) {
                    $uptry++;
                    if ( range_overlap( $query, $database->[$uptry]{range} )
                         && ( $strand eq $database->[$uptry]{strand} ) )
                    {
                        $try = $uptry;
                        goto RANGE_FOUND;
                    }
                    $downtry--;
                    if ( range_overlap( $query, $database->[$downtry]{range} )
                         && ( $strand eq $database->[$downtry]{strand} ) )
                    {
                        $try = $downtry;
                        goto RANGE_FOUND;
                    }
                }

            }

          RANGE_FOUND:
            # if we make it here, we've found a source $try which overlaps with the
            # target $query

            my ( $down, $up ) = ($try) x 2;
            my %seen = ();

            # create an iterator which, on every call, returns an overlapping
            # source range, starting with $query and then going down or up the
            # list.
            my $brs_iterator = sub {

                # if $query overlaps with $up + 1, and it's new, return it
                if (     $up + 1 <= $maxtry
                     and !exists $seen{ $up + 1 }
                     and range_overlap( $database->[ $up + 1 ]{range}, $query )
                     && ( $strand eq $database->[ $up + 1 ]{strand} ) )
                {
                    $seen{ $up + 1 } = undef;
                    return $database->[ ++$up ];
                }

                # if $query overlaps with $down - 1, and it's new, return it
                elsif (     $down - 1 >= 0
                        and !exists $seen{ $down - 1 }
                        and range_overlap( $database->[ $down - 1 ]{range}, $query )
                        && ( $strand eq $database->[ $down - 1 ]{strand} ) )
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

meRanAnnotate is a tool to annotate genomic m5C positions from meRanCall/meRanCompare results.
It can use either ensembl GTF or NCBI GFF3 files to annotate m5Cs using selected features
like 'mRNA', 'gene', 'ncRNA' and so on.
Besides meRanTK result files it can also take BED files and intersect the ranges with selected
GFF3 features.

USAGE: $self [options] ][-h] [--version]

Required options any of:
    -tab|-t               : meRanCall/meRanCompare result file to intersect with gff.
                            Format:
                            <chr><tab><position><tab><strand>[<tab><field>]...
     OR

    -bed|-b               : BED file to intersect with gff.
                            Format:
                            <chr><tab><start><tab><end><tab><name><tab><score><tab><strand>[<tab><field>...]

    -gff|-g               : Sorted or unsorted GFF3/GTF file.

Options:

    -ensGTF|-gtf          : Annotation file is a GTF file
                           (default: not set, assuming it is GFF3)

    -feautre|-f           : GFF3 features you want to intersect with your meRanTK or BED file,
                            e.g. if you are interested in mRNA and ncRNA use

                                 -f 'mRNA|ncRNA'

                            (default: 'gene|mRNA|transcript|ncRNA' for NCBI GFF3
                                      'gene|transcript' for Ensembl GTF)

    -outfile|-o           : Result file where to store the intersecting/overlapping
                            features.
                            (default: STDOUT)

    -parallel|-p          : Number of CPUs to use. Running in parallel mode is way faster but
                            it requires the MCE Perl module to be installed.
                            (default: 1, no parallel processing)

    -chrPrefix|-cp        : Prepend this string to chromosome/sequence name from the meRanCall/meRanCompare
                            result of BED file. The chromosome/sequence name has to
                            match the one used in the GFF3/GTF file.

                            e.g. 
                            If you have m5C calls or BED ranges that have chromosome names of the type 
                            '1, 2, 3, ...' and want to use a GFF3/GTF file that has chromosome names
                            of the type 'chr1, chr2, chr3, ...' then you might use -cp 'chr'

                           (default: not set, assuming chromosome names are matching)

    -reportDist|-rd       : Calculate distances to features 5' and 3' ends.
                            The following distances are calculated:
                            
                                  query 5' to feature 5'
                                  query 3' to feature 3'
                                  query 5' to feature 3'
                                  query 3' to feature 5'
                                  query center to feature 5'
                                  query center to feature 3'
                                  
                            All distances are reported in stranded mode.

                           (default: not set, not distances are reported)

    -relativeDist|-reld   : Calculate relative distances. Distances will be reported as
                            %-feature length, that is genomic feature 3' - genomic feature 5' end coordinate.

                           (default: not set, distances are reported in # of bp)

    -expandResults|-er    : Expand results. Report each GFF-feature match as a separate line.

                           (default: not set, forced if -rd is set)

    --version             : Print the program version and exit.

    -h|help               : Print the program help information.


Acknowledgements:
    Thanks go to Pedro Silva <psilva (at) pedrosilva (dot) pt> for his binary range search algorithm,
    which was a useful inspiration for this program.

EOF
}
