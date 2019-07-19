#!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.24.0/bin/perl -w
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

use List::Util qw( min );
use Math::CDF qw( pnorm qnorm pchisq);
use Text::NSP::Measures::2D::Fisher::twotailed;

use constant { LOG10 => log(10),
               LOG2  => log(2),
              };

my $aFilesCL     = "";
my $bFilesCL     = "";
my $minRep       = 2;
my $sig          = 0.01;
my $FDR          = 0.05;

my $minFC;
my $log2_minFC   = -1;

my $condNameA    = "conditionA";
my $condNameB    = "conditionB";
my $sizeFactorLA = "";
my $sizeFactorLB = "";
my $azaMode;

my $help    = undef;
my $man     = undef;
my $version = undef;
my $VERSION = "1.2.1b";

my $DEBUG = 0;
my $debug;

GetOptions(
            'condition-files-a|fa=s' => \$aFilesCL,
            'condition-files-b|fb=s' => \$bFilesCL,
            'condition-name-a|na=s'  => \$condNameA,
            'condition-name-b|nb=s'  => \$condNameB,
            'size-factors-a|sfa=s'   => \$sizeFactorLA,
            'size-factors-b|sfb=s'   => \$sizeFactorLB,
            'minRep|mr=i'            => \$minRep,
            'minFC|fc=f'             => \$minFC,
            'sig|s=f'                => \$sig,
            'fdr=f'                  => \$FDR,
            'azaMode|aza'            => \$azaMode,
            'debug|d'                => \$DEBUG,
            'help|h'                 => \$help,
            'man|m'                  => \$man,
            'version'                => \$version,
            );

if ($DEBUG) {
    eval "use Data::Dumper";
    $debug = sub { print STDERR "\n" . join( " ", @_ ) . "\n"; }
}
else {
    $debug = sub { return; };
}

if ($minFC) {
    $log2_minFC = log2($minFC);
}

usage() and exit(0) if ($help);
pod2usage( -verbose => 3 ) if $man;

say $VERSION and exit(0) if ($version);

my $resultLine;
my %m5Cdata;
my $m5C   = \%m5Cdata;
my $i     = 0;
my $resNr = 0;
my @fdrData;

my ( @aFiles, @bFiles );
if ( $aFilesCL && $bFilesCL ) {
    @aFiles = split( ',', $aFilesCL );
    @bFiles = split( ',', $bFilesCL );
}
else {
    say STDERR "Need at least one bs call file from condition A and condition B";
    exit(1);
}

my @sizeFactorsA;
my @sizeFactorsB;
if ( $sizeFactorLA && $sizeFactorLB ) {
    @sizeFactorsA = split( ',', $sizeFactorLA );
    @sizeFactorsB = split( ',', $sizeFactorLB );
}
else {
    map { $sizeFactorsA[$_] = 1; } 0 .. ( @aFiles - 1 );
    map { $sizeFactorsB[$_] = 1; } 0 .. ( @bFiles - 1 );
}

if ( $azaMode && ( ! ($sizeFactorLA && $sizeFactorLB) ) ) {
    say STDERR "Warning: no library size factors specified - library size normalization will be ommitted and enrichment (fold change) over IP may be incorrect.";
}

# Read BS data for condition A
my $fileNrA = 0;
for my $bsDataF (@aFiles) {
    my $fh = IO::File->new( $bsDataF, O_RDONLY ) || die( $bsDataF . ": " . $! );
    readBScallData( $fh, 'a', $sizeFactorsA[$fileNrA], $m5C );
    $fh->close();
    $fileNrA++;
}

# Read BS data for condition B
my $fileNrB = 0;
for my $bsDataF (@bFiles) {
    my $fh = IO::File->new( $bsDataF, O_RDONLY ) || die( $bsDataF . ": " . $! );
    readBScallData( $fh, 'b', $sizeFactorsB[$fileNrB], $m5C );
    $fh->close();
    $fileNrB++;
}

# use Fisher's exact test if no replicates
my $useFET = 0;
if ( ($fileNrA + $fileNrB) == 2 ) {
    $useFET = 1;
    $minRep = 1;
}

# Define output file names
my $condAuniq_bed = "uniq_" . $condNameA . ".bed";
my $condAuniq_res = "uniq_" . $condNameA . ".txt";
my $condBuniq_bed = "uniq_" . $condNameB . ".bed";
my $condBuniq_res = "uniq_" . $condNameB . ".txt";

my $condABintersect_bed = "intersect_" . $condNameA . "_" . $condNameB . ".bed";
my $condABintersect_res = "intersect_" . $condNameA . "_" . $condNameB . ".txt";

# open file handles
my $uAbedFH = IO::File->new( $condAuniq_bed, O_RDWR | O_CREAT | O_TRUNC ) || die( $condAuniq_bed . ": " . $! );
my $uAresFH = IO::File->new( $condAuniq_res, O_RDWR | O_CREAT | O_TRUNC ) || die( $condAuniq_res . ": " . $! );
my $uBbedFH = IO::File->new( $condBuniq_bed, O_RDWR | O_CREAT | O_TRUNC ) || die( $condBuniq_bed . ": " . $! );
my $uBresFH = IO::File->new( $condBuniq_res, O_RDWR | O_CREAT | O_TRUNC ) || die( $condBuniq_res . ": " . $! );

my $iABbedFH = IO::File->new( $condABintersect_bed, O_RDWR | O_CREAT | O_TRUNC )
  || die( $condABintersect_bed . ": " . $! );
my $iABresFH = IO::File->new( $condABintersect_res, O_RDWR | O_CREAT | O_TRUNC )
  || die( $condABintersect_res . ": " . $! );

my $uniqHL = "";
my $intersectHL = "";
if ($useFET) {
    $uniqHL = join("\t", ('#SeqID',
                          'refPos',   'refStrand',  'Name',            'cov',
                          'C_count',  'norm(cov)',  'norm(C_count)',   'methRate',
                          'p-value_mState',         'p-value_mRate',   'rep' . "\n") );

    $intersectHL = join( "\t", ('#SeqID',
                                'refPos',
                                'refStrand',
                                'Name',
                                'p-value',
                                'cov A',
                                'C_count A',
                                'norm(cov) A',
                                'norm(C_count)',
                                'methRate A',
                                'p-value_mState A',
                                'p-value_mRate A',
                                'rep A',
                                'cov B',
                                'C_count B',
                                'norm(cov)',
                                'norm(C_count)',
                                'methRate B',
                                'p-value_mState B',
                                'p-value_mRate B',
                                'rep B',
                                ($azaMode) ? ('log2(IP/CTRL)' . "\n") : ('log2(mRate-FC)' . "\n") ) );
}
else {
    $uniqHL = join("\t", ('#SeqID',
                          'refPos',       'refStrand',      'Name',                'avg/cov',
                          'avg/C_count',  'norm(avg/cov)',  'norm(avg/C_count)',   'avg/methRate',
                          'sem/methRate', 'stdev/methRate', 'comb/p-value_mState', 'comb/p-value_mRate',
                          'rep' . "\n") );

    $intersectHL = join( "\t", ('#SeqID',
                                'refPos',
                                'refStrand',
                                'Name',
                                'p-value',
                                'avg/cov A',
                                'avg/C_count A',
                                'norm(avg/cov) A',
                                'norm(avg/C_count) A',
                                'avg/methRate A',
                                'sem/methRate A',
                                'stdev/methRate A',
                                'comb/p-value_mState A',
                                'comb/p-value_mRate A',
                                'rep A',
                                'avg/cov B',
                                'avg/C_count B',
                                'norm(avg/cov) B',
                                'norm(avg/C_count) B',
                                'avg/methRate B',
                                'sem/methRate B',
                                'stdev/methRate B',
                                'comb/p-value_mState B',
                                'comb/p-value_mRate B',
                                'rep B',
                                ($azaMode) ? ('log2(IP/CTRL)' . "\n") : ('log2(mRate-FC)' . "\n") ) );
}
$uAresFH->print($uniqHL);
$uBresFH->print($uniqHL);
$iABresFH->print($intersectHL);

while ( my ( $key, $value ) = each %$m5C ) {

    #    last if $i == 2000;
    my @t;
    my $r = 0;

    my @m5Ccoords = split( '###', $key );
    my $m5Cname = 'm5C_' . $i;

    if ( ( !defined( $m5C->{$key}->{b}->{r} ) ) && ( $m5C->{$key}->{a}->{r} >= $minRep ) ) {
        my @stats = ($useFET) ? getStatsNoRep( $m5C->{$key}->{a} ) : getStats( $m5C->{$key}->{a} );
        $uAbedFH->print(
              join( "\t", ( $m5Ccoords[0], $m5Ccoords[1] - 1, $m5Ccoords[1], $m5Cname, 1000, $m5Ccoords[2] ) ) . "\n" );
        $uAresFH->print( join( "\t", ( $m5Ccoords[0], $m5Ccoords[1], $m5Ccoords[2], $m5Cname, @stats ) ) . "\n" );
        $i++;
        next;
    }
    elsif ( ( !defined( $m5C->{$key}->{a}->{r} ) ) && ( $m5C->{$key}->{b}->{r} >= $minRep ) ) {
        my @stats = ($useFET) ? getStatsNoRep( $m5C->{$key}->{b} ) : getStats( $m5C->{$key}->{b} );
        $uBbedFH->print(
              join( "\t", ( $m5Ccoords[0], $m5Ccoords[1] - 1, $m5Ccoords[1], $m5Cname, 1000, $m5Ccoords[2] ) ) . "\n" );
        $uBresFH->print( join( "\t", ( $m5Ccoords[0], $m5Ccoords[1], $m5Ccoords[2], $m5Cname, @stats ) ) . "\n" );
        $i++;
        next;
    }
    elsif ( ( !defined( $m5C->{$key}->{b}->{r} ) ) && ( $m5C->{$key}->{a}->{r} < $minRep ) ) {
        next;
    }
    elsif ( ( !defined( $m5C->{$key}->{a}->{r} ) ) && ( $m5C->{$key}->{b}->{r} < $minRep ) ) {
        next;
    }
    elsif ( ( $m5C->{$key}->{a}->{r} < $minRep ) && ( $m5C->{$key}->{b}->{r} < $minRep ) ) {
        next;
    }

    #
    #         m       n
    #   A   [0][0]  [0][1]
    #   B   [1][0]  [1][1]

    my $aReps  = $m5C->{$key}->{a}->{r};
    my $bReps  = $m5C->{$key}->{b}->{r};
    my $maxRep = ( ( $bReps > $aReps ) ? $bReps : $aReps ) - 1;
    for ( 0 .. $maxRep ) {
        my @st = undef;
        $st[0][0] = $m5C->{$key}->{a}->{m}->[$r] // 0;
        $st[0][1] = $m5C->{$key}->{a}->{n}->[$r] // 0;

        $st[1][0] = $m5C->{$key}->{b}->{m}->[$r] // 0;
        $st[1][1] = $m5C->{$key}->{b}->{n}->[$r] // 0;

        # avoid emtpy rep
        next if ( sum( $st[0][0], $st[0][1], $st[1][0], $st[1][1] ) == 0 );
        $t[$r] = [@st];
        $r++;
    }

    my $pv = 1;
    if ($useFET) {
        $pv = fisherExactTest(@t);
    }
    else {
        $pv = cmhTest(@t);
    }

    my $score = ( $pv != 0 ) ? 100 * -log10($pv) : 1000;
    $score = ( $score > 1000 ) ? 1000 : int($score);

    if ( $pv < $sig ) {
        my @statsA = ($useFET) ? getStatsNoRep( $m5C->{$key}->{a} ) : getStats( $m5C->{$key}->{a} );
        my @statsB = ($useFET) ? getStatsNoRep( $m5C->{$key}->{b} ) : getStats( $m5C->{$key}->{b} );
        
        my $mRateFC = ($azaMode)
                        ? ( ($sizeFactorLA && $sizeFactorLB)
                              ? ( log2($statsA[2] / $statsB[2]) )
                              : ( log2($statsA[0] / $statsB[0]) ) ) 
                        : ( log2($statsA[4] / $statsB[4]) );
        

        if ( abs($mRateFC) >= $log2_minFC ) {
            $iABbedFH->print(
                join( "\t",
                       ( $m5Ccoords[0], $m5Ccoords[1] - 1, $m5Ccoords[1], $m5Cname, $score, $m5Ccoords[2] )
                   )
                   . "\n" 
                   );
            
            $resultLine = join( "\t",
                                 ( $m5Ccoords[0], $m5Ccoords[1], $m5Ccoords[2], $m5Cname, sprintf( "%01.6e", $pv ), @statsA, @statsB, sprintf( "%01.3f", $mRateFC) )
                              )
                              . "\n";
                     
            $iABresFH->print( $resultLine );
            $fdrData[$resNr] = [ $pv, $resultLine ];
            $resNr++;
        }
    }
    $i++;
}
$uAbedFH->close();
$uAresFH->close();
$uBbedFH->close();
$uBresFH->close();
$iABbedFH->close();
$iABresFH->close();

if ($FDR) {
    my $condABintersect_res_fdr = "intersect_" . $condNameA . "_" . $condNameB . "_FDR_" . $FDR . ".txt";
    my $iABresFDRFH = IO::File->new( $condABintersect_res_fdr, O_RDWR | O_CREAT | O_TRUNC )
      || die( $condABintersect_res_fdr . ": " . $! );

    my $condABintersect_bed_fdr = "intersect_" . $condNameA . "_" . $condNameB . "_FDR_" . $FDR . ".bed";
    my $iABbedFDRFH = IO::File->new( $condABintersect_bed_fdr, O_RDWR | O_CREAT | O_TRUNC )
      || die( $condABintersect_bed_fdr . ": " . $! );

    my @HF = split(/\t/, $intersectHL);
    splice( @HF, 5, 0, "adjusted p-val" );
    $iABresFDRFH->print( join("\t", @HF) );

    my $fdrResults = FDRc( \@fdrData, $FDR );

    if ( defined($fdrResults) ) {
        foreach my $record (@$fdrResults) {

            my @resFields = split( '\t', $record->[2] );
            my $pv_adj_nice = sprintf( "%01.6e", $record->[1]);
            splice( @resFields, 5, 0, $pv_adj_nice );

            my $score = ( $record->[0] != 0 ) ? 100 * -log10($record->[0]) : 1000;
            $score = ( $score > 1000 ) ? 1000 : int($score);

            my $resLine = join( "\t", @resFields );
            my $bedLine = join( "\t", $resFields[0], ($resFields[1] - 1), $resFields[1], $resFields[3], $score, $resFields[2]) . "\n";

            $iABresFDRFH->print($resLine);
            $iABbedFDRFH->print($bedLine);
        }
    }
    else {
        $iABresFDRFH->print( "No significant m5Cs found at FDR: " . $FDR . "\n" );
        $iABbedFDRFH->print( "No significant m5Cs found at FDR: " . $FDR . "\n" );
    }
    
    $iABresFDRFH->close();
    $iABbedFDRFH->close();
}

exit;

############################## Subroutines ############################
sub readBScallData {
    my $dataFH    = $_[0];
    my $condition = $_[1];
    my $ls        = $_[2];
    my $m5C       = $_[3];

    my $dataLine;

    while ( $dataLine = $dataFH->getline() ) {
        next if $dataLine =~ /^\#|^No significant/;
        chomp($dataLine);
        my @fields = split( '\t', $dataLine );
        my ( $mCount, $nCount, $mRate, $pvState, $pvRate ) =
          ( $fields[5], ( $fields[4] - $fields[5] ), $fields[6], $fields[14], $fields[15] );
        my $key = join( "###", @fields[ 0 .. 2 ] );

        push( @{ $m5C->{$key}->{$condition}->{m} },     $mCount );
        push( @{ $m5C->{$key}->{$condition}->{n} },     $nCount );
        push( @{ $m5C->{$key}->{$condition}->{mNorm} }, sprintf( '%01.3f', ($mCount / $ls) ) );
        push( @{ $m5C->{$key}->{$condition}->{nNorm} }, sprintf( '%01.3f', ($nCount / $ls) ) );
        push( @{ $m5C->{$key}->{$condition}->{mr} },    $mRate );
        push( @{ $m5C->{$key}->{$condition}->{ps} },    $pvState );
        push( @{ $m5C->{$key}->{$condition}->{pr} },    $pvRate );
        $m5C->{$key}->{$condition}->{r} += 1;
    }
}

sub log10 {
    return ( log( $_[0] ) / LOG10 );
}

sub log2 {
    return ( log( $_[0] ) / LOG2 );
}

sub getStats {
    my $data = $_[0];

    my @stats;

    my @cov     = map { $data->{m}->[$_] + $data->{n}->[$_] } 0 .. $#{ $data->{m} };
    my @covNorm = map { $data->{mNorm}->[$_] + $data->{nNorm}->[$_] } 0 .. $#{ $data->{mNorm} };

    my $avgCov        = avg( \@cov,          1, 0 );
    my $avgCcount     = avg( $data->{m},     1, 0 );
    my $avgCovNorm    = avg( \@covNorm,      1, 0 );
    my $avgCcountNorm = avg( $data->{mNorm}, 1, 0 );
    my $avgMrate      = avg( $data->{mr},    0, '%01.3f' );
    my $semMrate    = sprintf( '%01.3f', sem( $data->{mr} ) );
    my $stdMrate    = sprintf( '%01.3f', stdev( $data->{mr} ) );
    my $combPVstate = (grep { /NaN/ } @{ $data->{ps} }) ? 'NaN' : sprintf( '%01.6e', fisherTest( @{ $data->{ps} } ) );
    my $combPVrate  = (grep { /NaN/ } @{ $data->{pr} }) ? 'NaN' : sprintf( '%01.6e', fisherTest( @{ $data->{pr} } ) );
    my $reps        = $data->{r};

    if ( $sizeFactorLA && $sizeFactorLB ) {
        @stats = (
                   $avgCov,   $avgCcount, $avgCovNorm,  $avgCcountNorm, $avgMrate,
                   $semMrate, $stdMrate,  $combPVstate, $combPVrate,    $reps
                   );
    }
    else {
        @stats =
          ( $avgCov, $avgCcount, "NaN", "NaN", $avgMrate, $semMrate, $stdMrate, $combPVstate, $combPVrate, $reps );
    }

    return (@stats);
}

sub getStatsNoRep {
    my $data = $_[0];

    my @stats;

    my $cov        = $data->{m}->[0] + $data->{n}->[0];
    my $Ccount     = $data->{m}->[0];
    my $covNorm    = $data->{mNorm}->[0] + $data->{nNorm}->[0];
    my $CcountNorm = $data->{mNorm}->[0];
    my $mRate      = $data->{mr}->[0];
    my $pVstate    = $data->{ps}->[0];
    my $pVrate     = $data->{pr}->[0];
    my $reps       = $data->{r};

    if ( $sizeFactorLA && $sizeFactorLB ) {
        @stats = (
                   $cov, $Ccount, $covNorm, $CcountNorm,
                   $mRate,  $pVstate,  $pVrate,  $reps
                   );
    }
    else {
        @stats =
          ( $cov, $Ccount, "NaN", "NaN", $mRate, $pVstate, $pVrate, $reps );
    }

    return (@stats);
}

sub stdev {
    my ($data) = @_;

    if ( @$data == 1 ) {
        return 0;
    }
    my $average = avg( $data, 0, 0 );
    my $sqtotal = 0;

    map { $sqtotal += ( $average - $_ )**2 } @$data;

    my $std = ( $sqtotal / ( scalar @$data ) )**0.5;
    return $std;
}

sub sem {
    my ($data) = @_;

    if ( @$data == 1 ) {
        return 0;
    }
    my $std = stdev($data);
    my $sem = $std / ( scalar @$data )**0.5;

    return $sem;
}

sub avg {
    my $valArryRef = $_[0];
    my $int        = $_[1] // 0;
    my $round      = $_[2] // 0;

    my $sum = sum(@$valArryRef);

    my $avg = $sum / ( scalar @$valArryRef );

    return ( ($int) ? int($avg) : ( ($round) ? sprintf( $round, $avg ) : $avg ) );
}

sub stoufferTest {    # not used  here
    my @pValues = @_;

    my $l = scalar @pValues;
    my @w = map { 1 / $l } 1 .. $l;

    my @Zi = map { qnorm( 1 - $_ ) } @pValues;

    my $num;
    map { $num += ( $w[$_] * $Zi[$_] ) } 0 .. $l - 1;
    my $denom;
    map { $denom += $w[$_] * $w[$_] } 0 .. $l - 1;

    my $Z = $num / sqrt($denom);

    my $pValue = 1 - pnorm($Z);

    return ($pValue);
}

# Fisher test for combinded p-value
sub fisherTest {
    my @pValues = @_;

    my $df = 2 * scalar @pValues;

    my $Xsq;
    for (@pValues) { my $pv = ( $_ != 0 ) ? $_ : 1e-10; $Xsq += log($pv) }
    $Xsq = -2 * $Xsq;

    my $pValue = 1 - pchisq( $Xsq, $df, 0 );

    return ($pValue);

}

# Fisher's exact test if we do not deal with replicates
sub fisherExactTest {
    my $table = $_[0];
    
    my $pValue = 1;

    #          word2   ~word2
    #  word1    n11      n12 | n1p
    # ~word1    n21      n22 | n2p
    #           --------------
    #           np1      np2   npp    
    
    my $npp = $table->[0]->[0] + $table->[0]->[1] + $table->[1]->[0] + $table->[1]->[1];
    my $n11 = $table->[0]->[0];
    my $n1p = $table->[0]->[0] + $table->[0]->[1];
    my $np1 = $table->[0]->[0] + $table->[1]->[0];

    $pValue = calculateStatistic( n11=>$n11,
                                  n1p=>$n1p,
                                  np1=>$np1,
                                  npp=>$npp
                                 );

    # TODO: Check this! Wont't work when packed into a standalone executable using pp
    #if ( (my $errorCode = getErrorCode()) ) {
    #    say STDERR $errorCode . " - " . getErrorMessage();
    #    $pValue = 1;
    #}
    
    return ($pValue);

}


# Cochran–Mantel–Haenszel test
sub cmhTest {
    my @table = @_;

    my @numerator;
    my @denominator;
    my $n;

    # Table
    # a  b
    # c  d

    foreach my $t (@table) {

        # n = a + b + c + d
        $n = $t->[0]->[0] + $t->[0]->[1] + $t->[1]->[0] + $t->[1]->[1];

        push( @numerator, calcNumerator_( $t, $n ) );
        push( @denominator, calcDenominator_( $t, $n ) );
    }

    my $sumDenominator = sum(@denominator);

    return (0) if $sumDenominator == 0;

    my $Xsq = ( abs( sum(@numerator) ) - 0.5 )**2 / sum(@denominator);
    my $pValue = 1 - pchisq( $Xsq, 1, 0 );

    return ($pValue);

}

sub sum {
    my $sum;
    map { $sum += $_ } @_;
    return ($sum);
}

sub calcNumerator_ {
    my $table = $_[0];
    my $n     = $_[1];

    return (

        # a-(a+b)(a+c)/n
        $table->[0]->[0] - ( $table->[0]->[0] + $table->[0]->[1] ) * ( $table->[0]->[0] + $table->[1]->[0] ) / $n
        );
}

sub calcDenominator_ {
    my $table = $_[0];
    my $n     = $_[1];

    return (

        # (a+b)(a+c)(b+d)(c+d)/(n^3-n^2)
        ( $table->[0]->[0] + $table->[0]->[1] ) *
          ( $table->[0]->[0] + $table->[1]->[0] ) *
          ( $table->[0]->[1] + $table->[1]->[1] ) *
          ( $table->[1]->[0] + $table->[1]->[1] ) /
          ( $n**3 - $n**2 )
          );
}

sub FDRc {
    my $data = shift;
    my $FDR  = shift;
    
    my @fdrResults;

    # sort p-values ascending
    my @rankedData = sort { $a->[0] <=> $b->[0] } @$data;
    my $tests      = scalar @rankedData;

    # Top down
    my $k = $tests - 1;
    my $i = $k;
    while ($k >= 0) {
        my $rank = $rankedData[$k];
        my $pvT = ($k + 1) / $tests * $FDR;
        if ( $rank->[0] <= $pvT ) {
            $i = $k;
            last;
        }
        $k--;
    }

    if ( $k == - 1 ) {    # No significant m5C found
        return(undef);
    }

    my $l           = 0;
    my $pv_adj      = 0;
    
    my @fdrCorrPV = map { $rankedData[$_][0] * $tests / ( $_ + 1 ) } 0..$i;

    while ( $l <= $i ) {
        $pv_adj = min(@fdrCorrPV[$l..$i]);

        # return FDR controlled data: 
        # [l][0]: raw p-value
        # [l][1]: FDR adjusted p-value
        # [l][2]: optional data if exists (untouched)
        $fdrResults[$l][0] = $rankedData[$l][0];
        $fdrResults[$l][1] = $pv_adj;
        $fdrResults[$l][2] = $rankedData[$l][1] if ( defined($rankedData[$l][1]) );

        $l++;
    }
    return \@fdrResults;
}

sub usage {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self [options] [-h] [-man] [--version]

Required options all of:
    -condition-files-a|-fa : meRanCall result files from condition A. All result
                             files from the first condition can be specified as
                             a comma separated list.
                             (default: not set)

    -condition-files-b|-fb : meRanCall result files from condition B. All result
                             files from the second condition can be specified as
                             a comma separated list.
                             (default: not set)

Options:

    -condition-name-a|-na  : Name of the condition A. This name is used in the
                             file names for the meRanCompare results.
                             (default: ConditionA)

    -condition-name-b|-nb  : Name of the condition B. This name is used in the
                             file names for the meRanCompare results.
                             (default: ConditionB)

    -size-factors-a|-sfa   : Library size factors for samples in condition A/B
    -size-factors-b|-sfb     specified as comma separated list. These size factors
                             are used to calculate normalized counts.
    
                             e.g. -sfa 0.6673,0.6609,0.7347
                                  -sfb 0.9559,1.4098,2.3802
                             
                             The ordering of the individual size factors has to match
                             the ordering of the meRanCall result files for the
                             corresponding conditions (see -fa and -fb option).
                             
                             The size factors can be calculated using the meRanTK
                             tool "estimatSizeFactors.pl". Alternatively one can
                             use htseq-count and DESeq2 "estimateSizeFactors"
                             
                             (default: not set, no normailzed counts are reported)

    -minRep|-mr            : number of replicates a m5C candidate has to be present
                             so that it is considered as high confidence call.
                             (default: 2)

    -sig|-s                : p-value below which the differential methylation will
                             be reported as significant.
                             (default: 0.01)

    -minFC|-fc             : minimum fold change above which the differential methylation
                             will be reported. In bisulfite mode (default) the foldchange
                             will be calculated as ratio of the methylation rate in
                             condition A and B. In "aza" mode this will be the ratio of
                             the (normalized) coverage at the specific positon in
                             condition A and B.
                             (default: not set, report all significant (see -sig) changes)

    -fdr                   : FDR, false discovery rate
                             (default: 0.01)

    -azaMode|-aza          : run meRanCompare in Aza-IP mode. In this mode the IP enrichment (A)
                             over control (B) is calculated for each candidate m5C which is compared.
                             (default: not set)


    --version              : Print the program version and exit.
    -h|help                : Print the program help information.
    -man                   : Print a detailed documentation.
    
    -debug|-d              : Print some debugging information.
    
EOF
}

__END__

=head1 NAME

meRanCompare - Find differntial RNA cytosine methylation

=head1 SYNOPSIS

=head2 compare cytosine methylation from two conditions and find differences

=over 2
 
### Example with 3 replicates

 meRanCompare \
 -fa meRanCallResult_A_rep1.txt,meRanCallResult_A_rep2.txt,meRanCallResult_A_rep3.txt \
 -fb meRanCallResult_B_rep1.txt,meRanCallResult_B_rep2.txt,meRanCallResult_B_rep3.txt \
 -na "myConditonA" \
 -nb "myConditonB" \
 -sfa 0.6673,0.6609,0.7347 \
 -sfb 0.9559,1.4098,2.3802 \
 -mr 2 \
 -sig 0.01


 The command above finds differentially methylated C's between two conditions using result
 files obtained from meRanCall. It considers m5C calls from all replicates of
 each condition and uses a Cochran–Mantel–Haenszel test or Fisher's exact test (no replicates)
 to find significant methylation differences of any m5C called. The significance threshold is
 specified by the "-sig" option.
 
 meRanCompare will generate 6 result files, containing the names specified by the "-na, -nb"
 options:
        - uniq_myConditionA.txt
        - uniq_myConditionA.bed
        - uniq_myConditionB.txt
        - uniq_myConditionB.bed
        - intersect_myConditonA_myCondition_B.txt
        - intersect_myConditonA_myCondition_B.bed
 
 The uniq_*.txt result files are in tab separated format and contain the following
 data-fields for each methylated C that was only found in one condition:

 1.  SeqID               : sequence ID from reference database (e.g. chr1)
 2.  refPos              : postion of the methylated C on the reference sequence
 3.  refStrand           : genomic strand
 4.  Name                : arbitrary methylation site name
 5.  avg/cov             : average coverage
 6.  avg/C_count         : average C count
 7.  norm(avg/cov)       : normalized average coverage (if size factors are provided)
 8.  norm(avg/C_count)   : normailzed average C count (if size factors are provided)
 9.  avg/methRate        : average methylation rate
 10. sem/methRate        : standard error of mean methylation rate
 11. stdev/methRate      : standard deviation methylation rate
 12. comb/p-value_mState : combined p-value of methylation state
 13. comb/p-value_mRate  : combined p-value of methylation rate
 14. rep                 : found in this many replicates

 For these uniq m5C also corresponding BED files will be generated.
 
 The intersect_*.txt result file is in tab separated format and contains the following
 data-fields for each methylated C that was differentially methylated between the two
 conditions:

 1.  SeqID                 : sequence ID from reference database (e.g. chr1)
 2.  refPos                : postion of the methylated C on the reference sequence
 3.  refStrand             : genomic strand
 4.  Name                  : arbitrary methylation site name
 5.  p-value               : p-value for differential methylation test
                            (Cochran–Mantel–Haenszel or Fishers exact test)
 6.  avg/cov A             : average coverage in conditon A
 7.  avg/C_count A         : average C count in conditon A
 8.  norm(avg/cov)         : normalized average coverage (if size factors are provided)
 9.  norm(avg/C_count)     : normailzed average C count (if size factors are provided)
 10. avg/methRate A        : average methylation rate in conditon A
 11. sem/methRate A        : standard error of mean methylation rate in conditon A
 12. stdev/methRate A      : standard deviation methylation rate in conditon A
 13. comb/p-value_mState A : combined p-value of methylation state in conditon A
 14. comb/p-value_mRate A  : combined p-value of methylation rate in conditon A
 15. rep A                 : found in this many replicates in conditon A
 16. avg/cov B             : average coverage in conditon B
 17. avg/C_count B         : average C count in conditon B
 18. norm(avg/cov)         : normalized average coverage (if size factors are provided)
 19. norm(avg/C_count)     : normailzed average C count (if size factors are provided)
 20. avg/methRate B        : average methylation rate in conditon B
 21. sem/methRate B        : standard error of mean methylation rate in conditon B
 22. stdev/methRate B      : standard deviation methylation rate in conditon B
 23. comb/p-value_mState B : combined p-value of methylation state in conditon B
 24. comb/p-value_mRate B  : combined p-value of methylation rate in conditon B
 25. rep B                 : found in this many replicates in conditon B


 For these differentially methylated m5Cs a corresponding BED file will be generated.
 
=back

=head1 DEPENDENCIES

 The program requires additional perl modules and depending on your perl installation
 you might need do add the following modules:

 Math::CDF
 
 These modules should be availble via CPAN or depending on your OS via the package
 manager.

=head1 TESTED WITH:

=over

=item *
Perl 5.18.1 (Centos 6.5)

=item *
Perl 5.18.2 (Centos 6.5)

=item *
Perl 5.18.2 (RHEL 6.5)

=item *
Perl 5.24.0 (RHEL 6.8)

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

  Copyright (C) 2016 D.Rieder

=head1 COPYRIGHT

  Copyright (C) 2016 Dietmar Rieder. All Rights Reserved.

=head1 AUTHOR

Dietmar Rieder

=head1 CONTACT

dietmar . rieder (at) i-med . ac . at

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.


exit;
