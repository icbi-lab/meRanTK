#!/usr/local/bioinf/perl/perl-5.24.0/bin/perl -w
##!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.18.2/bin/perl -w
#
#  meRanCall.pl
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

use Pod::Usage;
use Bio::DB::Sam;
use List::Util qw( reduce sum min max );
use IO::File;
use File::Basename;

use Parallel::ForkManager;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Fcntl qw(:flock SEEK_END);
use Math::CDF qw( pnorm pbinom);
use Text::NSP::Measures::2D::Fisher::twotailed;

use constant {
               LOG10          => log(10),
               REGIONSIZE     => 2_000_000,
               Z95            => 1.96,
               Z95SQ          => 1.96**2,
               MAX_PILEUP_CNT => 500_000,
               };

my $DEBUG = 0;

my $VERSION = '1.2.1b';
my $version;
my $help;
my $man;

my $samInFile       = "";
my $resOutFile      = "";
my $fastaInFile     = "";
my $procs           = 1;
my $fskip5          = 0;
my $fskip3          = 0;
my $rskip5          = 0;
my $rskip3          = 0;
my $oriReadLen      = 100;
my $readDir         = "fr";
my $minMethRate     = 0.2;
my $minMutRate      = 0.8;
my $minBaseQ        = 30;
my $qsOffset        = 33;
my $minCov          = 10;
my $maxDuplicates   = 0;
my $errorInterval   = 0;
my $conversionRate;
my $reportUP        = "";    # report mutation but unmethylated
my $genomeDBref     = "";    # we handle mappings to a nondirected DB (e.g. genomic)
my $transcriptDBref = "";    # we handle mappings to a directed DB (e.g. transcriptome)
my $haveZGsamTag    = "";    # we have a sam/bam file that holds ZG tag with gene name associated with the read.
my $seqContext;
my $narrowPeak;
my $bed63;
my $bedFileIN;
my $chrPrefix;
my $FDR;
my $FDRrate;
my $calculateConvR;
my @controlSeqIDs;
my @exTargets;
my $azaMode;
my $fisherTest = 0;

GetOptions(
            'fasta|f=s'              => \$fastaInFile,
            'sam|bam|s=s'            => \$samInFile,
            'procs|p=i'              => \$procs,
            'result|o=s'             => \$resOutFile,
            'regions|bi=s'           => \$bedFileIN,
            'chrRrefix|cpf=s'        => \$chrPrefix,
            'fskip5|fs5=i'           => \$fskip5,
            'fskip3|fs3=i'           => \$fskip3,
            'rskip5|rs5=i'           => \$rskip5,
            'rskip3|rs3=i'           => \$rskip3,
            'QSoffset|qso=i'         => \$qsOffset,
            'readLength|rl=i'        => \$oriReadLen,
            'readDir|rd=s'           => \$readDir,
            'minMethR|mr=f'          => \$minMethRate,
            'minMutR|mutR=f'         => \$minMutRate,
            'minBaseQ|mBQ=i'         => \$minBaseQ,
            'minCov|mcov=i'          => \$minCov,
            'maxDup|md=i'            => \$maxDuplicates,
            'fisherTest|fet'         => \$fisherTest,
            'fdr=f'                  => \$FDR,
            'fdrRate'                => \$FDRrate,
            'conversionRate|cr=f'    => \$conversionRate,
            'calcConvRate|ccr'       => \$calculateConvR,
            'controlSeqID|cSeqID=s'  => \@controlSeqIDs,
            'excludeSeqID|exSeqID=s' => \@exTargets,
            'errorInterval|ei=f'     => \$errorInterval,
            'reportUP|rUP'           => \$reportUP,
            'genomeDBref|gref'       => \$genomeDBref,
            'transcriptDBref|tref'   => \$transcriptDBref,
            'narrowPeak|np'          => \$narrowPeak,
            'bed63'                  => \$bed63,
            'seqContext|sc=i'        => \$seqContext,
            'haveZG|zg'              => \$haveZGsamTag,
            'azaMode|aza'            => \$azaMode,
            'help|h'                 => \$help,
            'man|m'                  => \$man,
            'version'                => \$version,
            'debug|d'                => \$DEBUG,
            );

my $debug;
if ($DEBUG) {
    eval "use Data::Dumper";
    $debug = sub { print STDERR "\n" . join( " ", @_ ) . "\n"; }
}
else {
    $debug = sub { return; };
}

my $callMethStatus;
if ($azaMode) {
    $callMethStatus = \&callMethStatusAZA;
    $conversionRate = ($conversionRate) ? $conversionRate : -1;
}
else {
    $callMethStatus = \&callMethStatus;
    $conversionRate = ($conversionRate) ? $conversionRate : -1;
}
if ( ($conversionRate == -1) && (!$calculateConvR) ) {
    if ($FDR) {
        say STDERR "WARNING: No conversion rate: p-value for methylation state will be calculate only based on the expected sequencing error";
    }
}

# globals
my %m5Cs;
my @strands;
my $bedwriter;
my $bed63File;
my $bednpFile;
my $bed63FH;
my $bednpFH;
my $resFH;
my $fdrRESoutFile;
my $FDRresFH;
my $FDRbed63File;
my $FDRbednpFile;
my $FDRbed63FH;
my $FDRbednpFH;
my %regionData;
my $regionCount = 0;

# set baseline sequencing error
my $BASELINE_ERR = 10**(-$minBaseQ/10);

my $meth_caller = \&meth_caller;

usage() and exit(0) if ($help);
pod2usage( -verbose => 3 ) if $man;

say $VERSION and exit(0) if ($version);

if ( ( $minMethRate > 1 ) ) {
    say STDERR "ERROR: min methlyation rate should be 0 < mr <= 1 ";
    exit(1);
}
if ( ( $minMutRate > 1 ) ) {
    say STDERR "ERROR: min mutation rate should be 0 < mutR <= 1 ";
    exit(1);
}
if ( ( $conversionRate > 1 ) ) {
    say STDERR "ERROR: conversion rate should be 0 < cr <= 1 ";
    exit(1);
}
if ( ( $errorInterval > 1 ) ) {
    say STDERR "ERROR: error interval should be 0 < mutR <= 1 ";
    exit(1);
}
if ( ( defined($FDR) ) && ( $FDR > 1 ) ) {
    say STDERR "ERROR: FDR should be 0 < mutR <= 1 ";
    exit(1);
}
if ( ( $minBaseQ > 41 ) ) {
    say STDERR "ERROR: min base quality score should be 0 < mBQ <= 41 ";
    exit(1);
}
if ( ( $readDir ne "fr" ) && ( $readDir ne "rf" ) ) {
    say STDERR "ERROR: readDir can only be fr or rf, got: " . $readDir;
    exit(1);
}

if ( !-r $samInFile ) {
    say STDERR "ERROR: can not read BAM/SAM file, please provide either a BAM or a SAM file";
    exit(1);
}

if ( !-r $fastaInFile ) {
    say STDERR
"ERROR: can not read FASTA reference file, please provide the FASTA reference that was used to create the BAM/SAM file";
    exit(1);
}

my $bamFile;
my ( $fname, $fpath, $fsuffix ) = fileparse( $samInFile, qr/\.[^.]*$/ );
if ( $fsuffix eq '.sam' ) {
    $bamFile = sam_to_bam($samInFile);
}
elsif ( $fsuffix eq '.bam' ) {
    $bamFile = $samInFile;
}
else {
    say STDERR "ERROR: invalid file format, please provide either a BAM or a SAM file";
    exit(1);
}

if ( $qsOffset == 64 ) {
    say STDOUT "Your sequencing base quality score offset is " . $qsOffset . " ... adjusting PHERD scores accordingly";
    $qsOffset = 31;
}
elsif ( $qsOffset == 33 ) {
    $qsOffset = 0;
}
else {
    say STDERR "Your sequencing base quality score offset " . $qsOffset . " is not supported. Exiting ...";
    exit(1);
}

my $skipPos = $fskip5 + $fskip3 + $rskip5 + $rskip3;

my $getGeneName;
if ($haveZGsamTag) {
    $getGeneName = sub {
        return ( $_[0]->aux_get('ZG') );
    };
}
else {
    $getGeneName = sub { return "-" };
}

# Create SAM/BAM object
my $sam = Bio::DB::Sam->new(
                             -bam           => $bamFile,
                             -fasta         => $fastaInFile,
                             -expand_flags  => 1,
                             -split_splices => 1,
                             -autoindex     => 1,
                             );

# Set max pileup count
$sam->max_pileup_cnt(MAX_PILEUP_CNT);

# get the targets
my @targets = $sam->seq_ids;
if ( ( scalar @targets ) == 0 ) {
    say STDERR "ERROR: No valid alignments found in the supplied BAM/SAM file: " . $samInFile;
    exit(1);
}

if ($calculateConvR) {
    $conversionRate = 1;
    $errorInterval  = 0;
    $FDR            = undef;
    $minMethRate    = 0.0;
    if ( !@controlSeqIDs ) {
        say STDERR
          "ERROR: sequenceID(s) of the control sequcence(s) not specified, use -cSeqID <seqID1> [-cSeqID <seqID2> ...]";
        exit(1);
    }
    my %seen;
    @targets = grep { $seen{$_}++ } @controlSeqIDs, @targets;
    if ( ( scalar @targets ) == 0 ) {
        say STDERR
          "ERROR: No valid alignment for the supplied control sequnece IDs found in the supplied BAM/SAM file: "
          . $samInFile;
        exit(1);
    }
    my %checkIDs;
    @checkIDs{@targets} = @targets;
    my @notFound = grep { !exists $checkIDs{$_} } @controlSeqIDs;
    foreach my $id (@notFound) {
        say STDERR "WARNING: No valid alignment for " . $id . " found in the supplied BAM/SAM file: " . $samInFile;
    }
}

my $faidxFile = $fastaInFile . ".fai";
my $fai;
if ( !-r $faidxFile ) {
    if ( checkWirteToDir($fpath) ) {
        $fai = $sam->fai;    # try to open the bam index, should be recreated if this fails.
        undef($fai);
    }
    else {
        say STDERR "ERROR: could neither read nor create the fasta index file " . $faidxFile . " : check permissions";
        exit(1);
    }
}

# Test if fasta matches SAM/BAM
$fai = $sam->fai;
if ( defined( $targets[0] ) ) {
    my $seqCheck = $targets[0] . ":1-2";
    if ( !$fai->fetch($seqCheck) ) {
        say STDERR "ERROR: the provided fasta file and SAM/BAM seem not to match:\n\n\t"
          . $fastaInFile . " <> "
          . $samInFile
          . "\n\nPlease provide the fasta file that was used for mapping the reads.";
        exit(1);
    }
    
    my $refLen = getRefLengths($faidxFile);

    foreach my $seqid (@targets) {
        if ( $sam->length($seqid) != $refLen->{$seqid} ) {
            say STDERR "ERROR: the provided fasta file and SAM/BAM seem not to match:\n\n\t"
              . $fastaInFile . " <> "
              . $samInFile
              . "\n\nPlease provide the fasta file that was used for mapping the reads.";
            exit(1);
        }
    }
}
undef($fai);

( $fname, $fpath, $fsuffix ) = fileparse( $bamFile, ".bam" );
my $baidxFile = $fpath . "/" . $fname . ".bai";
my $bai;
my $header;
if ( !-r $baidxFile ) {
    if ( checkWirteToDir($fpath) ) {
        $bai = $sam->bam_index;    # try to open the bam index, should be recreated if this fails.
        undef($bai);
    }
    else {
        say STDERR "ERROR: could neither read nor create the BAM index file " . $baidxFile . " : check permissions";
        exit(1);
    }
}

if ( $bedFileIN ) {
    readBED( $bedFileIN, \%regionData, $regionCount, $chrPrefix );
}

# result file handle
if ( ( !$resOutFile ) && ( !$calculateConvR ) ) {
    say STDERR "ERROR: No result file specified";
    exit(1);
}
else {
    my ( $fname, $fpath, $fsuffix ) = fileparse( $resOutFile, qw/\.[^.]*$/ );
    checkDir($fpath);
}
$resFH = IO::File->new( $resOutFile, O_RDWR | O_CREAT | O_TRUNC )
  || die( "Can not create the result file: " . $resOutFile . " : " . $! )
  unless ($calculateConvR);

if ($FDR) {
    ( $fname, $fpath, $fsuffix ) = fileparse( $resOutFile, qr/\.[^.]*$/ );
    $fdrRESoutFile = nicePath($fpath . "/" . $fname . "_FDR_" . $FDR . ".txt");
    $FDRresFH = IO::File->new( $fdrRESoutFile, O_RDWR | O_CREAT | O_TRUNC )
      || die( "Can not create the FDR contolled result file: " . $fdrRESoutFile . " : " . $! );
}

if ( !$genomeDBref && !$transcriptDBref ) {
    say STDERR "ERROR: please specify the type of reference your reads in the SAM|BAM file were aligned to:";
    say STDERR "\t\ttranscriptome reference? use -tref";
    say STDERR "\t\tgenome reference? use -gref";
    exit(1);
}

if ( $genomeDBref && $transcriptDBref ) {
    say STDERR
      "ERROR: please specify only one type of reference your reads in the SAM|BAM file were aligned to: -gref|-tref";
    exit(1);
}

if ($genomeDBref) {
    $m5Cs{'+'} = {};
    $m5Cs{'-'} = {};
    $m5Cs{dup} = {};

    @strands = ( "+", "-" );

    unless ($calculateConvR) {

        my ( $fname, $fpath, $fsuffix ) = fileparse( $resOutFile, qr/\.[^.]*$/ );
        if ($bed63) {
            $bed63File = nicePath($fpath . "/" . $fname . ".bed");
            $bed63FH = IO::File->new( $bed63File, O_RDWR | O_CREAT | O_TRUNC )
              || die( "Can not create the BED6+3 file: " . $bed63File . ": " . $! );

            if ($FDR) {
                $FDRbed63File = nicePath($fpath . "/" . $fname . "_FDR_" . $FDR . ".bed");
                $FDRbed63FH = IO::File->new( $FDRbed63File, O_RDWR | O_CREAT | O_TRUNC )
                  || die( "Can not create the FDR contolled BED6+3 file: " . $FDRbed63File . " : " . $! );
            }
        }
        if ($narrowPeak) {
            $bednpFile = nicePath($fpath . "/" . $fname . "_narrowPeak.bed");
            $bednpFH = IO::File->new( $bednpFile, O_RDWR | O_CREAT | O_TRUNC )
              || die( "Can not create the narrowPeak BED file: " . $bednpFile . ": " . $! );

            if ($FDR) {
                $FDRbednpFile = nicePath($fpath . "/" . $fname . "_FDR_" . $FDR . "_narrowPeak.bed");
                $FDRbednpFH = IO::File->new( $FDRbednpFile, O_RDWR | O_CREAT | O_TRUNC )
                  || die( "Can not create the FDR contolled narrowPeak BED file: " . $FDRbednpFile . " : " . $! );
            }
        }

        $bedwriter = \&bedWriter;
    }
}
else {
    $m5Cs{'+'} = {};
    $m5Cs{dup} = {};

    @strands = ("+");

    $bedwriter = sub { return; };

    if ( $bed63 || $narrowPeak ) {
        say STDERR "WARNING: not creating BED files since reads were mapped to transcripts";
    }
}

####### Let the fun begin #######
my $MAX_PROCESSES = $procs;
my $pm = new Parallel::ForkManager( $MAX_PROCESSES, "/dev/shm" );

my $count              = 0;
my $totalMethRefCs     = 0;
my $totalMutMethCs     = 0;
my $totalRefCs         = 0;
my $totalRefCsanalyzed = 0;
my $methCT             = 0;
my $methCG             = 0;
my $methC              = 0;
my $methCTanalyzed     = 0;
my $methCanalyzed      = 0;
my $totalFilteredDups  = 0;
my $m5Cnr              = 0;
my $nrOftargets        = ( $bedFileIN ) ? $regionCount : scalar @targets;
my $totpctDone         = 0;
my $seqpctDone         = 0;
my $end                = 0;
my $start              = -99999;
my $tlength            = 0;
my $halfminCov         = ( $minCov * 0.5 );
my $fetMinSuccess      = $minMethRate * $minCov;
my $fetMaxFail         = 10 - $fetMinSuccess;
my $resNr              = 0;
my @fdrData;

print "Working on: " . $samInFile . "\n";
if ( ! $bedFileIN ) {
    print "No region BED file specified: calling m5Cs on entire alignment file: " . $samInFile . " ...\n";
}
else {
    print "Region BED file specified: calling m5Cs on regions in file: " . $bedFileIN . " ...\n";
}
print "Starting to process: " . $nrOftargets . " targets on " . $procs . " CPUs ...\n\n";

# print Header line
my $headerLine = "\#"
  . join( "\t",
          qw(SeqID refPos refStrand refBase cov C_count methRate mut_count mutRate CalledBase CB_count state 95_CI_lower 95_CI_upper p-value_mState p-value_mRate score seqContext geneName candidateName)
          ) . "\n";
$resFH->print($headerLine) unless ($calculateConvR);

if ($FDR) {
    if ( !$FDRrate ) {
        $headerLine = "\#"
          . join( "\t",
                  qw(SeqID refPos refStrand refBase cov C_count methRate mut_count mutRate CalledBase CB_count state 95_CI_lower 95_CI_upper p-value_mState p-value_mRate p-value_mState_adj score seqContext geneName candidateName)
                  ) . "\n";
    }
    else {
        $headerLine = "\#"
          . join( "\t",
                  qw(SeqID refPos refStrand refBase cov C_count methRate mut_count mutRate CalledBase CB_count state 95_CI_lower 95_CI_upper p-value_mState p-value_mRate p-value_mRate_adj score seqContext geneName candidateName)
                  ) . "\n";
    }

    $FDRresFH->print($headerLine);
}

# run parallel jobs
$pm->run_on_finish( \&resultCollector );
$pm->run_on_wait( \&tellStatus, 0.5 )  unless ( $bedFileIN );

my $wSid;
my $REGIONSIZE;

if ( ! $bedFileIN ) {
    foreach my $seqid (@targets) {
    
        $wSid = $seqid;

        next if (grep /^$seqid$/, @exTargets);
        
        printf( "processing %i sequences: \[%s - %02.2f%%\] [overall - %02.2f%%] done ...\r",
                $nrOftargets, $wSid, $seqpctDone, $totpctDone );
    
        $tlength = $sam->length($seqid);
        $REGIONSIZE = ( $tlength <= REGIONSIZE ) ? $tlength : REGIONSIZE;
    
        $end   = 0;
        $start = -1 * $REGIONSIZE + 1;
    
        while ( $end < $tlength ) {
    
            $start += $REGIONSIZE;
            $end = ( ( $end + $REGIONSIZE ) > $tlength ) ? $tlength : ( $end + $REGIONSIZE );
            my $region = $seqid . ":" . $start . "-" . $end;

            my $pid = $pm->start($seqid) and next;
            $sam->clone;

            $sam->fast_pileup( $region, $meth_caller );

            $pm->finish( 0, \%m5Cs );
    
        }
        $totpctDone = ++$count * 100 / $nrOftargets;
    }
}
else {
 
    foreach my $seqid ( keys(%regionData) ) {
        next if ( ! grep( /^$seqid$/, @targets ) );
        next if (grep /^$seqid$/, @exTargets);

        $wSid = $seqid;
                
        foreach my $range ( @{ $regionData{$seqid}->{range} } ) {
            
            my $region = $seqid . ":" . $range->[0] . "-" . $range->[1];
            &$debug($region);

            printf( "\nWorking on: %s ", $region);
 
            $tlength = $range->[1] - $range->[0];
            $REGIONSIZE = ( $tlength <= REGIONSIZE ) ? $tlength : REGIONSIZE;

            $end   = $range->[0];
            $start = $range->[0] - $REGIONSIZE;
    
            while ( $end < $range->[1] ) {
    
                $start += $REGIONSIZE;
                $end = ( ( $end + $REGIONSIZE ) > $range->[1] ) ? $range->[1] : ( $end + $REGIONSIZE );
                my $sregion = $seqid . ":" . $start . "-" . $end;
                printf( "...\n\t... subregion: %s ", $sregion) unless ( ($start == $range->[0]) && ($end == $range->[1]) );

                my $pid = $pm->start($seqid) and next;
                $sam->clone;
        
                $sam->fast_pileup( $sregion, $meth_caller );
        
                $pm->finish( 0, \%m5Cs );    
            }

            
        }
        
        $totpctDone = ++$count * 100 / $nrOftargets;
    }    
}
$pm->wait_all_children;
$resFH->close() unless ($calculateConvR);
$bed63FH->close() if ( $genomeDBref && $bed63 );
$bednpFH->close() if ( $genomeDBref && $narrowPeak );
undef($resFH);
undef($bed63FH);
undef($bednpFH);

if ($FDR) {

    my $fdrResults = FDRc( \@fdrData, $FDR );

    if ( defined($fdrResults) ) {
        foreach my $record (@$fdrResults) {

            my @resFields = split( '\t', $record->[2] );
            my $pv_adj_nice = sprintf( "%01.6e", $record->[1]);
            splice( @resFields, 16, 0, $pv_adj_nice );

            my $resLine = join( "\t", @resFields );
            $FDRresFH->print($resLine);

            # remove newline from m5C name, coming form orignal result line
            chomp($resFields[20]);
            
            if ($narrowPeak) {
                &$bedwriter( $resFields[0],  $resFields[1], $resFields[20], $resFields[2], $resFields[6],
                             $resFields[14], $resFields[4], $resFields[17], 64,            $FDRbednpFH );
            }
            if ($bed63) {
                &$bedwriter( $resFields[0],  $resFields[1], $resFields[20], $resFields[2], $resFields[6],
                             $resFields[14], $resFields[4], $resFields[17], 63,            $FDRbed63FH );
            }
        }
    }
    else {
        $FDRresFH->print( "No significant m5Cs found at FDR: " . $FDR . "\n" );
        $FDRbed63FH->print( "No significant m5Cs found at FDR: " . $FDR . "\n" ) if ( $genomeDBref && $bed63 );
        $FDRbednpFH->print( "No significant m5Cs found at FDR: " . $FDR . "\n" ) if ( $genomeDBref && $narrowPeak );
    }

    $FDRresFH->close();
    $FDRbed63FH->close() if ( $genomeDBref && $bed63 );
    $FDRbednpFH->close() if ( $genomeDBref && $narrowPeak );
    undef($FDRresFH);
    undef($FDRbed63FH);
    undef($FDRbednpFH);
}

# done methylation calling
print "\nDone...\n";

my $totalConvRate_analyzed =
  ( $methCTanalyzed > 0 ) ? sprintf( "%01.4f", 1 - ( $methCanalyzed / $methCTanalyzed ) ) : 0;
my $totalConvRate_estimated = 0;
if ($azaMode) {
    $totalConvRate_estimated = ( $methCG > 0 ) ? sprintf( "%01.4f", 1 - ( $methC / $methCG ) ) : 0;
}
else {
    $totalConvRate_estimated = ( $methCT > 0 ) ? sprintf( "%01.4f", 1 - ( $methC / $methCT ) ) : 0;
}

### Print out some stats #####
print STDOUT "

Analyzed: " . $samInFile . "

Summary:
Analysis of Cs with minimum coverage of " . $minCov . "

Total duplicate reads filtered:\t" . $totalFilteredDups . "

Total analyzed Cs on reference:\t" . $totalRefCsanalyzed . "
Total analyzed methylated Cs ("
  . (($azaMode) ? (">= " . ($minMethRate * 100)) : (( "<= " . ( 1 - $minMethRate ) * 100 ))) . "% " . (($azaMode) ? ("C->G") : ("C->T")) . " conversion) on reference:\t" . $totalMethRefCs . "
Total analyzed " . (($azaMode) ? "converted" : "unconverted") . " Cs from queries:\t" . $methCanalyzed . "
Total analyzed " . (($azaMode) ? "converted" : "unconverted") . " Cs from mutation:\t" . $totalMutMethCs . "
Total analyzed C to " . (($azaMode) ? "G" : "T") . " conversion rate:\t" . $totalConvRate_analyzed . "


Summary over all Cs:

Total Cs on reference covered by seq data:\t" . $totalRefCs . "
Total Cs from queries " . (($azaMode) ? "converted" : "unconverted") . ":\t" . $methC . "
Total C to " . (($azaMode) ? "G" : "T") . " conversion rate estimated:\t" . $totalConvRate_estimated . "

Result file: " . $resOutFile . "\n";

if ($FDR) {
    print STDOUT "Result file at FDR " . $FDR . ": " . $fdrRESoutFile . "\n";
}
if ($bed63) {
    print STDOUT "Results in BED6+3: " . $bed63File . "\n";
    if ($FDR) {
        print STDOUT "Results at FDR " . $FDR . " in BED6+3: " . $FDRbed63File . "\n";
    }
}
if ($narrowPeak) {
    print STDOUT "Results in narrowPeak BED: " . $bednpFile . "\n";
    if ($FDR) {
        print STDOUT "Results at FDR " . $FDR . " in narrowPeak BED: " . $FDRbednpFile . "\n";
    }
}
print STDOUT "\n";

#### Subroutines #####

sub meth_caller {
    my ( $seqid, $pos, $p ) = @_;
    my $refbase;
    my $context;
    my $geneName;
    my %baseCount_;    #= (A => 0, C => 0, T => 0, G => 0, N => 0);
    my $baseCount = \%baseCount_;
    my %baseStrandCount;
    my %readMapStart = ();

    return unless $pos >= $start && $pos <= $end;

    my $total = 0;

    if (@$p) {
        $total = 0;

        ( $refbase, $context ) = getRefSeq( $sam, $seqid, $pos, $seqContext );

        $baseCount->{'+'} = { A => 0, C => 0, T => 0, G => 0, N => 0 };
        $baseCount->{'-'} = { A => 0, C => 0, T => 0, G => 0, N => 0 };

        for my $pileup (@$p) {


            next if $pileup->indel or $pileup->is_refskip;    # don't deal with these ;-)

            my $a = $pileup->alignment;

            # $readMapStart{$a->start} += 1;
            my $dupMapID = $a->start . '_' . $a->cigar_str . '_' . $a->flag;
            $readMapStart{$dupMapID} += 1;

            # filter Duplicates
            # if ( ($readMapStart{$a->start} > $maxDuplicates) && ($maxDuplicates > 0) ) {
            if ( ( $readMapStart{$dupMapID} > $maxDuplicates ) && ( $maxDuplicates > 0 ) ) {
                $m5Cs{dup}->{ $a->qname } = 1;
                &$debug( "Found read duplicate: ", $a->start, $a->length, $a->qname );
                next;
            }

            # print $a->qname . " : " . $a->cigar_str . " : ". length($a->qseq) .  " : " . $pileup->pos . " : " . $pos . "\n";
            my $reversed = $a->reversed;

            my $paired    = $a->paired;
            my $firstMate = 1;            # set initially true also for SE reads
            if ($paired) {
                $firstMate = ( ( $a->flag & 0x40 ) != 0 ) ? 1 : 0;
            }

            # Get the read associated gene Name
            $geneName = &$getGeneName($a);

            unless ( $skipPos == 0 ) {

                my $posOnread = $pileup->qpos;    # 0-based read position

                my $qLen            = length( $a->qseq );
                my $trimmedBases    = $oriReadLen - $qLen;
                my $corr            = $fskip3 - $trimmedBases;
                my $fskip3corrected = ( $corr >= 0 ) ? $corr : 0;

               # Unfortuantely we can not infer a correction at 5', we assume the reads did not get trimed at 5', FIX me
                my $fskip5corrected = $fskip5;

                if ($paired) {                    # paired end
                    if ($firstMate) {             # First mate
                        if ($reversed) {          # reversed
                            next if ( $posOnread > ( $qLen - $fskip5corrected ) );
                            next if ( $posOnread < $fskip3corrected );
                        }
                        else {                    # forward
                            next if ( $posOnread < $fskip5corrected );
                            next if ( $posOnread > ( $qLen - $fskip3corrected ) );
                        }
                    }
                    else {                        # Second mate
                        my $corr            = $rskip3 - $trimmedBases;
                        my $rskip3corrected = ( $corr >= 0 ) ? $corr : 0;
                        my $rskip5corrected = $rskip5;

                        if ($reversed) {          # reversed
                            next if ( $posOnread > ( $qLen - $rskip5corrected ) );
                            next if ( $posOnread < $rskip3corrected );
                        }
                        else {                    # forward
                            next if ( $posOnread < $rskip5corrected );
                            next if ( $posOnread > ( $qLen - $rskip3corrected ) );
                        }
                    }
                }
                else {                            # single end
                    if ($reversed) {
                        next if ( $posOnread > ( $qLen - $fskip5corrected ) );
                        next if ( $posOnread < $fskip3corrected );
                    }
                    else {
                        next if ( $posOnread < $fskip5corrected );
                        next if ( $posOnread > ( $qLen - $fskip3corrected ) );
                    }
                }
            }    # end skip

            my $qbase = uc( substr( $a->qseq, $pileup->qpos, 1 ) );
            # print "RB: " . $refbase . " QB: " . $qbase . "\n";
            next if ( ( $qbase eq 'N' ) || ( $refbase eq 'N' ) );

            my $strand = "+";

            ### genomeDBref && strand
            if ($genomeDBref) {
              if($readDir eq "fr") {
                if ( ($reversed) && ($firstMate) ) {    # reverse the mapped base
                    $qbase =~ tr/[ACGT]/[TGCA]/;
                    $strand = "-";
                }
                elsif ( ($reversed) && ( !$firstMate ) ) {
                    $strand = "+";
                }
                elsif ( ( !$reversed ) && ( !$firstMate ) ) {    # reverse the mapped base
                    $qbase =~ tr/[ACGT]/[TGCA]/;
                    $strand = "-";
                }
                elsif ( ( !$reversed ) && ($firstMate) ) {       # reverse the mapped base
                    $strand = "+";
                }
              }
              else {
                if ( ($reversed) && ($firstMate) ) {
                   $strand = "+";
                }
                elsif ( ($reversed) && ( !$firstMate ) ) {    # reverse the mapped base
                    $qbase =~ tr/[ACGT]/[TGCA]/;
                    $strand = "-";
                }
                elsif ( ( !$reversed ) && ( !$firstMate ) ) {    # reverse the mapped base
                    $qbase =~ tr/[ACGT]/[TGCA]/;
                    $strand = "+";
                }
                elsif ( ( !$reversed ) && ($firstMate) ) {       # reverse the mapped base
                    $strand = "-";
                }
              }
            }
            else {
                if ( ($reversed) && ($firstMate) && ( $readDir ne "rf") ) {
                    warn("First mate read mapped in wrong direction, skipping read...\n");
                    next;
                }
                elsif ( ( !$reversed ) && ( !$firstMate ) && ( $readDir ne "rf")  ) {
                    warn("Second mate read mapped in wrong direction, skipping read...\n");
                    next;
                }

                $strand = "+";

            }

            my $qscore = $a->qscore->[ $pileup->qpos ];
            &$debug( "QS: ", $qscore );
            next unless ( defined($qscore) && ( ( $qscore - $qsOffset ) >= $minBaseQ ) );

            $baseStrandCount{$strand}++;
            $baseCount->{$strand}->{$qbase}++;

        }    # end for @$p

        my $refbaseRev = $refbase;
        $refbaseRev    =~ tr/[ACGT]/[TGCA]/;
        my $contextRev = reverse($context);
        $contextRev =~ tr/[ACGTacgt]/[TGCAtgca]/;
        my %strandRefbase = ( '+' => $refbase, '-' => $refbaseRev );
        my %strandContext = ( '+' => $context, '-' => $contextRev );
        foreach my $s ( keys(%baseStrandCount) ) {

            # debug("strand ", $s);
            if ( $baseStrandCount{$s} > 0 ) {
                &$callMethStatus( $baseCount->{$s}, $s, $strandRefbase{$s}, $strandContext{$s}, $seqid, $pos,
                                  $baseStrandCount{$s}, $geneName );
            }
        }
    }
}

sub callMethStatus {
    my $baseCount = shift;
    my $strand    = shift;
    my $refbase   = shift;
    my $context   = shift;
    my $seqid     = shift;
    my $pos       = shift;
    my $total     = shift;
    my $geneName  = shift;

    my $m5Cs__ = \%m5Cs;
    my $m5Cs_  = $m5Cs__->{$strand};

    if ( $refbase eq 'C' ) {
        $m5Cs_->{$seqid}->{refC}   += 1;
        $m5Cs_->{$seqid}->{methCT} += $baseCount->{C} + $baseCount->{T};
        $m5Cs_->{$seqid}->{methC}  += $baseCount->{C};
    }

    my $convRate = 0;
    my $methNr   = 0;
    my $C2T      = 0;
    my $mutNr    = 0;
    my $mutRate  = 0;
    my $methRate = 0;
    my $coverage = 0;
    my $reportMe = 0;

    if ( $total >= $minCov ) {
        my $maxCalledBase = reduce { $baseCount->{$a} > $baseCount->{$b} ? $a : $b } keys %{$baseCount};
        my @equalCalledBases = grep { $baseCount->{$_} eq $baseCount->{$maxCalledBase} } keys %{$baseCount};
        $maxCalledBase = join( '', ( sort(@equalCalledBases) ) );

        if ( $refbase eq 'C' ) {
            $m5Cs_->{$seqid}->{refCanalyzed} += 1;
            $methNr   = $baseCount->{C};
            $C2T      = $baseCount->{T};
            $mutNr    = $total - $methNr - $C2T;
            $mutRate  = $mutNr / $total;
            $coverage = $total - $mutNr;
            $methRate = ( $coverage > 0 ) ? ( $methNr / $coverage ) : 0;
            $convRate = 1 - $methRate;
            $m5Cs_->{$seqid}->{ConvRate} += $convRate;

            $m5Cs_->{$seqid}->{methCTanalyzed} += $C2T;
            $m5Cs_->{$seqid}->{methCanalyzed}  += $methNr;

            if ( ( $mutRate >= $minMutRate ) && ( $methRate < $minMethRate ) && ( $coverage >= $minCov ) ) {
                if ($reportUP) { $reportMe = 1; }
                $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'V';
            }
            else {
                if ( ( $methRate >= $minMethRate ) && ( $coverage >= $minCov ) ) {
                    $reportMe = 1;
                    $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'M';

                    # get statistics
                    my ( $lci, $uci, $rBinomP, $binomP, $score ) = getSignificance( $methRate, $coverage, $methNr );

                    # Wilson score interval
                    @{ $m5Cs_->{$seqid}->{pos}{$pos}{CI95} } = ( $lci, $uci );

                    # bionomial p-val for methylation rate
                    $m5Cs_->{$seqid}->{pos}{$pos}{rpv} = $rBinomP;

                    # binomial p-value Lister
                    $m5Cs_->{$seqid}->{pos}{$pos}{pv} = $binomP;

                    # methylation Score
                    $m5Cs_->{$seqid}->{pos}{$pos}{sc} = $score;
                }
            }
        }
        elsif ( ( $refbase ne 'T' ) && ( $refbase ne 'C' ) ) {    # refbase not C or T
            $methNr   = $baseCount->{C};
            $C2T      = $baseCount->{T};
            $mutNr    = $methNr + $C2T;
            $coverage = $mutNr;                                      # Bases mutated to C/T
            $mutRate  = $mutNr / $total;
            $methRate = ( $mutNr > 0 ) ? ( $methNr / $mutNr ) : 0;

            if ( ( $mutRate >= $minMutRate ) && ( $methRate >= $minMethRate ) && ( $mutNr >= $minCov ) ) {
                $reportMe = 1;
                $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'MV';

                # get statistics
                my ( $lci, $uci, $rBinomP, $binomP, $score ) = getSignificance( $methRate, $coverage, $methNr );

                # Wilson score interval
                @{ $m5Cs_->{$seqid}->{pos}{$pos}{CI95} } = ( $lci, $uci );

                # bionomial p-val for methylation rate
                $m5Cs_->{$seqid}->{pos}{$pos}{rpv} = $rBinomP;

                # binomial p-value Lister
                $m5Cs_->{$seqid}->{pos}{$pos}{pv} = $binomP;

                # methylation Score
                $m5Cs_->{$seqid}->{pos}{$pos}{sc} = $score;

            }
            else {
                if ( ( $mutRate >= $minMutRate ) && ( $methRate < $minMethRate ) && ( $mutNr >= $minCov ) ) {
                    if ($reportUP) { $reportMe = 1; }
                    $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'UV';
                }
            }
        }
        else {    # refbase = T
            $methNr   = $baseCount->{C};
            $methRate = $methNr / $total;
            $mutNr    = ( sum values %{$baseCount} ) - $baseCount->{$refbase};
            $coverage = $total;    # conservative since we can not say how many T mutated to C

            # $mutRate = ($mutNr > 0) ? $baseCount{$refbase} / $mutNr : 0;
            $mutRate = -1;         # we don't really know if T is from reference or if it is a converted C
            if ( ( $methRate >= $minMethRate ) && ( $methNr >= $halfminCov ) ) {
                $reportMe = 1;
                $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'MV';
            }
            elsif ( !grep( $refbase, @equalCalledBases ) ) {
                if ($reportUP) { $reportMe = 1; }
                $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'V';
            }
        }

        if ( $reportMe > 0 ) {
            $m5Cs_->{$seqid}->{pos}{$pos}{geneName} = $geneName;
            $m5Cs_->{$seqid}->{pos}{$pos}{mrate}    = $methRate;
            $m5Cs_->{$seqid}->{pos}{$pos}{mutRate}  = $mutRate;
            $m5Cs_->{$seqid}->{pos}{$pos}{cov}      = $total;
            $m5Cs_->{$seqid}->{pos}{$pos}{methNr}   = $methNr;
            $m5Cs_->{$seqid}->{pos}{$pos}{mutNr}    = $mutNr;
            $m5Cs_->{$seqid}->{pos}{$pos}{refbase}  = $refbase;
            $m5Cs_->{$seqid}->{pos}{$pos}{context}  = $context;
            $m5Cs_->{$seqid}->{pos}{$pos}{mCBase}   = $maxCalledBase;
            $m5Cs_->{$seqid}->{pos}{$pos}{mCBaseC}  = $baseCount->{ substr( $maxCalledBase, 0, 1 ) };
        }
        else {
            delete( $m5Cs_->{$seqid}->{pos}{$pos} );
        }
        $reportMe = 0;
    }

}

sub callMethStatusAZA {
    my $baseCount = shift;
    my $strand    = shift;
    my $refbase   = shift;
    my $context   = shift;
    my $seqid     = shift;
    my $pos       = shift;
    my $total     = shift;
    my $geneName  = shift;

    my $m5Cs__ = \%m5Cs;
    my $m5Cs_  = $m5Cs__->{$strand};

    if ( $refbase eq 'C' ) {
        $m5Cs_->{$seqid}->{refC}   += 1;
        $m5Cs_->{$seqid}->{methCG} += $baseCount->{C} + $baseCount->{G};
        $m5Cs_->{$seqid}->{methC}  += $baseCount->{G};
    }

    my $convRate = 0;
    my $methNr   = 0;
    my $C        = 0;
    my $mutNr    = 0;
    my $mutRate  = 0;
    my $methRate = 0;
    my $coverage = 0;
    my $reportMe = 0;

    if ( $total >= $minCov ) {
        my $maxCalledBase = reduce { $baseCount->{$a} > $baseCount->{$b} ? $a : $b } keys %{$baseCount};
        my @equalCalledBases = grep { $baseCount->{$_} eq $baseCount->{$maxCalledBase} } keys %{$baseCount};
        $maxCalledBase = join( '', ( sort(@equalCalledBases) ) );

        if ( $refbase eq 'C' ) {
            $m5Cs_->{$seqid}->{refCanalyzed} += 1;
            $methNr   = $baseCount->{G};
            $C        = $baseCount->{C};
            $mutNr    = $total - $methNr - $C;
            $mutRate  = $mutNr / $total;
            $coverage = $total - $mutNr;
            $methRate = ( $coverage > 0 ) ? ( $methNr / $coverage ) : 0;
            $convRate = $methRate;
            $m5Cs_->{$seqid}->{ConvRate} += $convRate;

            $m5Cs_->{$seqid}->{methCTanalyzed} += $C;    # TODO: change methCTanalyzed to something meaningful
            $m5Cs_->{$seqid}->{methCanalyzed}  += $methNr;

            if ( ( $mutRate >= $minMutRate ) && ( $methRate < $minMethRate ) && ( $coverage >= $minCov ) ) {
                if ($reportUP) { $reportMe = 1; }
                $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'V';
            }
            else {
                if ( ( $methRate >= $minMethRate ) && ( $coverage >= $minCov ) ) {
                    $reportMe = 1;
                    $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'M';

                    # get statistics
                    my ( $lci, $uci, $rBinomP, $binomP, $score ) = getSignificance( $methRate, $coverage, $methNr );

                    # Wilson score interval
                    @{ $m5Cs_->{$seqid}->{pos}{$pos}{CI95} } = ( $lci, $uci );

                    # bionomial p-val for methylation rate
                    $m5Cs_->{$seqid}->{pos}{$pos}{rpv} = $rBinomP;

                    # binomial p-value Lister
                    $m5Cs_->{$seqid}->{pos}{$pos}{pv} = $binomP;

                    # methylation Score
                    $m5Cs_->{$seqid}->{pos}{$pos}{sc} = $score;

                }
            }
        }
        elsif ( ( $refbase ne 'G' ) && ( $refbase ne 'C' ) ) {    # refbase not C or G
            $methNr   = $baseCount->{G};
            $C        = $baseCount->{C};
            $mutNr    = $methNr + $C;
            $coverage = $mutNr;                                      # Bases mutated to C/G
            $mutRate  = $mutNr / $total;
            $methRate = ( $mutNr > 0 ) ? ( $methNr / $mutNr ) : 0;

            if ( ( $mutRate >= $minMutRate ) && ( $methRate >= $minMethRate ) && ( $mutNr >= $minCov ) ) {
                $reportMe = 1;
                $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'MV';

                # get statistics
                my ( $lci, $uci, $rBinomP, $binomP, $score ) = getSignificance( $methRate, $coverage, $methNr );

                # Wilson score interval
                @{ $m5Cs_->{$seqid}->{pos}{$pos}{CI95} } = ( $lci, $uci );

                # bionomial p-val for methylation rate
                $m5Cs_->{$seqid}->{pos}{$pos}{rpv} = $rBinomP;

                # binomial p-value Lister
                $m5Cs_->{$seqid}->{pos}{$pos}{pv} = $binomP;

                # methylation Score
                $m5Cs_->{$seqid}->{pos}{$pos}{sc} = $score;

            }
            else {
                if ( ( $mutRate >= $minMutRate ) && ( $methRate < $minMethRate ) && ( $mutNr >= $minCov ) ) {
                    if ($reportUP) { $reportMe = 1; }
                    $m5Cs_->{$seqid}->{pos}{$pos}{sate} = 'UV';
                }
            }
        }

        if ( $reportMe > 0 ) {
            $m5Cs_->{$seqid}->{pos}{$pos}{geneName} = $geneName;
            $m5Cs_->{$seqid}->{pos}{$pos}{mrate}    = $methRate;
            $m5Cs_->{$seqid}->{pos}{$pos}{mutRate}  = $mutRate;
            $m5Cs_->{$seqid}->{pos}{$pos}{cov}      = $total;
            $m5Cs_->{$seqid}->{pos}{$pos}{methNr}   = $methNr;
            $m5Cs_->{$seqid}->{pos}{$pos}{mutNr}    = $mutNr;
            $m5Cs_->{$seqid}->{pos}{$pos}{refbase}  = $refbase;
            $m5Cs_->{$seqid}->{pos}{$pos}{context}  = $context;
            $m5Cs_->{$seqid}->{pos}{$pos}{mCBase}   = $maxCalledBase;
            $m5Cs_->{$seqid}->{pos}{$pos}{mCBaseC}  = $baseCount->{ substr( $maxCalledBase, 0, 1 ) };
        }
        else {
            delete( $m5Cs_->{$seqid}->{pos}{$pos} );
        }
        $reportMe = 0;
    }

}

sub getSignificance {
    my ( $methRate, $coverage, $methNr ) = @_;
    
    # shortcut, do not caclulate stats when estimating the C->T conversion rate
    return (0, 0, 0, 0, 0) if ( $calculateConvR );

    # Wilson score interval
    my ( $lci, $uci ) = calculateCi95( $methRate, $coverage );

    # caclulate p-values using a Fisher's exact test when no C->T conversion rate is known
    # for this we compare the methylation rate with the max baseline sequencing error provided
    # by the minimum Q-value (default Q30)
    if ( $conversionRate == -1 ) {
        my $fisherP = 1;
        
        # BASELINE_ERR is the probability that the basecalls are not correct
        # = 10^(-Q/10) but we assume to only expect 25% of the misscalls to
        # be a C or G in aza mode
        my $expMeth = int( $coverage * $BASELINE_ERR * 0.25 );
        my $expConv = $coverage - $expMeth;

        # Setup contingency table for Fisher's exact test
        #         m       n
        #   A   [0][0]  [0][1]
        #   B   [1][0]  [1][1]

        my @t;
        $t[0][0] = $methNr;
        $t[0][1] = $coverage - $methNr;

        $t[1][0] = $expMeth;
        $t[1][1] = $expConv;
        
        $fisherP = fisherExactTest([@t]);
        
        return ($lci, $uci, undef, $fisherP, ($methNr * $lci * 10));

    }

    # Barturen et al. 2013
    # maximum false C calls at the given error interval & conversion rate
    my $mFCs = maxFalseCs( $methRate, $coverage, $methNr, $errorInterval );
    my $rBinomP = 1 - pbinom( $mFCs, $methNr, ( 1 - $conversionRate ) );
    $rBinomP = ( $rBinomP < 0 ) ? 0 : $rBinomP;    # return

    if ( $fisherTest == 1 ) {
        my $fisherP = 1;
        
        my $expMeth = int( $coverage * (1 - $conversionRate) );
        my $expConv = $coverage - $expMeth;

        # Setup contingency table for Fisher's exact test
        #         m       n
        #   A   [0][0]  [0][1]
        #   B   [1][0]  [1][1]

        my @t;
        $t[0][0] = $methNr;
        $t[0][1] = $coverage - $methNr;

        $t[1][0] = $expMeth;
        $t[1][1] = $expConv;

        $fisherP = fisherExactTest([@t]);

        my $fisherP_ = ( $fisherP == 0 ) ? 1e-20 : $fisherP;
        my $score = $methNr * $lci * -log10($fisherP_);

        return ( $lci, $uci, $rBinomP, $fisherP, $score );

    }
    else {

        # Lister et al. 2009
        my $binomP = 1 - pbinom( $methNr, $coverage, ( 1 - $conversionRate ) );

        my $binomP_ = ( $binomP == 0 ) ? 1e-20 : $binomP;
        my $score = $methNr * $lci * -log10($binomP_);

        return ( $lci, $uci, $rBinomP, $binomP, $score );
    }
}

sub log10 {
    return ( log( $_[0] ) / LOG10 );
}

sub calculateCi95 {
    my ( $p, $n ) = @_;    # ($methRate, $n)

    # Wilson score interval
    my $ic = $p + Z95SQ / ( 2 * $n );
    my $sd  = Z95 * sqrt( ( $p * ( 1 - $p ) ) / $n + Z95SQ / ( 4 * $n**2 ) );
    my $nom = 1 + Z95SQ / $n;
    my $uCI = ( $ic + $sd ) / $nom;
    my $lCI = ( $ic - $sd ) / $nom;

    return ( $lCI, $uCI );
}

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

    # TODO: Check this! Won't work when packed into a standalone binary via pp
    #if ( (my $errorCode = getErrorCode()) ) {
    #    say STDERR $errorCode . " - " . getErrorMessage();
    #    $pValue = 1;
    #}
    
    return ($pValue);

}

sub calculatePV {          # NOT used
    my ( $p, $n, $nullHyp ) = @_;    # ($methRate, $n)

    my $u = ( $p - $nullHyp ) / sqrt( $nullHyp * ( 1 - $nullHyp ) / $n );
    my $pval = pnorm( -abs($u) );

    return ($pval);
}

sub binomCoeff {          # NOT used
    my ( $n, $k ) = @_;

    my ( $u, $l ) = ( 0, 0 );

    if ( $k > $n ) {
        return (-1);
    }
    elsif ( $k == 0 ) {
        return (0);
    }
    elsif ( $k < ( $n - $k ) ) {
        my $i = $n;
        $u += log($i), $i-- while $i > $n - $k;
        $l += log($k), $k-- while $k > 0;
    }
    else {
        my $i = 1;
        $u += log($n), $n-- while $n > $k;
        $l += log($i), $i++ while $i > $n - $k;
    }

    return ( $u - $l );

}

# maxFalseCs($methRate, $coverage, $methNr, $errorInterval)
sub maxFalseCs {
    my $methRate      = shift;
    my $coverage      = shift;
    my $methNr        = shift;
    my $errorInterval = shift // 0.2;

    my $mFCs = 0;
    for ( my $i = 0 ; $i <= $methNr ; $i++ ) {
        my $corr = ( $methNr - $i ) / $coverage;
        if ( ( $methRate - $corr ) <= $errorInterval ) {
            $mFCs = $i;
        }
        else {
            return $mFCs;
        }
    }

    return $mFCs;
}

sub getRefSeq {
    my $sam        = shift;
    my $seqid      = shift;
    my $pos        = shift;
    my $seqContext = shift;

    my ( $refbase, $context );

    if ( !$seqContext ) {
        $refbase = uc( $sam->segment( $seqid, $pos, $pos )->dna );
        return ( $refbase, "" );
    }

    my ( $seqContextL, $seqContextR ) = ( abs($seqContext), abs($seqContext) );
    my $sLength = $sam->length($seqid);
    
    if ( ( $pos > $sLength ) || ( $pos < 1 ) ) {
        warn("Base Position " . $pos . " out of reference sequence length: " . $sLength);
        return ( "N", ("N")x(2*$seqContextR) );
    }
    
    my $start   = $pos - $seqContextL;
    my $end     = $pos + $seqContextR;

    ( $start, $seqContextL ) = ( $start > 0 )        ? ( $start, $seqContextL ) : ( 1, 0 );
    ( $end,   $seqContextR ) = ( $end > $sLength )   ? ( $sLength, $end - $sLength ) : ( $end, $seqContextR );

    $context = $sam->segment( $seqid, $start, $end )->dna;
    $refbase = uc( substr( $context, $seqContextL, 1 ) );
    $context =~ s/(.{$seqContextL})(.)(.{$seqContextR})/\L$1\U$2\L$3/;

    return ( $refbase, $context );
}

sub resultCollector {
    my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $m5CsAllData ) = @_;

    if ( defined($m5CsAllData) ) {

        $totalFilteredDups += scalar( keys %{ $m5CsAllData->{dup} } ) if ( $maxDuplicates > 0 );

        foreach my $s (@strands) {    # strands $s
            my %m5CsAll = %{ $m5CsAllData->{$s} };
            my $seqid   = $ident;
            my $m5Cdata = $m5CsAll{$seqid};

            $totalRefCs         += ( defined( $m5Cdata->{refC} ) )         ? $m5Cdata->{refC}         : 0;
            $totalRefCsanalyzed += ( defined( $m5Cdata->{refCanalyzed} ) ) ? $m5Cdata->{refCanalyzed} : 0;

            if ($azaMode) {
                $methCG += ( defined( $m5Cdata->{methCG} ) ) ? $m5Cdata->{methCG} : 0;
            }
            else {
                $methCT += ( defined( $m5Cdata->{methCT} ) ) ? $m5Cdata->{methCT} : 0;
            }
            $methC  += ( defined( $m5Cdata->{methC} ) )  ? $m5Cdata->{methC}  : 0;

            $methCTanalyzed += ( defined( $m5Cdata->{methCTanalyzed} ) ) ? $m5Cdata->{methCTanalyzed} : 0;
            $methCanalyzed  += ( defined( $m5Cdata->{methCanalyzed} ) )  ? $m5Cdata->{methCanalyzed}  : 0;

            my $m5CseqData = $m5Cdata->{pos};
            if ( scalar keys( %{$m5CseqData} ) > 0 ) {

                foreach my $pos ( keys( %{$m5CseqData} ) ) {
                    my $m5CposData = $m5CseqData->{$pos};
                    my $mrate      = sprintf( "%01.3f", $m5CposData->{mrate} );
                    my $mutRate    = sprintf( "%01.3f", $m5CposData->{mutRate} );
                    my ( $lCI, $uCI ) =
                        ( defined( $m5CposData->{CI95} ) )
                      ? ( sprintf( "%01.3f", $m5CposData->{CI95}->[0] ), sprintf( "%01.3f", $m5CposData->{CI95}->[1] ) )
                      : ( "NaN", "NaN" );
                    my $pvalNice =
                      ( defined( $m5CposData->{pv} ) )
                      ? sprintf( "%01.6e", $m5CposData->{pv} )
                      : "NaN";
                    my $rpvalNice =
                      ( defined( $m5CposData->{rpv} ) )
                      ? sprintf( "%01.6e", $m5CposData->{rpv} )
                      : "NaN";
                    my $scoreNice =
                      ( defined( $m5CposData->{sc} ) )
                      ? sprintf( "%01.5f", $m5CposData->{sc} )
                      : "NaN";
                    my $context =
                      ( defined( $m5CposData->{context} ) )
                      ? $m5CposData->{context}
                      : $m5CposData->{refbase};
                    my $geneName = $m5CposData->{geneName};

                    my $candidateName = "m5C_" . ++$m5Cnr;
                    my $resultLine = join(
                                           "\t",
                                           (
                                              $seqid,                 $pos,                   $s,
                                              $m5CposData->{refbase}, $m5CposData->{cov},     $m5CposData->{methNr},
                                              $mrate,                 $m5CposData->{mutNr},   $mutRate,
                                              $m5CposData->{mCBase},  $m5CposData->{mCBaseC}, $m5CposData->{sate},
                                              $lCI,                   $uCI,                   $pvalNice,
                                              $rpvalNice,             $scoreNice,             $context,
                                              $geneName,              $candidateName
                                             )
                                             ) . "\n";

                    &resultWirter($resultLine) unless ($calculateConvR);

                    $totalMethRefCs += 1 if ( $m5CposData->{sate} eq 'M' );
                    $totalMutMethCs += 1 if ( $m5CposData->{sate} eq 'MV' );

                    if ( defined( $m5CposData->{pv} ) && $FDR && ( !$FDRrate ) ) {
                        $fdrData[$resNr] = [ $m5CposData->{pv}, $resultLine ];
                        $resNr++;
                    }
                    if ( defined( $m5CposData->{rpv} ) && $FDR && ($FDRrate) ) {
                        $fdrData[$resNr] = [ $m5CposData->{rpv}, $resultLine ];
                        $resNr++;
                    }

                    # pass data to bed writer
                    if ($narrowPeak) {
                        &$bedwriter( $seqid, $pos, $candidateName, $s, $mrate, $m5CposData->{pv}, $m5CposData->{cov},
                                     $m5CposData->{sc}, 64, $bednpFH );
                    }
                    if ($bed63) {
                        &$bedwriter( $seqid, $pos, $candidateName, $s, $mrate, $m5CposData->{pv}, $m5CposData->{cov},
                                     $m5CposData->{sc}, 63, $bed63FH );
                    }
                }
            }
        }
    }
    else {
        print STDERR "No result received from child process" . $pid . "!\n";
    }

    $seqpctDone = $end * 100 / $tlength;
    printf( "processing %i sequences: \[%s - %02.2f%%\] [overall - %02.2f%%] done ...\r",
            $nrOftargets, $wSid, $seqpctDone, $totpctDone ) unless ( $bedFileIN );

}

sub tellStatus {
    printf( "processing %i sequences: \[%s - %02.2f%%\] [overall - %02.2f%%] done ...\r",
            $nrOftargets, $wSid, $seqpctDone, $totpctDone );
}

sub lockF {
    my ($fh) = @_;
    flock( $fh, LOCK_EX ) or die "Cannot lock resultfile - $!\n";
    seek( $fh, 0, SEEK_END ) or die "Cannot seek - $!\n";
}

sub unlockF {
    my ($fh) = @_;
    flock( $fh, LOCK_UN ) or die "Cannot unlock resultfile - $!\n";
}

sub resultWirter {
    my $resultLine = shift;

    lockF($resFH);
    $resFH->print($resultLine);
    $resFH->flush();
    unlockF($resFH);
}

sub bedWriter {
    my ( $chrom, $m5Cpos, $candidateName, $strand, $mrate, $pv, $cov, $score, $bedType, $bedFH ) = @_;

    my $chromStart = $m5Cpos - 1;                    # BED format is 0 based!
    my $chromEnd   = $chromStart + 1;
    my $pVal_      = ( !defined($pv) ) ? -1 : $pv;
    my $pVal =
      ( $bedType != 63 ) ? sprintf( "%01.6e", $pVal_ )
      : (
          ( $pVal_ <= 0 ) ? -1
          : sprintf( "%01.6e", -log10($pVal_) )
          );

    $score =
        ( !defined($score) ) ? 0
      : ( ( 10 * $score ) > 1000 ) ? 1000
      :                              int( 10 * $score );

    my $bedLine;
    if ( $bedType != 63 ) {
        $bedLine =
          join( "\t", ( $chrom, $chromStart, $chromEnd, $candidateName, $score, $strand, $mrate, $pVal, $cov ) ) . "\n";
    }
    else {
        $bedLine = join(
                         "\t",
                         (
                            $chrom, $chromStart, $chromEnd, $candidateName, $score, $strand, ( $mrate * 100 ),
                            $pVal, -1, 0
                           )
                           ) . "\n";
    }

    lockF($bedFH);
    $bedFH->print($bedLine);
    $bedFH->flush();
    unlockF($bedFH);

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

sub sam_to_bam {
    my $sam = shift;

    my ( $fname, $fpath, $fsuffix ) = fileparse( $sam, ".sam" );
    my $bam = $fpath . "/" . $fname . ".bam";

    print STDERR "Converting SAM to BAM...\n";

    my $tam = Bio::DB::Tam->open($sam) || die("Could not open SAM file for reading: $!");
    my $header = $tam->header_read();

    my $out = Bio::DB::Bam->open( $bam, 'w' ) || die("Could not open BAM file for writing: $!");
    $out->header_write($header);

    my $alignment = Bio::DB::Bam::Alignment->new();
    my $lines     = 0;

    while ( $tam->read1( $header, $alignment ) > 0 ) {
        $out->write1($alignment);
        print STDERR "converted " . $lines . " lines...\n" if ( ++$lines % 100000 == 0 );
    }
    undef $tam;
    undef $out;

    print STDERR "converted " . $lines . " lines\n";

    my $sorted_bam = sort_bam($bam);
    index_bam($sorted_bam);

    return $sorted_bam;
}

sub sort_bam {
    my $bam = shift;

    my ( $fname, $fpath, $fsuffix ) = fileparse( $bam, ".bam" );

    print STDERR "sorting " . $fname . ".bam...\n";

    my $sorted = $fpath . "/" . $fname . "_sorted";

    Bio::DB::Bam->sort_core( 0, $bam, $sorted, 1000000000 );
    unlink($bam);    # we don't keep the unsorted bam file
    return ( $sorted . ".bam" );
}

sub index_bam {
    my $bam = shift;
    print STDERR "Indexing BAM...\n";
    Bio::DB::Bam->index_build($bam);
    return;
}

sub getRefLengths {
    my $faidxFile = $_[0];
    
    my $faiFH = IO::File->new( $faidxFile, O_RDONLY ) 
      || die("ERROR: could not read the fasta index file " . $faidxFile . " : check permissions");

    my $refLen;
    while ( my $l = $faiFH->getline() ) {
        my ($seqID, $len) = split("\t", $l, 3);
        $refLen->{$seqID} = $len;
    }
    return ($refLen);
}

sub readBED {
    my $bedFile   = $_[0];
    my $bedData   = $_[1];
    my $recCount  = $_[2];
    my $chrPrefix = $_[3] // "";
    
    my $bedFH = IO::File->new( $bedFile, O_RDONLY | O_EXCL ) || die( $bedFile . ": " . $! );

    while ( defined( my $bedLine = $bedFH->getline() ) ) {
        chomp($bedLine);

        if ( $bedLine =~ /^#/ ) {
            next;
        }

        my ( $chrom, $chromStart, $chromEnd, $name) = split( '\t', $bedLine );
        if ($chrPrefix) {
            $chrom = $chrPrefix . $chrom;
        }
        if ( ($chromStart !~ /\d/) || ($chromEnd !~ /\d/) || ($chromEnd < $chromStart) ) {
            say STDERR "This does not look like a BED file: " . $bedFile;
            exit(1);
        }

        push( @{ $bedData->{$chrom}->{range} }, [ ( $chromStart ), $chromEnd ] );
        push( @{ $bedData->{$chrom}->{line} },   $bedLine );
        
        $recCount++;
    }
    $bedFH->close();
    undef($bedFH);
    
    $_[2] = $recCount;
}

sub checkWirteToDir {
    my $dir = shift;

    if ( -d $dir ) {
        my $test_tmpF = $dir . "/_DirTest_" . $$;
        open( TEST, ">$test_tmpF" ) || ( return (0) );
        close(TEST);
        unlink($test_tmpF);
    }
    return (1);
}

sub checkDir {
    my $dir = shift;

    if ( -d $dir ) {
        my $test_tmpF = $dir . "/_DirTest_" . $$;
        open( TEST, ">$test_tmpF" ) || ( warn( $dir . " exists but is not writable. Exiting" ) and exit(1) );
        close(TEST);
        unlink($test_tmpF);
    }
    elsif ( -e $dir ) {
        warn( $dir . " exists but is not a directory. Exiting" );
        exit(1);
    }
    else {
        mkdir( $dir, 0755 ) || die( $dir . ": " . $! );
    }
}

sub nicePath {
    my $rawPath = shift;
    
    $rawPath =~ s/\/+/\//g;
    return($rawPath);
}

sub usage {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self [options] [-h] [-man] [--version]

Required options any of:
    -fasta|-f              : Reference sequence FASTA file.

    -sam|-bam|-s           : Sequence read alignment file in SAM or BAM format.

    -result|-o             : Result file where to store the metylation calls.

    -genomeDBref|-gref     : SAM/BAM file was genereated by aligning bs-reads to
                             a genome reference (DNA database): e.g. mouse
                             genome using meRanG.
                             If set, a BED6 + 3 file will be created in
                             addition to the standard result file.
                             (default: not set)

    -transcriptDBref|-tref : SAM/BAM file was genereated by aligning bs-reads to
                             a transcript reference (Transcript database): e.g. mouse
                             refSeqRNA using meRanT.
                             (default: not set)


Options:

    -procs|-p              : Number or processors (CPUs) to use in parallel.
                             Setting this option significantly reduces the
                             processing time. E.g. when set to "-p 16" 16 sequences
                             (e.g. chromosomes) will be processed in parallel.
                             (default: 1)

    -regions|-bi           : BED file with regions to scan for m5Cs.
                             If specified meRanCall will only call m5Cs in the regions
                             present in the BED file.
                             (default: not set, scan entire SAM/BAM file)

    -fskip5|-fs5           : number or bases to ignore on the 5' end of a forward
                             read. This helps to aviod biased results. See m-bias
                             plot from meRanGs or meRanT output to get
                             an estimate for this number.
                             (defualt: 0)

    -fskip3|-fs3           : number or bases to ignore on the 3' end of a forward
                             read. This helps to aviod biased results. See m-bias
                             plot from meRanG or meRanT output to get
                             an estimate for this number.
                             (defualt: 0)

    -rskip5|-rs5           : number or bases to ignore on the 5' end of a reverse
                             read. This helps to aviod biased results. See m-bias
                             plot from meRanG or meRanT output to get
                              an estimate for this number.
                             (defualt: 0)

    -rskip3|-rs3           : number or bases to ignore on the 3' end of a reverse
                             read. This helps to aviod biased results. See m-bias
                             plot from meRanG or meRanT output to get
                             an estimate for this number.
                             (defualt: 0)

    -readLength|-rl        : If set to the original read length, then the 3' end
                             skipping will be adjusted for 3' trimming. In other
                             words: if you trimmed some of your reads before
                             mapping, than the number of trimmed bases on the 3'
                             end will be treated as already skipped.
                             This has no effect if fskip5, fskip3, rskip5 or rskip3
                             is 0.
                             (default: 100)

    -readDir|-rd           : read direction for single end reads: may be fr or rf.
                             fr: standard Illumina forward reads
                             rf: reverse reads
                             (default: fr)

    -minMethR|-mr          : Minimum methylation ratio of a single C, that is needed
                             to consider this C as potentially methylated
                             (default: 0.2)

    -minMutR|-mutR         : Minimum ratio (bases on reads at a given reference
                             position different from reference base) above which
                             a base will be considered as mutated in respect to
                             the base on the reference sequence.
                             (default: 0.8)

    -minBaseQ|-mBQ         : Minimum read base quality (phred score) to condsider
                             for methylation calling.
                             (default: 30)

    -minCov|-mcov          : Minimum coverage at a given reference position above
                             which methylation calling will be performed.
                             (default: 10)

    -maxDup|-md            : Maximum number of read duplicates covering a given
                             position. Read duplicates have the same start positon
                             on the reference and map to the same sequence.
                             (default: 0, do not filter duplicates)

    -fisherTest|-fet       : Use Fisher exact Test to statistically test the methylation
                             state against background non conversion levels at the given
                             coverage.
                             (default: not set, default in future)

    -conversionRate|-cr    : C->T Conversion rate (0 < cr < 1)
                             (default: 1)

    -errorInterval|-ei     : Error interval for methylation rate p-value calculation
                             (default: 0)

    -fdr                   : Control the false discovery rate of methylated cytosines
                             at the specified FDR (0 < fdr < 1).
                             (default: not set)
                            
    -fdrRate               : Use the probability that the real methylation level or rate
                             instead of the methylation state p-value to control the
                             false discovery rate at -fdr FDR (0 < fdr < 1).
                             (default: not set)
                            
    -calcConvRate|-ccr     : Caluclate the C->T conversion rate from an unmehtylated
                             control sequence.
                             (default: not set)

    -controlSeqID|-cSeqID  : Control sequence ID(s) for C->T conversion rate caluclation
                             Can be specified multiple times for multiple control
                             sequences.
                             (default: not set)

    -excludeSeqID|-exSeqID : Sequence ID(s) to exclude from methylation calling.
                             Can be specified multiple times for multiple control
                             sequences. E.g. -exSeqID chr1 -exSeqID chrUn_gl000220
                             (default: not set)

    -reportUP|-rUP         : report unmethylated mutated bases?
                            (default: not set)
                            
    -bed63                 : Generate a BED6 + 3 file - only relevant for genome mapped
                             data!
                            (default: not set)

    -narrowPeak|-np        : Generate a narrowPeak BED file - only relevant for genome
                             mapped data!
                            (default: not set)
                            
    -seqContext|-sc        : If set to a number, this number of bases 5' and 3' of
                             the methylated C will be displayed in the result file.
                             (default: not set)

    -havZG|-zg             : If set, the methylation caller will look for the "ZG"
                             custom SAM tag and use it a gene name associated with
                             the methylated positon in the result file.
                           
                             meRanT adds this tag to the SAM entries.
                           
                             meRanG does not, however you can use the BED6 + 3
                             file and run the "meRanAnnotate" tool from meRanTK to
                             associate methylated C's  with gene (transcript) names.
                           
                             (default: not set)

    -azaMode|-aza          : If set, the methylation caller will run in the Aza-IP mode
                             and enables methylation calling from Aza-IP data by looking
                             for C->G conversions, which are characteristic for Aza-IP data.
                           
                             (default: not set)

    --version              : Print the program version and exit.
    -h|help                : Print the program help information.
    -man                   : Print a detailed documentation.
    
    -debug|-d              : Print some debugging information.
    
EOF
}

__END__

=head1 NAME

meRanCall - RNA cytosine methylation caller

=head1 SYNOPSIS

=head2 methylation calling from RNA-BSseq short reads mapped to a transcriptome

=over 2
 
### Single End reads

 meRanCall \
 -p 32 \
 -o ./meRanCallResult.txt \
 -bam ./meRanT_mapping_sorted.bam \
 -f ./mm10.refSeqRNA.fa \
 -fs5 6 \
 -rl 51 \
 -sc 5 \
 -zg \
 -ei 0.1 \
 -cr 0.99 \
 -tref

 The command above calls methylated C's from reads in "meRanT_mapping_sorted.bam" mapped to
 the reference transcriptome database (transcripts) in “mm10.refSeqRNA.fa”. The methylation
 calling process will use (-p) 32 CPUs in parallel.
 The mapping was created using meRanT, therefore the "-zg" option is added in the example. This
 way, the gene names associated with the individual transcripts will be extracted from the BAM
 file and reported in the methylation calling result file. Since the type of reference for the
 alignment was a transcript database the "-tref" option has to be used.
 Let’s assume that the reads are C biased at the first 6 base positions on the 5' end (this could
 be estimated from the m-Bias plot produced by meRanT). We tell the meRanCall program to ignore
 these biased positions by specifying the option "-fs5 6" (forward read skip on 5’ end 6 bases).
 The original reads had a length of 100 base pairs before any trimming, we tell this by setting
 “–rl 100”, this way meRanCall ignores only up to 6 bases from the 5’ end of reads that were
 longer than 93 bps after read trimming (that you did in your QC before aligning the reads).
 We want also to get the sequence context 10 bps around the methylated C's (-sc 10). We set the
 error interval (-ei 0.1) for calculating the methylation level p-value to 0.1, that means that we
 calculate the probability that the real methylation level lies within that interval
 (Barturen et al., 2013). Our C→T conversion rate (-cr 0.99) is 0.99 as we determined from a
 un-methylated in-vitro transcribed control RNA that was spiked into our sample.

 The result meRanCallResult.txt is a tab separated file and contains the following
 data-fields for each potentially methylated C:
 
 1.  SeqID          : sequence ID from reference database
 2.  refPos         : postion of the methylated C on the reference sequence
 3.  refStrand      : strand (will always be '+' when using a reference transcriptome)
 4.  refBase        : base on the reference sequence
 5.  cov            : coverage (# of reads covering this position)
 6.  C_count        : # of C's counted at this position
 7.  methRate       : methylation rate
 8.  mut_count      : # of non-reference bases at the position
 9.  mutRate        : mutation rate (#non reference bases / coverage)
 10. CalledBase     : prevailing base(s) at the position
 11. CB_count       : CalledBase count
 12. state          : methylation status (M|MV|UV|V)
                      M : methylated C, C on reference
                      MV: methylated C, NO C on reference (mutated)
                      UV: unmethylated C, NO C on reference (mutated)
                      V : mutated base
 13. 95_CI_lower    : lower bound of the 95% confidence interval (Wilson score interval)
 14. 95_CI_upper    : upper bound of the 95% confidence interval (Wilson score interval)
 15. p-value_mState : p-value of the methylation State (Lister et al. 2009)
 16. p-value_mRate  : p-value of the methylation Rate (Barturen et al. 2013)
 17. Score          : methylation call score
 18. seqContext     : sequence Context arround the mehtylated C
 19. geneName       : gene name associated with the methylated C
 20. candidateName  : name assigned to the methylated C candidate
 
 The methylation calling process will use (-p) 32 cpus in parallel. 

  
### Paired End reads

 For paired end reads a command with analogous options as for single ends can be
 used. In addition, if you have 3’ or 5’ “C” biased reverse read ends that you want
 to be ignored by meRanCall, you can specify this for the reverse reads using the
 “-rsikp5” and/or “-rsikp3” option.

=back

=head2 methylation calling from RNA-BSseq short reads mapped to a genome

=over 2
 
### Single End reads

 meRanCall \
 -p 23 \
 -o ./meRanCallResultGenome.txt \
 -bam ./meRanG_mapping_sorted.bam \
 -f ./mm10.allchr.fa \
 -fs5 6 \
 -rl 51 \
 -sc 5 \
 -gref \
 -ei 0.1 \
 -cr 0.99

 The command above calls methylated C's from reads in "meRanG_mapping_sorted.bam"
 mapped to a reference genome database (DNA), e.g.: mouse refseq databases
 in "mm10.allchr.fa".
 The type of reference for the alignment (performed by meRanG) was a genome
 database therefore the "-gref" option has to be used.
 The reads are C biased at the first 6 base position at the 5' end (this
 could be estimated from the m-bias plot produced by meRanG), therefore
 in the example the option "-fs5 6" is set. The original reads head a length of 51
 base pairs before trimming. We want also to get the sequence context 5 bps arround the
 methylated C's. We set the error interval for calculating the methylation level p-value to
 0.1, that means that we calcultae the probability that the real methylation level lies
 within that interval. Our C->t conversion rate is 0.99 as we determined from a unmethylated
 invitro transcrbed control RNA that was spiked in.
 
 The methylation calling process will use (-p) 23 cpus in parallel.

 The result file meRanCallResultGenome.txt is a tab separated and contains the same
 datafields as described above, with the exception, that the 'geneName' field (19)
 will be empty '-', since the meRanG mappers will not add the "ZG" custom tag to
 the SAM lines. 

 However a BED6 + 3 file 'meRanCallResultGenome.bed' will be created in addition to the
 'meRanCallResultGenome.txt' result file. By using the
 
               "meRanAnnotate"

 tool from meRanTK you can associate methylated C's in these files (txt/BED) with gene
 (transcript) names. Or use the bed file with any other tool that can intersect BED with
 annotaion files.
 
 
 
### Paired End reads

 For paired end reads a command with analogous options can be used and in addition
 the biased end ignoring can be specified for the reverse reads using: '-rsikp5 -rsikp3'

=back

### Region based methylation calls

 If one is only interested in methylation calls for specific regions, one can use the
 "-regions" option and supply a BED file with the regions of interest. meRanCall will
 then only call methylated C's in these regions.

=back

=head1 DEPENDENCIES

 The program requires additional perl modules and depending on your perl installation
 you might need do add the following modules:

 Math::CDF
 Bio::DB::Sam
 Parallel::ForkManager
 
 These modules should be availble via CPAN or depending on your OS via the package
 manager.
 
 Bio::DB:Sam requires the samtools libraries (version 0.1.10 or higher, version 1.0
 is not compatible with Bio::DB::Sam, yet) and header files in order to compile
 sucessfully.

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
