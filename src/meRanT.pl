#!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.24.0/bin/perl -w
#
#  meRanT.pl
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
use Config;
$Config{useithreads} or die('Recompile Perl with threads to run this program.');
use IO::File;
use Fcntl qw(SEEK_END SEEK_SET SEEK_CUR);
use File::Basename;
use Cwd qw(abs_path);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use POSIX qw(mkfifo);
use threads;
use threads::shared;
use Thread::Queue;
use List::Util qw(max min sum);

use Parallel::ForkManager;
use Bio::DB::Sam;

my $DEBUG = 0;

my $VERSION = '1.2.1b';

my $version;
my $help;
my $man;
my $self               = basename( $0, () );
my $instDir            = dirname( abs_path($0) );
my $arch               = $Config{archname};
my $extUtilDir         = $instDir . '/extutil/' . $arch;
my $useShippedExtUtils = ( -x $extUtilDir ) ? 1 : 0;

my $max_threads = 1;
my $fastQfwd    = "";
my $fastQrev    = "";
my $firstNreads = -1;
my $mbQSt       = 30;

my $fixDirectionaly     = 0;
my $forceDirectionality = 0;
my $useIlluminaQC       = 0;

my $bowtie2_cmd   = ($useShippedExtUtils) ? ( $extUtilDir . '/bowtie2/bowtie2' )       : 'bowtie2';
my $bowtie2_build = ($useShippedExtUtils) ? ( $extUtilDir . '/bowtie2/bowtie2-build' ) : 'bowtie2-build';
my $bowtie2_idx = $ENV{BS_BWT2IDX} // "";
my $bowtie2_N   = 0;
my $bowtie2_L   = 20;
my $bowtie2_D   = 30;
my $bowtie2_R   = 2;
my $bowtie2_I   = 0;
my $bowtie2_X   = 1000;
my $bowtie2_min_score      = "";            # defaults to 'G,20,8' in local mode and to 'L,-0.4,-0.4' in end-to-end mode
my $bowtie2_k              = 10;
my $bowtie2_un             = "";
my $reportAmbiguosMultiMap = "";
my $bowtie2_mode           = "end-to-end";

my @invariableBowtie2Settings = qw(--reorder -q);

my $fasta    = "";
my $bsIdxDir = "";

my $max_edit_dist = 2;
my $max_mm_ratio  = 0.05;
my $samout        = "";
my $aMM_fname     = "";     # TXT filename for ambiguos alignments
my $samoutMM      = "";     # SAM filename for ambiguos alignments
my $ommitBAM      = 0;
my $deleteSAM     = 0;
my $saveMultiMapperSeparately;

my $id2genMapFile = "";

my @fqFiles;
my %fwdFL;
my %revFL;
my $qsOffset;

my %fqHash;
my %idMap;

my $mBiasPlot     = 0;
my $mBiasPlotOutF = "";
my $mBiasPlotOutR = "";
my %mBiasData;
my $mbq;
my $mBiasThr;

my $fMO;
my $hardClipMateOverlap;

my %mappingStatsHash = (
                         'reads'                   => 0,
                         'uniqGeneMultiTranscript' => 0,
                         'uniqGeneUniqTranscript'  => 0,
                         'mtpMapperFilter'         => 0,
                         'mgMapperFilter'          => 0,
                         'flagFilter'              => 0,
                         'editDistFilter'          => 0,
                         'bowtie2unal'             => 0,
                         'filteredReads'           => undef,
                         'editDistFilter_ID'       => undef,
                         );

my $mappingStats = \%mappingStatsHash;

our $outDir   = "";
our $unoutDir = "";

our @suffixlist = qw(.fastq .fq .fastq.gz .fq.gz .fastq.gzip .fq.gzip);

my @SAMhdr_cmd_line = @ARGV;
GetOptions(
            'fastqF|f=s'           => \$fastQfwd,
            'fastqR|r=s'           => \$fastQrev,
            'illuminaQC|iqc'       => \$useIlluminaQC,
            'forceDir|fDir'        => \$forceDirectionality,
            'fasta|fa=s'           => \$fasta,
            'first|fn=i'           => \$firstNreads,
            'QSoffset|qso=i'       => \$qsOffset,
            'outdir|o=s'           => \$outDir,
            'sam|S=s'              => \$samout,
            'threads|t=i'          => \$max_threads,
            'bowtie2cmd|bwt2=s'    => \$bowtie2_cmd,
            'bowtie2build|bwt2b=s' => \$bowtie2_build,
            'bsidxdir|id=s'        => \$bsIdxDir,
            'bsidx|x=s'            => \$bowtie2_idx,
            'bowtie2N|N=i'         => \$bowtie2_N,
            'bowtie2L|L=i'         => \$bowtie2_L,
            'bowtie2D|D=i'         => \$bowtie2_D,
            'bowtie2R|R=i'         => \$bowtie2_R,
            'bowtie2I|I=i'         => \$bowtie2_I,
            'bowtie2X|X=i'         => \$bowtie2_X,
            'min-score=s'          => \$bowtie2_min_score,
            'bowtie2k|k=i'         => \$bowtie2_k,
            'bowtie2un|un'         => \$bowtie2_un,
            'unalDir|ud=s'         => \$unoutDir,
            'samMM|MM'             => \$saveMultiMapperSeparately,
            'ommitBAM|ob'          => \$ommitBAM,
            'deleteSAM|ds'         => \$deleteSAM,
            'reportAM|ra'          => \$reportAmbiguosMultiMap,
            'bowtie2mode|m=s'      => \$bowtie2_mode,
            'max-edit-dist|e=i'    => \$max_edit_dist,
            'max-mm-rate|mmr=f'    => \$max_mm_ratio,
            'id2gene|i2g=s'        => \$id2genMapFile,
            'mbiasplot|mbp'        => \$mBiasPlot,
            'mbiasQS|mbQS=i'       => \$mbQSt,
            'fixMateOverlap|fmo'   => \$fMO,
            'hardClipMO|hcmo'      => \$hardClipMateOverlap,
            'help|h'               => \$help,
            'man'                  => \$man,
            'version'              => \$version,
            'debug|d'              => \$DEBUG,
            );

my $debug;
if ($DEBUG) {
    eval "use Data::Dumper";
    $debug = sub { print STDERR "\n" . join( " ", @_ ) . "\n"; }
}
else {
    $debug = sub { return; };
}

# List of supported Bowtie2 versions
my %supportedBowtie2Versions = ( '2.2.9' => 1
                                );
my $bowtie2Version           = "";

# see how we are called
my $runMode = shift;

usage() and exit(0) if ( !$runMode );
usage() and exit(0) if ( $help && !$runMode );
pod2usage( -verbose => 2 ) if $man;
say $VERSION and exit(0) if ($version);

if ( !defined($runMode) ) {
    say STDERR "ERROR: No runMode specified!";
    usage();
    exit(1);
}
if ( $runMode eq 'mkbsidx' ) {
    if ($help) {
        usage_mkbsidx();
        exit(0);
    }
    say STDOUT "Starting BS index generation mode...";
    if ( !$bsIdxDir || !$bowtie2_build || !$fasta ) {
        say STDERR "ERROR: invalid/insufficient arguments";
        usage_mkbsidx();
        exit(1);
    }
    checkBowtie2();
    mkbsidx();
    exit(0);
}
elsif ( $runMode eq 'align' ) {
    if ($help) {
        usage_align();
        exit(0);
    }
    checkBowtie2();
    say STDOUT "Starting BS mapping mode...";

    # continue with main
}
else {
    say STDERR "ERROR: Invalid runMode specified!";
    usage();
    exit(1);
}

#### align mode invoced, so lets continue the main program ####

# SAM outfile
if ( !$samout ) {
    my ( $sec, $min, $hour, $mday, $mon, $year ) = localtime(time);
    my $timeStr = $year . $mon . $mday . "-" . $hour . $min . $sec;
    $samout = "meRanT_" . $timeStr . ".sam";
    warn( "No output filename for SAM alignment file specified, using: " . $samout . "\n" );
}

if ($outDir) {
    my ( $fname, $fpath, $fsuffix ) = fileparse( $samout, qw(\.sam) );
    $samout = $outDir . "/" . $fname . $fsuffix;
    checkDir($outDir);
    if ($outDir ne $fpath) {
        warn( "Output directory specified, changing SAM file location to: " . $samout . "\n" );
    }
}
else {
    my ( $fname, $fpath, $fsuffix ) = fileparse( $samout, qw(\.sam) );
    $outDir = $fpath;
    checkDir($outDir);
    warn( "No output directory specified, using: " . $outDir . "\n" );
}

if ($saveMultiMapperSeparately) {
    my ( $fname, $fpath, $fsuffix ) = fileparse( $samout, qw(\.sam) );
    $samoutMM = $outDir . "/" . $fname . "_multimappers.sam";
}

if ($reportAmbiguosMultiMap) {
    my ( $fname, $fpath, $fsuffix ) = fileparse( $samout, qw(\.sam) );
    if ($unoutDir) {
        checkDir($unoutDir);
        $aMM_fname = $unoutDir . "/" . $fname . "_ambiguos.txt";
    }
    else {
        $aMM_fname = $outDir . "/" . $fname . "_ambiguos.txt";
    }
}

my $bowtie2_fwd_convfq;
my $bowtie2_fwd_convfq_FH;
my @FWDfiles;
if ($fastQfwd) {
    @FWDfiles = split( ',', $fastQfwd );
    $bowtie2_fwd_convfq = $outDir . "/meRanT_fwd_reads_C2Tconv_" . $$ . ".fastq";
    unlink($bowtie2_fwd_convfq);
    $bowtie2_fwd_convfq_FH = IO::File->new( $bowtie2_fwd_convfq, O_RDWR | O_CREAT | O_TRUNC )
      || die( $bowtie2_fwd_convfq . ": " . $! );
}

my $bowtie2_rev_convfq;
my $bowtie2_rev_convfq_FH;
my @REVfiles;
if ($fastQrev) {
    @REVfiles = split( ',', $fastQrev );
    $bowtie2_rev_convfq = $outDir . "/meRanT_rev_reads_G2Aconv_" . $$ . ".fastq";
    unlink($bowtie2_rev_convfq);
    $bowtie2_rev_convfq_FH = IO::File->new( $bowtie2_rev_convfq, O_RDWR | O_CREAT | O_TRUNC )
      || die( $bowtie2_rev_convfq . ": " . $! );
}
push( @fqFiles, @FWDfiles, @REVfiles );

my $bowtie2_p    = $max_threads;
my $singleEnd    = 1;

# got both fwd and rev reads
if ( ($fastQfwd) && ($fastQrev) ) {
    $singleEnd = 0;

    if ( $#FWDfiles != $#REVfiles ) {
        say STDERR "Unequal number of forward and reverse read files in PE align mode";
        usage_align();
        exit(1);
    }

}

my $readCountFactor = scalar @fqFiles;
if ( $readCountFactor < 1 ) {
    say STDERR "Please specify at least one fastq file";
    usage_align();
    exit(1);
}

if ( scalar @fqFiles < 1 ) {
    say STDERR "Please specify at least one fastq file";
    usage_align();
    exit(1);
}

if ( !$bowtie2_cmd ) {
    say STDERR "bowtie2 command not sepcified";
    usage_align();
    exit(1);
}
if ( !$bowtie2_idx ) {
    say STDERR "No bowtie2 BS index specified";
    usage_align();
    exit(1);
}

if ( !$id2genMapFile ) {
    say STDERR
      "Please specify a seqID to Gene mapping file (e.g. refSeqIDs to Gene Symbols, created by mkRefSeq2GeneMap.pl)
         Format: [refseqID\tGenesymbol\tsequence length])";
    usage_align();
    exit(1);
}
else {
    readIDtoGeneMap( \%idMap, $id2genMapFile );
}

if ( ( $max_mm_ratio < 0 ) || ( $max_mm_ratio >= 1 ) ) {
    say STDERR "mismatch ratio should be >= 0 and  < 1";
    exit(1);
}

if ( !$qsOffset ) {
    say STDOUT "No quality score offset specified, trying to determin from fastq file....";
    $qsOffset = getQualityScoreOffset( $fqFiles[0] );
    if ( $qsOffset > 0 ) {
        say STDOUT "Your fastq files seem to have a quality score offset of " . $qsOffset . " ... using this value";
    }
    else {
        say STDERR "Could not determine quality score offset.";
        exit(1);
    }
    if ( $qsOffset == 33 ) {
        push( @invariableBowtie2Settings, "--phred33" );
    }
    if ( $qsOffset == 64 ) {
        push( @invariableBowtie2Settings, "--phred64" );
    }
}

# Full read BS conversion
if ($singleEnd) {
    if ($fastQfwd) {
        runSEconv( \@FWDfiles, 'C2T', $firstNreads );
    }
    if ($fastQrev) {
        runSEconv( \@REVfiles, 'G2A', $firstNreads );
    }
}
else {
    runPEconv( \@FWDfiles, \@REVfiles, $firstNreads );
}

if ($bowtie2_fwd_convfq_FH) {
    close($bowtie2_fwd_convfq_FH);
}
if ($bowtie2_rev_convfq_FH) {
    close($bowtie2_rev_convfq_FH);
}

# correct for multiple fq files
my $firstNreadsTotal = $firstNreads * $readCountFactor;

my @bowtie2args;

if ( $bowtie2_mode eq "local" ) {
    my $arg = "--local";
    push( @bowtie2args, $arg );
    if ( !$bowtie2_min_score ) {
        $bowtie2_min_score = 'G,20,8';
    }
}
elsif ( $bowtie2_mode eq "end-to-end" ) {
    my $arg = "--end-to-end";
    push( @bowtie2args, $arg );
    if ( !$bowtie2_min_score ) {
        $bowtie2_min_score = 'L,-0.4,-0.4';
    }
}
else {
    say STDERR "Unknown bowtie2 alignment mode! Can be either local or end-to-end";
    usage_align();
    exit(1);
}
if ($bowtie2_N) {
    my $arg = "-N " . $bowtie2_N;
    push( @bowtie2args, $arg );
}
if ($bowtie2_L) {
    my $arg = "-L " . $bowtie2_L;
    push( @bowtie2args, $arg );
}
if ($bowtie2_D) {
    my $arg = "-D " . $bowtie2_D;
    push( @bowtie2args, $arg );
}
if ($bowtie2_R) {
    my $arg = "-R " . $bowtie2_R;
    push( @bowtie2args, $arg );
}
if ($bowtie2_min_score) {
    my $arg = "--min-score " . $bowtie2_min_score;
    push( @bowtie2args, $arg );
}
if ($bowtie2_k) {
    my $arg = "-k " . $bowtie2_k;
    push( @bowtie2args, $arg );
}
if ($bowtie2_p) {
    my $arg = "-p " . $bowtie2_p;
    push( @bowtie2args, $arg );
}
if ( ( !$bowtie2_rev_convfq ) && ($bowtie2_fwd_convfq) ) {    # single end fwd
    $singleEnd = 1;
    my $arg = "-U " . $bowtie2_fwd_convfq;
    push( @invariableBowtie2Settings, "--norc" );
    push( @bowtie2args,               $arg );
}
elsif ( ( !$bowtie2_fwd_convfq ) && ($bowtie2_rev_convfq) ) {    # single end rev
    $singleEnd = 1;
    my $arg = "-U " . $bowtie2_rev_convfq;
    push( @invariableBowtie2Settings, "--nofw" );
    push( @bowtie2args,               $arg );
}
else {                                                           # paired end
    $singleEnd = 0;
    # attention due to a change in the bowtie2 wrapper script introduced in
    # version 2.2.5 we have to specify the fwd and reverse read without space
    # after "-1" and "-2" parameter respectively.
    my $arg = "-1" . $bowtie2_fwd_convfq;
    push( @bowtie2args, $arg );
    $arg = "-2" . $bowtie2_rev_convfq;
    push( @bowtie2args, $arg );
    $arg = "-I" . $bowtie2_I;
    push( @bowtie2args, $arg );
    $arg = "-X" . $bowtie2_X;
    push( @bowtie2args,               $arg );
    push( @invariableBowtie2Settings, "--no-discordant" );
    push( @invariableBowtie2Settings, "--no-mixed" );
    push( @invariableBowtie2Settings, "--norc" );
}
if ($bowtie2_idx) {
    my $arg = "-x " . $bowtie2_idx;
    push( @bowtie2args, $arg );
}
else {
    say STDERR "No bowtie2 index specified";
    usage_align();
    exit(1);
}

my $fixMateOverlap;
my $getMateOverlap;
if ( ($fMO) && ( !$singleEnd ) ) {
    $fixMateOverlap = \&fixMateOverlap;
    $getMateOverlap = \&getMateOverlap;
}
else {
    $fixMateOverlap = sub { return (0); };
    $getMateOverlap = sub { return (0); };
}

my @bowtie2cmd = ( $bowtie2_cmd, @invariableBowtie2Settings, @bowtie2args );

# Start thread for m-bias plot data generation
if ($mBiasPlot) {

    # Load GD modules if we want to create a m-bias plot
    eval("use GD;");
    eval("use GD::Text;");
    eval("use GD::Text::Align;");
    eval("use GD::Graph::lines;");

    my ( $fname, $fpath, $fsuffix ) = fileparse( $samout, qw(\.sam) );
    $mBiasPlotOutF = $outDir . "/" . $fname . "_m-bias_fwd.png";
    $mBiasPlotOutR = $outDir . "/" . $fname . "_m-bias_rev.png";

    $mbq = Thread::Queue->new();
    ($mBiasThr) = threads->create( \&mBiasCounter );
}
else {
    sub enqueue { }
    my %dummy;
    $mbq = \%dummy;
    bless $mbq;
}

# Start the Bowtie mapping process
# using the parameters specified above
runBowtie2( \@bowtie2cmd );
say STDOUT "DONE: BS mapping";

# Stop thread for m-bias plot data generation and generate the plot
if ($mBiasPlot) {
    $mbq->end();
    (%mBiasData) = $mBiasThr->join();

    plot_mBias( $mBiasData{mBiasDataF},
                $mBiasData{mBiasDataFhq},
                $mBiasData{mBiasReadDataF},
                $mBiasPlotOutF, 'forward' );

    if ( $singleEnd == 0 ) {
        plot_mBias( $mBiasData{mBiasDataR},
                    $mBiasData{mBiasDataRhq},
                    $mBiasData{mBiasReadDataR},
                    $mBiasPlotOutR, 'reverse' );
    }
}

###################################### Print out some stats ######################################
if ( $mappingStats->{reads} == 0 ) {
    say STDERR "Not enough data for printing statistics! " . __LINE__;
}
else {
    my $uGmTpct = sprintf( "%.2f", $mappingStats->{uniqGeneMultiTranscript} * 100 / $mappingStats->{reads} );
    my $uGuTpct = sprintf( "%.2f", $mappingStats->{uniqGeneUniqTranscript} * 100 / $mappingStats->{reads} );
    my $totalMapped               = $mappingStats->{uniqGeneMultiTranscript} + $mappingStats->{uniqGeneUniqTranscript};
    my $totalMappedpct            = sprintf( "%.2f", $totalMapped * 100 / $mappingStats->{reads} );
    my $totalBWTunalpct           = sprintf( "%.2f", $mappingStats->{bowtie2unal} * 100 / $mappingStats->{reads} );
    my $totalDiscard              = ( $mappingStats->{reads} - $totalMapped );
    my $totalEditDistFiltered_a   = $mappingStats->{editDistFilter};
    my $totalFlagFiltered_a       = $mappingStats->{flagFilter};
    my $totalMGMapperFiltered     = $mappingStats->{mgMapperFilter};
    my $totalMGMapperFilteredpct  = sprintf( "%.2f", $totalMGMapperFiltered * 100 / $mappingStats->{reads} );
    my $totalMTpMapperFiltered    = $mappingStats->{mtpMapperFilter};
    my $totalMTpMapperFilteredpct = sprintf( "%.2f", $totalMTpMapperFiltered * 100 / $mappingStats->{reads} );

    my $col2width = length( $mappingStats->{reads} );
    my $f_str     = "%*i\t\(%6.2f%%\)";
    my $f_str_b   = "%*i\t\  ------ \ ";
    my $sl1       = sprintf( $f_str, $col2width, $mappingStats->{reads}, 100 );
    my $sl2       = sprintf( $f_str, $col2width, $totalMapped, $totalMappedpct );
    my $sl2b      = sprintf( $f_str, $col2width, $mappingStats->{bowtie2unal}, $totalBWTunalpct );
    my $sl3       = sprintf( $f_str, $col2width, $totalDiscard, ( 100 - $totalMappedpct ) );
    my $sl4       = sprintf( $f_str_b, $col2width, $totalFlagFiltered_a );
    my $sl5       = sprintf( $f_str_b, $col2width, $totalEditDistFiltered_a );
    my $sl6       = sprintf( $f_str, $col2width, $totalMGMapperFiltered, $totalMGMapperFilteredpct );
    my $sl7       = sprintf( $f_str, $col2width, $totalMTpMapperFiltered, $totalMTpMapperFilteredpct );
    my $sl8       = sprintf( $f_str, $col2width, $mappingStats->{uniqGeneUniqTranscript}, $uGuTpct );
    my $sl9       = sprintf( $f_str, $col2width, $mappingStats->{uniqGeneMultiTranscript}, $uGmTpct );

    print STDERR "

FASTQ input files:\n\t\t" . join( "\n\t\t", @fqFiles ) . "

Total # of reads:                                    " . $sl1 . "
Total # of mapped reads:                             " . $sl2 . "
Total # of unmapped reads:                           " . $sl2b . "
Total # of unmapped + filtered reads:                " . $sl3 . "
Alignments filtered with incorrect mapping flag:     " . $sl4 . "
Alignments filtered with edit distance > " . $max_edit_dist . ":          " . $sl5 . "
Reads filtered mapping to multiple genes:            " . $sl6 . "
Reads filtered mapping to multiple places on a gene: " . $sl7 . " 
Reads mapped to uniq genes and transcripts:          " . $sl8 . "
Reads mapped multiple transcripts of a uniq gene:    " . $sl9 . " 

";
}
###################################### END stats ######################################

# Cleanup temprary files created sofar
cleanup();

# Generate a sorted and indexed BAM file
my $bamFile;
my $bamFileMM;
if ( $ommitBAM == 0 ) {
    $bamFile = sam_to_bam($samout) unless ( ($mappingStats->{uniqGeneMultiTranscript} + $mappingStats->{uniqGeneUniqTranscript}) == 0 );
    if ($saveMultiMapperSeparately) {
        $bamFileMM = sam_to_bam($samoutMM) unless ( ($mappingStats->{mgMapperFilter} + $mappingStats->{mtpMapperFilter}) == 0 );
    }
    if ($deleteSAM) {
        unlink($samout);
        unlink($samoutMM) if ($saveMultiMapperSeparately);
    }

}

# DONE our job exit now and have a rest
exit(0);

###################################### SUBROUTINES  ######################################
sub cleanup {
    unlink($bowtie2_fwd_convfq) if ($bowtie2_fwd_convfq);
    unlink($bowtie2_rev_convfq) if ($bowtie2_rev_convfq);
}

sub mkbsidx {
    my $bowtie2buildcmd;

    checkDir($bsIdxDir);

    my @FastafilesCmdLine = glob $fasta;
    my @Fastafiles;

    if ( ( scalar(@FastafilesCmdLine) == 1 ) && ( $FastafilesCmdLine[0] =~ /\,/ ) ) {
        @Fastafiles = split( ',', $FastafilesCmdLine[0] );
    }
    elsif ( scalar(@FastafilesCmdLine) == 0 ) {
        die("Could not parse the fasta file argument");
    }
    else {
        @Fastafiles = @FastafilesCmdLine;
    }

    # Check
    checkFastaFiles( \@Fastafiles );

    my @FastaConvFiles;
    my $convertor = sub { $_[0] =~ tr/Cc/Tt/; };

    my $procs = ( $max_threads < ( ( scalar @Fastafiles ) * 2 ) ) ? $max_threads : ( ( scalar @Fastafiles ) * 2 );
    print STDOUT "Running $procs processes for BS index geneartion\n\n";

    my $pm = Parallel::ForkManager->new($max_threads);

    $pm->run_on_finish( \&checkChild );

    # loop through the fasta files
    for my $child ( 0 .. $#Fastafiles ) {

        my ( $fname, $fpath, $fsuffix ) = fileparse( $Fastafiles[$child], qr/\.[^.]*$/ );
        my $fastaConv = abs_path($bsIdxDir . "/" . $fname . "_C2T" . $fsuffix);

        push( @FastaConvFiles, $fastaConv );

        say STDOUT "BS converting reference: " . $fasta . " -> " . $fastaConv;

        my $pid = $pm->start( $Fastafiles[$child] ) and next;

        ### in child ###
        # do the job
        my $returnCode = bsconvertFA( $Fastafiles[$child], $fastaConv, $convertor );

        # finished with this chromosome
        $pm->finish($returnCode);
    }
    $pm->wait_all_children;

    print STDOUT "\n\n";

    my $idxName = getIdxName( \@FastaConvFiles );
    my $fastaC2T = join( ',', @FastaConvFiles );

    my $buildThreads = ($max_threads > 1) ? $max_threads : 1;

    $bowtie2buildcmd = $bowtie2_build . " --threads " . $buildThreads . " -q -o 3 -f " . $fastaC2T . " " . $idxName . " 2>&1";

    chdir($bsIdxDir) || die( $bsIdxDir . ": " . $! );

    # fork the BOWTIE2 index builder processes
    my $pid = open( my $IDXBUILD, "-|", $bowtie2buildcmd ) || die($!);

    say STDOUT "Indexing ... it may take a while to create: " . $idxName;

    my $bwt2b_error = 0;
    
    $IDXBUILD->autoflush();
    while ( defined( my $idxBuildOut = $IDXBUILD->getline() ) ) {
        print STDOUT $idxBuildOut;
        $bwt2b_error = 1 if ($idxBuildOut =~ /^Error:/);
    }
    waitpid( $pid, 0 );
    close($IDXBUILD);
    
    if ($bwt2b_error) {
        say STDERR "Could not generate the BS index...";
        exit(1);
    }

    print STDOUT "\n\n";
    say STDOUT "Finished building the BS index for: " . $fasta;
    say STDOUT "BS index is stored in: " . $bsIdxDir;
    say STDOUT "BS index name: " . $idxName . "\n";
    say STDOUT "Use: '" . $idxName
      . "' as '-x' option or assign it to the 'BS_BWT2IDX' environment variable in the 'align' runMode";
    print STDOUT "\n\n";
}

sub checkChild {
    my ( $pid, $exit_code, $ident ) = @_;
    print "** " . $ident . " just got finished conversion with PID " . $pid . " and exit code: " . $exit_code . "\n";
    die($exit_code) if ( $exit_code ne 0 );
    return ($exit_code);
}

sub getIdxName {
    my $fastaFiles = shift;

    my $idxName;

    my @fnames = stripPathSuffix($fastaFiles);
    my @sortedFnames = sort { length($a) <=> length($b) } @fnames;

    my @charArray = split( //, $sortedFnames[0] );
    my $pos       = 0;
    my $match     = 1;

    foreach my $char (@charArray) {
        foreach my $file ( @sortedFnames[ 0 .. $#sortedFnames ] ) {
            if ( $char ne substr( $file, $pos, 1 ) ) {
                $match = 0;
                last;
            }
        }
        last if ( $match == 0 );
        $pos++;
    }
    if ( $pos > 0 ) {
        return ( abs_path($bsIdxDir) . "/" . substr( $sortedFnames[0], 0, $pos ) );
    }
    else {
        return ( abs_path($bsIdxDir) . "/meRanT_C2T" );
    }
}

sub stripPathSuffix {
    my $fileList = shift;

    my @fn;
    foreach ( @{$fileList} ) {
        my ( $fname, $fpath, $fsuffix ) = fileparse( $_, qr/\.[^.]*$/ );
        push( @fn, $fname );
    }

    return (@fn);
}

sub checkFastaFiles {
    my $fileList = shift;

    foreach my $file ( @{$fileList} ) {
        die( "Could not read: " . $file ) if ( !-r $file );
    }
    return;
}

sub bsconvertFA {
    my $fastaFile = shift;
    my $fastaConv = shift;
    my $convertor = shift;

    unlink($fastaConv);

    my $fastaConv_FH = IO::File->new( $fastaConv, O_RDWR | O_CREAT | O_TRUNC ) || return ( $fastaConv . ": " . $! );
    my $fasta_FH     = IO::File->new( $fastaFile, O_RDONLY | O_EXCL )          || return ( $fastaFile . ": " . $! );

    my $fastaLine;
    while ( defined( $fastaLine = $fasta_FH->getline() ) ) {
        if ( $fastaLine =~ /^\>/ ) {
            $fastaConv_FH->print($fastaLine);
        }
        else {
            $convertor->($fastaLine);
            $fastaConv_FH->print($fastaLine);
        }
    }
    $fasta_FH->close();
    undef($fasta_FH);
    $fastaConv_FH->close();
    undef($fastaConv_FH);
}

sub runBowtie2 {
    my $bowtie2cmd = shift;

    my %mappedSeqs;
    my $mappedSeqsRef = \%mappedSeqs;

    my @samRawHrd;
    my $samRawHrdRef = \@samRawHrd;

    my %lastFQrec = ( readID => "", seq => "", readID2 => "", QS => "", readCounter => 0 );
    my $lastFQrecRef = \%lastFQrec;

    my @samL;
    my $samLref = \@samL;

    my $samout_tmp = $samout . "_" . $$;
    unlink($samout_tmp);
    my $sam_tmp_FH = IO::File->new( $samout_tmp, O_RDWR | O_CREAT | O_TRUNC ) || die( $samout_tmp . ": " . $! );

    my $aMM_FH = undef;
    if ($reportAmbiguosMultiMap) {
        unlink($aMM_fname);
        $aMM_FH = IO::File->new( $aMM_fname, O_RDWR | O_CREAT | O_TRUNC ) || die( $aMM_fname . ": " . $! );
    }

    my $sam_MM_tmp_FH = undef;
    my $samout_MM_tmp = $samoutMM . "_" . $$;
    if ($saveMultiMapperSeparately) {
        unlink($samout_MM_tmp);
        $sam_MM_tmp_FH = IO::File->new( $samout_MM_tmp, O_RDWR | O_CREAT | O_TRUNC )
          || die( $samout_MM_tmp . ": " . $! );
    }
    else {
        $sam_MM_tmp_FH = undef;
    }

    # fork the BOWTIE2 mapper processes
    &$debug(@$bowtie2cmd);
    my $pid = open( my $BOWTIE2, "-|", @$bowtie2cmd ) || die($!);

    $BOWTIE2->autoflush();

    if ( $singleEnd == 1 ) {    # single end

        my $fqFIFO = $outDir . "/meRanT_FQ_spooler_" . $$;
        unlink($fqFIFO);
        mkfifo( $fqFIFO, 0700 ) || die("can't mkfifo $fqFIFO: $!");
        my $fq_FH = IO::File->new( $fqFIFO, O_RDWR ) || die( $fqFIFO . ": " . $! );

        # Create and detach a fastq multifile spooler fifo

        # Create and detach a fastq multifile spooler fifo
        my $fqFilesRef = \@FWDfiles;
        my $readDirection = 0;
        if ($fastQrev) {
            $fqFilesRef = \@REVfiles;
            $readDirection = 1;
        }
        my $fqSpoolerThr = threads->create( \&fqSpooler, $fq_FH, $fqFilesRef, $firstNreads );

        my @unmappedReadsFHs;    # file handles for unmapped reads
        foreach my $fqf (@fqFiles) {
            push( @unmappedReadsFHs, getUnmappedFH($fqf) ) if ($bowtie2_un);
        }

        my ($bowtie2_line);
        while (1) {              # process SAM header
            $bowtie2_line = <$BOWTIE2>;
            if ( $bowtie2_line =~ /^@/ ) {
                processbowtie2samHdr( $bowtie2_line, $samRawHrdRef );
            }
            else {
                last;
            }
        }

        # mTypes: -1: unmapped, 0: uniquely mapped, 1: multimapper
        while (1) {
            processbowtie2samSE( $fq_FH, $lastFQrecRef, $bowtie2_line, $mappedSeqsRef, $sam_tmp_FH,
                                 \@unmappedReadsFHs, $readDirection );
            if ( $lastFQrecRef->{mType} == 1 ) {
                my @multiMappers;
                while (1) {    # advance sam output
                    push( @multiMappers, $bowtie2_line );
                    last unless ( defined( $bowtie2_line = <$BOWTIE2> ) );
                    my @samF = split( '\t', $bowtie2_line );
                    last if ( $samF[0] ne $lastFQrecRef->{readID} );
                }
                process_multimappersSE( \@multiMappers, $lastFQrecRef,  $mappedSeqsRef,
                                        $sam_tmp_FH,    $sam_MM_tmp_FH, $aMM_FH, $readDirection );
                last unless ( defined($bowtie2_line) );
            }
            else {
                last unless ( defined( $bowtie2_line = <$BOWTIE2> ) );
            }
        }

        $fqSpoolerThr->join();

        $fq_FH->close();
        undef($fq_FH);
        unlink($fqFIFO);
        foreach my $fqfh (@unmappedReadsFHs) {
            $fqfh->close();
            undef($fqfh);
        }
        $mappingStats->{reads} = $lastFQrecRef->{readCounter};

    }
    else {    # paired end

        %lastFQrec = (
                       FWD => {
                                readID      => "",
                                seq         => "",
                                readID2     => "",
                                QS          => "",
                                readCounter => 0
                         },
                       REV => {
                                readID      => "",
                                seq         => "",
                                readID2     => "",
                                QS          => "",
                                readCounter => 0
                         },
                         );

        # mk fifo for forward reads
        my $fqFIFO_fwd = $outDir . "/meRanT_FQ_spooler_fwd_" . $$;
        unlink($fqFIFO_fwd);
        mkfifo( $fqFIFO_fwd, 0700 ) || die("can't mkfifo $fqFIFO_fwd: $!");
        my $fq_FH_fwd = IO::File->new( $fqFIFO_fwd, O_RDWR ) || die( $fqFIFO_fwd . ": " . $! );

        # mk fifo for reverse reads
        my $fqFIFO_rev = $outDir . "/meRanT_FQ_spooler_rev_" . $$;
        unlink($fqFIFO_rev);
        mkfifo( $fqFIFO_rev, 0700 ) || die("can't mkfifo $fqFIFO_rev: $!");
        my $fq_FH_rev = IO::File->new( $fqFIFO_rev, O_RDWR ) || die( $fqFIFO_fwd . ": " . $! );

        my %fq_FH = ( FWD => $fq_FH_fwd, REV => $fq_FH_rev );
        my $fq_FHref = \%fq_FH;

        # Create and detach a fastq multifile spooler fifo
        my $fqFWDspoolerThr = threads->create( \&fqSpooler, $fq_FH_fwd, \@FWDfiles, $firstNreads );
        my $fqREVspoolerThr = threads->create( \&fqSpooler, $fq_FH_rev, \@REVfiles, $firstNreads );

        my @unmappedFWDReadsFHs;
        foreach my $fqf (@FWDfiles) {
            push( @unmappedFWDReadsFHs, getUnmappedFH($fqf) ) if ($bowtie2_un);
        }
        my @unmappedREVReadsFHs;
        foreach my $fqf (@REVfiles) {
            push( @unmappedREVReadsFHs, getUnmappedFH($fqf) ) if ($bowtie2_un);
        }
        my %unmappedReadsFHs = ( FWD => [@unmappedFWDReadsFHs], REV => [@unmappedREVReadsFHs] );
        my $unmappedReadsFHsRef = \%unmappedReadsFHs;

        my ( $bowtie2_line_1, $bowtie2_line_2 );

        while (1) {    # Start Header
            last unless ( defined( $bowtie2_line_1 = <$BOWTIE2> ) && defined( $bowtie2_line_2 = <$BOWTIE2> ) );
            if ( ( $bowtie2_line_1 =~ /^\@/ ) && ( $bowtie2_line_2 =~ /^\@/ ) ) {
                processbowtie2samHdr( $bowtie2_line_1, $samRawHrdRef );
                processbowtie2samHdr( $bowtie2_line_2, $samRawHrdRef );
                next;
            }
            elsif ( ( $bowtie2_line_1 =~ /^\@/ ) && ( $bowtie2_line_2 !~ /^\@/ ) ) {
                processbowtie2samHdr( $bowtie2_line_1, $samRawHrdRef );
                $bowtie2_line_1 = $bowtie2_line_2;
                $bowtie2_line_2 = <$BOWTIE2>;
            }
            else {
                last;
            }
        }    # End Header

        # Fill read pair array with first read pair
        fillReadPairArray( $samLref, $bowtie2_line_1, $bowtie2_line_2 );

        while (1) {    # Start process
            processbowtie2samPE( $fq_FHref, $lastFQrecRef, $samLref, $mappedSeqsRef, $sam_tmp_FH,
                                 $unmappedReadsFHsRef );
            if ( $lastFQrecRef->{mType} == 1 ) {
                my @multiMappers;
                my $i = 0;
                while (1) {    # advance sam output
                    $multiMappers[$i]->[0] = $bowtie2_line_1;
                    $multiMappers[$i]->[1] = $bowtie2_line_2;
                    $i++;
                    last
                      unless (    ( defined( $bowtie2_line_1 = <$BOWTIE2> ) )
                               && ( defined( $bowtie2_line_2 = <$BOWTIE2> ) ) );
                    my @samF = split( '\t', $bowtie2_line_1 );
                    last if ( $samF[0] ne $lastFQrecRef->{FWD}->{readID} );
                }
                process_multimappersPE( \@multiMappers, $lastFQrecRef,  $mappedSeqsRef,
                                        $sam_tmp_FH,    $sam_MM_tmp_FH, $aMM_FH );
                fillReadPairArray( $samLref, $bowtie2_line_1, $bowtie2_line_2 );
                last unless ( defined($bowtie2_line_1) && defined($bowtie2_line_2) );
            }
            else {
                last unless ( defined( $bowtie2_line_1 = <$BOWTIE2> ) && defined( $bowtie2_line_2 = <$BOWTIE2> ) );
                fillReadPairArray( $samLref, $bowtie2_line_1, $bowtie2_line_2 );
            }
        }

        $fqFWDspoolerThr->join();
        $fqREVspoolerThr->join();

        $fq_FH_fwd->close();
        $fq_FH_rev->close();
        undef($fq_FH_fwd);
        undef($fq_FH_rev);
        unlink($fqFIFO_fwd);
        unlink($fqFIFO_rev);
        foreach my $fqfh ( @unmappedFWDReadsFHs, @unmappedREVReadsFHs ) {
            $fqfh->close();
            undef($fqfh);
        }
        $mappingStats->{reads} = $lastFQrecRef->{FWD}->{readCounter} + $lastFQrecRef->{REV}->{readCounter};
    }

    waitpid( $pid, 0 );
    close($BOWTIE2);
    close($sam_tmp_FH);
    $sam_tmp_FH = undef;

    if ($saveMultiMapperSeparately) {
        close($sam_MM_tmp_FH);
        $sam_MM_tmp_FH = undef;
    }

    if ($reportAmbiguosMultiMap) {
        close($aMM_FH);
        $aMM_FH = undef;
    }

    my $data_buffer;
    unlink($samout);

    my $samFH = IO::File->new( $samout, O_RDWR | O_CREAT | O_TRUNC ) || die( $samout . ": " . $! );
    $sam_tmp_FH = IO::File->new( $samout_tmp, O_RDONLY | O_EXCL ) || die( $samout_tmp . ": " . $! );

    my $headerT = 'u';
    if ( !$saveMultiMapperSeparately ) {
        $headerT = 'um';
    }
    writeFinalSAM( $samFH, $sam_tmp_FH, \@samRawHrd, $mappedSeqsRef, $headerT );

    close($sam_tmp_FH);
    close($samFH);
    $sam_tmp_FH = undef;
    $samFH      = undef;
    unlink($samout_tmp);

    if ($saveMultiMapperSeparately) {
        unlink($samoutMM);
        my $samMMFH = IO::File->new( $samoutMM, O_RDWR | O_CREAT | O_TRUNC ) || die( $samoutMM . ": " . $! );
        $sam_MM_tmp_FH = IO::File->new( $samout_MM_tmp, O_RDONLY | O_EXCL ) || die( $samout_MM_tmp . ": " . $! );

        writeFinalSAM( $samMMFH, $sam_MM_tmp_FH, \@samRawHrd, $mappedSeqsRef, 'm' );

        close($sam_MM_tmp_FH);
        close($samMMFH);
        $sam_MM_tmp_FH = undef;
        $samMMFH       = undef;
        unlink($samout_MM_tmp);
    }

}

sub fillReadPairArray {
    my $samLref        = shift;
    my $bowtie2_line_1 = shift;
    my $bowtie2_line_2 = shift;

    # Fill read pair array with first read pair
    $samLref->[0] = $bowtie2_line_1;
    $samLref->[1] = $bowtie2_line_2;
}

sub writeFinalSAM {
    my $samFH         = shift;
    my $sam_tmp_FH    = shift;
    my $samRawHrd     = shift;
    my $mappedSeqsRef = shift;
    my $type          = shift;

    my $data_buffer;

    foreach my $rawHrdLine ( @{$samRawHrd} ) {
        if ( $rawHrdLine !~ /^\@SQ/ ) {
            $samFH->syswrite($rawHrdLine);
            next;
        }
        my @samHf = split( '\t', $rawHrdLine );
        my $hdrSeqID = ( split( ':', $samHf[1] ) )[1];
        if ( $type eq 'um' ) {
            if ( defined( $mappedSeqsRef->{$hdrSeqID} ) ) {
                $samFH->syswrite($rawHrdLine);
            }
        }
        else {
            if ( defined( $mappedSeqsRef->{$hdrSeqID}->{$type} ) ) {
                $samFH->syswrite($rawHrdLine);
            }
        }
    }

    my $SAMpg_line =
      "\@PG\tID:meRanT\tPN:meRanT\tVN:" . $VERSION . "\tCL:\"" . $0 . " " . join( ' ', @SAMhdr_cmd_line ) . "\"\n";
    $samFH->syswrite($SAMpg_line);

    $sam_tmp_FH->sysseek( 0, SEEK_SET );
    $samFH->syswrite($data_buffer) while $sam_tmp_FH->sysread( $data_buffer, 8192 );
}

sub processbowtie2samHdr {
    my $sam_rec = shift;
    my $rawHrd  = shift;

    if ( $sam_rec =~ /^@/ ) {

        # print $samFH $sam_rec;
        $rawHrd->[ $#$rawHrd + 1 ] = $sam_rec;
        return;
    }
}

sub processbowtie2samSE {
    my $fq_FH            = shift;
    my $lastFQrec        = shift;
    my $sam_rec          = shift;
    my $mappedSeqs       = shift;
    my $samFH            = shift;
    my $unmappedReadsFHs = shift;
    my $readDirection    = shift;

    my $fq_rec_id     = "";
    my $fq_rec_id_sam = "";
    my $unmappedReadsFH;

    my @samFields = split( '\t', $sam_rec );
    my $readID = $samFields[0];

    # orignal FQ record
    unless ( $readID eq $lastFQrec->{readID} ) {

        ( $fq_rec_id, $fq_rec_id_sam, $unmappedReadsFH ) =
          getOrigFQrec( $readID, $fq_FH, $lastFQrec, $unmappedReadsFHs );

    }

    my $mapped = ( $samFields[1] == 4 ) ? 0 : 1;
    if ( !$mapped ) {
        if ( $bowtie2_un ) {
            print $unmappedReadsFH $fq_rec_id . $lastFQrec->{seq} . "\n" . $lastFQrec->{readID2} . $lastFQrec->{QS};
            $mappingStats->{bowtie2unal} += 1;
        }
        $lastFQrec->{mType} = -1;
        return;
    }

    my %auxTags = samGet_auxTags(@samFields);
    if ( exists( $auxTags{XS} ) ) {    # Alternative alignments found
        $lastFQrec->{mType} = 1;
        return;
    }

    my $editDist   = $auxTags{NM}{v};
    my $readLength = length( $samFields[9] );
    if ( checkEditDist( $editDist, $readLength, $readID ) > 0 ) {
        $lastFQrec->{mType} = -1;
        return;
    }

    my $flag = $samFields[1];
    my $read = $lastFQrec->{seq};
    if ( ( $flag & 0x10 ) != 0 ) {
        $samFields[9] = reverseComp($read);
    }
    else {
        $samFields[9] = $read;
    }

    # print  join("\t", @samFields);

    $lastFQrec->{mType} = 0;
    $mappedSeqs->{ $samFields[2] }->{u} = 1;

    # add GeneTag to Sam record
    my $auxTagGeneName = (defined($idMap{ $samFields[2] }->{'geneName'})) ? $idMap{ $samFields[2] }->{'geneName'} : "unknown";
    splice( @samFields, 11, 0, ( "ZG:Z:" . $auxTagGeneName ) );

    my $softClippedBases = getSoftClippedSE(\@samFields);
    
    # add to m-bias plot process queue
    if ( ( $samFields[1] & 0x10 ) == 0 ) {
        $mbq->enqueue( $read, $lastFQrec->{QS}, $softClippedBases->[0], $softClippedBases->[1], $readDirection );
    }
    else {
        $mbq->enqueue( $read, $lastFQrec->{QS}, $softClippedBases->[1], $softClippedBases->[0], $readDirection );
    }

    # write sam line
    print $samFH join( "\t", @samFields );
    $mappingStats->{uniqGeneUniqTranscript}++;
}

# process_multimappersSE(\@multiMappers, $lastFQrecRef, $mappedSeqsRef, $sam_tmp_FH, $sam_MM_tmp_FH, $aMM_FH);
sub process_multimappersSE {
    my $multiMappers  = shift;    # arrayRef
    my $lastFQrec     = shift;
    my $mappedSeqs    = shift;
    my $samFH         = shift;
    my $samMMFH       = shift;
    my $aMMFH         = shift;
    my $readDirection = shift;

    my ( $mgM, $mtpM, $uGmt, $bl ) = ( 0, 0, 0, 0 );

    my %bestAlignment = (
                          geneName  => '',
                          ltr       => '',
                          t_length  => 0,
                          AS        => -1000,
                          samFields => [],
                          samLine   => '',
                          read      => '',
                          QS        => '',
                          );

    my %ambiguosMultiMappers;
    my $i = 0;

    foreach my $samLine ( @{$multiMappers} ) {
        my @samFields  = split( '\t', $samLine );
        my $seqID      = $samFields[2];
        my $geneName   = $idMap{$seqID}->{'geneName'};
        my $t_length   = $idMap{$seqID}->{'t_length'};
        my %auxTags    = samGet_auxTags(@samFields);
        my $editDist   = $auxTags{NM}->{v};
        my $readLength = length( $samFields[9] );

        if ( checkEditDist( $editDist, $readLength, $samFields[0] ) > 0 ) {
            next;
        }

        # add GeneTag to Sam record
        my $auxTagGeneName = (defined($idMap{$seqID}->{'geneName'})) ? $idMap{$seqID}->{'geneName'} : "unknown";
        splice( @samFields, 11, 0, ( "ZG:Z:" . $auxTagGeneName ) );

        my %thisAlignment = (
                              geneName  => $geneName,
                              ltr       => $seqID,
                              t_length  => $t_length,
                              AS        => $auxTags{AS}->{v},
                              samFields => [@samFields],
                              samLine   => $samLine,
                              read      => '',
                              QS        => '',
                              );

        # get the BS seq
        my $flag = $samFields[1];
        my $read = $lastFQrec->{seq};
        if ( ( $flag & 0x10 ) != 0 ) {
            $thisAlignment{samFields}->[9] = reverseComp($read);
        }
        else {
            $thisAlignment{samFields}->[9] = $read;
        }
        $thisAlignment{read} = $read;
        $thisAlignment{QS}   = $lastFQrec->{QS};

        if ( $i == 0 ) {    # first alignment
            %bestAlignment = %thisAlignment;
        }
        elsif ( $auxTags{AS}->{v} > $bestAlignment{AS} ) {
            %bestAlignment = %thisAlignment;
            $bl            = 0;
        }
        elsif ( ( $bestAlignment{geneName} ne $geneName ) && ( $auxTags{AS}->{v} == $bestAlignment{AS} ) ) {

            # multi gene mapper
            $ambiguosMultiMappers{mgM}->{$i}->{ambData} =
              [ "MG", $lastFQrec->{readID}, $samFields[2], $geneName, $bestAlignment{geneName} ];
            $ambiguosMultiMappers{mgM}->{$i}->{samFields} = $thisAlignment{samFields};
            $bl = 1;
        }
        elsif ( ( $bestAlignment{ltr} eq $seqID ) && ( $auxTags{AS}->{v} == $bestAlignment{AS} ) ) {

            # multi mapper on Transcript
            $ambiguosMultiMappers{mtpM}->{$i}->{ambData} =
              [ "MP", $lastFQrec->{readID}, $samFields[2], $samFields[3], $bestAlignment{samFields}->[3], $geneName ];
            $ambiguosMultiMappers{mtpM}->{$i}->{samFields} = $thisAlignment{samFields};
            $bl = 1;
        }
        else {
            if ( ( $t_length > $bestAlignment{'t_length'} ) && ( $auxTags{AS}->{v} >= $bestAlignment{AS} ) ) {

                # update if maps to same gene longer transcript with same or better AS
                %bestAlignment = %thisAlignment;
                $bl            = 0;
            }
        }
        $i++;
    }

    # No vaild alignment found
    return if ( $i == 0 );

    if ( $bl == 1 ) {
        foreach my $idx ( keys %{ $ambiguosMultiMappers{mgM} } ) {
            if ($aMMFH) {
                print $aMMFH join( "\t", @{ $ambiguosMultiMappers{mgM}->{$idx}->{ambData} } ) . "\n";
            }
            if ($saveMultiMapperSeparately) {
                print $samMMFH join( "\t", @{ $ambiguosMultiMappers{mgM}->{$idx}->{samFields} } );
                $mappedSeqs->{ $ambiguosMultiMappers{mgM}->{$idx}->{samFields}->[2] }->{m} = 1;
            }
            $mgM = 1;
        }
        foreach my $idx ( keys %{ $ambiguosMultiMappers{mtpM} } ) {
            if ($aMMFH) {
                print $aMMFH join( "\t", @{ $ambiguosMultiMappers{mtpM}->{$idx}->{ambData} } ) . "\n";
            }
            if ($saveMultiMapperSeparately) {
                print $samMMFH join( "\t", @{ $ambiguosMultiMappers{mtpM}->{$idx}->{samFields} } );
                $mappedSeqs->{ $ambiguosMultiMappers{mtpM}->{$idx}->{samFields}->[2] }->{m} = 1;
            }
            $mtpM = 1;
        }
    }
    else {

        my $softClippedBases = getSoftClippedSE($bestAlignment{samFields});
        
        # add to m-bias plot process queue
        if ( ( $bestAlignment{samFields}->[1] & 0x10 ) == 0 ) {
            $mbq->enqueue( $bestAlignment{read}, $bestAlignment{QS}, $softClippedBases->[0], $softClippedBases->[1], $readDirection );
        }
        else {
            $mbq->enqueue( $bestAlignment{read}, $bestAlignment{QS}, $softClippedBases->[1], $softClippedBases->[0], $readDirection );
        }

        # make primary alignment
        $bestAlignment{samFields}->[1] -= 256 if ( ( $bestAlignment{samFields}->[1] & 0x100 ) == 256 );

        # write sam line
        print $samFH join( "\t", @{ $bestAlignment{samFields} } );
        $mappedSeqs->{ $bestAlignment{samFields}->[2] }->{u} = 1;
        $uGmt = 1;
        ( $mgM, $mtpM ) = ( 0, 0 );
    }

    $mappingStats->{mgMapperFilter}          += $mgM;
    $mappingStats->{mtpMapperFilter}         += $mtpM;
    $mappingStats->{uniqGeneMultiTranscript} += $uGmt;

    $lastFQrec->{mType} = -1;

    return;
}

sub checkEditDist {
    my $editDist   = shift;
    my $readLength = shift;
    my $readID     = shift;

    my $misMatchRatio = $editDist / $readLength;
    if ( ( $editDist > $max_edit_dist ) || ( $misMatchRatio > $max_mm_ratio ) ) {
        $mappingStats->{editDistFilter}++;
        return ($editDist);
    }

    return (0);
}

sub sortPE {
    my $sam_rec = shift;

    # Make sure that fwd read is idx 0 and rev read is idx 1
    my @mates     = ();
    my @mates_tmp = ();

    $mates_tmp[0]{samFields} = [ ( split( '\t', $sam_rec->[0] ) ) ];
    $mates_tmp[0]{readID} = $mates_tmp[0]->{samFields}->[0];

    $mates_tmp[1]{samFields} = [ ( split( '\t', $sam_rec->[1] ) ) ];
    $mates_tmp[1]{readID} = $mates_tmp[1]->{samFields}->[0];

    if ( $mates_tmp[0]{readID} ne $mates_tmp[1]{readID} ) {
        print STDERR "Not a paired alignment: " . $mates_tmp[0]{readID} . " : " . $mates_tmp[1]{readID} . "\n";
        die();
    }    # Not a pair

    # sort mates
    my $fqMateNr_0 = ( ( $mates_tmp[0]->{samFields}->[1] & 0x40 ) != 0 ) ? 1 : 2;
    my $fqMateNr_1 = ( ( $mates_tmp[1]->{samFields}->[1] & 0x80 ) != 0 ) ? 2 : 1;

    $mates[0]{samFields} = ( $fqMateNr_0 == 1 ) ? $mates_tmp[0]{samFields} : $mates_tmp[1]{samFields};
    $mates[1]{samFields} = ( $fqMateNr_1 == 2 ) ? $mates_tmp[1]{samFields} : $mates_tmp[0]{samFields};

    my $readID_core = $mates_tmp[0]{readID};
    $mates[0]{readID}      = $readID_core;    # . "/1";
    $mates[0]{readID_core} = $readID_core;
    $mates[1]{readID}      = $readID_core;    # . "/2";
    $mates[1]{readID_core} = $readID_core;

    $mates[0]{samFields}->[0] = $mates[0]{readID};
    $mates[1]{samFields}->[0] = $mates[1]{readID};

    # done mate sorting

    return (@mates);
}

# processbowtie2samPE($fq_FHref, $lastFQrecRef, $samLref, $mappedSeqsRef, $sam_tmp_FH, $unmappedReadsFHsRef);
sub processbowtie2samPE {
    my $fq_FH            = shift;
    my $lastFQrec        = shift;
    my $sam_rec          = shift;
    my $mappedSeqs       = shift;
    my $samFH            = shift;
    my $unmappedReadsFHs = shift;

    # Make sure that fwd read is idx 0 and rev read is idx 1
    my @mates = sortPE($sam_rec);

    my $readID_core = $mates[0]{readID_core};

    my $fwd_fq_rec_id     = "";
    my $fwd_fq_rec_id_sam = "";
    my $rev_fq_rec_id     = "";
    my $rev_fq_rec_id_sam = "";
    my $unmappedFWDReadsFH;
    my $unmappedREVReadsFH;

    if ( $mates[0]{readID_core} ne $mates[1]{readID_core} ) {
        print STDERR "Not a paired alignment: " . $mates[0]{readID} . " : " . $mates[1]{readID} . "\n";
        return;
    }    # Not a pair

    # if new pair get the orignal FQ record
    unless (    ( ($readID_core) eq $lastFQrec->{FWD}->{readID} )
             && ( ($readID_core) eq $lastFQrec->{REV}->{readID} ) )
    {

        ( $fwd_fq_rec_id, $fwd_fq_rec_id_sam, $unmappedFWDReadsFH ) =
          getOrigFQrec( $mates[0]{readID}, $fq_FH->{FWD}, $lastFQrec->{FWD}, $unmappedReadsFHs->{FWD} );
        ( $rev_fq_rec_id, $rev_fq_rec_id_sam, $unmappedREVReadsFH ) =
          getOrigFQrec( $mates[1]{readID}, $fq_FH->{REV}, $lastFQrec->{REV}, $unmappedReadsFHs->{REV} );

    }


    if (    ( ( $mates[0]->{samFields}->[1] & 0x4 ) != 0 )
            || ( ( $mates[1]->{samFields}->[1] & 0x8 ) != 0 )
            || ( ( $mates[0]->{samFields}->[1] & 0x8 ) != 0 )
            || ( ( $mates[1]->{samFields}->[1] & 0x4 ) != 0 ) )
        {   # unmapped or discordant (which should not happen here)
        if ($bowtie2_un) {
            print $unmappedFWDReadsFH $fwd_fq_rec_id
              . $lastFQrec->{FWD}->{seq} . "\n"
              . $lastFQrec->{FWD}->{readID2}
              . $lastFQrec->{FWD}->{QS};
            print $unmappedREVReadsFH $rev_fq_rec_id
              . $lastFQrec->{REV}->{seq} . "\n"
              . $lastFQrec->{REV}->{readID2}
              . $lastFQrec->{REV}->{QS};
            $mappingStats->{bowtie2unal} += 2;
        }
        $lastFQrec->{mType} = -1;
        return;
    }

    $mates[0]{auxTags} = { samGet_auxTags( @{ $mates[0]->{samFields} } ) };
    $mates[1]{auxTags} = { samGet_auxTags( @{ $mates[1]->{samFields} } ) };

    if ( exists( $mates[0]->{auxTags}->{XS} ) || exists( $mates[1]->{auxTags}->{XS} ) )
    {    # Alternative alignments found, handle multimappers
        $lastFQrec->{mType} = 1;
        return;
    }

    # Filter the alignments
    my $validMates = filterPEalignments( \@mates );

    if ( $validMates < 1 ) {    # no mate passed the filters
        $lastFQrec->{mType} = -1;
        return;
    }

    # lets get the BS seq
    getBSseqPE( \@mates, $lastFQrec );

    # Paranoia: both mates have map to the same transcript
    if ( $validMates == 2 ) {
        if ( $mates[0]->{samFields}->[2] ne $mates[1]->{samFields}->[2] ) {
            print STDERR "Mate pair aligns discordantly to differnt transcripts: " . $readID_core . "\n";
            $lastFQrec->{mType} = -1;
            return;
        }
    }

    if ( $validMates == 1 ) {
        my $seIdx = ( grep { defined( $mates[$_] ) } 0 .. 1 )[0];
        mkReadSE( $mates[$seIdx]->{samFields}, $mates[$seIdx]->{auxTags} );
    }

    if ( $validMates > 0 ) {

        my $mateOverlap = 0;
        my $softClippedBases = [ [ 0, 0 ], [ 0, 0 ] ];
        my $samLineFinal;

        if ( $validMates == 2 ) {
            my $mappedPair = [ $mates[0]->{samFields}, $mates[1]->{samFields} ];

            if ( $bowtie2_mode eq 'local' ) {
                $softClippedBases = getSoftClippedPE($mappedPair);
            }

            my $fwdMateMappedLength =
              length( $mates[0]->{samFields}->[9] ) - $softClippedBases->[0]->[0] - $softClippedBases->[0]->[1];
            my $revMateMappedLength =
              length( $mates[1]->{samFields}->[9] ) - $softClippedBases->[1]->[0] - $softClippedBases->[1]->[1];

            $mateOverlap = $getMateOverlap->( $mappedPair, $fwdMateMappedLength, $revMateMappedLength );
        }

        foreach my $i ( 0 .. 1 ) {
            if ( defined( $mates[$i] ) ) {
                $lastFQrec->{mType} = 0;
                $mappedSeqs->{ $mates[$i]->{samFields}->[2] }->{u} = 1;

                if ( $mateOverlap != 0 ) {
                    $fixMateOverlap->(
                                       $mates[$i]->{samFields},
                                       $mates[$i]->{auxTags},
                                       $mateOverlap, $softClippedBases->[$i], $i
                                       );
                }

                # add GeneTag to Sam record
                my $auxTagGeneName = (defined($idMap{ $mates[$i]->{samFields}->[2] }->{'geneName'})) ?
                                      $idMap{ $mates[$i]->{samFields}->[2] }->{'geneName'} : "unknown";
                splice( @{ $mates[$i]->{samFields} }, 11, 0, ( "ZG:Z:" . $auxTagGeneName ) );



                # add to m-bias plot process queue
                if ( ( $mates[$i]->{samFields}->[1] & 0x10 ) == 0 ) {
                    $mbq->enqueue( $mates[$i]->{read}, $mates[$i]->{QS}, $softClippedBases->[$i]->[0], $softClippedBases->[$i]->[1], $i );
                }
                else {
                    $mbq->enqueue( $mates[$i]->{read}, $mates[$i]->{QS}, $softClippedBases->[$i]->[1], $softClippedBases->[$i]->[0], $i );
                }

                # print $samFH join( "\t", @{ $mates[$i]->{samFields} } );
                $samLineFinal->[$i] = [ @{ $mates[$i]->{samFields} } ];

            }
        }

        if ( $mateOverlap != 0 ) {
            $samLineFinal->[0]->[7] = $samLineFinal->[1]->[3];
            $samLineFinal->[1]->[7] = $samLineFinal->[0]->[3];
        }

        foreach my $line (@$samLineFinal) {
            if ( defined($line) ) {
                print $samFH join( "\t", @$line );
                $mappingStats->{uniqGeneUniqTranscript}++;
            }
        }

    }    #end valid mate
}

sub filterPEalignments {
    my $mates = shift;

    my $validMates = 2;

    my $i = 0;
    for my $m ( @{$mates} ) {
        if ($m) {
            if ( !checkMappingFLAG( $m->{samFields}->[1] ) ) {
                &$debug( $m->{readID} . " checkFLAG filtered" );
                $mappingStats->{flagFilter}++;
                $mates->[$i] = undef;
                $validMates--;
            }
        }
        $i++;
    }

    $i = 0;
    foreach my $m ( @{$mates} ) {
        if ($m) {
            my $editDist   = $m->{auxTags}->{NM}->{v};
            my $readLength = length( $m->{samFields}->[9] );
            if ( checkEditDist( $editDist, $readLength, $m->{readID} ) > 0 ) {
                &$debug( $m->{readID} . " checkEditDist filtered" );
                $mates->[$i] = undef;
                $validMates--;
            }
        }
        $i++;
    }

    return ($validMates);
}

sub getBSseqPE {
    my $mates     = shift;
    my $lastFQrec = shift;

    my $i = 0;
    foreach my $m ( @{$mates} ) {
        if ($m) {
            my $readDir = ( $i == 0 ) ? "FWD" : "REV";
            my $flag    = $m->{samFields}->[1];
            my $read    = $lastFQrec->{$readDir}->{seq};
            if ( ( $flag & 0x10 ) != 0 ) {
                $m->{samFields}->[9] = reverseComp($read);
            }
            else {
                $m->{samFields}->[9] = $read;
            }
            $m->{read} = $read;
            $m->{QS}   = $lastFQrec->{$readDir}->{QS};
        }
        $i++;
    }
}

# process_multimappersPE(\@multiMappers, $lastFQrecRef, $mappedSeqsRef, $sam_tmp_FH, $sam_MM_tmp_FH, $aMM_FH);
sub process_multimappersPE {
    my $multiMappers = shift;    # arrayRef
    my $lastFQrec    = shift;
    my $mappedSeqs   = shift;
    my $samFH        = shift;
    my $samMMFH      = shift;
    my $aMMFH        = shift;

    my ( $mgM, $mtpM, $uGmt, $bl ) = ( 0, 0, 0, 0 );

    my @bestAlignment;
    $bestAlignment[0] = {
                          geneName  => '',
                          ltr       => '',
                          t_length  => 0,
                          AS        => -1000,
                          ASsum     => -2000,
                          samFields => [],
                          auxTags   => {},
                          read      => '',
                          QS        => '',
                          };
    $bestAlignment[1] = {
                          geneName  => '',
                          ltr       => '',
                          t_length  => 0,
                          AS        => -1000,
                          ASsum     => -2000,
                          samFields => [],
                          auxTags   => {},
                          read      => '',
                          QS        => '',
                          };

    my %ambiguosMultiMappers;
    my $i = 0;

    foreach my $samLinePair ( @{$multiMappers} ) {

        # sort the mates
        my @mates = sortPE($samLinePair);

        # lets get the BS seq
        getBSseqPE( \@mates, $lastFQrec );

        my $samFields_0 = $mates[0]->{samFields};
        my $samFields_1 = $mates[1]->{samFields};
        my $seqID;
        my $seqID_0 = $samFields_0->[2];
        my $seqID_1 = $samFields_1->[2];
        if ( $seqID_0 ne $seqID_1 ) {

            # discordant: read pair maps to different transcripts
            next;
        }
        else {
            $seqID = $seqID_0;
        }

        my $geneName = (defined($idMap{$seqID}->{'geneName'})) ? $idMap{$seqID}->{'geneName'} : "unknown";
        my $t_length = $idMap{$seqID}->{'t_length'};

        # add GeneTag to Sam record
        splice( @{$samFields_0}, 11, 0, ( "ZG:Z:" . $geneName ) );
        splice( @{$samFields_1}, 11, 0, ( "ZG:Z:" . $geneName ) );

        # get aux tags
        my %auxTags_0 = samGet_auxTags( @{ $mates[0]->{samFields} } );
        my %auxTags_1 = samGet_auxTags( @{ $mates[1]->{samFields} } );
        $mates[0]{auxTags} = {%auxTags_0};
        $mates[1]{auxTags} = {%auxTags_1};

        my $ASsum = $auxTags_0{AS}{v} + $auxTags_1{AS}{v};

        my @thisAlignment;
        $thisAlignment[0] = {
                              geneName  => $geneName,
                              ltr       => $seqID,
                              t_length  => $t_length,
                              AS        => $auxTags_0{AS}{v},
                              ASsum     => $ASsum,
                              samFields => [ @{$samFields_0} ],
                              auxTags   => {%auxTags_0},
                              read      => $mates[0]->{read},
                              QS        => $mates[0]->{QS},
                              };

        $thisAlignment[1] = {
                              geneName  => $geneName,
                              ltr       => $seqID,
                              t_length  => $t_length,
                              AS        => $auxTags_1{AS}{v},
                              ASsum     => $ASsum,
                              samFields => [ @{$samFields_1} ],
                              auxTags   => {%auxTags_1},
                              read      => $mates[1]->{read},
                              QS        => $mates[1]->{QS},
                              };

        # Filter the alignments
        my $validMates = filterPEalignments( \@mates );

        if ( $validMates < 1 ) {    # no mate passed the filters
            next;
        }

        my @validMatesIdx = grep { defined( $mates[$_] ) } 0 .. 1;
        my $vTA = $thisAlignment[ $validMatesIdx[0] ];

        my @invalidMatesIdx = grep { !defined( $mates[$_] ) } 0 .. 1;
        if ( defined( $invalidMatesIdx[0] ) ) {
            $thisAlignment[ $invalidMatesIdx[0] ] = undef;
            mkReadSE( $vTA->{samFields}, $vTA->{auxTags} );
        }

        my @validBAidx = grep { defined( $bestAlignment[$_] ) } 0 .. 1;
        my $vBA = $bestAlignment[ $validBAidx[0] ];

        if ( $samLinePair == 0 ) {    # first alignment
            @bestAlignment = @thisAlignment;
        }
        elsif ( $vTA->{ASsum} > $vBA->{ASsum} ) {
            @bestAlignment = @thisAlignment;
            $bl            = 0;
        }
        elsif ( ( $vBA->{geneName} ne $vTA->{geneName} ) && ( $vTA->{ASsum} == $vBA->{ASsum} ) ) {

            # multi gene mapper
            if ( defined( $mates[0] ) ) {
                $ambiguosMultiMappers{mgM}->{$samLinePair}->[0]->{ambData} =
                  [ "MG", $lastFQrec->{FWD}->{readID}, $samFields_0->[2], $vTA->{geneName}, $vBA->{geneName} ];
                $ambiguosMultiMappers{mgM}->{$samLinePair}->[0]->{samFields} = $thisAlignment[0]->{samFields};
            }
            if ( defined( $mates[1] ) ) {
                $ambiguosMultiMappers{mgM}->{$samLinePair}->[1]->{ambData} =
                  [ "MG", $lastFQrec->{REV}->{readID}, $samFields_1->[2], $vTA->{geneName}, $vBA->{geneName} ];
                $ambiguosMultiMappers{mgM}->{$samLinePair}->[1]->{samFields} = $thisAlignment[1]->{samFields};
            }
            $bl = 1;
        }
        elsif ( ( $vBA->{ltr} eq $vTA->{ltr} ) && ( $vTA->{ASsum} == $vBA->{ASsum} ) ) {

            # multi Mapper on Transcript
            if ( defined( $mates[0] ) ) {
                $ambiguosMultiMappers{mtpM}->{$samLinePair}->[0]->{ambData} = [
                                                      "MP",              $lastFQrec->{FWD}->{readID}, $samFields_0->[2],
                                                      $samFields_0->[3], $vBA->{samFields}->[3],      $vTA->{geneName}
                                                      ];
                $ambiguosMultiMappers{mtpM}->{$samLinePair}->[0]->{samFields} = $thisAlignment[0]->{samFields};
            }
            if ( defined( $mates[1] ) ) {
                $ambiguosMultiMappers{mtpM}->{$samLinePair}->[1]->{ambData} = [
                                                      "MP",              $lastFQrec->{REV}->{readID}, $samFields_1->[2],
                                                      $samFields_1->[3], $vBA->{samFields}->[3],      $vTA->{geneName}
                                                      ];
                $ambiguosMultiMappers{mtpM}->{$samLinePair}->[1]->{samFields} = $thisAlignment[1]->{samFields};
            }
            $bl = 1;
        }
        else {
            if ( ( $vTA->{t_length} > $vBA->{t_length} ) && ( $vTA->{ASsum} >= $vBA->{ASsum} ) ) {

                # update if maps to same gene longer transcript with same or better AS
                @bestAlignment = @thisAlignment;
                $bl            = 0;
            }
        }
        $i++;
    }

    # No vaild alignment found
    return if ( $i == 0 );

    if ( $bl == 1 ) {
        foreach my $idx ( keys %{ $ambiguosMultiMappers{mgM} } ) {
            if ($aMMFH) {
                if ( defined( $ambiguosMultiMappers{mgM}->{$idx}->[0] ) ) {
                    print $aMMFH join( "\t", @{ $ambiguosMultiMappers{mgM}->{$idx}->[0]->{ambData} } ) . "\n";
                }
                if ( defined( $ambiguosMultiMappers{mgM}->{$idx}->[1] ) ) {
                    print $aMMFH join( "\t", @{ $ambiguosMultiMappers{mgM}->{$idx}->[1]->{ambData} } ) . "\n";
                }
            }
            if ($saveMultiMapperSeparately) {
                if ( defined( $ambiguosMultiMappers{mgM}->{$idx}->[0] ) ) {
                    print $samMMFH join( "\t", @{ $ambiguosMultiMappers{mgM}->{$idx}->[0]->{samFields} } );
                    $mappedSeqs->{ $ambiguosMultiMappers{mgM}->{$idx}->[0]->{samFields}->[2] }->{m} = 1;
                }
                if ( defined( $ambiguosMultiMappers{mgM}->{$idx}->[1] ) ) {
                    print $samMMFH join( "\t", @{ $ambiguosMultiMappers{mgM}->{$idx}->[1]->{samFields} } );
                    $mappedSeqs->{ $ambiguosMultiMappers{mgM}->{$idx}->[1]->{samFields}->[2] }->{m} = 1;
                }

            }

            # Clean me up, could be nicer :-)
            my $mgM_ = ( defined( $ambiguosMultiMappers{mgM}->{$idx}->[0] ) ) ? 1 : 0;
            $mgM_ += ( defined( $ambiguosMultiMappers{mgM}->{$idx}->[1] ) ) ? 1 : 0;
            $mgM = ( $mgM_ > $mgM ) ? $mgM_ : $mgM;
        }
        foreach my $idx ( keys %{ $ambiguosMultiMappers{mtpM} } ) {
            if ($aMMFH) {
                if ( defined( $ambiguosMultiMappers{mtpM}->{$idx}->[0] ) ) {
                    print $aMMFH join( "\t", @{ $ambiguosMultiMappers{mtpM}->{$idx}->[0]->{ambData} } ) . "\n";
                }
                if ( defined( $ambiguosMultiMappers{mtpM}->{$idx}->[1] ) ) {
                    print $aMMFH join( "\t", @{ $ambiguosMultiMappers{mtpM}->{$idx}->[1]->{ambData} } ) . "\n";
                }
            }
            if ($saveMultiMapperSeparately) {
                if ( defined( $ambiguosMultiMappers{mtpM}->{$idx}->[0] ) ) {
                    print $samMMFH join( "\t", @{ $ambiguosMultiMappers{mtpM}->{$idx}->[0]->{samFields} } );
                    $mappedSeqs->{ $ambiguosMultiMappers{mtpM}->{$idx}->[0]->{samFields}->[2] }->{m} = 1;
                }
                if ( defined( $ambiguosMultiMappers{mtpM}->{$idx}->[1] ) ) {
                    print $samMMFH join( "\t", @{ $ambiguosMultiMappers{mtpM}->{$idx}->[1]->{samFields} } );
                    $mappedSeqs->{ $ambiguosMultiMappers{mtpM}->{$idx}->[1]->{samFields}->[2] }->{m} = 1;
                }

            }
            unless ( $mgM > 0 ) {

                # Clean me up, could be nicer :-)
                my $mtpM_ = ( defined( $ambiguosMultiMappers{mtpM}->{$idx}->[0] ) ) ? 1 : 0;
                $mtpM_ += ( defined( $ambiguosMultiMappers{mtpM}->{$idx}->[1] ) ) ? 1 : 0;
                $mtpM = ( $mtpM_ > $mtpM ) ? $mtpM_ : $mtpM;
            }
        }
    }
    else {

        my $mateOverlap = 0;
        my $softClippedBases = [ [ 0, 0 ], [ 0, 0 ] ];
        my $samLineFinal;

        if ( defined( $bestAlignment[0] ) && defined( $bestAlignment[1] ) ) {
            my $mappedPair = [ $bestAlignment[0]->{samFields}, $bestAlignment[1]->{samFields} ];

            if ( $bowtie2_mode eq 'local' ) {
                $softClippedBases = getSoftClippedPE($mappedPair);
            }

            my $fwdMateMappedLength =
              length( $bestAlignment[0]->{samFields}->[9] ) - $softClippedBases->[0]->[0] - $softClippedBases->[0]->[1];
            my $revMateMappedLength =
              length( $bestAlignment[1]->{samFields}->[9] ) - $softClippedBases->[1]->[0] - $softClippedBases->[1]->[1];

            $mateOverlap = $getMateOverlap->( $mappedPair, $fwdMateMappedLength, $revMateMappedLength );
        }

        if ( defined( $bestAlignment[0] ) ) {

            if ( $mateOverlap != 0 ) {
                $fixMateOverlap->(
                                   $bestAlignment[0]->{samFields},
                                   $bestAlignment[0]->{auxTags},
                                   $mateOverlap, $softClippedBases->[0], 0
                                   );
            }

            # add to m-bias plot process queue
            if ( ( $bestAlignment[0]->{samFields}->[1] & 0x10 ) == 0 ) {
                $mbq->enqueue( $bestAlignment[0]->{read}, $bestAlignment[0]->{QS}, $softClippedBases->[0]->[0], $softClippedBases->[0]->[1], 0 );
            }
            else {
                $mbq->enqueue( $bestAlignment[0]->{read}, $bestAlignment[0]->{QS}, $softClippedBases->[0]->[1], $softClippedBases->[0]->[0], 0 );
            }

            # make primary alignment
            $bestAlignment[0]->{samFields}->[1] -= 256 if ( ( $bestAlignment[0]->{samFields}->[1] & 0x100 ) == 256 );

            # print $samFH join( "\t", @{ $bestAlignment[0]->{samFields} } );
            $samLineFinal->[0] = [ @{ $bestAlignment[0]->{samFields} } ];

            $mappedSeqs->{ $bestAlignment[0]->{samFields}->[2] }->{u} = 1;
        }
        if ( defined( $bestAlignment[1] ) ) {

            if ( $mateOverlap != 0 ) {
                $fixMateOverlap->(
                                   $bestAlignment[1]->{samFields},
                                   $bestAlignment[1]->{auxTags},
                                   $mateOverlap, $softClippedBases->[1], 1
                                   );
            }

            # add to m-bias plot process queue
            if ( ( $bestAlignment[1]->{samFields}->[1] & 0x10 ) == 0 ) {
                $mbq->enqueue( $bestAlignment[1]->{read}, $bestAlignment[1]->{QS}, $softClippedBases->[1]->[0], $softClippedBases->[1]->[1], 1 );
            }
            else {
                $mbq->enqueue( $bestAlignment[1]->{read}, $bestAlignment[1]->{QS}, $softClippedBases->[1]->[1], $softClippedBases->[1]->[0], 1 );
            }

            # make primary alignment
            $bestAlignment[1]->{samFields}->[1] -= 256 if ( ( $bestAlignment[1]->{samFields}->[1] & 0x100 ) == 256 );

            #print $samFH join( "\t", @{ $bestAlignment[1]->{samFields} } );
            $samLineFinal->[1] = [ @{ $bestAlignment[1]->{samFields} } ];

            $mappedSeqs->{ $bestAlignment[1]->{samFields}->[2] }->{u} = 1;
        }

        if ( $mateOverlap != 0 ) {
            $samLineFinal->[0]->[7] = $samLineFinal->[1]->[3];
            $samLineFinal->[1]->[7] = $samLineFinal->[0]->[3];
        }

        foreach my $line (@$samLineFinal) {
            if ( defined($line) ) {
                print $samFH join( "\t", @$line );
                $mappingStats->{uniqGeneMultiTranscript}++;
            }
        }

        $uGmt = 1;
        ( $mgM, $mtpM ) = ( 0, 0 );
    }

    $mappingStats->{mgMapperFilter}  += $mgM;
    $mappingStats->{mtpMapperFilter} += $mtpM;

    $lastFQrec->{mType} = -1;

    return;
}

sub mkReadSE {
    my $samFields = shift;
    my $aux       = shift;

    if ( ( $samFields->[1] & 0x1 ) == 1 ) {
        $samFields->[1] -= 2 if ( ( $samFields->[1] & 0x2 ) == 2 );    # read NOT mapped in proper pair
        $samFields->[1] += 8 if ( ( $samFields->[1] & 0x8 ) != 8 );    # mate unmapped
        $samFields->[1] -= 32 if ( ( $samFields->[1] & 0x20 ) == 32 ); # remove mate reversed
    }
    my $YSf = $aux->{YS}->{pos};
    splice( @{$samFields}, $YSf, 1 ) if ( defined($YSf) );
    $samFields->[6] = "*";                                             # edit RNEXT
    $samFields->[7] = 0;                                               # edit PNEXT
    $samFields->[8] = 0;                                               # edit TLEN

}

sub checkMappingFLAG {
    my $flag = shift;

    my $ok = 0;

    if ( ( $flag & 0x1 ) != 0 ) {                                      # PE
        if ( ( ( $flag & 0x10 ) == 0 ) && ( ( $flag & 0x40 ) != 0 ) ) {    # mate1 not reversed
            $ok = 1;
        }
        if ( ( ( $flag & 0x10 ) != 0 ) && ( ( $flag & 0x80 ) != 0 ) ) {    # mate2 reversed
            $ok = 1;
        }
    }
    elsif ( ( $flag & 0x1 ) == 0 ) {                                       # SE
        if ( ( ( $flag & 0x10 ) == 0 ) ) { $ok = 1; }
    }
    return $ok;
}

sub getSoftClippedPE {
    my $alignedPair = shift;

    my $fwdMate = $alignedPair->[0];
    my $revMate = $alignedPair->[1];

    my $fwdMateSoftClipped = getSoftClippedSE($fwdMate);
    my $revMateSoftClipped = getSoftClippedSE($revMate);

    return ( [ $fwdMateSoftClipped, $revMateSoftClipped ] );
}

sub getSoftClippedSE {
    my $alignment = shift;

    my $softClipped;
    $softClipped->[0] = 0;
    $softClipped->[1] = 0;
        
    if ( $bowtie2_mode eq "local" ) {
        my @tmpC = $alignment->[5] =~ /^(\d+H)?((\d+)S)?((\d+)[MIDNP=X])*((\d+)S)?(\d+H)?$/;
        if ( defined( $tmpC[2] ) ) {
            $softClipped->[0] = $tmpC[2];
        }
        if ( defined( $tmpC[6] ) ) {
            $softClipped->[1] = $tmpC[6];
        }
    }
    return( $softClipped );
}

sub getMateOverlap {
    my $alignedPair         = shift;
    my $fwdMateMappedLength = shift;
    my $revMateMappedLength = shift;

    my $fwdMate = $alignedPair->[0];
    my $revMate = $alignedPair->[1];

    my $overlap = 0;

    # C2T strand
    if ( $fwdMate->[3] < $revMate->[3] ) {
        $overlap = $fwdMate->[3] + $fwdMateMappedLength - $revMate->[3];
    }

    # G2A strand
    elsif ( $fwdMate->[3] > $revMate->[3] ) {
        $overlap = $revMate->[3] + $revMateMappedLength - $fwdMate->[3];
    }

    else {
        $overlap = min( $fwdMateMappedLength, $revMateMappedLength );
    }

    $overlap = ( $overlap > 0 ) ? $overlap : 0;

    return ($overlap);

}

sub fixMateOverlap {
    my $alignment        = $_[0];
    my $auxTags          = $_[1];
    my $mateOverlap      = $_[2];
    my $softClippedBases = $_[3];
    my $dir              = $_[4];

    my $trimL  = int( $mateOverlap / 2 );
    my $trimR  = $mateOverlap - $trimL;
    my $trimIL = 0;                         # trimmed leading I in CIGAR after planned trimL
    my $trimIR = 0;                         # trimmed leading I in CIGAR after planned trimR

    # forward mapped mates
    if ( ( $trimR > 0 ) && ( $dir == 0 ) ) {

        if ($hardClipMateOverlap) {
            $alignment->[9]  = substr( $alignment->[9],  0, -( $trimR + $softClippedBases->[1] ) );
            $alignment->[10] = substr( $alignment->[10], 0, -( $trimR + $softClippedBases->[1] ) );

            my $alignedSeqLengthTrimmed = length( $alignment->[9] ) - $softClippedBases->[0];

            $auxTags->{MD}->{v} = updateMDtag( $auxTags->{MD}->{v}, $alignment->[5], 0, $trimR );
            $auxTags->{NM}->{v} = updateNMtag( $auxTags->{MD}->{v}, $alignedSeqLengthTrimmed );

            ( $alignment->[5], $alignment->[3], $trimIL, $trimIR ) =
              hardClipCigar( $alignment->[5], 0, $trimR, $alignment->[3] );

        }
        else {
            my $alignedSeqLengthSoftClipped = length( $alignment->[9] ) - $softClippedBases->[0] - $softClippedBases->[1] - $trimR;

            $auxTags->{MD}->{v} = updateMDtag( $auxTags->{MD}->{v}, $alignment->[5], 0, $trimR );
            $auxTags->{NM}->{v} = updateNMtag( $auxTags->{MD}->{v}, $alignedSeqLengthSoftClipped );

            ( $alignment->[5], $alignment->[3], $trimIL, $trimIR ) =
              softClipCigar( $alignment->[5], 0, $trimR, $alignment->[3] );

        }

        $softClippedBases->[1] += ( $trimR + $trimIR );
    }

    # reversed mapped mates
    if ( ( $trimL > 0 ) && ( $dir == 1 ) ) {

        if ($hardClipMateOverlap) {
            $alignment->[9]  = substr( $alignment->[9],  ( $trimL + $softClippedBases->[0] ) );
            $alignment->[10] = substr( $alignment->[10], ( $trimL + $softClippedBases->[0] ) );

            my $alignedSeqLengthTrimmed = length( $alignment->[9] ) - $softClippedBases->[1];

            $auxTags->{MD}->{v} = updateMDtag( $auxTags->{MD}->{v}, $alignment->[5], $trimL, 0 );
            $auxTags->{NM}->{v} = updateNMtag( $auxTags->{MD}->{v}, $alignedSeqLengthTrimmed );

            ( $alignment->[5], $alignment->[3], $trimIL, $trimIR ) =
              hardClipCigar( $alignment->[5], $trimL, 0, $alignment->[3] );
        }
        else {
            my $alignedSeqLengthSoftClipped = length( $alignment->[9] ) - $softClippedBases->[0] - $softClippedBases->[1] - $trimL;

            $auxTags->{MD}->{v} = updateMDtag( $auxTags->{MD}->{v}, $alignment->[5], $trimL, 0 );
            $auxTags->{NM}->{v} = updateNMtag( $auxTags->{MD}->{v}, $alignedSeqLengthSoftClipped );

            ( $alignment->[5], $alignment->[3], $trimIL, $trimIR ) =
              softClipCigar( $alignment->[5], $trimL, 0, $alignment->[3] );
        }

        $softClippedBases->[0] += ( $trimL + $trimIL );
    }

    $_[0] = $alignment;
    $_[1] = $auxTags;
    $_[3] = $softClippedBases;

}

sub hardClipCigar {
    my $cigar  = shift;
    my $trimL  = shift;    # Left trim length
    my $trimR  = shift;    # Right trim length
    my $mapPos = shift;

    my $new_cigar;
    my @operationLength_new;
    my @operation_new;

    my $D_count  = 0;
    my $I_count  = 0;
    my $trimIL   = 0;
    my $trimIR   = 0;
    my $trimLori = $trimL;
    my $trimRori = $trimR;

    my $opType;

    my ( $operationLength, $operation ) = tokenizeCigar($cigar);

    #&$debug(Dumper($operation), "Line:", __LINE__);
    #&$debug(Dumper($operationLength), "Line:", __LINE__);

    if ( $#$operationLength != $#$operation ) {
        &$debug( "Invalid CIGAR string:", $cigar, "Line:", __LINE__ );
        return ( $cigar, $mapPos, 0, 0 );
    }

    if ( $trimL > 0 ) {

        foreach my $opIdx ( 0 .. $#$operation ) {
            $opType = $operation->[$opIdx];
            if ( $opIdx == 0 ) {
                push( @operation_new, 'H' );
                $operationLength_new[0] = ( $opType =~ /[HS]/ ) ? ( $operationLength->[0] + $trimL ) : $trimL;
                if ( defined( $operation->[1] ) ) {
                    $operationLength_new[0] =
                        ( $operation->[1] eq 'S' )
                      ? ( $operationLength_new[0] + $operationLength->[1] )
                      : $operationLength_new[0];
                }
            }

            if ( $opType eq 'M' ) {
                $operationLength->[$opIdx] -= $trimL;
                if ( $operationLength->[$opIdx] > 0 ) {
                    push( @operation_new,       $opType );
                    push( @operationLength_new, $operationLength->[$opIdx] );
                    push( @operation_new,       @$operation[ $opIdx + 1 .. $#$operation ] );
                    push( @operationLength_new, @$operationLength[ $opIdx + 1 .. $#$operation ] );
                    last;
                }
                elsif ( $operationLength->[$opIdx] < 0 ) {
                    $trimL = abs( $operationLength->[$opIdx] );
                }
                else {
                    my $opOffset = 1;
                    $trimL = 0;

                    next if ( $operation->[ $opIdx + $opOffset ] eq 'I' );

                    while ( $operation->[ $opIdx + $opOffset ] =~ /[ND]/ ) {
                        $D_count += $operationLength->[ $opIdx + $opOffset ];
                        $opOffset++;
                    }

                    push( @operation_new,       @$operation[ $opIdx + $opOffset .. $#$operation ] );
                    push( @operationLength_new, @$operationLength[ $opIdx + $opOffset .. $#$operation ] );
                    last;
                }

            }
            elsif ( $opType eq 'I' ) {
                $trimIL += $operationLength->[$opIdx];
                $operationLength_new[0] += $operationLength->[$opIdx];
            }
            elsif ( $opType =~ /[DN]/ ) {
                $D_count += $operationLength->[$opIdx];
            }
            elsif ( $opType !~ /[SH]/ ) {
                &$debug( "Unsupported CIGAR operation:", $opType, "Line:", __LINE__ );
                return ( $cigar, $mapPos, $trimIL, $trimIR );
            }
        }

        foreach my $opIdx ( 0 .. $#operation_new ) {
            $new_cigar .= $operationLength_new[$opIdx] . $operation_new[$opIdx];
        }
        $cigar = $new_cigar;

        # adjust the mapping start position
        $mapPos += $trimLori + $D_count - $I_count;

        #&$debug(Dumper(\@operation_new), "Line:", __LINE__);
        #&$debug(Dumper(\@operationLength_new), "Line:", __LINE__);
    }

    # RTRIM
    if ( $trimR > 0 ) {
        ( $operationLength, $operation ) = tokenizeCigar($cigar);

        @operationLength_new = ();
        @operation_new       = ();
        $new_cigar           = "";

        foreach my $opIdx ( reverse 0 .. $#$operation ) {
            $opType = $operation->[$opIdx];
            if ( $opIdx == $#$operation ) {
                unshift( @operation_new, 'H' );
                $operationLength_new[0] =
                  ( $opType =~ /[HS]/ ) ? ( $operationLength->[$#$operation] + $trimR ) : $trimR;
                if ( defined( $operation->[ $#$operation - 2 ] ) ) {
                    $operationLength_new[0] =
                        ( $operation->[ $#$operation - 2 ] eq 'S' )
                      ? ( $operationLength_new[0] + $operationLength->[ $#$operation - 2 ] )
                      : $operationLength_new[0];
                }
            }

            if ( $opType eq 'M' ) {
                $operationLength->[$opIdx] -= $trimR;
                if ( $operationLength->[$opIdx] > 0 ) {
                    unshift( @operation_new,       $opType );
                    unshift( @operationLength_new, $operationLength->[$opIdx] );
                    unshift( @operation_new,       @$operation[ 0 .. $opIdx - 1 ] );
                    unshift( @operationLength_new, @$operationLength[ 0 .. $opIdx - 1 ] );
                    last;
                }
                elsif ( $operationLength->[$opIdx] < 0 ) {
                    $trimR = abs( $operationLength->[$opIdx] );
                }
                else {
                    my $opOffset = 1;
                    $trimR = 0;

                    next if ( $operation->[ $opIdx - $opOffset ] eq 'I' );

                    while ( $operation->[ $opIdx - $opOffset ] =~ /[ND]/ ) {
                        $D_count += $operationLength->[ $opIdx - $opOffset ];
                        $opOffset++;
                    }

                    unshift( @operation_new,       @$operation[ 0 .. $opIdx - $opOffset ] );
                    unshift( @operationLength_new, @$operationLength[ 0 .. $opIdx - $opOffset ] );

                    last;
                }

            }
            elsif ( $opType eq 'I' ) {
                $trimIR += $operationLength->[$opIdx];
                $operationLength_new[0] += $operationLength->[$opIdx];
            }
            elsif ( $opType =~ /[DN]/ ) {
                $D_count += $operationLength->[$opIdx];
            }
            elsif ( $opType !~ /[SH]/ ) {
                &$debug( "Unsupported CIGAR operation:", $opType, "Line:", __LINE__ );
                return ( $cigar, $mapPos, $trimIL, $trimIR );
            }
        }

        foreach my $opIdx ( 0 .. $#operation_new ) {
            $new_cigar .= $operationLength_new[$opIdx] . $operation_new[$opIdx];
        }
        $cigar = $new_cigar;

        #&$debug(Dumper(\@operation_new));
        #&$debug(Dumper(\@operationLength_new));
    }

    return ( $cigar, $mapPos, $trimIL, $trimIR );
}

sub tokenizeCigar {
    my $cigar = shift;

    my @operationLength = split( /\D+/, $cigar );
    my @operation       = split( /\d+/, $cigar );

    shift(@operation);

    return ( [@operationLength], [@operation] );
}

sub softClipCigar {
    my $cigar  = shift;
    my $trimL  = shift;    # Left trim length
    my $trimR  = shift;    # Right trim length
    my $mapPos = shift;

    my $new_cigar;
    my @operationLength_new;
    my @operation_new;

    my $D_count  = 0;
    my $trimIL   = 0;
    my $trimIR   = 0;
    my $trimLori = $trimL;
    my $trimRori = $trimR;

    my $opType;

    my ( $operationLength, $operation ) = tokenizeCigar($cigar);

    #&$debug(Dumper($operation), "Line:", __LINE__);
    #&$debug(Dumper($operationLength), "Line:", __LINE__);

    if ( $#$operationLength != $#$operation ) {
        &$debug( "Invalid CIGAR string:", $cigar, "Line:", __LINE__ );
        return ( $cigar, $mapPos, 0, 0 );
    }

    if ( $trimL > 0 ) {

        my $softClip_opIdx = 0;
        foreach my $opIdx ( 0 .. $#$operation ) {
            $opType = $operation->[$opIdx];

            if ( $opIdx == 0 ) {
                if ( $opType ne 'H' ) {
                    push( @operation_new, 'S' );
                    if ( $opType eq 'S' ) {
                        push( @operationLength_new, ( $operationLength->[$opIdx] + $trimL ) );
                    }
                    else {
                        push( @operationLength_new, $trimL );
                    }
                }
                else {
                    $softClip_opIdx = 1;
                    push( @operation_new, ( 'H', 'S' ) );
                    push( @operationLength_new, ( $operationLength->[$opIdx], $trimL ) );
                    next;
                }
            }

            if ( $opType eq 'M' ) {
                $operationLength->[$opIdx] -= $trimL;
                if ( $operationLength->[$opIdx] > 0 ) {
                    push( @operation_new,       $opType );
                    push( @operationLength_new, $operationLength->[$opIdx] );
                    push( @operation_new,       @$operation[ $opIdx + 1 .. $#$operation ] );
                    push( @operationLength_new, @$operationLength[ $opIdx + 1 .. $#$operation ] );
                    last;
                }
                elsif ( $operationLength->[$opIdx] < 0 ) {
                    $trimL = abs( $operationLength->[$opIdx] );
                }
                else {
                    my $opOffset = 1;
                    $trimL = 0;

                    next if ( $operation->[ $opIdx + $opOffset ] eq 'I' );

                    while ( $operation->[ $opIdx + $opOffset ] =~ /[ND]/ ) {
                        $D_count += $operationLength->[ $opIdx + $opOffset ];
                        $opOffset++;
                    }

                    push( @operation_new,       @$operation[ $opIdx + $opOffset .. $#$operation ] );
                    push( @operationLength_new, @$operationLength[ $opIdx + $opOffset .. $#$operation ] );
                    last;
                }

            }
            elsif ( $opType eq 'I' ) {
                $trimIL += $operationLength->[$opIdx];
                $operationLength_new[$softClip_opIdx] += $operationLength->[$opIdx];
            }
            elsif ( $opType =~ /[DN]/ ) {
                $D_count += $operationLength->[$opIdx];
            }
            elsif ( $opType !~ /[SH]/ ) {
                &$debug( "Unsupported CIGAR operation:", $opType, "Line:", __LINE__ );
                return ( $cigar, $mapPos, $trimIL, $trimIR );
            }
        }

        foreach my $opIdx ( 0 .. $#operation_new ) {
            $new_cigar .= $operationLength_new[$opIdx] . $operation_new[$opIdx];
        }
        $cigar = $new_cigar;

        # adjust the mapping start position
        $mapPos += $trimLori + $D_count;

        #&$debug(Dumper(\@operation_new), "Line:", __LINE__);
        #&$debug(Dumper(\@operationLength_new), "Line:", __LINE__);
    }

    # RTRIM
    if ( $trimR > 0 ) {
        ( $operationLength, $operation ) = tokenizeCigar($cigar);

        @operationLength_new = ();
        @operation_new       = ();
        $new_cigar           = "";

        my $softClip_opIdx = 0;
        foreach my $opIdx ( reverse 0 .. $#$operation ) {
            $opType = $operation->[$opIdx];

            if ( $opIdx == $#$operation ) {
                if ( $opType ne 'H' ) {
                    unshift( @operation_new, 'S' );
                    if ( $opType eq 'S' ) {
                        unshift( @operationLength_new, ( $operationLength->[$opIdx] + $trimR ) );
                    }
                    else {
                        unshift( @operationLength_new, $trimR );
                    }
                }
                else {
                    $softClip_opIdx = 1;
                    unshift( @operation_new, ( 'S', 'H' ) );
                    unshift( @operationLength_new, ( $trimR, $operationLength->[$opIdx] ) );
                    next;
                }
            }

            if ( $opType eq 'M' ) {
                $operationLength->[$opIdx] -= $trimR;
                if ( $operationLength->[$opIdx] > 0 ) {
                    unshift( @operation_new,       $opType );
                    unshift( @operationLength_new, $operationLength->[$opIdx] );
                    unshift( @operation_new,       @$operation[ 0 .. $opIdx - 1 ] );
                    unshift( @operationLength_new, @$operationLength[ 0 .. $opIdx - 1 ] );
                    last;
                }
                elsif ( $operationLength->[$opIdx] < 0 ) {
                    $trimR = abs( $operationLength->[$opIdx] );
                }
                else {
                    my $opOffset = 1;
                    $trimR = 0;

                    next if ( $operation->[ $opIdx - $opOffset ] eq 'I' );

                    while ( $operation->[ $opIdx - $opOffset ] =~ /[ND]/ ) {
                        $D_count += $operationLength->[ $opIdx - $opOffset ];
                        $opOffset++;
                    }

                    unshift( @operation_new,       @$operation[ 0 .. $opIdx - $opOffset ] );
                    unshift( @operationLength_new, @$operationLength[ 0 .. $opIdx - $opOffset ] );
                    last;
                }

            }
            elsif ( $opType eq 'I' ) {
                $trimIR += $operationLength->[$opIdx];
                $operationLength_new[$softClip_opIdx] += $operationLength->[$opIdx];
            }
            elsif ( $opType =~ /[DN]/ ) {
                $D_count += $operationLength->[$opIdx];
            }
            elsif ( $opType !~ /[SH]/ ) {
                &$debug( "Unsupported CIGAR operation:", $opType, "Line:", __LINE__ );
                return ( $cigar, $mapPos, $trimIL, $trimIR );
            }
        }

        foreach my $opIdx ( 0 .. $#operation_new ) {
            $new_cigar .= $operationLength_new[$opIdx] . $operation_new[$opIdx];
        }
        $cigar = $new_cigar;

        #&$debug(Dumper(\@operation_new));
        #&$debug(Dumper(\@operationLength_new));
    }

    return ( $cigar, $mapPos, $trimIL, $trimIR );
}

# updateNMtag($MDtag, $seqLength)
sub updateNMtag {
    my $MDs       = shift;
    my $seqLength = shift;

    my @matches = split( /\D/, $MDs );
    my $delOPS = 0;

    @matches = removeEmpty(@matches);
    my $sumMatches = sum(@matches);

    if ( $MDs =~ /\^/ ) {
        my @deletions    = ();
        my $deletedBases = "";

        push( @deletions, $1 ) while ( $MDs =~ /\^([A-Z]+)/g );
        $deletedBases = join( "", @deletions );
        $delOPS = length($deletedBases);
    }
    my $new_NM = ( $seqLength - $sumMatches + $delOPS );

    return $new_NM;
}

sub removeEmpty {
    my @array = @_;
    @array = grep length($_), @array;
    return @array;
}

# updateMDtag($MDtag, $oriCigar, $trimL, $trimR);
sub updateMDtag {
    my $MDs      = shift;
    my $oriCigar = shift;
    my $trimL    = shift;
    my $trimR    = shift;

    # match count
    my @mcount = split( /\D+/, $MDs );

    # mismatches
    my @mms = split( /\d+/, $MDs );
    shift(@mms);    # delete the first empty field

    # expand the MDs, mismatches toupper, deletions to lower, matches M
    my $maskedExpandedMD = "";
    my $mm;
    for ( my $i = 0 ; $i <= $#mcount ; $i++ ) {
        if ( defined( $mms[$i] ) ) {
            if ( $mms[$i] =~ /^\^/ ) {
                $mm = lc( substr( $mms[$i], 1 ) );
            }
            else {
                $mm = uc( $mms[$i] );
            }
            $maskedExpandedMD .= 'M' x $mcount[$i] . $mm;
        }
        else {
            $maskedExpandedMD .= 'M' x $mcount[$i];
        }
    }

    # split expanded md to an array
    my @comp_maskedExpandedMD = split( //, $maskedExpandedMD );

    # parse the cigar string (original not trimmed)
    my ( $operationLength, $operation ) = tokenizeCigar($oriCigar);

    if ( $#$operationLength != $#$operation ) {
        &$debug( "Invalid CIGAR string:", $oriCigar, "Line:", __LINE__ );
        return ($MDs);
    }

    # MDs trimming on the expanded MDs
    my $opType;
    if ( $trimL > 0 ) {
        foreach my $opIdx ( 0 .. $#$operation ) {
            $opType = $operation->[$opIdx];

            # if ( $opType !~ /[DNSHI]/ ) {
            if ( $opType eq 'M' ) {
                $operationLength->[$opIdx] -= $trimL;
                if ( $operationLength->[$opIdx] > 0 ) {
                    splice( @comp_maskedExpandedMD, 0, $trimL );
                    last;
                }
                else {
                    my $trimRpart = $trimL + $operationLength->[$opIdx];
                    splice( @comp_maskedExpandedMD, 0, $trimRpart );
                    $trimL = abs( $operationLength->[$opIdx] );
                }
            }
        }
    }
    if ( $trimR > 0 ) {
        foreach my $opIdx ( reverse 0 .. $#$operation ) {
            $opType = $operation->[$opIdx];

            #if ( $opType !~ /[DNSHI]/ ) {
            if ( $opType eq 'M' ) {
                $operationLength->[$opIdx] -= $trimR;
                if ( $operationLength->[$opIdx] >= 0 ) {
                    splice( @comp_maskedExpandedMD, -$trimR, $trimR );
                    last;
                }
                else {
                    my $trimRpart = $trimR + $operationLength->[$opIdx];
                    splice( @comp_maskedExpandedMD, -$trimRpart, $trimRpart );
                    $trimR = abs( $operationLength->[$opIdx] );
                }
            }
        }
    }

    # Repair MDs
    # fill the first position with 0 if not a match
    if ( $comp_maskedExpandedMD[0] ne 'M' ) {
        unshift( @comp_maskedExpandedMD, 0 );
    }

    # fill the last position with 0 if not a match
    if ( $comp_maskedExpandedMD[-1] ne 'M' ) {
        push( @comp_maskedExpandedMD, 0 );
    }

    # reconstruct the trimmed MDs
    my $mcount  = 0;
    my $new_MDs = "";
    my $last_op = "";

    foreach my $p (@comp_maskedExpandedMD) {
        if ( $p eq 'M' ) {
            $last_op = 'M';
            ++$mcount;
        }
        elsif ( $p =~ /[acgtn]/ ) {
            if ( $last_op eq 'M' ) {
                $new_MDs .= $mcount . '^' . uc($p);
                $mcount = 0;
            }
            elsif ( $last_op eq 'mm' ) {
                $new_MDs .= '0^' . uc($p);
            }
            else {
                $new_MDs .= uc($p);
            }
            $last_op = "D";
        }
        elsif ( $p =~ /[ACGTN]/ ) {
            if ( $last_op eq 'M' ) {
                $new_MDs .= $mcount . uc($p);
                $mcount = 0;
            }
            else {
                $new_MDs .= '0' . uc($p);
            }
            $last_op = 'mm';
        }
    }

    # Add last match count if last operation was match
    if ( $last_op eq 'M' ) { $new_MDs .= $mcount; $mcount = 0; }

    return ($new_MDs);
}

sub getOrigFQrec {
    my $readID           = shift;
    my $fq_FH            = shift;
    my $lastFQrec        = shift;
    my $unmappedReadsFHs = shift;

    # my $mate             = shift // "";

    my $fq_rec_id     = "";
    my $fq_rec_id_sam = "";

    my %fq_rec;

    my $unmappedReadsFH = $unmappedReadsFHs->[0];

  READ_LOOP:
    while (1) {
        %fq_rec = getFQrec( $fq_FH, 0 );

        last READ_LOOP
          unless ( ( $fq_rec{id} && $fq_rec{seq} && $fq_rec{id2} && $fq_rec{qs} ) );

        $fq_rec_id = $fq_rec_id_sam = $fq_rec{id};
        $lastFQrec->{seq}     = $fq_rec{seq};
        $lastFQrec->{readID2} = $fq_rec{id2};
        $lastFQrec->{QS}      = $fq_rec{qs};

        chomp( $fq_rec_id_sam, $lastFQrec->{seq} );

        if ( $fq_rec_id_sam =~ /_meRanT_FQ_FILE_SWITCH_/ ) {
            shift( @{$unmappedReadsFHs} );
            $unmappedReadsFH = $unmappedReadsFHs->[0];
        }
        else {
            $lastFQrec->{readCounter}++;
        }

        # $fq_rec_id_sam =~ s/([ \t]+.+)?$/$mate/;
        # $fq_rec_id_sam = substr( $fq_rec_id_sam, 1 );

        $fq_rec_id_sam =~ s/([ \t]+.+)?$//;
        $fq_rec_id_sam =~ s/(\/1|\/2)?$//;
        $fq_rec_id_sam = substr( $fq_rec_id_sam, 1 );

        chomp($fq_rec_id_sam);

        if ( $fq_rec_id_sam eq $readID ) {
            last;
        }
        else {
            if ( ( $fq_FH->eof ) || ( ( $lastFQrec->{readCounter} >= $firstNreadsTotal ) && ( $firstNreads != -1 ) ) ) {
                die( "Fastq EOF reached! Orphan read: " . $readID . " : " . $fq_rec_id_sam . "\n" );
            }
        }
    }
    $lastFQrec->{readID} = $fq_rec_id_sam;
    return ( $fq_rec_id, $fq_rec_id_sam, $unmappedReadsFH );
}

sub samGet_auxTags {
    my @samFields = @_;

    my %sam_auxTags;
    my $c = 11;
    foreach my $sf ( @samFields[ 11 .. $#samFields ] ) {
        chomp($sf);

        # TAG:TYPE:VALUE
        my @f = split( ':', $sf );
        $sam_auxTags{ $f[0] }{t}   = $f[1];
        $sam_auxTags{ $f[0] }{v}   = $f[2];
        $sam_auxTags{ $f[0] }{pos} = $c;
        $c++;
    }
    return (%sam_auxTags);
}

sub runSEconv {
    my ( $fqL, $conv, $nrOfReads ) = @_;

    my $i = 0;
    while ( $i <= $#$fqL ) {
        bsconvertFQse( $fqL->[$i], $conv, $nrOfReads );
        $i++;
    }
}

sub runPEconv {
    my ( $FWDfqL, $REVfqL, $nrOfReads ) = @_;

    my $i = 0;
    while ( $i <= $#$FWDfqL ) {
        bsconvertFQpe( $FWDfqL->[$i], $REVfqL->[$i], $nrOfReads );
        $i++;
    }
}

sub bsconvertFQse {
    my $fq       = shift;
    my $conv     = shift;
    my $maxReads = shift;

    my $conversion = undef;
    my $readNr     = 0;

    my $IQCfiltered = 0;
    my $DIRfiltered = 0;

    my $fqFH      = getFQFH($fq);
    my $convfq_FH = undef;

    my %fq_rec;

    if ( $conv eq 'C2T' ) {
        $conversion = sub { $_[0] =~ tr/C/T/; };
        $convfq_FH = $bowtie2_fwd_convfq_FH;
    }
    elsif ( $conv eq 'G2A' ) {
        $conversion = sub { $_[0] =~ tr/G/A/; };
        $convfq_FH = $bowtie2_rev_convfq_FH;
    }
    else {
        die("Unknown conversion, this should not happen!");
    }

    while (1) {

        %fq_rec = getFQrec( $fqFH, $fq );

        last
          unless (    ( $fq_rec{id} && $fq_rec{seq} && $fq_rec{id2} && $fq_rec{qs} )
                   && ( ( $readNr < $maxReads ) || ( $maxReads == -1 ) ) );

        chomp( $fq_rec{id} );

        # not passing quality control
        if ($useIlluminaQC) {
            if ( !checkIlluminaFilterFlag( $fq_rec{id} ) ) {
                $readNr++;
                $IQCfiltered++;
                next;
            }
        }

        #$fq_rec{id} =~ s/(([ \t]+.+)?$)//;
        $fq_rec{id} =~ s/(([ \t]+.+)?$)|(\/(1|2)$)//;

        if ($forceDirectionality) {
            if ( !checkDirectionality( $fq_rec{seq}, 0 ) ) {
                &$debug(
                    "Skipping read: " . substr( $fq_rec{id}, 1, -2 ) . " - directionality unassured: " . $fq_rec{seq} );
                $readNr++;
                $DIRfiltered++;
                next;
            }
        }

        $conversion->( $fq_rec{seq} );

        $convfq_FH->print( $fq_rec{id} . "\n", $fq_rec{seq}, $fq_rec{id2}, $fq_rec{qs} );

        $readNr++;
    }

    $convfq_FH->flush();
    close($fqFH);

    say STDOUT "Skipped " . $DIRfiltered . " reads due to unassured directionality";
    say STDOUT "Filtered " . $IQCfiltered . " reads not passing Illuminal QC" if ($useIlluminaQC);
}

sub bsconvertFQpe {
    my $FWDfq    = shift;
    my $REVfq    = shift;
    my $maxReads = shift;

    my $readNr = 0;

    my $IQCfiltered = 0;
    my $DIRfiltered = 0;

    my $FWDfqFH = getFQFH($FWDfq);
    my $REVfqFH = getFQFH($REVfq);

    my %FWDfq_rec;
    my %REVfq_rec;

    my $C2T = sub { $_[0] =~ tr/C/T/; };
    my $G2A = sub { $_[0] =~ tr/G/A/; };

    while (1) {

        %FWDfq_rec = getFQrec( $FWDfqFH, $FWDfq );
        %REVfq_rec = getFQrec( $REVfqFH, $REVfq );

        last
          unless (
                   (
                        ( $FWDfq_rec{id} && $FWDfq_rec{seq} && $FWDfq_rec{id2} && $FWDfq_rec{qs} )
                     && ( $REVfq_rec{id} && $REVfq_rec{seq} && $REVfq_rec{id2} && $REVfq_rec{qs} )
                   )
                   && ( ( $readNr < $maxReads ) || ( $maxReads == -1 ) )
                   );

        chomp( $FWDfq_rec{id}, $REVfq_rec{id} );

        # not passing quality control
        if ($useIlluminaQC) {
            if ( ( ( !checkIlluminaFilterFlag( $FWDfq_rec{id} ) ) || ( !checkIlluminaFilterFlag( $REVfq_rec{id} ) ) ) )
            {
                $IQCfiltered++;
                $readNr++;
                next;
            }
        }

        # $FWDfq_rec{id} =~ s/(([ \t]+.+)?$)|(\/1$)/\/1/;
        # $REVfq_rec{id} =~ s/(([ \t]+.+)?$)|(\/2$)/\/2/;
        $FWDfq_rec{id} =~ s/(([ \t]+.+)?$)|(\/1$)//;
        $REVfq_rec{id} =~ s/(([ \t]+.+)?$)|(\/2$)//;

        # check if reads are properly paired
        if ( substr( $FWDfq_rec{id}, 0, -2 ) ne substr( $REVfq_rec{id}, 0, -2 ) ) {
            die(   $FWDfq . " and " 
                 . $REVfq
                 . " are not properly paired, you may use pairfq \(S. Evan Staton\) tool to pair and sort the mates" );
        }

        if ($fixDirectionaly) {    # EXPERIMENTAL: fixing directionality
            if ( ( !checkDirectionality( $FWDfq_rec{seq}, 0 ) ) || ( !checkDirectionality( $REVfq_rec{seq}, 1 ) ) ) {
                &$debug( "Fixing read pair: " . substr( $FWDfq_rec{id}, 1, -2 ) . " - directionality problem?" );
                ( %FWDfq_rec, %REVfq_rec ) = ( %REVfq_rec, %FWDfq_rec );
                $DIRfiltered++;
            }
        }
        else {
            if ($forceDirectionality) {
                if ( ( !checkDirectionality( $FWDfq_rec{seq}, 0 ) ) || ( !checkDirectionality( $REVfq_rec{seq}, 1 ) ) )
                {
                    &$debug(   "Skipping read pair: "
                             . substr( $FWDfq_rec{id}, 1, -2 )
                             . " - directionality unassured: "
                             . $FWDfq_rec{seq} . " : "
                             . $REVfq_rec{seq} );
                    $DIRfiltered++;
                    $readNr++;
                    next;
                }
            }
        }

        $C2T->( $FWDfq_rec{seq} );
        $G2A->( $REVfq_rec{seq} );

        $bowtie2_fwd_convfq_FH->print( $FWDfq_rec{id} . "\n", $FWDfq_rec{seq}, $FWDfq_rec{id2}, $FWDfq_rec{qs} );
        $bowtie2_rev_convfq_FH->print( $REVfq_rec{id} . "\n", $REVfq_rec{seq}, $REVfq_rec{id2}, $REVfq_rec{qs} );

        $readNr++;
    }

    $bowtie2_fwd_convfq_FH->flush();
    $bowtie2_rev_convfq_FH->flush();

    close($FWDfqFH);
    close($REVfqFH);

    if ($fixDirectionaly) {    # EXPERIMENTAL: fixing directionality
        say STDOUT "Fixed directionality of " . $DIRfiltered . " read pairs";
    }
    else {
        say STDOUT "Skipped " . $DIRfiltered . " read pairs due to unassured directionality";
    }
    say STDOUT "Filtered " . $IQCfiltered . " read pairs not passing Illuminal QC" if ($useIlluminaQC);

}

sub checkIlluminaFilterFlag {
    my $readID = shift;

    $readID =~ /(.+[ \t])([1|2]:)([Y|N])(:.+)?$/;
    my $res = ( $3 eq "Y" ) ? 0 : 1;
    return ($res);
}

sub checkDirectionality {
    my $seq        = shift;
    my $assumedDir = shift;    # 0 = FWD, 1 = REV

    my $C_count = ( $seq =~ tr/C// );
    my $G_count = ( $seq =~ tr/G// );
    my $A_count = ( $seq =~ tr/A// );
    my $T_count = ( $seq =~ tr/T// );

    if ( $assumedDir == 0 ) {
        if ( ( $C_count > $G_count ) && ( $A_count > $G_count ) && ( $C_count > $T_count ) ) {
            return (0);        # Sequence in wrong direction
        }
    }
    else {
        if ( ( $G_count > $C_count ) && ( $T_count > $C_count ) && ( $G_count > $A_count ) ) {
            return (0);        # Sequence in wrong direction
        }
    }
    return (1);
}

sub getQualityScoreOffset {
    my $fq = shift;

    my $fqFH = getFQFH($fq);
    my $qo   = 0;

  FQRECORDS:
    while ( my %fq_rec = getFQrec( $fqFH, $fq ) ) {

        my @qvals = unpack( '(W)*', $fq_rec{qs} );

        foreach my $qv (@qvals) {
            if ( $qv > 74 ) {
                $qo = 64;
                last;
            }
            elsif ( $qv < 59 ) {
                $qo = 33;
                last;
            }
        }
        last FQRECORDS if $qo;
    }

    $fqFH->close();
    undef($fqFH);

    return $qo;
}

sub getFQrec {
    my $fqFH = shift;
    my $fq = shift // "Fastqfile";

    # READ_4_FQLINES:
    my $readID  = $fqFH->getline();
    my $seq     = $fqFH->getline();
    my $readID2 = $fqFH->getline();
    my $QS      = $fqFH->getline();

    if ( !defined($readID) ) {
        return (
                 id  => undef,
                 seq => undef,
                 id2 => undef,
                 qs  => undef,
                 );
    }
    elsif ( ( $readID =~ /\@_meRanT_FQ_FILE_END_1/ ) ) {
        return (
                 id  => undef,
                 seq => undef,
                 id2 => undef,
                 qs  => undef,
                 );
    }
    elsif ( ( $readID !~ /^\@/ ) || ( $readID2 !~ /^\+/ ) ) { die( $fq . ": Not in proper fastq format " ) }

    return (
             id  => $readID,
             seq => uc($seq),
             id2 => $readID2,
             qs  => $QS,
             );
}

sub getFQFH {
    my $fq = shift;

    my $gzip;

    my ( $fname, $fpath, $fsuffix ) = fileparse( $fq, @suffixlist );
    if ( $fsuffix eq "" ) {
        die(   $fq
             . ": Unknown filetype!\nPlease specify a supported fastq file to convert\n[supported filetypes: "
             . join( ", ", @suffixlist )
             . "]\n" );
    }

    if ( $fsuffix =~ /gz$|gzip$/ ) {
        $gzip = 1;
    }

    my $fqFH_ = IO::File->new;

    if ($gzip) {
        open( $fqFH_, "-|", "zcat $fq" ) || die( $fq . ": " . $! );
    }
    else {
        open( $fqFH_, "<", $fq ) || die( $fq . ": " . $! );
    }

    return ($fqFH_);
}

sub getUnmappedFH {

    my $fq = shift;

    my ( $gzip, $unmapped_fq );

    my ( $fname, $fpath, $fsuffix ) = fileparse( $fq, @suffixlist );

    if ( $fsuffix =~ /gz$|gzip$/ ) {
        $gzip = 1;
    }

    if ($unoutDir) {
        checkDir($unoutDir);
        $unmapped_fq = $unoutDir . "/" . $fname . "_unmapped" . $fsuffix;
    }
    else {
        $unmapped_fq = $outDir . "/" . $fname . "_unmapped" . $fsuffix;
    }

    my $fqFH          = IO::File->new;
    my $unmapped_fqFH = IO::File->new;

    if ($gzip) {
        open( $unmapped_fqFH, "|-", "gzip -c - > $unmapped_fq" ) || die( $unmapped_fq . ": " . $! );
    }
    else {
        open( $unmapped_fqFH, ">", $unmapped_fq ) || die( $unmapped_fq . ": " . $! );
    }

    return ($unmapped_fqFH);
}

# readIDtoGeneMap(\%idMap, $mapFile);
sub readIDtoGeneMap {
    my $idMap   = shift;
    my $mapFile = shift;

    my $seqID;
    my $geneName;

    my $mapfh = IO::File->new;
    open( $mapfh, "<", $mapFile ) || die( $mapFile . ": " . $! );
    while (<$mapfh>) {
        chomp;
        my ( $seqID, $geneName, $len ) = split( '\t', $_ );
        $idMap->{$seqID}->{'geneName'} = $geneName;
        $idMap->{$seqID}->{'t_length'} = $len;
    }
    close($mapfh);
    undef $mapfh;
}

sub fqSpooler {
    my $fifoout     = shift;
    my $fqFs        = shift;
    my $firstNreads = shift;

    my $filesToSpool = scalar @{$fqFs};
    my $lineCounter  = 0;
    my $fileCounter  = 0;
    my $maxLines     = 4 * $firstNreads;

    $fifoout->autoflush(1);
    foreach my $fqF ( @{$fqFs} ) {
        $lineCounter = 0;

        my $gzip = 0;
        my $fq;

        my ( $fname, $fpath, $fsuffix ) = fileparse( $fqF, @suffixlist );

        if ( $fsuffix =~ /gz$|gzip$/ ) {
            $gzip = 1;
        }

        if ($gzip) {
            open( $fq, "-|", "zcat $fqF" ) || die( $fqF . ": " . $! );
        }
        else {
            open( $fq, "<$fqF" ) || die($!);
        }
                
        while (1) {
            if (    ( ( $filesToSpool > 1 ) && ( $fileCounter < $filesToSpool ) )
                 && ( ( $fq->eof() ) || ( ( $lineCounter >= $maxLines ) && ( $firstNreads != -1 ) ) ) )
            {
                print $fifoout "\@_meRanT_FQ_FILE_SWITCH_1\n";
                print $fifoout "_meRanT_FQ_FILE_SWITCH_2\n";
                print $fifoout "\+_meRanT_FQ_FILE_SWITCH_3\n";
                print $fifoout "_meRanT_FQ_FILE_SWITCH_4\n";
                $fileCounter++;
                last;
            }
            elsif ( $fq->eof() ) {
                last;
            }
            elsif ( ( $lineCounter >= $maxLines ) && ( $firstNreads != -1 ) ) {
                last;
            }
            else {
                $lineCounter++;
                print $fifoout $fq->getline();
            }
        }
        $fq->close();
    }

    print $fifoout "\@_meRanT_FQ_FILE_END_1\n";
    print $fifoout "_meRanT_FQ_FILE_END_2\n";
    print $fifoout "\+_meRanT_FQ_FILE_END_3\n";
    print $fifoout "_meRanT_FQ_FILE_END_4\n";

    $fifoout->close();
    undef($fifoout);
    
    return;
}

sub mBiasCounter {

    my %mBiasData = (
                      mBiasReadDataF => [],
                      mBiasDataF     => [],
                      mBiasDataFhq   => [],
                      mBiasReadDataR => [],
                      mBiasDataR     => [],
                      mBiasDataRhq   => [],
                      );

    my @mBiasReadDataF;
    my @mBiasDataF;
    my @mBiasDataFhq;
    my @mBiasReadDataR;
    my @mBiasDataR;
    my @mBiasDataRhq;

    # get read and direction from m-bias process queue
    while ( my ( $read, $qs, $leftClip, $rightClip, $mateNr ) = $mbq->dequeue(5) ) {
        # stop the process we are done
        last if ( !defined($read) || !defined($mateNr) );

        my $pos;
        my $idxLength  = length($read) - 1;
        my $alignStart = $leftClip;
        my $alignEnd   = $idxLength - $rightClip;
        my $offset     = $alignStart;

        my @qualValues = map { $_ - $qsOffset } unpack( '(W)*', $qs );

        if ( $mateNr == 0 ) {
            map { $mBiasReadDataF[$_] += 0 } 0 .. $idxLength;      # avoid undefined values;
            map { $mBiasDataF[$_]     += 0 } 0 .. $idxLength;      # avoid undefined values;
            map { $mBiasDataFhq[$_]   += 0 } 0 .. $idxLength;      # avoid undefined values;
            map { $mBiasReadDataF[$_] += 1 } $alignStart .. $alignEnd;

            while (1) {
                $pos = index( $read, 'C', $offset );
                last if ( ($pos < 0)  || ($pos > $alignEnd) );
                ++$mBiasDataF[$pos];
                ++$mBiasDataFhq[$pos] unless ( $qualValues[$pos] < $mbQSt );
                $offset = $pos + 1;
                last if ( $offset >= $alignEnd );
            }
        }
        else {
            map { $mBiasReadDataR[$_] += 0 } 0 .. $idxLength;      # avoid undefined values;
            map { $mBiasDataR[$_]     += 0 } 0 .. $idxLength;      # avoid undefined values;
            map { $mBiasDataRhq[$_]   += 0 } 0 .. $idxLength;      # avoid undefined values;
            map { $mBiasReadDataR[$_] += 1 } $alignStart .. $alignEnd;
            
            while (1) {
                $pos = index( $read, 'G', $offset );
                last if ( ($pos < 0)  || ($pos > $alignEnd) );
                ++$mBiasDataR[$pos];
                ++$mBiasDataRhq[$pos] unless ( $qualValues[$pos] < $mbQSt );
                $offset = $pos + 1;
                last if ( $offset >= $alignEnd );
            }
        }
    }
    
    %mBiasData = (
                   mBiasReadDataF => [@mBiasReadDataF],
                   mBiasDataF     => [@mBiasDataF],
                   mBiasDataFhq   => [@mBiasDataFhq],
                   mBiasReadDataR => [@mBiasReadDataR],
                   mBiasDataR     => [@mBiasDataR],
                   mBiasDataRhq   => [@mBiasDataRhq],
                   );

    return (%mBiasData);
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

sub reverseComp {
    my $seq = shift;

    my $revcomp = reverse($seq);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;

    return ($revcomp);
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

sub by_number {
    if    ( $a < $b ) { -1 }
    elsif ( $a > $b ) { 1 }
    else              { 0 }
}

sub median {
    my @data = @_;

    my $median;

    my $n   = scalar @data;
    my $mid = int $n / 2;

    my @sorted_values = sort by_number @data;
    if ( @data % 2 ) {
        $median = $sorted_values[$mid];
    }
    else {
        $median = ( $sorted_values[ $mid - 1 ] + $sorted_values[$mid] ) / 2;
    }

    return ($median);
}

# borrowed from List::MoreUtils
sub firstidx (&@) {
    my $f = shift;
    foreach my $i ( 0 .. $#_ ) {
        local *_ = \$_[$i];
        return $i if $f->();
    }
    return -1;
}

sub plot_mBias {
    my $mBiasData     = shift;
    my $mBiasDataHQ   = shift;
    my $mBiasReadData = shift;
    my $mBiasPlotOut  = shift;
    my $direction     = shift;

    my @mCratio =
      map { sprintf( "%.1f", ( ( $mBiasData->[$_] / (($mBiasReadData->[$_] != 0) ? $mBiasReadData->[$_] : 1) ) * 100 ) ) } 0 .. $#$mBiasReadData;

    my @mCratioHQ =
      map { sprintf( "%.1f", ( ( $mBiasDataHQ->[$_] / (($mBiasReadData->[$_] != 0) ? $mBiasReadData->[$_] : 1) ) * 100 ) ) } 0 .. $#$mBiasReadData;


    my $n = scalar @mCratio;
    if ( $n < 1 ) {
        warn("No sufficient data for generating m-bias plot");
        return (0);
    }
    my $median = median(@mCratio);

    # calculate MAD (median absolut deviation)
    my @ad;
    foreach my $xi (@mCratio) {
        my ($adi) = abs( $xi - $median );
        push( @ad, $adi );
    }
    my $mad = median(@ad) + 0.0;

    my @madData = ( $median + 2 * $mad ) x $n;

    my @data = ( [ 1 .. $#$mBiasReadData + 1 ], [@mCratio], [@mCratioHQ], [@madData] );

    # print Dumper(\@data);

    my $maxY = max(@mCratio);

    my $graph = GD::Graph::lines->new( 650, 300 );
    my $gd;

    my $titleStr = 'Cytosine 5 methylation bias plot: ' . $direction . ' reads';

    $graph->set(
        title             => $titleStr,
        x_label           => 'Position in read',
        y_label           => 'methylation ratio [%]',
        x_label_position  => 0.5,
        x_label_skip      => 5,
        x_tick_offset     => 4,
        y_label_skip      => 0.5,
        y_max_value       => $maxY,
        y_min_value       => 0,
        x_min_value       => 1,
        line_width        => 2,
        box_axis          => 1,
        x_ticks           => 1,
        x_last_label_skip => 1,
        long_ticks        => 1,
        fgclr             => '#bbbbbb',
        axislabelclr      => '#333333',
        labelclr          => '#333333',
        textclr           => '#333333',
        legendclr         => '#333333',
        dclrs             => [ '#6daa3a', '#4d7296', '#ff782a', '#3366ff', '#6666ff', '#fc6c19', '#26875d', '#6b56b0' ],
        line_types       => [ 1, 1, 3 ],    # 1=solid, 2=dashed, 3=dotted, 4=dot-dashed
        legend_placement => 'RC',
        transparent      => 0,

        ) or warn $graph->error;

    # Try some common font paths
    GD::Text->font_path('/usr/share/fonts:/usr/share/fonts/liberation:/usr/fonts');
    my $font = [ 'LiberationSans-Regular', 'verdana', 'arial', GD::Font->MediumBold ];

    $graph->set_legend_font( $font, 10 );
    $graph->set_title_font( $font, 14 );
    $graph->set_x_label_font( $font, 12 );
    $graph->set_y_label_font( $font, 12 );
    $graph->set_x_axis_font( $font, 10 );
    $graph->set_y_axis_font( $font, 10 );

    my @r_mCratio    = reverse(@mCratio);
    my $ignoreThresh = $median + 2 * $mad;
    my $fivePignore  = firstidx { $_ <= $ignoreThresh } @mCratio;
    my $threePignore = firstidx { $_ <= $ignoreThresh } @r_mCratio;

    my @legend = ( "m5C ratio", "m5C ratio HQ \(Q >= " . $mbQSt . "\)", "threshold: " . $ignoreThresh );

    $graph->set_legend(@legend);

    $gd = $graph->plot( \@data );
    my $grey = $gd->colorAllocate( 155, 155, 155 );

    my $axes = $graph->get_feature_coordinates('axes');

    my $align = GD::Text::Align->new($gd);

    $align->set_font( $font, 8 );
    $align->set( color => $grey );
    $align->set_align( 'base', 'left' );

    $align->set_text( "3\' ignore: " . $threePignore );
    my @bb = $align->bounding_box( 0, 0, 0 );

    $align->draw( $axes->[3] + 5, $axes->[2], 0 );
    $align->set_text( "5\' ignore: " . $fivePignore );
    $align->draw( $axes->[3] + 5, $axes->[2] + $bb[5] - 2, 0 );

    open( PNG, ">$mBiasPlotOut" ) || die( "Could not write to " . $mBiasPlotOut . " " . $! );
    binmode(PNG);
    print PNG $gd->png();
    close PNG;
    print "\nm-bias plot generation completed.\n";

}

sub getBowtie2version {

    my $versionCheckCmd = $bowtie2_cmd . " --version";
    my $pid = open( my $VERSIONCHECK, "-|", $versionCheckCmd )
      || die( "\n\b\aERROR could not determine the Bowtie2 version: ", $! );

    $VERSIONCHECK->autoflush();
    my $versionStr = $VERSIONCHECK->getline();
    waitpid( $pid, 0 );
    close($VERSIONCHECK);

    if ( !defined($versionStr) ) {
        die( "Could not get Bowtie2 version by running: " . $versionCheckCmd );
    }
    chomp($versionStr);

    $versionStr =~ s/.+version\ //;
    return ($versionStr);
}

sub checkBowtie2 {
    $bowtie2Version = getBowtie2version();
    if ( !exists( $supportedBowtie2Versions{$bowtie2Version} ) ) {
        say STDERR "\nERROR: Bowtie2 version "
          . $bowtie2Version
          . " not supported please install one of the following versions:\n";
        foreach my $sv ( keys(%supportedBowtie2Versions) ) {
            say STDERR $sv;
        }
        exit(1);
    }
    return (1);
}

sub usage {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self <runmode> [-h] [-man] [--version]

Required <runmode> any of:
    mkbsidx     :   Generate the Bowtie2 BS index.
    align       :   Align bs reads to transcripts.

Options:
    --version   :   Print the program version and exit.
    -h|help     :   Print the program help information.
    -man        :   Print a detailed documentation.
    
EOF
}

sub usage_mkbsidx {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self mkbsidx [-fa] [-id] [-h] [-man]

Required all of :
    -fa|fasta           : Fasta file to use for BS index generation.
    -id|bsidxdir        : Directory where to store the BS index.

Options:
    -bwt2b|bowtie2build : Path to Bowtie2 indexer "bowtie2-build".
                          (default: bowtie2_build from the meRanTK installation or
                          your system PATH)
    -t|threads          : number of CPUs/threads to run
    --version           : Print the program version and exit.
    -h|help             : Print the program help information.
    -man                : Print a detailed documentation.
    
EOF
}

sub usage_align {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self align [-f] [-x] [-h] [-man]

Required all of :
    -fastqF|-f            : Fastq file with forward reads (required if no -r)
                            This file must contain the reads that align to the
                            5' end of the RNA, which is the left-most end of the
                            sequenced fragment (in transcript coordinates).
    -fastqR|-r            : Fastq file with reverse reads (required if no -f)
                            This file must contain the reads that align to the
                            opposite strand on the 3' end of the RNA, which is
                            the right-most end of the sequenced fragment (in
                            transcript coordinates).
    -id2gene|-i2g         : Transcript to gene mapping file.
                            This mapping file must in in the following tab
                            delimited format:

                            #seqID  Genesymbol  sequencelength

Options:
    -illuminaQC|-iqc      : Filter reads that did not pass the Illumina QC.
                            Only relevant if you have Illumina 1.8+ reads.
                            (default: not set)

    -forceDir|-fDir       : Filter reads that did not pass did not pass the internal
                            directionality check:
                                FWDreads: #C > #G && #C > #T && #A > #G)
                                REVreads: #G > #C && #T > #C && #G > #A)
                            (default: not set)

    -first|-fn            : Process only this many reads/pairs
                            (default: process all reads/pairs)

    -outdir|-o            : Directory where results get stored
                            (default: current directory)

    -sam|-S               : Name of the SAM file for uniq and resolved alignments
                            (default: meRanT_[timestamp].sam )

    -unalDir|-ud          : Directory where unaligned reads get stored
                            (default: outdir)

                            Note: if -bowtie2un|-un is not set, unaligned reads will
                            not get stored

    -threads|-t           : Use max. this many CPUs to process data
                            (default: 1)

    -bowtie2cmd|-bwt2     : Path to bowtie2
                            (default: bowtie2 from the meRanTK installation or your
                             system PATH)

    -bsidx|-x             : Name of bsindex created in mkbsidx runMode
                            (default: use BS_BWT2IDX environment variable)

    -samMM|-MM            : Save multimappers? If set multimappers will be stored
                            in SAM format '\$sam_multimappers.sam'
                            (default: not set)

    -ommitBAM|-ob         : Do not create an sorted and indexed BAM file
                            (default: not set)

    -deleteSAM|-ds        : Delete the SAM files after conversion to BAM format
                            (default: not set)

    -reportAM|-ra         : Report ambiguos mappings? If set ambiguos mappings
                            will be stored in '\$unalDir/\$sam_ambiguos.txt'
                            (default: not set)

    -bowtie2mode|-m       : Alignment mode. Can either be 'local' or 'end-to-end'
                            See Bowtie2 documentation for more information.
                            (default: end-to-end)

    -max-edit-dist|-e     : Maximum edit distance to allow for a valid alignment
                            (default: 2)

    -max-mm-rate|-mmr     : Maximum mismatch ratio (mismatches over read length)
                            [0 <= mmr < 1]
                            (default: 0.05)

    -mbiasplot|-mbp       : Create an m-bias plot, that shows potentially biased
                            read positions
                            (default: not set)

    -mbiasQS|-mbQS        : Quality score for a high quality m-bias plot. This plot
                            considers only basecalls with a quality score equal or
                            higher than specified by this option.
                            (default: 30)

    -fixMateOverlap|-fmo  : The sequenced fragment and read lengths might be such that
                            alignments for the two mates from a pair overlap each other.
                            If '-fmo' is set, deduplicate alignment subregions that are
                            covered by both, forward and reverse, reads of the same
                            read pair. Only relevant for paired end reads.
                            (default: not set)

    -hardClipMO|-hcmo     : If '-fmo' is set, hardclip instead of softclip the overlaping
                            sequence parts.
                            (default: not set = softclip)

    -bowtie2N|-N          : see Bowtie2 -N option (default: 0)
    -bowtie2L|-L          : see Bowtie2 -L option (default: 20)
    -bowtie2D|-D          : see Bowtie2 -D option (default: 30)
    -bowtie2R|-R          : see Bowtie2 -R option (default: 2)
    
    -bowtie2I|-I          : Minimum fragment length for valid paired-end alignments.
                            see Bowtie2 -I option (default: 0)

    -bowtie2X|-X          : Maximum fragment length for valid paired-end alignments. 
                            see Bowtie2 -X option (default: 1000)

    -min-score            : see Bowtie2 -score-min option
                            (default: 'G,20,8' local, 'L,-0.4,-0.4' end-to-end)

    -bowtie2k|-k          : Max. number of valid alignment to consider in mapping.
                            From these the programs will then choose the one with
                            the best score on the longest transcript of the gene
                            to which it maps unambiguosly.
                            see also Bowtie2 -k option
                            (default: 10)

    -bowtie2un|-un        : report unaligned reads. See also -unalDir|-ud
                           (default: not set)

    --version             : Print the program version and exit.
    -h|-help              : Print the program help information.
    -man                  : Print a detailed documentation.
    
    -debug|-d             : Print some debugging information.
        
EOF
}

__END__

=head1 NAME

meRanT - RNA bisulfite short read mapping to transcriptome

=head1 SYNOPSIS

=head2 Index generation

## Generate a transcriptome database BS index for the aligner

 meRanT mkbsidx -fa mm10.refSeqRNA.fa -id /export/data/mm10/BStranscriptomeIDX

 Generates an index for bisulfite mapping strand specific RNA-BSseq reads to a 
 transcriptome database provided as fasta file (e.g. RNA fasta file from RefSeq: 
 ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.rna.fna.gz)

 The example above assumes that the Bowtie2 index builder command "bowtie2-build"
 is found in the system path "$PATH" or “bowtie2-build” from the meRanTK shipped
 third party programs is used. Alternatively, the path to "bowtie2-build" can be
 specified by using the commandline option "-bwt2b".
 

=head2 Align directed/strand specific RNA-BSseq short reads to a transcriptome

=over 2
 
### Single End reads

 meRanT align \
 -o ./meRanTResult \
 -f ./FastqDir/01.fastq,./FastqDir/02.fastq,./FastqDir/03.fastq \
 -t 32 \
 -k 10 \
 -S RNA-BSseq.sam \
 -un \
 -ud ./meRanTunaligned \
 -ra \
 -MM \
 -i2g ./mm10.refSeqRNA2GeneName.map \
 -x /data/mm10/BSrefSeqIDX/mm10.refSeqRNA.C2T
 
 The command above maps the reads from three different fastq files, separated by
 commas, to all transcript sequences of the mouse refseq databases in "mm10.refSeqRNA.fa",
 using the index created as indicated in the previous section.
 The process for selecting the best alignment to a transcript representing a gene
 requires a transcript to gene mapping file (-i2g) "mm10.refSeqRNA2GeneName.map".
 This mapping file must be in the following tab delimited format:
    
    #seqID  Genesymbol  sequencelength
    [...]
    gi|568933834|ref|XR_376799.1|   Mpv17   1474
    gi|568933835|ref|XR_376800.1|   Mpv17   1301
    gi|568933836|ref|XR_376801.1|   Mpv17   1840
    gi|120444911|ref|NM_011960.2|   Parg    4391
    gi|58331157|ref|NM_017373.3|    Nfil3   2019
    gi|115298679|ref|NM_172673.3|   Frmd5   4218
    [...]
    
 Where each transcript in the transcript database (fasta) is mapped to a Genesymbol.
 

 The mapping process will use (-t) 32 CPUs and search for max. (-k) 10 valid alignments,
 from which the best one will be stored in the (-S) "RNA-BSseq.sam" result file. The
 program will save the unaligned reads (-un) in (-ud) the directory "meRanTunaligned"
 and it will also report ambiguous alignments in a separate tab delimited text file.
 The alignments of multi mapper reads (-MM) will additionally be stored in a separate
 SAM file.
  
 The example above assumes that the Bowtie2 aligner command “bowtie2” is found
 in the systems path “$PATH” or “bowtie2” from the meRanTK shipped third party programs
 is used. Alternatively, the path to “bowtie2” can be specified using command line
 option “-bwt2”.
 
### Paired End reads

 meRanT align \
 -o ./meRanTResult \
 -f ./FastqDir/fwd01-paired.fastq,./FastqDir/fwd02-paired.fastq \
 -r ./FastqDir/rev01-paired.fastq,./FastqDir/rev02-paired.fastq \
 -t 32 \
 -k 10 \
 -S RNA-BSseq.sam \
 -un \
 -ud ./meRanTunaligned \
 -ra \
 -MM \
 -i2g ./mm10.refSeqRNA2GeneName.map

 When using paired end reads, one can specify the forward and reverse reads using
 the commandline options "-f" and "-r" respectively. Multiple files for each read
 direction files can be specified separated by commas. Not only the sort order of
 the forward and reverse reads has to be the same within the fastq files but also
 the order in which one specifies the forward and reverse read fastq files (see
 example above). Note: The paired fastq files may not have unpaired reads. If
 this is the case, one can use for example the "pairfq" (S. Evan Staton) tool to
 pair and sort the mates.

=back

=head1 DEPENDENCIES

 The program requires additional perl modules and depending on your perl installation
 you might need do add the following modules:

 Bio::DB::Sam
 Parallel::ForkManager
 
 These modules should be availble via CPAN or depending on your OS via the package
 manager.
 
 Bio::DB:Sam requires the samtools libraries (version 0.1.10 or higher, version 1.0
 is not compatible with Bio::DB::Sam, yet) and header files in order to compile
 successfully.

 Optional modules for generating m-bias plots:
 GD
 GD::Text::Align
 GD::Graph::lines

 In addition to these Perl modules a working installation of Bowtie2 (>= v.2.2.3)
 is required.

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

=head1 REQUIRED ARGUMENTS

=over 2

=item The runMode must be either 'mkbsidx' or 'align'.
 
  mkbsidx : generate the transcriptome database BS index for the aligner
  align   : align RNA-BSseq reads to the transcriptome database.

=back

=head1 OPTIONS

=head2 general:

=over 2

=item -version

 Print the program version and exit.

=item -h|-help

 Print the program help information.

=item -man

 Print a detailed documentation.

=back

=head2 mkbsidx:

=over 2

=item -fa|-fasta

 Transcript fasta file to use for BS index generation.

=item -id|-bsidxdir

 Directory where to store the BS index.

=item -bowtie2build|-bwt2b

 Path to bowtie2-build (default: use BOWTIE2_BUILD environment variable)

=back

=head2 align:

=over 2

=item -fastqF|-f

 Fastq file with forward reads (required if no -r)
 
=item -fastqR|-r

 Fastq file with reverse reads (required if no -f)

=item -illuminaQC|-iqc

 Filter reads that did not pass the Illumina QC. Only relevant if you have
 Illumina 1.8+ reads. (default: not set)

=item -first|-fn

 Process only this many reads/pairs (default: process all reads/pairs)

=item -outdir|-o

 Directory where results get stored (default: current directory)
 
=item -sam|-S

 Name of the SAM file for uniq and resolved alignments (default: meRanT_[timestamp].sam )

=item -unalDir|-ud

 Directory where unaligned reads get stored (default: outdir)
 Note: if -bowtie2un|-un is set unaligned reads will not get stored

=item -threads|-t

 Use max. this many CPUs to process data (default: 1)

=item -bowtie2cmd|-bwt2

 Path to bowtie2 (default: use BOWTIE2 environment variable)

=item -bsidx|-x

 Path to bsindex created in mkbsidx runMode (default: use BS_BWT2IDX environment variable)

=item -samMM|-MM

 Save multimappers? (default: not set)
 If set multimappers will be stored in SAM format "$sam_multimappers.sam"

=item -ommitBAM|-ob

 Do not create an sorted and indexed BAM file (default: not set)
 
=item -deleteSAM|-ds

 Delete the SAM files after conversion to BAM format (default: not set)

=item -reportAM|-ra

 Report ambiguos mappings (default: not set)
 If set ambiguos mappings will be stored in text format "$unalDir/$sam_ambiguos.txt"

=item -bowtie2mode|-m

 Alignment mode. Can either be 'local' or 'end-to-end' (default: end-to-end)
 See Bowtie2 documentation for more information.
 
=item -max-edit-dist|-e

 Maximum edit distance to allow for a valid alignment (default: 2)

=item -max-mm-rate|-mmr

 Maximum mismatch ratio (mismatches over read length [0 <= mmr < 1]) (default: 0.05)

=item -id2gene|-i2g (required)

 Transcript to gene mapping file.
 This mapping file must in in the following tab delimited format:
    
    #seqID  Genesymbol  sequencelength

=item -mbiasplot|-mbp

 Create an m-bias plot, that shows potentially biased read positions (default: not set)

=item -bowtie2N|-N

 see Bowtie2 -N option (default: 0)

=item -bowtie2L|-L

 see Bowtie2 -L option (default: 20)

=item -bowtie2D|-D

 see Bowtie2 -D option (default: 30)

=item -bowtie2R|-R

 see Bowtie2 -R option (default: 2)

=item -bowtie2I|-I

 The minimum fragment length for valid paired-end alignments.
 see Bowtie2 -I option (default: 0)

=item -bowtie2X|-X

 The maximum fragment length for valid paired-end alignments.
 see Bowtie2 -X option (default: 1000)

=item -min-score

 see Bowtie2 -score-min option (default: 'G,20,8' in local, 'L,-0.4,-0.4' in end-to-end mode)

=item -bowtie2k|-k

 Max. number of valid alignment to consider in mapping. From these the programs
 will then choose the one with the best score on the longest transcript of the
 gene to which it maps unambiguosly.
 
 see also Bowtie2 -k option (default: 10)
 
=item -bowtie2un|-un

 do not report unaligned reads. (default: not set)
 
 See also -unalDir|-ud
            
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
