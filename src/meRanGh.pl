#!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.24.0/bin/perl -w
#
#  meRanGh
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
use IPC::Open3;
use Fcntl qw(SEEK_END SEEK_SET SEEK_CUR :flock);
use File::Basename;
use File::Path qw(remove_tree);
use Cwd qw(abs_path);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use POSIX qw(mkfifo);
use threads;
use Parallel::ForkManager;
use List::Util qw(max min sum);
use Bio::DB::Sam;

my $VERSION = '1.2.1b';
my $version;
my $help;
my $man;
my $self               = basename( $0, () );
my $instDir            = dirname( abs_path($0) );
my $arch               = $Config{archname};
my $extUtilDir         = $instDir . '/extutil/' . $arch;
my $useShippedExtUtils = ( -x $extUtilDir ) ? 1 : 0;

my $DEBUG = 0;

use constant {
               LOG2          => log(2),
               LOG10         => log(10),
               DUMP_INTERVAL => 1_000_000,
               READ_CACHE_SIZE_SE => 40_000,
               READ_CACHE_SIZE_PE => 15_000,
               SAM_CHUNK_SIZE     => 10_000,
               SAM_CHUNK_SIZE_MM  => 1_000,               
               };

# List of supported HISAT2 versions
my %supportedHiSat2Versions = ( 
                                 '2.0.4'  => 1,
                                );
my %notSupportedOpts         = ( 
                                 '2.0.4'  => [],
                                );
my $hisat2Version           = "";

####### Start default options ########
my $ommitBAM      = 0;
my $deleteSAM     = 0;
my $deleteBAMus   = 0;
my $mkBG          = 0;
my $strandedBG    = 1;
my $bgScale       = "";
my $nrOfBGTargets = 0;
my $bgCount       = 0;
my $bgPCTdone     = 0;
my $minBGcov      = 1;

my $max_threads = 1;
my $fastQfwd    = "";
my $fastQrev    = "";
my $firstNreads = -1;
my $mbQSt       = 30;
my $qsOffset;

my $fixDirectionaly     = 0;
my $forceDirectionality = 0;
my $useIlluminaQC       = 0;

my $max_edit_dist = 10;
my $max_mm_ratio  = 0.05;

my $mBiasPlot     = 0;
my $mBiasPlotOutF = "";
my $mBiasPlotOutR = "";
my $mBiasData;
my $mbq;
my $mBiasThr;

my $mappedSeqs = {};

my $fMO;
my $hardClipMateOverlap;
my $allowDoveTail = undef;

# set path to hisat2 for hisat2
$ENV{'PATH'} = ($useShippedExtUtils) ? ( $extUtilDir . '/hisat2:' . $ENV{'PATH'} ) : $ENV{'PATH'};

my $hisat2_cmd = ($useShippedExtUtils) ? $extUtilDir . '/hisat2/hisat2' : 'hisat2';
my $bsIdxDir         = $ENV{BS_HISAT2_IDX} // "";
my $hisat2_C2T_idx  = "";
my $hisat2_G2A_idx  = "";
my $hisat2_C2T_tidx = "";
my $hisat2_G2A_tidx = "";
my $GTF             = "";
my $spliceSitesFile = "";
my $fasta           = "";


my $hisat2_build = ($useShippedExtUtils) ? $extUtilDir . '/hisat2/hisat2-build' : 'hisat2-build';

my $hisat2_un;

our $outDir   = "";
our $lockDir  = ( -x "/tmp" ) ? "/tmp" : \$outDir;
our $unoutDir = "";
our $tmpDir   = ""; # hidden option
my $samout = "";
my $saveMultiMapperSeparately;

# Run + HiSat2 Options
my %run_opts = (
    'fastqF'        => \$fastQfwd,
    'fastqR'        => \$fastQrev,
    'first'         => \$firstNreads,
    'outdir'        => \$outDir,
    'forceDir'      => \$forceDirectionality,
    'sam'           => \$samout,
    'samMM'         => \$saveMultiMapperSeparately,
    'threads'       => \$max_threads,
    'hisat2cmd'     => \$hisat2_cmd,
    'hisat2build'   => \$hisat2_build,
    'bsidxW'        => \$hisat2_C2T_idx,
    'bsidxC'        => \$hisat2_G2A_idx,
    'hisat2un'      => \$hisat2_un,
    'unalDir'       => \$unoutDir,
    'max-edit-dist' => \$max_edit_dist,
    'max-mm-rate'   => \$max_mm_ratio,

#
    'hisat2_trim5' => 0,
    'hisat2_trim3' => 0,
    'hisat2_n-ceil' => undef, # "L,0,0.15",
    'hisat2_ignore-quals' => undef,
    'hisat2_nofw' => undef,
    'hisat2_norc' => undef,
    'hisat2_pen-cansplice' => undef, # 0,
    'hisat2_pen-noncansplice' => undef, # 12,
    'hisat2_pen-canintronlen' => undef, # "G,-8,1",
    'hisat2_pen-noncanintronlen' => undef, # "G,-8,1",
    'hisat2_min-intronlen' => undef, # 20,
    'hisat2_max-intronlen' => undef, # 500000,
    'hisat2_novel-splicesite-outfile' => undef,
    'hisat2_novel-splicesite-infile' => undef,
    'hisat2_no-temp-splicesite' => undef,
    'hisat2_no-spliced-alignment' => undef,
    # 'hisat2_rna-strandness' => "FR", # determined later
    'hisat2_transcriptome-mapping-only' => undef,
    
    'hisat2_remove-chrname' => undef,
    'hisat2_add-chrname' => undef,
    'hisat2_omit-sec-seq' => undef,
    'hisat2_offrate' => undef,
    
    
    'hisat2_ma' => undef,
    'hisat2_mp' => undef, # "6,2",
    'hisat2_sp' => undef, # "1,2",
    'hisat2_np' => undef, # 1,
    'hisat2_rdg' => undef, # "5,3",
    'hisat2_rfg' => undef, # "5,3",
    'hisat2_score-min' => undef, # "L,0.0,-0.2",
    'hisat2_k' => undef,
    
    'hisat2_minins' => undef, # 0,
    'hisat2_maxins' => undef, # 500,
    
    
    'hisat2_dovetail' => \$allowDoveTail,
    
    'hisat2_met-file' => undef

    );
####### End default options ########

my @offONopts = (
                  'dovetail',         'ignore-quals',
                  'nofw', 'norc', 'no-temp-splicesite',
                  'no-spliced-alignment', 'transcriptome-mapping-only',
                  'remove-chrname', 'add-chrname',
                  'omit-seq-seq'
                  );

my @invariableHiSat2Settings = ("--reorder", "--quiet");

my $samoutMM = "";    # filename for ambiguos alignments

my @fqFiles;
my %fwdFL;
my %revFL;

my $fqSpoolerThr;
my $fqFWDspoolerThr;
my $fqREVspoolerThr;

my %mappingStats = (
                     'reads'       => 0,
                     'hisat2unal' => 0,
                     'filteredAL'  => 0,
                     'uniqMapper'  => 0
                     );

my %mappingStatsC = %mappingStats;

my %bg_scaling = (
                   'log2'  => 1,
                   'log10' => 1,
                   );

our @suffixlist = qw(.fastq .fq .fastq.gz .fq.gz .fastq.gzip .fq.gzip);

####### Get commandline options ########
my @SAMhdr_cmd_line = @ARGV;
GetOptions(
    \%run_opts,

    'fastqF|f=s'              => \$fastQfwd,
    'fastqR|r=s'              => \$fastQrev,
    'illuminaQC|iqc'          => \$useIlluminaQC,
    'forceDir|fDir'           => \$forceDirectionality,
    'first|fn=i'              => \$firstNreads,
    'QSoffset|qso=i'          => \$qsOffset,
    'outdir|o=s'              => \$outDir,
    'sam|S=s'                 => \$samout,
    'threads|t=i'             => \$max_threads,
    'hisat2cmd|hisat2=s'      => \$hisat2_cmd,
    'hisat2build|hs2b=s'      => \$hisat2_build,
    'hisat2un|un'             => \$hisat2_un,
    'unalDir|ud=s'            => \$unoutDir,
    'samMM|MM'                => \$saveMultiMapperSeparately,
    'ommitBAM|ob'             => \$ommitBAM,
    'deleteSAM|ds'            => \$deleteSAM,
    'deleteBAMus|dbus'        => \$deleteBAMus,
    'mkbg|bg'                 => \$mkBG,
    'minbgCov|mbgc=i'         => \$minBGcov,
    'bgScale|bgs=s'           => \$bgScale,
    'GTF=s'                   => \$GTF,
    'fasta|fa=s'              => \$fasta,
    'bsidxdir|id=s'           => \$bsIdxDir,
    'mbiasplot|mbp'           => \$mBiasPlot,
    'mbiasQS|mbQS=i'          => \$mbQSt,
    'fixMateOverlap|fmo'      => \$fMO,
    'hardClipMO|hcmo'         => \$hardClipMateOverlap,
    'dovetail|dt'             => \$allowDoveTail,
    'max-edit-dist|e=i'       => \$max_edit_dist,
    'max-mm-rate|mmr=f'       => \$max_mm_ratio,
    'help|h'                  => \$help,
    'man|m'                   => \$man,
    'version'                 => \$version,
    'debug|d'                 => \$DEBUG,

    'hisat2_trim5=i',
    'hisat2_trim3=i',
    'hisat2_n-ceil=s',
    'hisat2_ignore-quals',
    'hisat2_nofw',
    'hisat2_norc',
    'hisat2_pen-cansplice=i',
    'hisat2_pen-noncansplice=i',
    'hisat2_pen-canintronlen=s',
    'hisat2_pen-noncanintronlen=s',
    'hisat2_min-intronlen=i',
    'hisat2_max-intronlen=i',
    'hisat2_novel-splicesite-outfile=s',
    'hisat2_novel-splicesite-infile=s',
    'hisat2_no-temp-splicesite',
    'hisat2_no-spliced-alignment',
    'hisat2_rna-strandness=s',
    'hisat2_transcriptome-mapping-only',
    
    'hisat2_remove-chrname',
    'hisat2_add-chrname',
    'hisat2_omit-sec-seq',
    'hisat2_offrate=i',
    
    
    'hisat2_ma=i',
    'hisat2_mp=s',
    'hisat2_sp=s',
    'hisat2_np=i',
    'hisat2_rdg=s',
    'hisat2_rfg=s',
    'hisat2_score-min=s',
    'hisat2_k=i',
    
    'hisat2_minins=i',
    'hisat2_maxins=i',
    
    'hisat2_met-file=s'

    );
####### End commandline options ########

my $debug;
if ($DEBUG) {
    eval "use Data::Dumper";
    $debug = sub { print STDERR "\n" . join( " ", @_ ) . "\n"; }
}
else {
    $debug = sub { return; };
}

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
    if ( !$bsIdxDir || !$hisat2_cmd || !$fasta ) {
        say STDERR "ERROR: invalid/insufficient arguments";
        usage_mkbsidx();
        exit(1);
    }
    checkHiSat2();
    mkbsidx();
    exit(0);
}
elsif ( $runMode eq 'align' ) {
    if ($help) {
        usage_align();
        exit(0);
    }
    checkHiSat2();
    say STDOUT "Starting BS mapping mode...";

    # continue with main
}
else {
    say STDERR "ERROR: Invalid runMode specified!";
    usage();
    exit(1);
}

#### align mode invoced, so lets continue the main program ####

if ( defined($run_opts{'hisat2_n-ceil'}) && ($run_opts{'hisat2_n-ceil'} !~ /[CLSG]\,\-?\d+(\.\d+)?,\-?\d+(\.\d+)?/) ) {
    warn 'Invalid hisat2_n-ceil: ' . $run_opts{'hisat2_n-ceil'} . ', using "default: L,0,0.15"';
    $run_opts{'hisat2_n-ceil'} = "L,0,0.15";
}

if ( defined($run_opts{'hisat2_pen-canintronlen'}) && ($run_opts{'hisat2_pen-canintronlen'} !~ /[CLSG]\,\-?\d+(\.\d+)?,\-?\d+(\.\d+)?/) ) {
    warn 'Invalid hisat2_pen-canintronlen: ' . $run_opts{'hisat2_pen-canintronlen'} . ', using "default: G,-8,1"';
    $run_opts{'hisat2_pen-canintronlen'} = "G,-8,1";
}

if ( defined($run_opts{'hisat2_pen-noncanintronlen'}) && ($run_opts{'hisat2_pen-noncanintronlen'} !~ /[CLSG]\,\-?\d+(\.\d+)?,\-?\d+(\.\d+)?/) ) {
    warn 'Invalid hisat2_pen-noncanintronlen: ' . $run_opts{'hisat2_pen-noncanintronlen'} . ', using "default: G,-8,1"';
    $run_opts{'hisat2_pen-noncanintronlen'} = "G,-8,1";
}

if ( defined($run_opts{'hisat2_score-min'}) && ($run_opts{'hisat2_score-min'} !~ /[CLSG]\,\-?\d+(\.\d+)?,\-?\d+(\.\d+)?/) ) {
    warn 'Invalid hisat2_score-min: ' . $run_opts{'hisat2_score-min'} . ', using "default: L,0.0,-0.2"';
    $run_opts{'hisat2_score-min'} = "L,0.0,-0.2";
}

if ( defined($run_opts{'hisat2_mp'}) && ($run_opts{'hisat2_mp'} !~ /\d+\,\d+/) ) {
    warn 'Invalid hisat2_mp: ' . $run_opts{'hisat2_mp'} . ', using "default: 2,6"';
    $run_opts{'hisat2_mp'} = "6,2";
}

if ( defined($run_opts{'hisat2_sp'}) && ($run_opts{'hisat2_sp'} !~ /\d+\,\d+/) ) {
    warn 'Invalid hisat2_sp: ' . $run_opts{'hisat2_sp'} . ', using "default: 1,2"';
    $run_opts{'hisat2_sp'} = "1,2";
}

if ( defined($run_opts{'hisat2_rdg'}) && ($run_opts{'hisat2_rdg'} !~ /\d+\,\d+/) ) {
    warn 'Invalid hisat2_rdg: ' . $run_opts{'hisat2_rdg'} . ', using "default: 5,3"';
    $run_opts{'hisat2_rdg'} = "5,3";
}

if ( defined($run_opts{'hisat2_rfg'}) && ($run_opts{'hisat2_rfg'} !~ /\d+\,\d+/) ) {
    warn 'Invalid hisat2_rfg: ' . $run_opts{'hisat2_rfg'} . ', using "default: 5,3"';
    $run_opts{'hisat2_rfg'} = "5,3";
}

if ( ( $max_mm_ratio < 0 ) || ( $max_mm_ratio >= 1 ) ) {
    say STDERR "mismatch ratio should be >= 0 and  < 1";
    exit(1);
}

if ( $allowDoveTail && $fMO ) {
    say STDERR "Allowing dovetail (-dt) reads and removing overlapping sequence from a pair (-fmo) are mutual exclusive options!";
    exit(1);
}

my $alignMode = "Local";

my $samAnalyzer_threads = 0;
if ( $max_threads > 1 ) {
    my $availble_threads = 0;
    my $hisat2_threads     = 1;
    if ( $max_threads >= 3 ) {
        $samAnalyzer_threads = 1;
        $availble_threads = $max_threads - $samAnalyzer_threads - ( $hisat2_threads * 2 );
    }
    $run_opts{'hisat2_threads'} = int( $availble_threads / 2 ) + $hisat2_threads;
}

# Default SE, but lets see later what we have
my $singleEnd = 1;

if ($bgScale) {
    if ( !exists( $bg_scaling{$bgScale} ) ) {
        warn 'Unknown BEDgraph scaling: ' . $bgScale
          . ', will not use scaling for BEDgraph for now (valid scalings: log2, log10)';
        $bgScale = "";
    }
}

if ( !$samout ) {
    my ( $sec, $min, $hour, $mday, $mon, $year ) = localtime(time);
    my $timeStr = $year . $mon . $mday . "-" . $hour . $min . $sec;
    $samout = "meRanGh_" . $timeStr . ".sam";
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

my $hisat2Dir_C2T = $outDir . "/meRanGh_C2T_" . $$ . "_hisat2_out";
my $hisat2Dir_G2A = $outDir . "/meRanGh_G2A_" . $$ . "_hisat2_out";

if ($saveMultiMapperSeparately) {
    my ( $fname, $fpath, $fsuffix ) = fileparse( $samout, qw(\.sam) );
    $samoutMM = $outDir . "/" . $fname . "_multimappers.sam";
}

my $hisat2_fwd_convfq;
my $hisat2_fwd_convfq_FH;
my @FWDfiles;
if ($fastQfwd) {
    @FWDfiles = split( ',', $fastQfwd );
    $hisat2_fwd_convfq = $outDir . "/meRanGh_fwd_reads_C2Tconv_" . $$ . ".fastq";
    unlink($hisat2_fwd_convfq);
    $hisat2_fwd_convfq_FH = IO::File->new( $hisat2_fwd_convfq, O_RDWR | O_CREAT | O_TRUNC )
      || die( $hisat2_fwd_convfq . ": " . $! );
}

my $hisat2_rev_convfq;
my $hisat2_rev_convfq_FH;
my @REVfiles;
if ($fastQrev) {
    @REVfiles = split( ',', $fastQrev );
    $hisat2_rev_convfq = $outDir . "/meRanGh_rev_reads_G2Aconv_" . $$ . ".fastq";
    unlink($hisat2_rev_convfq);
    $hisat2_rev_convfq_FH = IO::File->new( $hisat2_rev_convfq, O_RDWR | O_CREAT | O_TRUNC )
      || die( $hisat2_rev_convfq . ": " . $! );
}
push( @fqFiles, @FWDfiles, @REVfiles );

# got both fwd and rev reads
if ( ($fastQfwd) && ($fastQrev) ) {
    $singleEnd = 0;

    if ( $#FWDfiles != $#REVfiles ) {
        say STDERR "Unequal number of forward and reverse read files in PE align mode";
        usage_align();
        exit(1);
    }

    push( @invariableHiSat2Settings, "--no-discordant" ) unless defined($allowDoveTail);
    push( @invariableHiSat2Settings, "--no-mixed" );
}

my $readCountFactor = scalar @fqFiles;
if ( $readCountFactor < 1 ) {
    say STDERR "Please specify at least one fastq file";
    usage_align();
    exit(1);
}

if ( !$hisat2_cmd ) {
    say STDERR "hisat2 command not sepcified";
    usage_align();
    exit(1);
}
if ($bsIdxDir) {
    $hisat2_C2T_idx = $bsIdxDir . "/C2T/meRanGh_C2T";
    $hisat2_G2A_idx = $bsIdxDir . "/G2A/meRanGh_G2A";
}
else {
    say STDERR "No hisat2 index dir specified";
    usage_align();
    exit(1);
}

if ( !$qsOffset ) {
    say STDOUT "No quality score offset specified, trying to determine from fastq file....";
    $qsOffset = getQualityScoreOffset( $fqFiles[0] );
    if ( $qsOffset > 0 ) {
        say STDOUT "Your fastq files seem to have a quality score offset of " . $qsOffset . " ... using this value";
    }
    else {
        say STDERR "Could not determine quality score offset.";
        exit(1);
    }
    if ( $qsOffset == 64 ) {
        push( @invariableHiSat2Settings, "--phred64" );
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

if ($hisat2_fwd_convfq_FH) {
    close($hisat2_fwd_convfq_FH);
}
if ($hisat2_rev_convfq_FH) {
    close($hisat2_rev_convfq_FH);
}

# correct for multiple fq files
my $firstNreadsTotal = $firstNreads * $readCountFactor;

my @hisat2args;
my $cmd_vals;
if ($GTF) {    
    my ( $fname, $fpath, $fsuffix ) = fileparse( $GTF, qw(\.gtf \.GTF) );
    $spliceSitesFile = $outDir . "/" . $fname . "_splice_sites_" . $$ . ".txt";
    extractSpliceSites();
    $run_opts{'hisat2_known-splicesite-infile'} = $spliceSitesFile;
}

while ( my ( $opt, $val ) = each(%run_opts) ) {
    $cmd_vals = undef;
    if ( $opt =~ /^hisat2_/ ) {
        my $s_opt = substr( $opt, 7 );
        if ( ref($val) eq "ARRAY" ) {
            $cmd_vals = join( " ", @{$val} );
        }
        elsif ( ref($val) eq "SCALAR" ) {
            $cmd_vals = ${$val} if ( defined($val) );
        }
        else {
            $cmd_vals = $val if ( defined($val) );
        }

        if ( grep { /^$s_opt$/ } @{ $notSupportedOpts{$hisat2Version} } ) {
            say STDERR "WARNING: your HISAT2 version ("
              . $hisat2Version
              . ") does not support the option: "
              . $s_opt
              . " ... ignoring it";
            next;
        }
        if ( ( grep { /^$s_opt$/ } @offONopts ) && ( defined($cmd_vals) ) ) {
            push( @hisat2args, "--" . $s_opt );
        }
        elsif ( defined($cmd_vals) && ( $cmd_vals ne "" ) && ( grep { !/^$s_opt$/ } @offONopts ) ) {
            my $switch = ($s_opt eq 'k') ? '-' : '--';
            push( @hisat2args, $switch . $s_opt . " " . $cmd_vals );
        }
    }
}

if ( ( !$hisat2_rev_convfq ) && ($hisat2_fwd_convfq) ) {    # single end fwd
    $singleEnd = 1;
    my $arg = "-U " . $hisat2_fwd_convfq;
    push( @hisat2args,               $arg );
    $arg = "--rna-strandness F";
    push( @hisat2args,               $arg );
}
elsif ( ( !$hisat2_fwd_convfq ) && ($hisat2_rev_convfq) ) {    # single end rev
    $singleEnd = 1;
    my $arg = "-U " . $hisat2_rev_convfq;
    push( @hisat2args,               $arg );
    $arg = "--rna-strandness R";
    push( @hisat2args,               $arg );

}
else {                                                           # paired end
    $singleEnd = 0;
    # attention due to a change in the bowtie2 wrapper script introduced in
    # version 2.2.5 we have to specify the fwd and reverse read without space
    # after "-1" and "-2" parameter respectively.
    my $arg = "-1 " . $hisat2_fwd_convfq;
    push( @hisat2args, $arg );
    $arg = "-2 " . $hisat2_rev_convfq;
    push( @hisat2args, $arg );
    $arg = "--rna-strandness FR";
    push( @hisat2args,               $arg );
    $arg = "--fr";
    push( @hisat2args,               $arg );

}

my @hisat2cmd = ( $hisat2_cmd, @invariableHiSat2Settings, @hisat2args );

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

###############################  sssss
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

    $mBiasData = {
                   mBiasReadDataF => [],
                   mBiasDataF     => [],
                   mBiasDataFhq   => [],
                   mBiasReadDataR => [],
                   mBiasDataR     => [],
                   mBiasDataRhq   => [],
                   };

    $mbq = \&mBiasCounter;

}
else {
    $mbq = sub { return };
}

# Start the HiSat2 mapping process
# using the parameters specified above
runHiSat2( \@hisat2cmd );

say STDOUT "DONE: BS mapping";

if ($mBiasPlot) {

    plot_mBias( $mBiasData->{mBiasDataF},
                $mBiasData->{mBiasDataFhq},
                $mBiasData->{mBiasReadDataF},
                $mBiasPlotOutF, 'forward' );

    if ( $singleEnd == 0 ) {
        plot_mBias( $mBiasData->{mBiasDataR},
                    $mBiasData->{mBiasDataRhq},
                    $mBiasData->{mBiasReadDataR},
                    $mBiasPlotOutR, 'reverse' );
    }

}

###################################### Print out some stats ######################################
if ( $mappingStats{reads} == 0 ) {
    say STDERR "Not enough data for printing statistics! (Line: " . __LINE__ . ")";
}
else {
    my $totalMultiMapped    = $mappingStats{reads} - $mappingStats{uniqMapper} - $mappingStats{hisat2unal};
    my $totalMultiMappedpct = sprintf( "%.2f", $totalMultiMapped * 100 / $mappingStats{reads} );
    my $totalMapped         = $totalMultiMapped + $mappingStats{uniqMapper};
    my $totalMappedpct      = sprintf( "%.2f", $totalMapped * 100 / $mappingStats{reads} );
    my $totalunalpct        = sprintf( "%.2f", $mappingStats{hisat2unal} * 100 / $mappingStats{reads} );
    my $totalDiscard        = ( $mappingStats{reads} - $mappingStats{uniqMapper} );
    my $uniqMapperpct       = sprintf( "%.2f", $mappingStats{uniqMapper} * 100 / $mappingStats{reads} );

    my $col2width = length( $mappingStats{reads} );
    my $f_str     = "%*i\t\(%6.2f%%\)";
    my $f_str_b   = "%*i\t\  ------ \ ";
    my $sl1       = sprintf( $f_str, $col2width, $mappingStats{reads}, 100 );
    my $sl2       = sprintf( $f_str, $col2width, $totalMapped, $totalMappedpct );
    my $sl3       = sprintf( $f_str, $col2width, $mappingStats{hisat2unal}, $totalunalpct );
    my $sl4       = sprintf( $f_str, $col2width, $totalDiscard, ( 100 - $uniqMapperpct ) );
    my $sl5       = sprintf( $f_str, $col2width, $totalMultiMapped, $totalMultiMappedpct );
    my $sl6       = sprintf( $f_str, $col2width, $mappingStats{uniqMapper}, $uniqMapperpct );
    my $sl7       = sprintf( $f_str_b, $col2width, $mappingStats{filteredAL} );
    my $sl8       = sprintf( $f_str_b, $col2width, $mappingStats{discordantAL} );

    print STDERR "

FASTQ input files:\n\t\t" . join( "\n\t\t", @fqFiles ) . "

Total # of reads:                                    " . $sl1 . "
Total # of mapped reads:                             " . $sl2 . "
Total # of unmapped reads:                           " . $sl3 . "
Total # of unmapped + multimapper reads:             " . $sl4 . "
Reads mapping to multiple places on reference:       " . $sl5 . " 
Reads mapped to uniq place on reference:             " . $sl6 . "
----------------------------------------------------------------------------

Total # of alignments filtered with wrong direction: " . $sl7 . "
Total # of discordant alignments filtered: " . $sl8 . "


";
}
###################################### END stats ######################################

cleanup();
if ( $mappingStats{reads} == 0 ) {
    say STDERR "No reads aligned, check input files (Line: " . __LINE__ . ")";
    exit(1);
}

if ($mkBG) {
    if ($ommitBAM) {
        warn("BEDgraph requsted: creation of indexed BAM forced...");
        $ommitBAM = 0;
    }
}

my $bamFile;
my $bamFileMM;
if ( $ommitBAM == 0 ) {
    if ( $mappingStats{uniqMapper} > 0 ) {
        $bamFile = sam_to_bam($samout);
    }
    my $totalMultiMapped = $mappingStats{reads} - $mappingStats{uniqMapper} - $mappingStats{hisat2unal};
    if ( $saveMultiMapperSeparately && ( $totalMultiMapped > 0 ) ) {
        $bamFileMM = sam_to_bam($samoutMM);
    }
    if ($deleteSAM) {
        unlink($samout);
        unlink($samoutMM) if ($saveMultiMapperSeparately);
    }
}

my $sam = undef;    # global SAM object;
if ( $mkBG && ( $mappingStats{uniqMapper} > 0 ) ) {

    # makeBedGraph($bamFile);
    makeBedGraphParallel($bamFile);

    # I guess it does not make much sense to run this for multimappers
    #if($saveMultiMapperSeparately) {
    #	# makeBedGraph($bamFileMM);
    #    makeBedGraphParallel($bamFileMM);
    #}
}
$sam = undef;       # global SAM object;

exit(0);

###################################### SUBROUTINES  ######################################
sub cleanup {
    unlink($hisat2_fwd_convfq) if ($hisat2_fwd_convfq);
    unlink($hisat2_rev_convfq) if ($hisat2_rev_convfq);

    my $hisat2_C2T_files;
    my $hisat2_G2A_files;
    my $hisat2_C2T_genomeFiles;
    my $hisat2_G2A_genomeFiles;
    
    if ($tmpDir) {
        $hisat2_C2T_files = $tmpDir . "/meRanGh_C2T_" . $$ . "_hisat2_";
        $hisat2_G2A_files = $tmpDir . "/meRanGh_G2A_" . $$ . "_hisat2_";
        if ($GTF ne "") {
            $hisat2_C2T_genomeFiles = $tmpDir . "/meRanGh_C2T_" . $$ . "_hisat2__hisat2genome";
            $hisat2_G2A_genomeFiles = $tmpDir . "/meRanGh_G2A_" . $$ . "_hisat2__hisat2genome";
        }
    }
    else {
        $hisat2_C2T_files = $outDir . "/meRanGh_C2T_" . $$ . "_hisat2_";
        $hisat2_G2A_files = $outDir . "/meRanGh_G2A_" . $$ . "_hisat2_";
        if ($GTF ne "") {
            $hisat2_C2T_genomeFiles = $outDir . "/meRanGh_C2T_" . $$ . "_hisat2__hisat2genome";
            $hisat2_G2A_genomeFiles = $outDir . "/meRanGh_G2A_" . $$ . "_hisat2__hisat2genome";
        }

    }

    unlink( glob("$hisat2_C2T_files*") ) unless $DEBUG;
    unlink( glob("$hisat2_G2A_files*") ) unless $DEBUG;

    if ($GTF ne "") {
        unless ($DEBUG) {
            remove_tree( $hisat2_C2T_genomeFiles, $hisat2_G2A_genomeFiles, { error => \my $err } );
            if (@$err) {
                for my $diag (@$err) {
                    my ( $file, $message ) = %$diag;
                    if ( $file eq '' ) {
                        print "general error: $message\n";
                    }
                    else {
                        print "problem unlinking $file: $message\n";
                    }
                }
            }
        }
    }

}

###############################  end sssss


sub mkbsidx {
    my $hisat2_fwd_convfq;
    my $hisat2_fwd_convfq_FH;

    my $bsC2TidxDir   = $bsIdxDir . "/C2T";
    my $bsC2TFastaDir = $bsC2TidxDir;
    my $bsG2AidxDir   = $bsIdxDir . "/G2A";
    my $bsG2AFastaDir = $bsG2AidxDir;

    checkDir($bsIdxDir);
    checkDir($bsC2TidxDir);
    checkDir($bsC2TFastaDir);
    checkDir($bsG2AidxDir);
    checkDir($bsG2AFastaDir);

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

    my %FastaConvFiles = (
                           C2T => [],
                           G2A => []
                           );

    my %bsFastaDirs = (
                        C2T => $bsC2TFastaDir,
                        G2A => $bsG2AFastaDir
                        );
    my %convertors = (
                       C2T => sub { $_[0] =~ tr/Cc/Tt/; },
                       G2A => sub { $_[0] =~ tr/Gg/Aa/; }
                       );

    my $procs = ( $max_threads < ( ( scalar @Fastafiles ) * 2 ) ) ? $max_threads : ( ( scalar @Fastafiles ) * 2 );
    print STDOUT "Running $procs processes for BS index geneartion\n\n";

    my $pm = Parallel::ForkManager->new($max_threads);

    $pm->run_on_finish( \&checkChild );

    # loop through the fasta files
    for my $child ( 0 .. $#Fastafiles ) {
        foreach my $conversion ( keys(%convertors) ) {

            my ( $fname, $fpath, $fsuffix ) = fileparse( $Fastafiles[$child], qr/\.[^.]*$/ );
            my $dotSubfname = $fname;
            $dotSubfname =~ tr/\./\_/;
            my $fastaConv = $bsFastaDirs{$conversion} . "/" . $dotSubfname . "_" . $conversion . ".fa";

            push( @{ $FastaConvFiles{$conversion} }, $fastaConv );

            say STDOUT "BS converting reference: " . $Fastafiles[$child] . " -> " . $fastaConv;

            $pm->start( "bsCONV_" . $conversion . "_" . $Fastafiles[$child] ) and next;

            ### in child ###
            # do the job
            my $returnCode = bsconvertFA( $Fastafiles[$child], $fastaConv, $convertors{$conversion} );

            # finished with this chromosome
            $pm->finish($returnCode);
        }
    }
    $pm->wait_all_children;

    # run the indexer
    for my $conversion ( keys(%FastaConvFiles) ) {
        say STDOUT "Indexing " . $conversion . " converted genome: this may take a while ...";
        $pm->start( "bsIDX_" . $conversion ) and next;
        my $finalFastaConv = concatBSfa( $FastaConvFiles{$conversion}, $conversion );
        _runHiSat2Indexer( $finalFastaConv, $conversion );
        $pm->finish;
    }
    $pm->wait_all_children;

    print STDOUT "\n\n";
    say STDOUT "Done with BS index generation. Your meRanGh BS index was stored under:  " . abs_path($bsIdxDir);
    say STDOUT
"You may use this directoy with the '-id' option in the 'align' runMode, or assign it to the 'BS_HISAT2_IDX' environment varialbe.\n";

    if ($GTF) {
        say STDOUT "You may use the '-ts' option together with '-id "
          . abs_path($bsIdxDir)
          . "' in the 'align' runMode to ailgn to known transcripts as well.\n\n";
    }

}

sub checkChild {
    my ( $pid, $exit_code, $ident ) = @_;
    if ( $ident =~ /^bsCONV_/ ) {
        print "** "
          . substr( $ident, 11 )
          . " just finished "
          . substr( $ident, 7, 3 )
          . " conversion with PID "
          . $pid
          . " and exit code: "
          . $exit_code . "\n";
    }
    elsif ( $ident =~ /^bsIDX_/ ) {
        print "** "
          . substr( $ident, 6 )
          . " index generation finished with PID "
          . $pid
          . " and exit code: "
          . $exit_code . "\n";
    }
    else {
        print "** caught unknown process, with ident: " 
          . $ident 
          . " PID: " 
          . $pid
          . " and exit code: "
          . $exit_code . "\n";
    }

    if ( $exit_code ne 0 ) {
        say STDERR "FAILED building meRanGh index, exit code: " . $exit_code;
        exit($exit_code);
    }
    return ($exit_code);
}

sub checkFastaFiles {
    my $fileList = shift;

    foreach my $file ( @{$fileList} ) {
        die( "Could not read: " . $file ) if ( !-r $file );
    }
    return;
}

sub concatBSfa {
    my $FastaConvFiles = shift;
    my $conversion     = shift;

    my $data_buffer;

    my $finalFastaConv = getFinalFastaName( $FastaConvFiles, $conversion ) . ".fa";
    my $finalFastaConvFH = IO::File->new( $finalFastaConv, O_RDWR | O_CREAT | O_TRUNC )
      || die( $finalFastaConv . ": " . $! );

    my @sortedFnames = sort @{$FastaConvFiles};

    foreach my $FastaConvFile (@sortedFnames) {
        my $FastaConvFileFH = IO::File->new( $FastaConvFile, O_RDONLY | O_EXCL ) || die( $FastaConvFile . ": " . $! );
        $finalFastaConvFH->syswrite($data_buffer) while $FastaConvFileFH->sysread( $data_buffer, 8192 );
        $FastaConvFileFH->close();
        undef($FastaConvFileFH);
        unlink($FastaConvFile);
    }
    $finalFastaConvFH->close();
    undef($finalFastaConvFH);

    return ($finalFastaConv);
}

sub bsconvertFA {
    my $fastaFile = shift;
    my $fastaConv = shift;
    my $convertor = shift;

    unlink($fastaConv);

    my $fastaConv_FH = IO::File->new( $fastaConv, O_RDWR | O_CREAT | O_TRUNC ) || die( $fastaConv . ": " . $! );
    my $fasta_FH     = IO::File->new( $fastaFile, O_RDONLY | O_EXCL )          || die( $fastaFile . ": " . $! );

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

sub _runHiSat2Indexer {
    my $FastaConvFile = shift;
    my $conversion    = shift;

    my $pid;
    my $idxfailed = 0;

    local(*HISAT2_OUT, *HISAT2_ERR);
    
    my $hisat2OUT = \*HISAT2_OUT;
    my $hisat2ERR = \*HISAT2_ERR;

    my $idxThreads = ( $max_threads > 1 ) ? int( $max_threads / 2 ) : 1;    # we run 2 instances

    my $idxoutdir = $bsIdxDir . "/" . $conversion;

    my ( $fname, $fpath, $fsuffix ) = fileparse( $FastaConvFile, ".fa" );
    my $idxName      = $fpath . "/" . $fname;
    my $idxNameShort = $fname;

    my $hisat2buildcmd = $hisat2_build . " -p " . $idxThreads . " -q -o 3 -f " . $FastaConvFile . " " . $idxName;

    &$debug($hisat2buildcmd);

    say STDOUT "Indexing ... it may take a while to create: " . $idxName;

    # fork the HISAT2 index builder processes
    $pid = open3( 0, $hisat2OUT, $hisat2ERR, $hisat2buildcmd ) || die($!);
    waitpid( $pid, 0 );

    while ( defined( my $idxBuildOut = $hisat2OUT->getline() ) ) {
        print STDOUT $idxBuildOut;
    }
    while ( defined( my $idxBuildErr = $hisat2ERR->getline() ) ) {
        print STDERR $idxBuildErr;
        if ( $idxBuildErr =~ /\[FAILED\]/ ) {
            $idxfailed = 1;
        }
    }
    $hisat2OUT->close();
    $hisat2ERR->close();

    print STDOUT "\n\n";
    if ( $idxfailed == 0 ) {
        say STDOUT "Finished building the BS index for " . $conversion . " conversion";
        say STDOUT "BS index dir/name: " . $bsIdxDir;
        say STDOUT "BS " . $conversion . " idx storage location: " . $idxName;
    }
    else {
        say STERR "FAILED to build the BS index for " . $conversion . " conversion";
        say STERR "Please check your input files/parameters and the logs in " . $idxoutdir . "/hisat_out/logs";
        exit(1);
    }
    print STDOUT "\n\n";

}

sub getFinalFastaName {
    my $fastaFiles = shift;
    my $conversion = shift;

    return ( $bsIdxDir . "/" . $conversion . "/meRanGh_" . $conversion );
}


############################## sssssssttttttttttttaaaaaaaaaaaaaaa
sub runHiSat2 {
    my $hisat2cmd = shift;

    my @samRawHrd;
    my $samRawHrdRef = \@samRawHrd;

    my $samout_tmp    = $samout . "_" . $$;
    my $samout_MM_tmp = $samoutMM . "_" . $$;

    my $samFHs;
    $samFHs->[0] = IO::File->new( $samout_tmp, O_WRONLY | O_CREAT | O_TRUNC ) || die( $samout_tmp . ": " . $! );

    $samFHs->[1] = undef;
    if ($saveMultiMapperSeparately) {
        $samFHs->[1] = IO::File->new( $samout_MM_tmp, O_WRONLY | O_CREAT | O_TRUNC )
          || die( $samout_MM_tmp . ": " . $! );
    }
    else {
        $samFHs->[1] = $samFHs->[0];
    }

    my $hisat2Dir_C2T;
    my $hisat2Dir_G2A;
    
    if ($tmpDir) {
        $hisat2Dir_C2T = $tmpDir . "/meRanGh_C2T_" . $$ . "_hisat2_";
        $hisat2Dir_G2A = $tmpDir . "/meRanGh_G2A_" . $$ . "_hisat2_";
    }
    else {
        $hisat2Dir_C2T = $outDir . "/meRanGh_C2T_" . $$ . "_hisat2_";
        $hisat2Dir_G2A = $outDir . "/meRanGh_G2A_" . $$ . "_hisat2_";
    }

    my @hisat2C2Tcmd  = ( @$hisat2cmd, "-x $hisat2_C2T_idx");
    my @hisat2G2Acmd  = ( @$hisat2cmd, "-x $hisat2_G2A_idx");

    # fork the HiSat2 mapper processes
    my $pid_C2T = open( my $HISAT2_C2T, "-|", @hisat2C2Tcmd ) || die($!);
    my $pid_G2A = open( my $HISAT2_G2A, "-|", @hisat2G2Acmd ) || die($!);

    $HISAT2_C2T->autoflush(1);
    $HISAT2_G2A->autoflush(1);

    if ( $singleEnd == 1 ) {    # single end

        my $fqFIFO = $outDir . "/meRanGh_FQ_spooler_" . $$;
        unlink($fqFIFO);
        mkfifo( $fqFIFO, 0700 ) || die("can't mkfifo $fqFIFO: $!");
        my $fq_FH = IO::File->new( $fqFIFO, O_RDWR ) || die( $fqFIFO . ": " . $! );

        my @unmappedReadsFHs;
        foreach my $fqf (@fqFiles) {
            push( @unmappedReadsFHs, getUnmappedFH($fqf) ) if ($hisat2_un);
        }

        # Create and detach a fastq multifile spooler fifo
        my $fqFilesRef    = \@FWDfiles;
        my $readDirection = 1;
        if ($fastQrev) {
            $fqFilesRef    = \@REVfiles;
            $readDirection = -1;
        }

        # Create and detach a fastq multifile spooler fifo
        $fqSpoolerThr = threads->create( \&fqSpooler, $fq_FH, $fqFilesRef, $firstNreads )->detach();

        my $readCounter =
          parseHISAT2se( $HISAT2_C2T, $HISAT2_G2A, $samFHs, $fq_FH, $readDirection, \@unmappedReadsFHs, $samRawHrdRef );

        $fq_FH->close();
        $fq_FH = undef;
        unlink($fqFIFO);

        foreach my $fqfh (@unmappedReadsFHs) {
            $fqfh->flush();
            $fqfh->close();
            undef($fqfh);
        }

        $mappingStats{reads} = $readCounter;

    }
    else {    # paired end

        # mk fifo for forward reads
        my $fqFIFO_fwd = $outDir . "/meRanGh_FQ_spooler_fwd_" . $$;
        unlink($fqFIFO_fwd);
        mkfifo( $fqFIFO_fwd, 0700 ) || die("can't mkfifo $fqFIFO_fwd: $!");
        my $fq_FH_fwd = IO::File->new( $fqFIFO_fwd, O_RDWR ) || die( $fqFIFO_fwd . ": " . $! );

        # mk fifo for reverse reads
        my $fqFIFO_rev = $outDir . "/meRanGh_FQ_spooler_rev_" . $$;
        unlink($fqFIFO_rev);
        mkfifo( $fqFIFO_rev, 0700 ) || die("can't mkfifo $fqFIFO_rev: $!");
        my $fq_FH_rev = IO::File->new( $fqFIFO_rev, O_RDWR ) || die( $fqFIFO_rev . ": " . $! );

        my %fq_FH = ( FWD => $fq_FH_fwd, REV => $fq_FH_rev );
        my $fq_FHref = \%fq_FH;

        my @unmappedFWDReadsFHs;
        foreach my $fqf (@FWDfiles) {
            push( @unmappedFWDReadsFHs, getUnmappedFH($fqf) ) if ($hisat2_un);
        }
        my @unmappedREVReadsFHs;
        foreach my $fqf (@REVfiles) {
            push( @unmappedREVReadsFHs, getUnmappedFH($fqf) ) if ($hisat2_un);
        }
        my %unmappedReadsFHs = ( FWD => [@unmappedFWDReadsFHs], REV => [@unmappedREVReadsFHs] );
        my $unmappedReadsFHsRef = \%unmappedReadsFHs;

        # Create and detach a fastq multifile spooler fifo
        $fqFWDspoolerThr = threads->create( \&fqSpooler, $fq_FH_fwd, \@FWDfiles, $firstNreads )->detach();
        $fqREVspoolerThr = threads->create( \&fqSpooler, $fq_FH_rev, \@REVfiles, $firstNreads )->detach();

        my $readCounter = parseHISAT2pe( $HISAT2_C2T, $HISAT2_G2A, $samFHs, $fq_FHref, $unmappedReadsFHsRef, $samRawHrdRef );

        $fq_FH_fwd->close();
        $fq_FH_rev->close();
        undef($fq_FH_fwd);
        undef($fq_FH_rev);

        unlink($fqFIFO_fwd);
        unlink($fqFIFO_rev);

        foreach my $fqfh ( @unmappedFWDReadsFHs, @unmappedREVReadsFHs ) {
            $fqfh->flush();
            $fqfh->close();
            undef($fqfh);
        }

        $mappingStats{reads} = $readCounter;
    }
    waitpid( $pid_C2T, 0 );
    waitpid( $pid_G2A, 0 );
    $HISAT2_C2T->close();
    $HISAT2_G2A->close();

    $samFHs->[0]->flush();
    $samFHs->[0]->close();
    undef( $samFHs->[0] );

    if ($saveMultiMapperSeparately) {
        $samFHs->[1]->flush();
        $samFHs->[1]->close();
        undef( $samFHs->[1] );
    }

    writeFinalSAM( $samout, $samout_tmp, $samoutMM, $samout_MM_tmp, $samRawHrdRef );

}

sub writeFinalSAM {

    my ( $samout, $samout_tmp, $samoutMM, $samout_MM_tmp, $samRawHrd ) = @_;

    my $data_buffer;
    unlink($samout);

    unlink($samout);
    my $samFH      = IO::File->new( $samout,     O_RDWR | O_CREAT | O_TRUNC ) || die( $samout . ": " . $! );
    my $sam_tmp_FH = IO::File->new( $samout_tmp, O_RDONLY | O_EXCL )          || die( $samout_tmp . ": " . $! );

    my $headerT = 'u';
    if ( !$saveMultiMapperSeparately ) {
        $headerT = 'um';
    }
    writeFinalSAM_( $samFH, $sam_tmp_FH, $samRawHrd, $headerT );

    close($sam_tmp_FH);
    close($samFH);
    $sam_tmp_FH = undef;
    $samFH      = undef;
    unlink($samout_tmp);

    if ($saveMultiMapperSeparately) {
        unlink($samoutMM);
        my $samMMFH = IO::File->new( $samoutMM, O_RDWR | O_CREAT | O_TRUNC ) || die( $samoutMM . ": " . $! );
        my $sam_MM_tmp_FH = IO::File->new( $samout_MM_tmp, O_RDONLY | O_EXCL ) || die( $samout_MM_tmp . ": " . $! );

        writeFinalSAM_( $samMMFH, $sam_MM_tmp_FH, $samRawHrd, 'm' );

        close($sam_MM_tmp_FH);
        close($samMMFH);
        $sam_MM_tmp_FH = undef;
        $samMMFH       = undef;
        unlink($samout_MM_tmp);
    }
}

sub writeFinalSAM_ {
    my $samFH      = shift;
    my $sam_tmp_FH = shift;
    my $samRawHrd  = shift;
    my $type       = shift;

    my $data_buffer;

    foreach my $rawHrdLine ( @{$samRawHrd} ) {
        if ( $rawHrdLine !~ /^\@SQ/ ) {
            $samFH->syswrite($rawHrdLine);
            next;
        }
        my @samHf = split( '\t', $rawHrdLine );
        my $hdrSeqID = ( split( ':', $samHf[1] ) )[1];
        if ( $type eq 'um' ) {
            if ( defined( $mappedSeqs->{$hdrSeqID} ) ) {
                $samFH->syswrite($rawHrdLine);
            }
        }
        else {
            if ( defined( $mappedSeqs->{$hdrSeqID}->{$type} ) ) {
                $samFH->syswrite($rawHrdLine);
            }
        }
    }

    my $SAMpg_line =
      "\@PG\tID:meRanGh\tPN:meRanGh\tVN:" . $VERSION . "\tCL:\"" . $0 . " " . join( ' ', @SAMhdr_cmd_line ) . "\"\n";
    $samFH->syswrite($SAMpg_line);

    $sam_tmp_FH->sysseek( 0, SEEK_SET );
    $samFH->syswrite($data_buffer) while $sam_tmp_FH->sysread( $data_buffer, 8192 );
}

sub processhisat2samHdr {
    my $sam_rec = shift;
    my $rawHrd  = shift;

    if ( $sam_rec =~ /^@/ ) {
        $rawHrd->[ $#$rawHrd + 1 ] = $sam_rec;
        return;
    }
}

sub parseHISAT2se {
    my $samC2Tfh         = shift;
    my $samG2Afh         = shift;
    my $samFHs           = shift;
    my $fq_FH            = shift;
    my $readDirection    = shift;
    my $unmappedReadsFHs = shift;
    my $samRawHrd        = shift;

    # get the RAW SAM header
    my $samheader_C2T;
    my $samheader_G2A;
    while (1) {
        $samheader_C2T = $samC2Tfh->getline();
        if ( $samheader_C2T =~ /^@/ ) {
            processhisat2samHdr( $samheader_C2T, $samRawHrd );
        }
        else {
            last;
        }
    }
    while (1) {
        $samheader_G2A = $samG2Afh->getline();
        if ( $samheader_G2A =~ /^@/ ) {
            next;
        }
        else {
            last;
        }
    }

    # To be checked if valid for all FQ ids
    my $cleanSAMIDre = qr/(\/1|\/2|[ \t]+.+)*$/;

    # get the first alignment to C2T genome
    my $samRecC2T;
    my @samFieldsC2T;
    my $sam_rec_readID_C2T;

    if ( defined( $samRecC2T = $samheader_C2T ) ) {
        chomp($samRecC2T);
        @samFieldsC2T = split( '\t', $samRecC2T );
        $sam_rec_readID_C2T = $samFieldsC2T[0];
        $sam_rec_readID_C2T =~ s/$cleanSAMIDre//;
        undef($samRecC2T);
    }

    # get the first alignment to G2A genome
    my $samRecG2A;
    my @samFieldsG2A;
    my $sam_rec_readID_G2A;

    if ( defined( $samRecG2A = $samheader_G2A ) ) {
        chomp($samRecG2A);
        @samFieldsG2A = split( '\t', $samRecG2A );
        $sam_rec_readID_G2A = $samFieldsG2A[0];
        $sam_rec_readID_G2A =~ s/$cleanSAMIDre//;
        undef($samRecG2A);
    }

    # No valid alignments ???
    if ( ( !$sam_rec_readID_C2T ) && ( !$sam_rec_readID_G2A ) ) {
        say STDERR "No valid alignments found in HISAT2 output";
        exit(1);
    }

    my $readNr        = 0;
    my $fq_rec_id     = "";
    my $fq_rec_id_sam = "";

    my $skipC2Ta  = 1;
    my $skipG2Aa  = 1;
    my $endC2Tsam = 0;
    my $endG2Asam = 0;
    my $endSAM    = 0;

    my $readConv = ( defined($hisat2_rev_convfq) ) ? "G2A" : "C2T";
    my $YR = "YR:Z:" . $readConv;

    my $maxReads = $firstNreads * scalar(@fqFiles);

    my $unmappedReadsFH = $unmappedReadsFHs->[0] if ($hisat2_un);

    my $j = 0;
    my $todo;
    my $ident = 0;

    my $samout_tmp_lockF    = $lockDir . "/meRanGhSM_" . $$ . ".lock";
    my $samout_MM_tmp_lockF = $lockDir . "/meRanGhMM_" . $$ . ".lock";
    my $unal_lockF          = $lockDir . "/meRanGhUA_" . $$ . ".lock";

    my $unalChunk_parent = "";

    my $pm = new Parallel::ForkManager( $samAnalyzer_threads, "/dev/shm" );

    $pm->run_on_finish( \&updateParentData );

    my $cleanIDre = qr/(\/1|\/2|[ \t]+.+)?$/;

  READ_LOOP:
    while (1) {

        my %fq_rec = getFQrec( $fq_FH, 0 );

        last READ_LOOP
          unless (    ( $fq_rec{id} && $fq_rec{seq} && $fq_rec{id2} && $fq_rec{qs} )
                   && ( ( $readNr < $maxReads ) || ( $firstNreads == -1 ) ) );

        $fq_rec_id = $fq_rec_id_sam = $fq_rec{id};

        $fq_rec_id_sam =~ s/$cleanIDre//;
        $fq_rec_id_sam = substr( $fq_rec_id_sam, 1 );

        chomp($fq_rec_id_sam);

        if ( $fq_rec_id_sam =~ /_meRanGh_FQ_FILE_SWITCH_/ ) {
            if ($hisat2_un) {
                $unmappedReadsFH->flush();
                $unmappedReadsFH->close();
                undef($unmappedReadsFH);
                shift(@$unmappedReadsFHs);
                $unmappedReadsFH = $unmappedReadsFHs->[0];
            }
            next READ_LOOP;
        }

        my $foundInC2TSAM = 0;
        my $foundInG2ASAM = 0;
        my $thisRead      = 1;

        my $alignments = undef;
        my ( $C2Ta, $G2Aa ) = ( 0, 0 );

      SAM_LOOP:
        while ( $thisRead && ( $endSAM != 1 ) ) {

            if ( $samC2Tfh->eof() ) { $endC2Tsam = 1; }
            if ( $samG2Afh->eof() ) { $endG2Asam = 1; }

            unless ( $skipC2Ta or $endC2Tsam ) {
                $samRecC2T = $samC2Tfh->getline();
                chomp($samRecC2T);
                @samFieldsC2T = split( '\t', $samRecC2T );
                $sam_rec_readID_C2T = $samFieldsC2T[0];
                $sam_rec_readID_C2T =~ s/$cleanSAMIDre//;
            }

            unless ( $skipG2Aa or $endG2Asam ) {
                $samRecG2A = $samG2Afh->getline();
                chomp($samRecG2A);
                @samFieldsG2A = split( '\t', $samRecG2A );
                $sam_rec_readID_G2A = $samFieldsG2A[0];
                $sam_rec_readID_G2A =~ s/$cleanSAMIDre//;
            }

            if ( defined($sam_rec_readID_C2T) && ( $fq_rec_id_sam eq $sam_rec_readID_C2T ) ) {
                $foundInC2TSAM                           = 1;
                $skipC2Ta                                = 0;
                $alignments->{C2T}->{samFields}->[$C2Ta] = [@samFieldsC2T];
                $C2Ta++;
                undef($sam_rec_readID_C2T);
            }
            else {
                $skipC2Ta = 1;
            }

            if ( defined($sam_rec_readID_G2A) && ( $fq_rec_id_sam eq $sam_rec_readID_G2A ) ) {
                $foundInG2ASAM                           = 1;
                $skipG2Aa                                = 0;
                $alignments->{G2A}->{samFields}->[$G2Aa] = [@samFieldsG2A];
                $G2Aa++;
                undef($sam_rec_readID_G2A);
            }
            else {
                $skipG2Aa = 1;
            }

            if ( $skipC2Ta + $skipG2Aa == 2 ) {
                $thisRead = 0;
            }

            if (    ( $endC2Tsam + $endG2Asam == 2 )
                 && ( !defined($sam_rec_readID_C2T) )
                 && ( !defined($sam_rec_readID_G2A) ) )
            {
                $endSAM   = 1;
                $thisRead = undef;
            }

            # &$debug("skipC2T", $skipC2Ta, "skipGA", $skipG2Aa, "endC2T", $endC2Tsam, "endGA", $endG2Asam, "endSAM", $endSAM, "ThisRead", $thisRead);
        }

        my $foundInSAM = $foundInC2TSAM + $foundInG2ASAM;

        if ( $foundInSAM != 2 ) {
            if ($hisat2_un) {
                $unalChunk_parent .= $fq_rec{id} . $fq_rec{seq} . $fq_rec{id2} . $fq_rec{qs};
            }
            $mappingStats{hisat2unal} += 1;
            &$debug( "Could not find:", $fq_rec_id_sam, "not in HISAT2 output (Illumina filtered?)", "Line:", __LINE__ );
        }
        elsif (    ( $alignments->{C2T}->{samFields}->[0]->[1] == 4 )
                && ( $alignments->{G2A}->{samFields}->[0]->[1] == 4 ) )
        {
            if ($hisat2_un) {
                $unalChunk_parent .= $fq_rec{id} . $fq_rec{seq} . $fq_rec{id2} . $fq_rec{qs};
            }
            $mappingStats{hisat2unal} += 1;
        }
        else {
            $alignments->{C2T}->{NH} = $C2Ta;
            $alignments->{G2A}->{NH} = $G2Aa;
            $alignments->{fq_rec}    = {%fq_rec};
            $alignments->{YR}        = $YR;

            # cache alignments for parallel batch processing
            $todo->[ $j++ ] = $alignments;
        }

        # need alignments for 40000 reads before forking a SAMANALYZER process
        if ( ( $j == READ_CACHE_SIZE_SE ) || ( ($endSAM) && ( $j != 0 ) ) ) {
            $ident += $j;
            $j = 0;
            my @childTasks = @{$todo};
            $todo = undef;

            %mappingStatsC = (
                               'reads'        => 0,
                               'hisat2unal'   => 0,
                               'uniqMapper'   => 0,
                               'filteredAL'   => 0,
                               'discordantAL' => 0,
                               );
            my $mBiasDataC = {
                               mBiasReadDataF => [],
                               mBiasDataF     => [],
                               mBiasDataFhq   => [],
                               mBiasReadDataR => [],
                               mBiasDataR     => [],
                               mBiasDataRhq   => [],
                               };

            my $mappedSeqsC = {};

            #  SAMANALYZER:
            {    # BEGIN child scope
                if ( $samAnalyzer_threads != 0 ) {

                    # increment counter since immediately jump to next read after forking
                    $readNr++;
                    $pm->start($ident) and next READ_LOOP;
                }

                my $samSMchunk = "";
                my $samMMchunk = "";

                my $chunkSMcount = 0;
                my $chunkMMcount = 0;

                my $unalChunk = "";

                foreach my $a (@childTasks) {
                    my ( $SM, $samChunk, $mappedSeqsC, $mBiasDataC ) =
                      meRanSAManalyzerSE( $a, $readDirection, $mappedSeqsC, $mBiasDataC );

                    # Single mapper SAM chunks
                    if ( defined($samChunk) && ( $SM == 1 ) ) {
                        $samSMchunk .= $samChunk;
                        $chunkSMcount++;

                        # Flush to SAM if more than 10000 chunks
                        if ( $chunkSMcount > SAM_CHUNK_SIZE ) {
                            &$debug( "Flushing SM - SAM end:",
                                     $endSAM, "chunk count:", $chunkSMcount, "- Line:", __LINE__ );
                            $chunkSMcount = 0;
                            flushSAMchunk( $samFHs->[0], $samSMchunk, $samout_tmp_lockF );
                            undef($samSMchunk);
                        }
                    }

                    # Multimapper SAM chunks
                    elsif ( defined($samChunk) && ( $SM == 0 ) ) {
                        $samMMchunk .= $samChunk;
                        $chunkMMcount++;

                        # Flush to SAM if more than 1000 chunks
                        if ( $chunkMMcount > SAM_CHUNK_SIZE_MM ) {
                            &$debug( "Flushing MM - SAM end:",
                                     $endSAM, "chunk count:", $chunkMMcount, "- Line:", __LINE__ );
                            $chunkMMcount = 0;
                            flushSAMchunk( $samFHs->[1], $samMMchunk, $samout_MM_tmp_lockF );
                            undef($samMMchunk);
                        }
                    }

                    # No valid alignment
                    else {
                        if ($hisat2_un) {
                            $unalChunk .=
                              $a->{fq_rec}->{id} . $a->{fq_rec}->{seq} . $a->{fq_rec}->{id2} . $a->{fq_rec}->{qs};
                        }
                        $mappingStatsC{hisat2unal} += 1;
                    }

                }    # END foreach alignment in this childTask

                # Flush leftovers to SAM and unaligned fastq if requested
                if ($samSMchunk) {
                    &$debug( "Flushing SM leftovers - SAM end:",
                             $endSAM, "chunk count:", $chunkSMcount, "- Line:", __LINE__ );
                    flushSAMchunk( $samFHs->[0], $samSMchunk, $samout_tmp_lockF );
                    undef($samSMchunk);
                    $chunkSMcount = 0;
                }
                if ($samMMchunk) {
                    &$debug( "Flushing MM leftovers - SAM end:",
                             $endSAM, "chunk count:", $chunkMMcount, "- Line:", __LINE__ );
                    flushSAMchunk( $samFHs->[1], $samMMchunk, $samout_MM_tmp_lockF );
                    undef($samMMchunk);
                    $chunkMMcount = 0;
                }
                if ( $unalChunk && $hisat2_un ) {
                    &$debug( "Flushing unaligned reads - SAM end:",
                             $endSAM,
                             "unaligned count:",
                             $mappingStatsC{hisat2unal},
                             "- Line:", __LINE__ );
                    flushFQchunk( $unmappedReadsFH, $unalChunk, $unal_lockF );
                    undef($unalChunk);
                }

                &$debug( "END SAMANALYZER child - SAM end:",
                         $endSAM, "chunk count SM:",
                         $chunkSMcount, "chunk count MM:",
                         $chunkMMcount, "- Line:", __LINE__ );

            }    # END SAMANALYZER child scope

            my $childData =
              { 'mappingStats' => \%mappingStatsC, 'mBiasData' => $mBiasDataC, 'mappedSeqs' => $mappedSeqsC };

            if ( $samAnalyzer_threads != 0 ) {
                $pm->finish( 0, $childData );
            }
            else {
                updateParentData( 0, 0, 0, 0, 0, $childData );
            }

        }    # END 40000 alignments

        # &$debug( "IN parent - SAM end:", $endSAM, "todo:", $#$todo, "- Line:", __LINE__ );

        if ( ( $unalChunk_parent && $hisat2_un ) || ( $endSAM && $hisat2_un ) ) {
            &$debug( "Flushing unaligned reads parent - SAM end:",
                     $endSAM, "unaligned count:",
                     $mappingStats{hisat2unal}, "- Line:", __LINE__ );
            flushFQchunk( $unmappedReadsFH, $unalChunk_parent, $unal_lockF );
            undef($unalChunk_parent);
        }

        # &$debug($readNr, $fq_rec_id_sam);

        $readNr++;

    }
    if ( $samAnalyzer_threads != 0 ) { $pm->wait_all_children; }

    # check for leftovers in HISAT2 sam
    # This should not produce any output, if it does then something went wrong
    while ( $samRecC2T = $samC2Tfh->getline() ) {
        print STDERR "Orphan alignment found in C2T HISAT2 BAM file:\n" . $samRecC2T;
    }
    while ( $samRecG2A = $samG2Afh->getline() ) {
        print STDERR "Orphan alignment found in G2A HISAT2 BAM file:\n" . $samRecG2A;
    }

    unlink($samout_tmp_lockF);
    unlink($samout_MM_tmp_lockF);
    unlink($unal_lockF);

    return ($readNr);

}

# ($SM, $samChunk, $mappedSeqsC, $mBiasDataC) = meRanSAManalyzerSE( $a, $readDirection, $mappedSeqsC, $mBiasDataC );
sub meRanSAManalyzerSE {
    my $alignments_   = $_[0];
    my $readDirection = $_[1];
    my $mappedSeqsC   = $_[2];
    my $mBiasDataC    = $_[3];

    # mTypes:
    #          0: C2T uniq mapper;
    #          1: G2A uniq mapper;
    #          2: C2T multi mapper;
    #          3: G2A multi mapper;
    #          4: C2T/G2A multi mapper;

    my $alignments = filterSEalignments( $alignments_, $readDirection );

    # no vaild alignments left after filtering
    if ( ( $alignments->{C2T}->{NH} == 0 ) && ( $alignments->{G2A}->{NH} == 0 ) ) {
        &$debug( "No valid alignments found for :", $alignments_->{fq_rec}->{id} );
        return ( undef, undef, undef, undef );
    }

    my $alignmentCounter = 0;

    my @finalSMalignments;
    my @finalMMalignments;
    my $finalAlignments;

    chomp( $alignments_->{fq_rec}->{seq} );
    my $Seq        = $alignments_->{fq_rec}->{seq};
    my $revcompSeq = reverseComp( $alignments_->{fq_rec}->{seq} );

    my $qs = $alignments_->{fq_rec}->{qs};

    my $seqs = { Seq => $Seq, revcompSeq => $revcompSeq };

    # C2T uniq/multi mapper
    if ( ( $alignments->{C2T}->{NH} > 0 ) && ( $alignments->{G2A}->{NH} == 0 ) ) {

        # 0: C2T uniq mapper
        if ( $alignments->{C2T}->{NH} == 1 ) {
            $finalAlignments = \@finalSMalignments;
        }

        # 2: C2T multi mapper
        else {
            $finalAlignments = \@finalMMalignments;
        }

        foreach my $samFields_C2T ( @{ $alignments->{C2T}->{samFields} } ) {
            $alignmentCounter +=
              registerC2Tsingle( $samFields_C2T, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
        }
    }

    # G2A uniq/multi mapper
    elsif ( ( $alignments->{C2T}->{NH} == 0 ) && ( $alignments->{G2A}->{NH} > 0 ) ) {

        # 1: G2A uniq mapper
        if ( $alignments->{G2A}->{NH} == 1 ) {
            $finalAlignments = \@finalSMalignments;
        }

        # 3: G2A multi mapper
        else {
            $finalAlignments = \@finalMMalignments;
        }

        foreach my $samFields_G2A ( @{ $alignments->{G2A}->{samFields} } ) {
            $alignmentCounter +=
              registerG2Asingle( $samFields_G2A, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
        }
    }

    # C2T/G2A uniq/uniq mapper
    elsif ( ( $alignments->{C2T}->{NH} == 1 ) && ( $alignments->{G2A}->{NH} == 1 ) ) {

        my $samFields_C2T = $alignments->{C2T}->{samFields}->[0];
        my $samFields_G2A = $alignments->{G2A}->{samFields}->[0];

        my $AS_C2T = getASsingle($samFields_C2T);
        my $AS_G2A = getASsingle($samFields_G2A);

        if ( $AS_C2T > $AS_G2A ) {
            $finalAlignments = \@finalSMalignments;
            $alignmentCounter +=
              registerC2Tsingle( $samFields_C2T, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection);
        }
        elsif ( $AS_G2A > $AS_C2T ) {
            $finalAlignments = \@finalSMalignments;
            $alignmentCounter +=
              registerG2Asingle( $samFields_G2A, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
        }
        else {
            $finalAlignments = \@finalMMalignments;
            $alignmentCounter +=
              registerC2Tsingle( $samFields_C2T, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
            $alignmentCounter +=
              registerG2Asingle( $samFields_G2A, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
        }
    }

    # C2T/G2A uniq/multi mapper
    elsif ( ( $alignments->{C2T}->{NH} == 1 ) && ( $alignments->{G2A}->{NH} > 1 ) ) {

        my $samFields_C2T = $alignments->{C2T}->{samFields}->[0];
        my $AS_C2T        = getASsingle($samFields_C2T);

        my $max_AS = 'C2T';
        foreach my $samFields_G2A ( @{ $alignments->{G2A}->{samFields} } ) {

            my $AS_G2A = getASsingle($samFields_G2A);

            if ( $AS_C2T > $AS_G2A ) {
                $max_AS = 'C2T';
            }
            elsif ( $AS_C2T <= $AS_G2A ) {
                $max_AS = ( $AS_C2T < $AS_G2A ) ? 'G2A' : 'C2TG2A';
                $finalAlignments = \@finalMMalignments;
                $alignmentCounter +=
                  registerG2Asingle( $samFields_G2A, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
            }
        }

        if ( $max_AS eq 'C2T' ) {
            $finalAlignments = \@finalSMalignments;
        }
        if ( $max_AS eq 'C2TG2A' ) {
            $finalAlignments = \@finalMMalignments;
        }
        $alignmentCounter +=
          registerC2Tsingle( $samFields_C2T, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
    }

    # C2T/G2A multi/uniq mapper
    elsif ( ( $alignments->{C2T}->{NH} > 1 ) && ( $alignments->{G2A}->{NH} == 1 ) ) {

        my $samFields_G2A = $alignments->{G2A}->{samFields}->[0];
        my $AS_G2A        = getASsingle($samFields_G2A);

        my $max_AS = 'G2A';
        foreach my $samFields_C2T ( @{ $alignments->{C2T}->{samFields} } ) {

            my $AS_C2T = getASsingle($samFields_C2T);

            if ( $AS_G2A > $AS_C2T ) {
                $max_AS = 'G2A';
            }
            elsif ( $AS_G2A <= $AS_C2T ) {
                $max_AS = ( $AS_G2A < $AS_C2T ) ? 'C2T' : 'C2TG2A';
                $finalAlignments = \@finalMMalignments;
                $alignmentCounter +=
                  registerC2Tsingle( $samFields_C2T, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
            }

        }

        if ( $max_AS eq 'G2A' ) {
            $finalAlignments = \@finalSMalignments;
        }
        if ( $max_AS eq 'C2TG2A' ) {
            $finalAlignments = \@finalMMalignments;
        }
        $alignmentCounter +=
          registerG2Asingle( $samFields_G2A, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
    }

    # C2T/G2A multi/multi mapper
    else {
        $finalAlignments = \@finalMMalignments;

        foreach my $samFields_C2T ( @{ $alignments->{C2T}->{samFields} } ) {
            $alignmentCounter +=
              registerC2Tsingle( $samFields_C2T, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
        }
        foreach my $samFields_G2A ( @{ $alignments->{G2A}->{samFields} } ) {
            $alignmentCounter +=
              registerG2Asingle( $samFields_G2A, $finalAlignments->[$alignmentCounter], $seqs, $alignments_->{YR}, $readDirection );
        }
    }

    # write final alignments
    my $HI = 1;
    my $NH = $alignmentCounter;
    my $SM = 0;                   # single mapper

    if ( exists( $finalSMalignments[0] ) ) {
        $SM              = 1;
        $finalAlignments = \@finalSMalignments;
    }
    elsif ( exists( $finalMMalignments[0] ) ) {
        $finalAlignments = \@finalMMalignments;
    }
    else {
        &$debug( "No valid alignments found for :", $alignments_->{fq_rec}->{id} );

        # print Dumper($alignments_);

        return ( undef, undef, undef, undef );
    }

    my $samChunk;

    foreach my $finalAlignment (@$finalAlignments) {

        my $softClippedBases = splice(@$finalAlignment, $#$finalAlignment - 2, 1);
        my %auxTags = samGet_auxTags(@$finalAlignment);
        
        $auxTags{NH}{t} = "i";
        $auxTags{NH}{v} = $NH;

        $auxTags{HI}{t} = "i";
        $auxTags{HI}{v} = $HI;

        if ( ( ( $finalAlignment->[1] & 0x100 ) == 0 ) && ( $HI > 1 ) ) {
            $finalAlignment->[1] += 256;
        }
        if ( ( ( $finalAlignment->[1] & 0x100 ) != 0 ) && ( ( $NH == 1 ) || ( $HI == 1 ) ) ) {
            $finalAlignment->[1] -= 256;
        }

        $samChunk .= join( "\t", ( @$finalAlignment[ 0 .. 10 ], ( get_auxTagsFinal( \%auxTags ) ) ) ) . "\n";

        if ($SM) {
            $mappingStatsC{uniqMapper}++;
            my $dir = ( $alignments_->{YR} eq 'YR:Z:C2T' ) ? 0 : 1;

            # my $softClippedBases = getSoftClippedSE($finalAlignment);

            if ( ( $finalAlignment->[1] & 0x10 ) == 0 ) {
                $mbq->( $mBiasDataC, $Seq, $qs, $softClippedBases->[0], $softClippedBases->[1], $dir );
            }
            else {
                $mbq->( $mBiasDataC, $Seq, $qs, $softClippedBases->[1], $softClippedBases->[0], $dir );
            }

            $mappedSeqsC->{ $finalAlignment->[2] }->{u} = 1;
        }
        else {
            $mappedSeqsC->{ $finalAlignment->[2] }->{m} = 1;
        }

        $HI++;
    }

    return ( $SM, $samChunk, $mappedSeqsC, $mBiasDataC );

}

sub registerC2Tsingle {

    my $alignment       = $_[0];
    my $finalAlignments = $_[1];
    my $seq             = $_[2];
    my $readConv        = $_[3];
    my $readDir         = $_[4];

    $alignment->[9] = ( $readDir == 1 ) ? $seq->{Seq} :  $seq->{revcompSeq};
    $finalAlignments = [ @{$alignment}, 'YG:Z:C2T', $readConv ];
    $_[1] = $finalAlignments;

    return (1);

}

sub registerG2Asingle {

    my $alignment       = $_[0];
    my $finalAlignments = $_[1];
    my $seq             = $_[2];
    my $readConv        = $_[3];
    my $readDir         = $_[4];

    $alignment->[9] = ( $readDir == 1 ) ? $seq->{revcompSeq} : $seq->{Seq};
    $finalAlignments = [ @{$alignment}, 'YG:Z:G2A', $readConv ];
    $_[1] = $finalAlignments;

    return (1);
}

sub getASsingle {
    my $alignment = shift;

    my $AS = undef;

    my %auxTags = samGet_auxTags( @{$alignment} );
    $AS += $auxTags{AS}{v};
    return ($AS);
}

sub filterSEalignments {
    my $alignments_   = shift;
    my $readDirection = shift;

    my $filteredAlignments;
    my $a;

    # filter alignments on C2T strand
    $a = 0;
    foreach my $alignment ( @{ $alignments_->{C2T}->{samFields} } ) {
        my %auxTags = samGet_auxTags( @{$alignment} );
        my $softClippedBases = getSoftClippedSE($alignment);
        my $mappedLength = length( $alignment->[9] ) - $softClippedBases->[0] - $softClippedBases->[1];

        my $flag = $alignment->[1];
        if ( $flag == 4 ) {
            $alignment = undef;
        }
        elsif ( ( ( $flag & 0x10 ) != 0 ) && ( $readDirection == 1 ) ) {    # read direction = FWD (1)
            $mappingStatsC{filteredAL}++;
            &$debug( $alignment->[0] . ": maps in wrong direction: FWDrC2Trev", "Line:", __LINE__ );
            $alignment = undef;
        }
        elsif ( ( ( $flag & 0x10 ) == 0 ) && ( $readDirection == -1 ) ) {    # read direction = REV (-1)
            $mappingStatsC{filteredAL}++;
            &$debug( $alignment->[0] . ": maps in wrong direction: REVrC2Tfwd", "Line:", __LINE__ );
            $alignment = undef;
        }
        elsif ( ( checkEditDist($auxTags{NM}{v}) || checkMisMatchRatio($auxTags{XM}{v}, $mappedLength) ) ) {
            &$debug( $alignment->[0] . ": to many missmatches or edit distance exceeded", "Line:", __LINE__ );
            $alignment = undef;
        }

        if ( defined($alignment) ) {
            $filteredAlignments->{C2T}->{samFields}->[$a] = [ @{$alignment},  $softClippedBases ];
            $a++;
        }
    }

    $filteredAlignments->{C2T}->{NH} = $a;

    # filter alignments on G2A strand
    $a = 0;
    foreach my $alignment ( @{ $alignments_->{G2A}->{samFields} } ) {
        my %auxTags = samGet_auxTags( @{$alignment} );
        my $softClippedBases = getSoftClippedSE($alignment);
        my $mappedLength = length( $alignment->[9] ) - $softClippedBases->[0] - $softClippedBases->[1];

        my $flag = $alignment->[1];
        if ( $flag == 4 ) {
            $alignment = undef;
        }
        elsif ( ( ( $flag & 0x10 ) == 0 ) && ( $readDirection == 1 ) ) {
            $mappingStatsC{filteredAL}++;
            &$debug( $alignment->[0] . ": maps in wrong direction: FWDrG2Afwd", "Line:", __LINE__ );
            $alignment = undef;
        }
        elsif ( ( ( $flag & 0x10 ) != 0 ) && ( $readDirection == -1 ) ) {
            $mappingStatsC{filteredAL}++;
            &$debug( $alignment->[0] . ": maps in wrong direction: REVrG2Arev", "Line:", __LINE__ );
            $alignment = undef;
        }
        elsif ( ( checkEditDist($auxTags{NM}{v}) || checkMisMatchRatio($auxTags{XM}{v}, $mappedLength) ) ) {
            &$debug( $alignment->[0] . ": to many missmatches or edit distance exceeded", "Line:", __LINE__ );
            $alignment = undef;
        }

        if ( defined($alignment) ) {
            $filteredAlignments->{G2A}->{samFields}->[$a] = [ @{$alignment}, $softClippedBases ];
            $a++;
        }
    }

    $filteredAlignments->{G2A}->{NH} = $a;

    return ($filteredAlignments);
}

# parseHISAT2pe($HISAT2_C2T, $HISAT2_G2A, $samFHs, $fq_FHref, $unmappedReadsFHsRef, $samRawHrdRef);
sub parseHISAT2pe {
    my $samC2Tfh         = shift;
    my $samG2Afh         = shift;
    my $samFHs           = shift;
    my $fq_FH            = shift;
    my $unmappedReadsFHs = shift;
    my $samRawHrd        = shift;

    # get the RAW SAM header
    my $samheader_C2T;
    my $samheader_G2A;
    while (1) {
        $samheader_C2T = $samC2Tfh->getline();
        if ( $samheader_C2T =~ /^@/ ) {
            processhisat2samHdr( $samheader_C2T, $samRawHrd );
        }
        else {
            last;
        }
    }
    while (1) {
        $samheader_G2A = $samG2Afh->getline();
        if ( $samheader_G2A =~ /^@/ ) {
            next;
        }
        else {
            last;
        }
    }

    # get the first alignment to C2T genome
    my $samRecC2T;
    my $samFieldsC2T;
    my $sam_rec_readID_C2T;

    if (    defined( $samRecC2T->[0] = $samheader_C2T )
         && defined( $samRecC2T->[1] = $samC2Tfh->getline() ) )
    {

        # 1st mate
        chomp( $samRecC2T->[0] );
        $samFieldsC2T->[0] = [ ( split( '\t', $samRecC2T->[0] ) ) ];
        $sam_rec_readID_C2T->[0] = $samFieldsC2T->[0]->[0];

        # 2nd mate
        chomp( $samRecC2T->[1] );
        $samFieldsC2T->[1] = [ ( split( '\t', $samRecC2T->[1] ) ) ];
        $sam_rec_readID_C2T->[1] = $samFieldsC2T->[1]->[0];

        undef($samRecC2T);
    }

    # get the first alignment to G2A genome
    my $samRecG2A;
    my $samFieldsG2A;
    my $sam_rec_readID_G2A;

    if (    defined( $samRecG2A->[0] = $samheader_G2A )
         && defined( $samRecG2A->[1] = $samG2Afh->getline() ) )
    {

        # 1st mate
        chomp( $samRecG2A->[0] );
        $samFieldsG2A->[0] = [ ( split( '\t', $samRecG2A->[0] ) ) ];
        $sam_rec_readID_G2A->[0] = $samFieldsG2A->[0]->[0];

        # 2nd mate
        chomp( $samRecG2A->[1] );
        $samFieldsG2A->[1] = [ ( split( '\t', $samRecG2A->[1] ) ) ];
        $sam_rec_readID_G2A->[1] = $samFieldsG2A->[1]->[0];

        undef($samRecG2A);
    }

    # No valid alignments ???
    if ( ( !$sam_rec_readID_C2T->[0] ) && ( !$sam_rec_readID_G2A->[0] ) ) {
        say STDERR "No valid alignments found in HISAT2 output (Line: " . __LINE__ . ")";
        exit(1);
    }

    my $pairNr            = 0;
    my $readNr            = 0;
    my $maxReads          = $firstNreads * 2;
    my $fwd_fq_rec_id     = "";
    my $rev_fq_rec_id     = "";
    my $fwd_fq_rec_id_sam = "";
    my $rev_fq_rec_id_sam = "";

    my $skipC2Ta  = 1;
    my $skipG2Aa  = 1;
    my $endC2Tsam = 0;
    my $endG2Asam = 0;
    my $endSAM    = 0;

    my $unmappedFWDReadsFH = $unmappedReadsFHs->{FWD}->[0] if ($hisat2_un);
    my $unmappedREVReadsFH = $unmappedReadsFHs->{REV}->[0] if ($hisat2_un);

    my $j = 0;
    my $todo;
    my $ident = 0;

    my $samout_tmp_lockF    = $lockDir . "/meRanGhSM_" . $$ . ".lock";
    my $samout_MM_tmp_lockF = $lockDir . "/meRanGhMM_" . $$ . ".lock";
    my $unal_lockF          = $lockDir . "/meRanGhUA_" . $$ . ".lock";

    my $fwd_unalChunk_parent = "";
    my $rev_unalChunk_parent = "";

    my $pm = new Parallel::ForkManager( $samAnalyzer_threads, "/dev/shm" );

    $pm->run_on_finish( \&updateParentData );

    my $cleanIDre = qr/(\/1|\/2|[ \t]+.+)?$/;

  READ_LOOP:
    while (1) {

        my %fwd_fq_rec = getFQrec( $fq_FH->{FWD}, 0 );
        my %rev_fq_rec = getFQrec( $fq_FH->{REV}, 0 );

        last READ_LOOP
          unless (    ( $fwd_fq_rec{id} && $fwd_fq_rec{seq} && $fwd_fq_rec{id2} && $fwd_fq_rec{qs} )
                   && ( $rev_fq_rec{id} && $rev_fq_rec{seq} && $rev_fq_rec{id2} && $rev_fq_rec{qs} )
                   && ( ( $pairNr <= $maxReads ) || ( $firstNreads == -1 ) ) );

        $fwd_fq_rec_id = $fwd_fq_rec_id_sam = $fwd_fq_rec{id};
        $rev_fq_rec_id = $rev_fq_rec_id_sam = $rev_fq_rec{id};

        $fwd_fq_rec_id_sam =~ s/$cleanIDre//;
        $fwd_fq_rec_id_sam = substr( $fwd_fq_rec_id_sam, 1 );

        $rev_fq_rec_id_sam =~ s/$cleanIDre//;
        $rev_fq_rec_id_sam = substr( $rev_fq_rec_id_sam, 1 );

        chomp( $fwd_fq_rec_id_sam, $rev_fq_rec_id_sam );

        if ( $fwd_fq_rec_id_sam ne $rev_fq_rec_id_sam ) {
            say STDERR "Reads not properly paired: "
              . $fwd_fq_rec_id_sam . " <-> "
              . $rev_fq_rec_id_sam
              . " (Line: "
              . __LINE__ . ")";
            exit(1);
        }

        # ( FWD => [@unmappedFWDReadsFHs], REV => [@unmappedREVReadsFHs] );
        if ( $fwd_fq_rec_id_sam =~ /_meRanGh_FQ_FILE_SWITCH_/ ) {
            if ($hisat2_un) {

                # switch FWD unmapped file
                $unmappedFWDReadsFH->flush();
                $unmappedFWDReadsFH->close();
                undef($unmappedFWDReadsFH);
                shift( @{$unmappedReadsFHs->{FWD}} );
                $unmappedFWDReadsFH = $unmappedReadsFHs->{FWD}->[0];

                # switch REV unmapped file
                $unmappedREVReadsFH->flush();
                $unmappedREVReadsFH->close();
                undef($unmappedREVReadsFH);
                shift( @{$unmappedReadsFHs->{REV}} );
                $unmappedREVReadsFH = $unmappedReadsFHs->{REV}->[0];
            }
            next READ_LOOP;
        }

        my $foundInC2TSAM = 0;
        my $foundInG2ASAM = 0;
        my $thisPair      = 1;

        my $alignments = undef;
        my ( $C2Ta, $G2Aa ) = ( 0, 0 );

      SAM_LOOP:
        while ( $thisPair && ( $endSAM != 1 ) ) {

            if ( $samC2Tfh->eof ) { $endC2Tsam = 1 }
            if ( $samG2Afh->eof ) { $endG2Asam = 1 }

            unless ( $skipC2Ta or $endC2Tsam ) {

                # 1st mate
                $samRecC2T->[0] = $samC2Tfh->getline();
                chomp( $samRecC2T->[0] );
                $samFieldsC2T->[0] = [ ( split( '\t', $samRecC2T->[0] ) ) ];
                $sam_rec_readID_C2T->[0] = $samFieldsC2T->[0]->[0];

                # 2nd mate
                $samRecC2T->[1] = $samC2Tfh->getline();
                chomp( $samRecC2T->[1] );
                $samFieldsC2T->[1] = [ ( split( '\t', $samRecC2T->[1] ) ) ];
                $sam_rec_readID_C2T->[1] = $samFieldsC2T->[1]->[0];

                if ( $sam_rec_readID_C2T->[0] ne $sam_rec_readID_C2T->[1] ) {
                    say STDERR "ERROR: unpaired alignment, this should not happen! (Line: " . __LINE__ . ")";
                    exit(1);
                }
            }

            unless ( $skipG2Aa or $endG2Asam ) {

                # 1st mate
                $samRecG2A->[0] = $samG2Afh->getline();
                chomp( $samRecG2A->[0] );
                $samFieldsG2A->[0] = [ ( split( '\t', $samRecG2A->[0] ) ) ];
                $sam_rec_readID_G2A->[0] = $samFieldsG2A->[0]->[0];

                # 2nd mate
                $samRecG2A->[1] = $samG2Afh->getline();
                chomp( $samRecG2A->[1] );
                $samFieldsG2A->[1] = [ ( split( '\t', $samRecG2A->[1] ) ) ];
                $sam_rec_readID_G2A->[1] = $samFieldsG2A->[1]->[0];

                if ( $sam_rec_readID_G2A->[0] ne $sam_rec_readID_G2A->[1] ) {
                    say STDERR "ERROR: unpaired alignment, this should not happen! (Line: " . __LINE__ . ")";
                    exit(1);
                }

            }

            if (    defined( $sam_rec_readID_C2T->[0] )
                 && defined( $sam_rec_readID_C2T->[1] )
                 && ( $fwd_fq_rec_id_sam eq $sam_rec_readID_C2T->[0] )
                 && ( $rev_fq_rec_id_sam eq $sam_rec_readID_C2T->[1] ) )
            {
                $foundInC2TSAM                                = 1;
                $skipC2Ta                                     = 0;
                $alignments->{C2T}->{samFields}->[$C2Ta]->[0] = [ @{ $samFieldsC2T->[0] } ];
                $alignments->{C2T}->{samFields}->[$C2Ta]->[1] = [ @{ $samFieldsC2T->[1] } ];
                $C2Ta++;
                undef($sam_rec_readID_C2T);
            }
            else {
                $skipC2Ta = 1;
            }

            if (    defined( $sam_rec_readID_G2A->[0] )
                 && defined( $sam_rec_readID_G2A->[1] )
                 && ( $fwd_fq_rec_id_sam eq $sam_rec_readID_G2A->[0] )
                 && ( $rev_fq_rec_id_sam eq $sam_rec_readID_G2A->[1] ) )
            {
                $foundInG2ASAM                                = 1;
                $skipG2Aa                                     = 0;
                $alignments->{G2A}->{samFields}->[$G2Aa]->[0] = [ @{ $samFieldsG2A->[0] } ];
                $alignments->{G2A}->{samFields}->[$G2Aa]->[1] = [ @{ $samFieldsG2A->[1] } ];
                $G2Aa++;
                undef($sam_rec_readID_G2A);
            }
            else {
                $skipG2Aa = 1;
            }

            if ( $skipC2Ta + $skipG2Aa == 2 ) {
                $thisPair = 0;
            }

            if (    ( $endC2Tsam + $endG2Asam == 2 )
                 && ( !defined( $sam_rec_readID_C2T->[0] ) )
                 && ( !defined( $sam_rec_readID_G2A->[0] ) ) )
            {
                $endSAM   = 1;
                $thisPair = undef;
            }
        }

        my $foundInSAM = $foundInC2TSAM + $foundInG2ASAM;

        if ( $foundInSAM != 2 ) {
            if ($hisat2_un) {
                $fwd_unalChunk_parent .= $fwd_fq_rec{id} . $fwd_fq_rec{seq} . $fwd_fq_rec{id2} . $fwd_fq_rec{qs};
                $rev_unalChunk_parent .= $rev_fq_rec{id} . $rev_fq_rec{seq} . $rev_fq_rec{id2} . $rev_fq_rec{qs};
            }
            $mappingStats{hisat2unal} += 2;
            &$debug( "Could not find:", $fwd_fq_rec_id_sam, "UNALIGNED", "Line:", __LINE__ );
        }
        elsif (    ( $alignments->{C2T}->{samFields}->[0]->[0]->[1] == 4 )
                && ( $alignments->{C2T}->{samFields}->[0]->[1]->[1] == 4 )
                && ( $alignments->{G2A}->{samFields}->[0]->[0]->[1] == 4 )
                && ( $alignments->{G2A}->{samFields}->[0]->[1]->[1] == 4 ) )
        {
            if ($hisat2_un) {
                $fwd_unalChunk_parent .= $fwd_fq_rec{id} . $fwd_fq_rec{seq} . $fwd_fq_rec{id2} . $fwd_fq_rec{qs};
                $rev_unalChunk_parent .= $rev_fq_rec{id} . $rev_fq_rec{seq} . $rev_fq_rec{id2} . $rev_fq_rec{qs};
            }
            $mappingStats{hisat2unal} += 2;
        }
        else {
            $alignments->{C2T}->{NH}  = $C2Ta;
            $alignments->{G2A}->{NH}  = $G2Aa;
            $alignments->{fwd_fq_rec} = {%fwd_fq_rec};
            $alignments->{rev_fq_rec} = {%rev_fq_rec};

            # cache alignments for parallel batch processing
            $todo->[ $j++ ] = $alignments;
        }

        # need alignments for 15000 reads before forking a SAMANALYZER process
        if ( ( $j == READ_CACHE_SIZE_PE ) || ( ($endSAM) && ( $j != 0 ) ) ) {
            $ident += $j;
            $j = 0;
            my @childTasks = @{$todo};
            $todo = undef;

            %mappingStatsC = (
                               'reads'        => 0,
                               'hisat2unal'   => 0,
                               'uniqMapper'   => 0,
                               'filteredAL'   => 0,
                               'discordantAL' => 0,
                               );
            my $mBiasDataC = {
                               mBiasReadDataF => [],
                               mBiasDataF     => [],
                               mBiasDataFhq   => [],
                               mBiasReadDataR => [],
                               mBiasDataR     => [],
                               mBiasDataRhq   => [],
                               };

            my $mappedSeqsC = {};

            #  SAMANALYZER:
            {    # BEGIN child scope
                if ( $samAnalyzer_threads != 0 ) {

                    # increment counters since immediately jump to next read after forking
                    $readNr += 2;
                    $pairNr++;
                    $pm->start($ident) and next READ_LOOP;
                }

                my $samSMchunk = "";
                my $samMMchunk = "";

                my $chunkSMcount = 0;
                my $chunkMMcount = 0;

                my $fwd_unalChunk = "";
                my $rev_unalChunk = "";

                foreach my $a (@childTasks) {
                    my ( $SM, $samChunk, $mappedSeqsC, $mBiasDataC ) =
                      meRanSAManalyzerPE( $a, $mappedSeqsC, $mBiasDataC );

                    # Single mapper SAM chunks
                    if ( defined($samChunk) && ( $SM == 1 ) ) {
                        $samSMchunk .= $samChunk;
                        $chunkSMcount++;

                        # Flush to SAM if more than 10000 chunks
                        if ( $chunkSMcount > SAM_CHUNK_SIZE ) {
                            &$debug( "Flushing SM - SAM end:",
                                     $endSAM, "chunk count:", $chunkSMcount, "- Line:", __LINE__ );
                            $chunkSMcount = 0;
                            flushSAMchunk( $samFHs->[0], $samSMchunk, $samout_tmp_lockF );
                            undef($samSMchunk);
                        }
                    }

                    # Multimapper SAM chunks
                    elsif ( defined($samChunk) && ( $SM == 0 ) ) {
                        $samMMchunk .= $samChunk;
                        $chunkMMcount++;

                        # Flush to SAM if more than 1000 chunks
                        if ( $chunkMMcount > SAM_CHUNK_SIZE_MM ) {
                            &$debug( "Flushing MM - SAM end:",
                                     $endSAM, "chunk count:", $chunkMMcount, "- Line:", __LINE__ );
                            $chunkMMcount = 0;
                            flushSAMchunk( $samFHs->[1], $samMMchunk, $samout_MM_tmp_lockF );
                            undef($samMMchunk);
                        }
                    }

                    # No valid alignment
                    else {
                        if ($hisat2_un) {
                            $fwd_unalChunk .=
                                $a->{fwd_fq_rec}->{id}
                              . $a->{fwd_fq_rec}->{seq}
                              . $a->{fwd_fq_rec}->{id2}
                              . $a->{fwd_fq_rec}->{qs};
                            $rev_unalChunk .=
                                $a->{rev_fq_rec}->{id}
                              . $a->{rev_fq_rec}->{seq}
                              . $a->{rev_fq_rec}->{id2}
                              . $a->{rev_fq_rec}->{qs};
                        }
                        $mappingStatsC{hisat2unal} += 2;
                    }

                }    # END foreach alignment in this childTask

                # Flush leftovers to SAM and unaligned fastq if requested
                if ($samSMchunk) {
                    &$debug( "Flushing SM leftovers - SAM end:",
                             $endSAM, "chunk count:", $chunkSMcount, "- Line:", __LINE__ );
                    flushSAMchunk( $samFHs->[0], $samSMchunk, $samout_tmp_lockF );
                    undef($samSMchunk);
                    $chunkSMcount = 0;
                }
                if ($samMMchunk) {
                    &$debug( "Flushing MM leftovers - SAM end:",
                             $endSAM, "chunk count:", $chunkMMcount, "- Line:", __LINE__ );
                    flushSAMchunk( $samFHs->[1], $samMMchunk, $samout_MM_tmp_lockF );
                    undef($samMMchunk);
                    $chunkMMcount = 0;
                }
                if ( $fwd_unalChunk && $rev_unalChunk && $hisat2_un ) {
                    &$debug( "Flushing unaligned reads - SAM end:",
                             $endSAM,
                             "unaligned count:",
                             $mappingStatsC{hisat2unal},
                             "- Line:", __LINE__ );
                    flushFQchunk( $unmappedFWDReadsFH, $fwd_unalChunk, $unal_lockF );
                    flushFQchunk( $unmappedREVReadsFH, $rev_unalChunk, $unal_lockF );
                    undef($fwd_unalChunk);
                    undef($rev_unalChunk);
                }

                # $todo = undef;

                &$debug( "END SAMANALYZER child - SAM end:",
                         $endSAM, "chunk count SM:",
                         $chunkSMcount, "chunk count MM:",
                         $chunkMMcount, "- Line:", __LINE__ );

            }    # END SAMANALYZER child scope

            my $childData =
              { 'mappingStats' => \%mappingStatsC, 'mBiasData' => $mBiasDataC, 'mappedSeqs' => $mappedSeqsC };

            if ( $samAnalyzer_threads != 0 ) {
                $pm->finish( 0, $childData );
            }
            else {
                updateParentData( 0, 0, 0, 0, 0, $childData );
            }

        }    # END 15000 alignments

        # &$debug( "IN parent - SAM end:", $endSAM, "todo:", $#$todo, "- Line:", __LINE__ );

        if ( ( $fwd_unalChunk_parent && $rev_unalChunk_parent && $hisat2_un ) || ( $endSAM && $hisat2_un ) ) {
            &$debug( "Flushing unaligned reads parent - SAM end:",
                     $endSAM, "unaligned count:",
                     $mappingStats{hisat2unal}, "- Line:", __LINE__ );
            flushFQchunk( $unmappedFWDReadsFH, $fwd_unalChunk_parent, $unal_lockF );
            flushFQchunk( $unmappedREVReadsFH, $rev_unalChunk_parent, $unal_lockF );
            undef($fwd_unalChunk_parent);
            undef($rev_unalChunk_parent);
        }

        $readNr += 2;
        $pairNr++;
    }
    if ( $samAnalyzer_threads != 0 ) { $pm->wait_all_children; }

    # check for leftovers in HISAT2 sam
    # This should not produce any output, if it does then something went wrong
    while ( $samRecC2T = $samC2Tfh->getline() ) {
        print STDERR "Orphan alignment found in C2T HISAT2 output:\n" . $samRecC2T;
    }
    while ( $samRecG2A = $samG2Afh->getline() ) {
        print STDERR "Orphan alignment found in G2A HISAT2 output:\n" . $samRecG2A;
    }

    unlink($samout_tmp_lockF);
    unlink($samout_MM_tmp_lockF);
    unlink($unal_lockF);

    return ($readNr);

}

BEGIN {
    *flushSAMchunk = \&flushChunk;
    *flushFQchunk  = \&flushChunk;
}

sub flushChunk {
    my $writeFH = $_[0];
    my $chunk   = $_[1];
    my $lockF   = $_[2];

    my $lockFH = IO::File->new( $lockF, O_WRONLY | O_CREAT ) || die( $lockF . ": " . $! );

    lockF($lockFH);
    $writeFH->seek( 0, SEEK_END );
    $writeFH->print($chunk);
    $writeFH->flush();
    unlockF($lockFH);
    $lockFH->close();
    $_[1] = undef;
}

sub updateParentData {
    my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $childData ) = @_;

    &$debug( "PID:",         $pid,         "Exit code:", $exit_code, "Ident:",  $ident,
             "Exit signal:", $exit_signal, "Core dump:", $core_dump, "- Line:", __LINE__ );

    updateMappingStats( $childData->{mappingStats} );
    updateMbiasData( $childData->{mBiasData} );
    updateMappedSeqs( $childData->{mappedSeqs} );

    &$debug( "DONE child PID:", $pid, "- Line:", __LINE__ );
}

sub updateMappingStats {
    my $mappingStatsData = $_[0];

    $mappingStats{reads}        += $mappingStatsData->{reads};
    $mappingStats{hisat2unal}   += $mappingStatsData->{hisat2unal};
    $mappingStats{uniqMapper}   += $mappingStatsData->{uniqMapper};
    $mappingStats{filteredAL}   += $mappingStatsData->{filteredAL};
    $mappingStats{discordantAL} += $mappingStatsData->{discordantAL};

}

sub updateMbiasData {
    my $mBiasDataC = $_[0];

    sumUpArray( $mBiasData->{mBiasReadDataF}, $mBiasDataC->{mBiasReadDataF} );
    sumUpArray( $mBiasData->{mBiasDataF},     $mBiasDataC->{mBiasDataF} );
    sumUpArray( $mBiasData->{mBiasDataFhq},   $mBiasDataC->{mBiasDataFhq} );
    sumUpArray( $mBiasData->{mBiasReadDataR}, $mBiasDataC->{mBiasReadDataR} );
    sumUpArray( $mBiasData->{mBiasDataR},     $mBiasDataC->{mBiasDataR} );
    sumUpArray( $mBiasData->{mBiasDataRhq},   $mBiasDataC->{mBiasDataRhq} );
}

sub updateMappedSeqs {
    my $mappedSeqsC = $_[0];

    $mappedSeqs = mergeMappedSeqHash( $mappedSeqs, $mappedSeqsC );

}

sub sumUpArray {
    my $a = $_[0];
    my $b = $_[1];

    my $i = 0;
    map { $a->[ $i++ ] += $_ } @{$b};

    $_[0] = $a;

}

sub meRanSAManalyzerPE {
    my $alignments_ = $_[0];
    my $mappedSeqsC = $_[1];
    my $mBiasDataC  = $_[2];

    # mTypes:
    #          0: C2T uniq mapper;
    #          1: G2A uniq mapper;
    #          2: C2T multi mapper;
    #          3: G2A multi mapper;
    #          4: C2T/G2A multi mapper;

    my $alignments = filterPEalignments($alignments_);

    # no vaild alignments left after filtering
    if ( ( $alignments->{C2T}->{NH} == 0 ) && ( $alignments->{G2A}->{NH} == 0 ) ) {
        return ( undef, undef, undef, undef );
    }

    my $alignmentCounter = 0;
    my $mateCount        = 0;
    my $a                = 0;

    my @finalSMalignments;
    my @finalMMalignments;
    my $finalAlignments;

    chomp( $alignments_->{fwd_fq_rec}->{seq} );
    chomp( $alignments_->{rev_fq_rec}->{seq} );
    my $fwd_Seq        = $alignments_->{fwd_fq_rec}->{seq};
    my $rev_Seq        = $alignments_->{rev_fq_rec}->{seq};
    my $fwd_revcompSeq = reverseComp( $alignments_->{fwd_fq_rec}->{seq} );
    my $rev_revcompSeq = reverseComp( $alignments_->{rev_fq_rec}->{seq} );

    my $read;
    $read->[0] = $fwd_Seq;
    $read->[1] = $rev_Seq;

    my $qs;
    $qs->[0] = $alignments_->{fwd_fq_rec}->{qs};
    $qs->[1] = $alignments_->{rev_fq_rec}->{qs};

    my $seqs = {
                 fwd_Seq        => $fwd_Seq,
                 rev_Seq        => $rev_Seq,
                 fwd_revcompSeq => $fwd_revcompSeq,
                 rev_revcompSeq => $rev_revcompSeq
                 };

    # C2T uniq/multi mapper
    if ( ( $alignments->{C2T}->{NH} > 0 ) && ( $alignments->{G2A}->{NH} == 0 ) ) {

        # 0: C2T uniq mapper
        if ( $alignments->{C2T}->{NH} == 1 ) {
            $finalAlignments = \@finalSMalignments;
        }

        # 2: C2T multi mapper
        else {
            $finalAlignments = \@finalMMalignments;
        }

        $a = 0;
        foreach my $pairedAlignment_C2T ( @{ $alignments->{C2T}->{samFields} } ) {
            $mateCount = registerC2Tpair( $pairedAlignment_C2T, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{C2T}->{mapData}->[ $a++ ] );
            $alignmentCounter++ if $mateCount;
        }
    }

    # G2A uniq/multi mapper
    elsif ( ( $alignments->{C2T}->{NH} == 0 ) && ( $alignments->{G2A}->{NH} > 0 ) ) {

        # 1: G2A uniq mapper
        if ( $alignments->{G2A}->{NH} == 1 ) {
            $finalAlignments = \@finalSMalignments;
        }

        # 3: G2A multi mapper
        else {
            $finalAlignments = \@finalMMalignments;
        }

        $a = 0;
        foreach my $pairedAlignment_G2A ( @{ $alignments->{G2A}->{samFields} } ) {
            $mateCount = registerG2Apair( $pairedAlignment_G2A, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{G2A}->{mapData}->[ $a++ ] );
            $alignmentCounter++ if $mateCount;
        }
    }

    # C2T/G2A uniq/uniq mapper
    elsif ( ( $alignments->{C2T}->{NH} == 1 ) && ( $alignments->{G2A}->{NH} == 1 ) ) {

        my $pairedAlignment_C2T = $alignments->{C2T}->{samFields}->[0];
        my $pairedAlignment_G2A = $alignments->{G2A}->{samFields}->[0];

        # C2T
        my $AS_C2T = getASpaired($pairedAlignment_C2T) // -1_000_000;

        # G2A
        my $AS_G2A = getASpaired($pairedAlignment_G2A) // -1_000_000;
        

        if ( $AS_C2T > $AS_G2A ) {
            $finalAlignments = \@finalSMalignments;
            $mateCount = registerC2Tpair( $pairedAlignment_C2T, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{C2T}->{mapData}->[0] );
            $alignmentCounter++ if $mateCount;
        }
        elsif ( $AS_G2A > $AS_C2T ) {
            $finalAlignments = \@finalSMalignments;
            $mateCount = registerG2Apair( $pairedAlignment_G2A, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{G2A}->{mapData}->[0] );
            $alignmentCounter++ if $mateCount;
        }
        else {
            $finalAlignments = \@finalMMalignments;

            $mateCount = registerC2Tpair( $pairedAlignment_C2T, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{C2T}->{mapData}->[0] );
            $alignmentCounter++ if $mateCount;

            $mateCount = registerG2Apair( $pairedAlignment_G2A, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{G2A}->{mapData}->[0] );
            $alignmentCounter++ if $mateCount;

        }
    }

    # C2T/G2A uniq/multi mapper
    elsif ( ( $alignments->{C2T}->{NH} == 1 ) && ( $alignments->{G2A}->{NH} > 1 ) ) {

        # C2T
        my $pairedAlignment_C2T = $alignments->{C2T}->{samFields}->[0];
        my $AS_C2T = getASpaired($pairedAlignment_C2T) // -1_000_000;

        my $max_AS = 'C2T';

        $a = 0;
        foreach my $pairedAlignment_G2A ( @{ $alignments->{G2A}->{samFields} } ) {

            # G2A
            my $AS_G2A = getASpaired($pairedAlignment_G2A) // -1_000_000;

            if ( $AS_C2T > $AS_G2A ) {
                $max_AS = 'C2T';
            }
            elsif ( $AS_C2T <= $AS_G2A ) {
                $max_AS = ( $AS_C2T < $AS_G2A ) ? 'G2A' : 'C2TG2A';
                $finalAlignments = \@finalMMalignments;
                $mateCount = registerG2Apair( $pairedAlignment_G2A, $finalAlignments->[$alignmentCounter],
                                              $seqs,                $alignments->{G2A}->{mapData}->[$a] );
                $alignmentCounter++ if $mateCount;
            }
            $a++;
        }

        if ( $max_AS eq 'C2T' ) {
            $finalAlignments = \@finalSMalignments;
        }
        if ( $max_AS eq 'C2TG2A' ) {
            $finalAlignments = \@finalMMalignments;
        }
        $mateCount = registerC2Tpair( $pairedAlignment_C2T, $finalAlignments->[$alignmentCounter],
                                      $seqs,                $alignments->{C2T}->{mapData}->[0] );
        $alignmentCounter++ if $mateCount;
    }

    # C2T/G2A multi/uniq mapper
    elsif ( ( $alignments->{C2T}->{NH} > 1 ) && ( $alignments->{G2A}->{NH} == 1 ) ) {

        # G2A
        my $pairedAlignment_G2A = $alignments->{G2A}->{samFields}->[0];
        my $AS_G2A = getASpaired($pairedAlignment_G2A) // -1_000_000;

        my $max_AS = 'G2A';

        $a = 0;
        foreach my $pairedAlignment_C2T ( @{ $alignments->{C2T}->{samFields} } ) {

            # C2T
            my $AS_C2T = getASpaired($pairedAlignment_C2T) // -1_000_000;

            if ( $AS_G2A > $AS_C2T ) {
                $max_AS = 'G2A';
            }
            elsif ( $AS_G2A <= $AS_C2T ) {
                $max_AS = ( $AS_G2A < $AS_C2T ) ? 'C2T' : 'C2TG2A';
                $finalAlignments = \@finalMMalignments;
                $mateCount = registerC2Tpair( $pairedAlignment_C2T, $finalAlignments->[$alignmentCounter],
                                              $seqs,                $alignments->{C2T}->{mapData}->[$a] );
                $alignmentCounter++ if $mateCount;
            }
            $a++;
        }

        if ( $max_AS eq 'G2A' ) {
            $finalAlignments = \@finalSMalignments;
        }
        if ( $max_AS eq 'C2TG2A' ) {
            $finalAlignments = \@finalMMalignments;
        }
        $mateCount = registerG2Apair( $pairedAlignment_G2A, $finalAlignments->[$alignmentCounter],
                                      $seqs,                $alignments->{G2A}->{mapData}->[0] );
        $alignmentCounter++ if $mateCount;
    }

    # C2T/G2A multi/multi mapper
    else {
        $finalAlignments = \@finalMMalignments;

        $a = 0;
        foreach my $pairedAlignment_C2T ( @{ $alignments->{C2T}->{samFields} } ) {
            $mateCount = registerC2Tpair( $pairedAlignment_C2T, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{C2T}->{mapData}->[ $a++ ] );
            $alignmentCounter++ if $mateCount;
        }

        $a = 0;
        foreach my $pairedAlignment_G2A ( @{ $alignments->{G2A}->{samFields} } ) {
            $mateCount = registerG2Apair( $pairedAlignment_G2A, $finalAlignments->[$alignmentCounter],
                                          $seqs,                $alignments->{G2A}->{mapData}->[ $a++ ] );
            $alignmentCounter++ if $mateCount;
        }
    }

    # write final alignments
    my $HI = 1;
    my $NH = $alignmentCounter;
    my $SM = 0;
    if ( exists( $finalSMalignments[0] ) ) {
        $SM              = 1;
        $finalAlignments = \@finalSMalignments;
    }
    elsif ( exists( $finalMMalignments[0] ) ) {
        $finalAlignments = \@finalMMalignments;
    }
    else {
        chomp( $alignments_->{fwd_fq_rec}->{id}, $alignments_->{rev_fq_rec}->{id} );
        &$debug( "No valid alignments found for read pair:",
                 $alignments_->{fwd_fq_rec}->{id},
                 "<->", $alignments_->{rev_fq_rec}->{id},
                 "Line:", __LINE__ );

        # print Dumper($alignments_);
        return ( undef, undef, undef, undef );
    }

    my $samChunk;

    foreach my $finalAlignmentPair (@$finalAlignments) {
        my $i = 0;
        my $samLineFinal;
        my $mateOverlap = 0;
        my $softClippedBases;
        my $fwdMateMappedLength;
        my $revMateMappedLength;

        if ( defined( $finalAlignmentPair->[0] ) && defined( $finalAlignmentPair->[1] ) ) {
            $softClippedBases    = $finalAlignmentPair->[2]->[0];
            $fwdMateMappedLength = $finalAlignmentPair->[2]->[1];
            $revMateMappedLength = $finalAlignmentPair->[2]->[2];
            $mateOverlap         = $getMateOverlap->( $finalAlignmentPair, $fwdMateMappedLength, $revMateMappedLength );
        }

        foreach my $finalAlignment ( @$finalAlignmentPair[ 0 .. 1 ] ) {
            if ( defined($finalAlignment) ) {
                my %auxTags = samGet_auxTags(@$finalAlignment);

                $auxTags{NH}{t} = "i";
                $auxTags{NH}{v} = $NH;

                $auxTags{HI}{t} = "i";
                $auxTags{HI}{v} = $HI;

                if ( ( ( $finalAlignment->[1] & 0x100 ) == 0 ) && ( $HI > 1 ) ) {
                    $finalAlignment->[1] += 256;
                }
                if ( ( ( $finalAlignment->[1] & 0x100 ) != 0 ) && ( ( $NH == 1 ) || ( $HI == 1 ) ) ) {
                    $finalAlignment->[1] -= 256;
                }

                if ( $mateOverlap != 0 ) {
                    $fixMateOverlap->( $finalAlignment, \%auxTags, $mateOverlap, $softClippedBases->[$i] );
                }

                $samLineFinal->[$i] = [ @$finalAlignment[ 0 .. 10 ], ( get_auxTagsFinal( \%auxTags ) ) ];

                if ($SM) {
                    $mappingStatsC{uniqMapper}++;

                    if ( ( $finalAlignment->[1] & 0x10 ) == 0 ) {
                        $mbq->(
                                $mBiasDataC, $read->[$i], $qs->[$i],
                                $softClippedBases->[$i]->[0],
                                $softClippedBases->[$i]->[1], $i
                                );
                    }
                    else {
                        $mbq->(
                                $mBiasDataC, $read->[$i], $qs->[$i],
                                $softClippedBases->[$i]->[1],
                                $softClippedBases->[$i]->[0], $i
                                );
                    }

                    $mappedSeqsC->{ $finalAlignment->[2] }->{u} = 1;
                }
                else {
                    $mappedSeqsC->{ $finalAlignment->[2] }->{m} = 1;
                }
            }
            $i++;
        }

        if ( $mateOverlap != 0 ) {
            $samLineFinal->[0]->[7] = $samLineFinal->[1]->[3];
            $samLineFinal->[1]->[7] = $samLineFinal->[0]->[3];
        }

        foreach my $line (@$samLineFinal) {
            if ( defined($line) ) {
                $samChunk .= join( "\t", @$line ) . "\n";
            }
        }

        $HI++;
    }

    return ( $SM, $samChunk, $mappedSeqsC, $mBiasDataC );

}

sub sortPE {
    my $sam_rec = shift;

    my $validMates = 2;

    # Make sure that fwd read is idx 0 and rev read is idx 1
    my $mates;
    my $mates_tmp;

    $mates_tmp->[0] = [ @{ $sam_rec->[0] } ];
    $mates_tmp->[1] = [ @{ $sam_rec->[1] } ];

    # sort mates
    my $fqMateNr_0 = ( ( $mates_tmp->[0]->[1] & 0x40 ) != 0 ) ? 1 : 2;
    my $fqMateNr_1 = ( ( $mates_tmp->[1]->[1] & 0x80 ) != 0 ) ? 2 : 1;

    $mates->[0] = ( $fqMateNr_0 == 1 ) ? [ @{ $mates_tmp->[0] } ] : [ @{ $mates_tmp->[1] } ];
    $mates->[1] = ( $fqMateNr_1 == 2 ) ? [ @{ $mates_tmp->[1] } ] : [ @{ $mates_tmp->[0] } ];

    # done mate sorting

    return ($mates);
}

sub registerC2Tpair {
    my $pairedAlignment = $_[0];
    my $finalAlignments = $_[1];
    my $seq             = $_[2];
    my $mapData         = $_[3];

    my $mateCount = 0;

    if ( defined( $pairedAlignment->[0] ) ) {
        $pairedAlignment->[0]->[9] = $seq->{fwd_Seq};
        $finalAlignments->[0] = [ @{ $pairedAlignment->[0] }, 'YG:Z:C2T', 'YR:Z:C2T' ];
        $mateCount++;
    }
    if ( defined( $pairedAlignment->[1] ) ) {
        $pairedAlignment->[1]->[9] = $seq->{rev_revcompSeq};
        $finalAlignments->[1] = [ @{ $pairedAlignment->[1] }, 'YG:Z:C2T', 'YR:Z:G2A' ];
        $mateCount++;
    }

    $finalAlignments->[2] = $mapData;

    $_[1] = $finalAlignments;

    return ($mateCount);
}

sub registerG2Apair {
    my $pairedAlignment = $_[0];
    my $finalAlignments = $_[1];
    my $seq             = $_[2];
    my $mapData         = $_[3];

    my $mateCount = 0;

    if ( defined( $pairedAlignment->[0] ) ) {
        $pairedAlignment->[0]->[9] = $seq->{fwd_revcompSeq};
        $finalAlignments->[0] = [ @{ $pairedAlignment->[0] }, 'YG:Z:G2A', 'YR:Z:C2T' ];
        $mateCount++;
    }
    if ( defined( $pairedAlignment->[1] ) ) {
        $pairedAlignment->[1]->[9] = $seq->{rev_Seq};
        $finalAlignments->[1] = [ @{ $pairedAlignment->[1] }, 'YG:Z:G2A', 'YR:Z:G2A' ];
        $mateCount++;
    }

    $finalAlignments->[2] = $mapData;

    $_[1] = $finalAlignments;

    return ($mateCount);
}

sub getASpaired {
    my $pairedAlignment = shift;

    my $AS = undef;

    if ( defined( $pairedAlignment->[0] ) ) {
        my %auxTags_fwd = samGet_auxTags( @{ $pairedAlignment->[0] } );
        $AS += $auxTags_fwd{AS}{v};
    }
    if ( defined( $pairedAlignment->[1] ) ) {
        my %auxTags_rev = samGet_auxTags( @{ $pairedAlignment->[1] } );
        $AS += $auxTags_rev{AS}{v};
    }

    return ($AS);
}

sub filterPEalignments {
    my $alignments_ = shift;

    my $filteredAlignments;
    my $a;

    # filter alignments on C2T strand
    $a = 0;
    foreach my $pairedAlignment ( @{ $alignments_->{C2T}->{samFields} } ) {
        my $sortedPairedAlignment = sortPE($pairedAlignment);

        my $softClippedBases = getSoftClippedPE($sortedPairedAlignment);
        my $fwdMateMappedLength =
          length( $sortedPairedAlignment->[0]->[9] ) - $softClippedBases->[0]->[0] - $softClippedBases->[0]->[1];
        my $revMateMappedLength =
          length( $sortedPairedAlignment->[1]->[9] ) - $softClippedBases->[1]->[0] - $softClippedBases->[1]->[1];

        if (
            ( ( $sortedPairedAlignment->[0]->[2] ne $sortedPairedAlignment->[1]->[2] ) )

            || (    ( ( $sortedPairedAlignment->[0]->[1] & 0x2 ) == 0 )
                 && ( ( $sortedPairedAlignment->[1]->[1] & 0x2 ) == 0 ) 
                 && ( !defined($allowDoveTail) ) )

            || (    ( ( $sortedPairedAlignment->[0]->[1] & 0x10 ) != 0 )
                 && ( ( $sortedPairedAlignment->[1]->[1] & 0x10 ) != 0 ) )

            || (    ( ( $sortedPairedAlignment->[0]->[1] & 0x10 ) == 0 )
                 && ( ( $sortedPairedAlignment->[1]->[1] & 0x10 ) == 0 ) )

            || (    ( $sortedPairedAlignment->[0]->[3] > $sortedPairedAlignment->[1]->[3] )
                 && ( !defined($allowDoveTail) ) )

            || (    ( $sortedPairedAlignment->[0]->[3] + $fwdMateMappedLength >
                      $sortedPairedAlignment->[1]->[3] + $revMateMappedLength )
                 && ( !defined($allowDoveTail) ) )
          )
        {
            &$debug( $sortedPairedAlignment->[0]->[0] . " : " . $sortedPairedAlignment->[1]->[0] . " map discordantly",
                     "Line:", __LINE__ );
            $mappingStatsC{discordantAL} += 2;
            next;
        }

        my $i = 0;    # walk through mates
        for ( 0 .. 1 ) {
            my $flag = $sortedPairedAlignment->[$i]->[1];
            my %auxTags = samGet_auxTags( @{$sortedPairedAlignment->[$i]} );
            if ( ( ( $flag & 0x10 ) != 0 ) && ( $i == 0 ) ) {
                $mappingStatsC{filteredAL}++;
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": maps in wrong direction: FWDrC2Trev",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }
            if ( ( ( $flag & 0x10 ) == 0 ) && ( $i == 1 ) ) {
                $mappingStatsC{filteredAL}++;
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": maps in wrong direction: REVrC2Tfwd",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }
            if ( defined($sortedPairedAlignment->[$i]) && ( checkEditDist($auxTags{NM}{v}) || checkMisMatchRatio($auxTags{XM}{v}, $fwdMateMappedLength) ) && ( $i == 0) ) {
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": to many missmatches or edit distance exceeded",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }
            if ( defined($sortedPairedAlignment->[$i]) && ( checkEditDist($auxTags{NM}{v}) || checkMisMatchRatio($auxTags{XM}{v}, $revMateMappedLength) ) && ( $i == 1) ) {
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": to many missmatches or edit distance exceeded",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }            
            $i++;
        }
        if ( defined( $sortedPairedAlignment->[0] ) || defined( $sortedPairedAlignment->[1] ) ) {
            checkPairedAlignment($sortedPairedAlignment);
            $filteredAlignments->{C2T}->{samFields}->[$a] = [ @{$sortedPairedAlignment} ];
            $filteredAlignments->{C2T}->{mapData}->[$a] =
              [ $softClippedBases, $fwdMateMappedLength, $revMateMappedLength ];
            $a++;
        }
    }
    $filteredAlignments->{C2T}->{NH} = $a;

    # filter alignments on G2A strand
    $a = 0;
    foreach my $pairedAlignment ( @{ $alignments_->{G2A}->{samFields} } ) {
        my $sortedPairedAlignment = sortPE($pairedAlignment);

        my $softClippedBases = getSoftClippedPE($sortedPairedAlignment);
        my $fwdMateMappedLength =
          length( $sortedPairedAlignment->[0]->[9] ) - $softClippedBases->[0]->[0] - $softClippedBases->[0]->[1];
        my $revMateMappedLength =
          length( $sortedPairedAlignment->[1]->[9] ) - $softClippedBases->[1]->[0] - $softClippedBases->[1]->[1];

        if (
            ( ( $sortedPairedAlignment->[0]->[2] ne $sortedPairedAlignment->[1]->[2] ) )

            || (    ( ( $sortedPairedAlignment->[0]->[1] & 0x2 ) == 0 )
                 && ( ( $sortedPairedAlignment->[1]->[1] & 0x2 ) == 0 )
                 && ( !defined($allowDoveTail) ) )

            || (    ( ( $sortedPairedAlignment->[0]->[1] & 0x10 ) != 0 )
                 && ( ( $sortedPairedAlignment->[1]->[1] & 0x10 ) != 0 ) )

            || (    ( ( $sortedPairedAlignment->[0]->[1] & 0x10 ) == 0 )
                 && ( ( $sortedPairedAlignment->[1]->[1] & 0x10 ) == 0 ) )

            || (    ( $sortedPairedAlignment->[0]->[3] < $sortedPairedAlignment->[1]->[3] )
                 && ( !defined($allowDoveTail) ) )

            || (    ( $sortedPairedAlignment->[0]->[3] + $fwdMateMappedLength <
                      $sortedPairedAlignment->[1]->[3] + $revMateMappedLength )
                 && ( !defined($allowDoveTail) ) )
          )
        {
            &$debug( $sortedPairedAlignment->[0]->[0] . " : " . $sortedPairedAlignment->[1]->[0] . " map discordantly",
                     "Line:", __LINE__ );
            $mappingStatsC{discordantAL} += 2;
            next;
        }

        my $i = 0;    # walk through mates
        for ( 0 .. 1 ) {
            my $flag = $sortedPairedAlignment->[$i]->[1];
            my %auxTags = samGet_auxTags( @{$sortedPairedAlignment->[$i]} );
            if ( ( ( $flag & 0x10 ) != 0 ) && ( $i == 1 ) ) {
                $mappingStatsC{filteredAL}++;
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": maps in wrong direction: REVrG2Arev",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }
            if ( ( ( $flag & 0x10 ) == 0 ) && ( $i == 0 ) ) {
                $mappingStatsC{filteredAL}++;
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": maps in wrong direction: FWDrG2Afwd",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }
            if ( defined($sortedPairedAlignment->[$i]) && ( checkEditDist($auxTags{NM}{v}) || checkMisMatchRatio($auxTags{XM}{v}, $fwdMateMappedLength) ) && ( $i == 0) ) {
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": to many missmatches or edit distance exceeded",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }
            if ( defined($sortedPairedAlignment->[$i]) && ( checkEditDist($auxTags{NM}{v}) || checkMisMatchRatio($auxTags{XM}{v}, $revMateMappedLength) ) && ( $i == 1) ) {
                &$debug( $sortedPairedAlignment->[$i]->[0] . ": to many missmatches or edit distance exceeded",
                         "Line:", __LINE__ );
                $sortedPairedAlignment->[$i] = undef;
            }
            $i++;
        }
        if ( defined( $sortedPairedAlignment->[0] ) || defined( $sortedPairedAlignment->[1] ) ) {
            checkPairedAlignment($sortedPairedAlignment);
            $filteredAlignments->{G2A}->{samFields}->[$a] = [ @{$sortedPairedAlignment} ];
            $filteredAlignments->{G2A}->{mapData}->[$a] =
              [ $softClippedBases, $fwdMateMappedLength, $revMateMappedLength ];
            $a++;
        }
    }
    $filteredAlignments->{G2A}->{NH} = $a;

    return ($filteredAlignments);
}

sub checkPairedAlignment {
    my $sortedPairedAlignment = shift;

    if ( !defined( $sortedPairedAlignment->[1] ) ) {
        unpairAlignment( $sortedPairedAlignment->[0] );
    }
    if ( !defined( $sortedPairedAlignment->[0] ) ) {
        unpairAlignment( $sortedPairedAlignment->[1] );
    }

}

sub unpairAlignment {
    my $samFields = shift;

    if ( ( $samFields->[1] & 0x1 ) == 1 ) {
        $samFields->[1] -= 2  if ( ( $samFields->[1] & 0x2  ) == 2  );    # read NOT mapped in proper pair
        $samFields->[1] += 8  if ( ( $samFields->[1] & 0x8  ) != 8  );    # mate unmapped
        $samFields->[1] -= 32 if ( ( $samFields->[1] & 0x20 ) == 32 );    # remove mate reversed
    }

    $samFields->[6] = "*";                                             # edit RNEXT
    $samFields->[7] = 0;                                               # edit PNEXT
    $samFields->[8] = 0;                                               # edit TLEN
}

sub checkMisMatchRatio {
    my $misMatches = $_[0];
    my $mapLength  = $_[1];

    my $misMatchRatio = $misMatches / $mapLength;
    if ( $misMatchRatio > $max_mm_ratio ) {
        return (1);
    }

    return (0);
}

sub checkEditDist {
    my $editDist  = $_[0];

   if ( $editDist > $max_edit_dist ) {
        return (1);
   }

    return (0);
}

sub get_auxTagsFinal {
    my $auxTags = shift;

    my @reportTags = qw(AS NM MD NH HI YG YR);
    my @finalTags;

    foreach my $tag (@reportTags) {
        if ( exists( $auxTags->{$tag} ) ) {
            my $finalTag = $tag . ":" . $auxTags->{$tag}->{t} . ":" . $auxTags->{$tag}->{v};
            push( @finalTags, $finalTag );
        }
    }

    return (@finalTags);
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

    if ( $alignMode ne 'EndToEnd' ) {
        my @tmpC = $alignment->[5] =~ /^(\d+H)?((\d+)S)?((\d+)[MIDNP=X])*((\d+)S)?(\d+H)?$/;
        if ( defined( $tmpC[2] ) ) {
            $softClipped->[0] = $tmpC[2];
        }
        if ( defined( $tmpC[6] ) ) {
            $softClipped->[1] = $tmpC[6];
        }
    }
    return ($softClipped);
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

    my $trimL  = int( $mateOverlap / 2 );
    my $trimR  = $mateOverlap - $trimL;
    my $trimIL = 0;                         # trimmed leading I in CIGAR after planned trimL
    my $trimIR = 0;                         # trimmed leading I in CIGAR after planned trimR

    # forward mapped mates
    if (
         ( $trimR > 0 )
         && (    ( ( $auxTags->{YG}->{v} eq 'C2T' ) && ( $auxTags->{YR}->{v} eq 'C2T' ) )
              || ( ( $auxTags->{YG}->{v} eq 'G2A' ) && ( $auxTags->{YR}->{v} eq 'G2A' ) ) )
      )
    {

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
            my $alignedSeqLengthSoftClipped =
              length( $alignment->[9] ) - $softClippedBases->[0] - $softClippedBases->[1] - $trimR;

            $auxTags->{MD}->{v} = updateMDtag( $auxTags->{MD}->{v}, $alignment->[5], 0, $trimR );
            $auxTags->{NM}->{v} = updateNMtag( $auxTags->{MD}->{v}, $alignedSeqLengthSoftClipped );

            ( $alignment->[5], $alignment->[3], $trimIL, $trimIR ) =
              softClipCigar( $alignment->[5], 0, $trimR, $alignment->[3] );

        }

        $softClippedBases->[1] += ( $trimR + $trimIR );
    }

    # reversed mapped mates
    if (
         ( $trimL > 0 )
         && (    ( ( $auxTags->{YG}->{v} eq 'C2T' ) && ( $auxTags->{YR}->{v} eq 'G2A' ) )
              || ( ( $auxTags->{YG}->{v} eq 'G2A' ) && ( $auxTags->{YR}->{v} eq 'C2T' ) ) )
      )
    {

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
            my $alignedSeqLengthSoftClipped =
              length( $alignment->[9] ) - $softClippedBases->[0] - $softClippedBases->[1] - $trimL;

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

        # &$debug("ORI:", $cigar, "L", $trimLori);
        $cigar = $new_cigar;

        # &$debug("NEW:", $cigar);

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

        # &$debug("ORI:", $cigar, "R", $trimRori);
        $cigar = $new_cigar;

        #&$debug("NEW:", $cigar);

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
        $convfq_FH = $hisat2_fwd_convfq_FH;
    }
    elsif ( $conv eq 'G2A' ) {
        $conversion = sub { $_[0] =~ tr/G/A/; };
        $convfq_FH = $hisat2_rev_convfq_FH;
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

        $fq_rec{id} =~ s/([ \t]+.+)?$/\/1/;

        if ($forceDirectionality) {
            if ( !checkDirectionality( $fq_rec{seq}, 0 ) ) {
                &$debug( "Skipping read: " . substr( $fq_rec{id}, 1, -2 ) . " - directionality unassured" );
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

    if ( ( $readNr - $IQCfiltered - $DIRfiltered ) < 1 ) {
        say STDOUT "No reads to map: all filtered!";
        exit(1);
    }

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

        $FWDfq_rec{id} =~ s/(([ \t]+.+)?$)|(\/1$)/\/1/;
        $REVfq_rec{id} =~ s/(([ \t]+.+)?$)|(\/2$)/\/2/;

        # check if reads are properly paired
        if ( substr( $FWDfq_rec{id}, 0, -2 ) ne substr( $REVfq_rec{id}, 0, -2 ) ) {
            &$debug( substr( $FWDfq_rec{id}, 0, -2 ), "<--->", substr( $REVfq_rec{id}, 0, -2 ), "Line:", __LINE__ );
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
                    &$debug( "Skipping read pair: " . substr( $FWDfq_rec{id}, 1, -2 ) . " - directionality unassured" );
                    $DIRfiltered++;
                    $readNr++;
                    next;
                }
            }
        }

        $C2T->( $FWDfq_rec{seq} );
        $G2A->( $REVfq_rec{seq} );

        $hisat2_fwd_convfq_FH->print( $FWDfq_rec{id} . "\n", $FWDfq_rec{seq}, $FWDfq_rec{id2}, $FWDfq_rec{qs} );
        $hisat2_rev_convfq_FH->print( $REVfq_rec{id} . "\n", $REVfq_rec{seq}, $REVfq_rec{id2}, $REVfq_rec{qs} );

        $readNr++;
    }

    $hisat2_fwd_convfq_FH->flush();
    $hisat2_rev_convfq_FH->flush();

    close($FWDfqFH);
    close($REVfqFH);

    if ($fixDirectionaly) {    # EXPERIMENTAL: fixing directionality
        say STDOUT "Fixed directionality of " . $DIRfiltered . " read pairs";
    }
    else {
        say STDOUT "Skipped " . $DIRfiltered . " read pairs due to unassured directionality";
    }
    say STDOUT "Filtered " . $IQCfiltered . " read pairs not passing Illuminal QC" if ($useIlluminaQC);

    if ( ( $readNr - $IQCfiltered - $DIRfiltered ) < 1 ) {
        say STDOUT "No reads to map: all filtered!";
        exit(1);
    }
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
    elsif ( ( $readID =~ /\@_meRanGh_FQ_FILE_END_1/ ) ) {
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

sub reverseComp {
    my $seq = shift;

    my $revcomp = reverse($seq);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;

    return $revcomp;
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
                print $fifoout "\@_meRanGh_FQ_FILE_SWITCH_1\n";
                print $fifoout "_meRanGh_FQ_FILE_SWITCH_2\n";
                print $fifoout "\+_meRanGh_FQ_FILE_SWITCH_3\n";
                print $fifoout "_meRanGh_FQ_FILE_SWITCH_4\n";
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

    print $fifoout "\@_meRanGh_FQ_FILE_END_1\n";
    print $fifoout "_meRanGh_FQ_FILE_END_2\n";
    print $fifoout "\+_meRanGh_FQ_FILE_END_3\n";
    print $fifoout "_meRanGh_FQ_FILE_END_4\n";

    $fifoout->close();
    undef($fifoout);

    return;
}

sub mBiasCounter {

    my ( $mBiasData, $read, $qs, $leftClip, $rightClip, $mateNr ) = @_;

    my $pos;
    my $idxLength  = length($read) - 1;
    my $alignStart = $leftClip;
    my $alignEnd   = $idxLength - $rightClip;
    my $offset     = $alignStart;

    my @qualValues = map { $_ - $qsOffset } unpack( '(W)*', $qs );

    if ( $mateNr == 0 ) {
        map { $mBiasData->{mBiasReadDataF}[$_] += 0 } 0 .. $idxLength;            # avoid undefined values;
        map { $mBiasData->{mBiasDataF}[$_]     += 0 } 0 .. $idxLength;            # avoid undefined values;
        map { $mBiasData->{mBiasDataFhq}[$_]   += 0 } 0 .. $idxLength;            # avoid undefined values;
        map { $mBiasData->{mBiasReadDataF}[$_] += 1 } $alignStart .. $alignEnd;

        while (1) {
            $pos = index( $read, 'C', $offset );
            last if ( ( $pos < 0 ) || ( $pos > $alignEnd ) );
            ++$mBiasData->{mBiasDataF}[$pos];
            ++$mBiasData->{mBiasDataFhq}[$pos] unless ( $qualValues[$pos] < $mbQSt );
            $offset = $pos + 1;
            last if ( $offset >= $alignEnd );
        }
    }
    else {
        map { $mBiasData->{mBiasReadDataR}[$_] += 0 } 0 .. $idxLength;            # avoid undefined values;
        map { $mBiasData->{mBiasDataR}[$_]     += 0 } 0 .. $idxLength;            # avoid undefined values;
        map { $mBiasData->{mBiasDataRhq}[$_]   += 0 } 0 .. $idxLength;            # avoid undefined values;
        map { $mBiasData->{mBiasReadDataR}[$_] += 1 } $alignStart .. $alignEnd;

        while (1) {
            $pos = index( $read, 'G', $offset );
            last if ( ( $pos < 0 ) || ( $pos > $alignEnd ) );
            ++$mBiasData->{mBiasDataR}[$pos];
            ++$mBiasData->{mBiasDataRhq}[$pos] unless ( $qualValues[$pos] < $mbQSt );
            $offset = $pos + 1;
            last if ( $offset >= $alignEnd );
        }
    }

    return;
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
    unlink $bam;    # we don't want the unsorted version!
    return $sorted . ".bam";
}

sub index_bam {
    my $bam = shift;
    print STDERR "Indexing BAM...\n";
    Bio::DB::Bam->index_build($bam);
    return;
}

sub makeBedGraph {
    my $bam = shift;

    my ( $fname, $fpath, $fsuffix ) = fileparse( $bam, "_sorted.bam" );

    my $bedGraph = $fpath . "/" . $fname . ".bg";
    unlink($bedGraph);

    my $bgFH = IO::File->new( $bedGraph, O_RDWR | O_CREAT | O_TRUNC ) || die( $bedGraph . ": " . $! );

    print STDERR "Creating BedGraph: " . $bedGraph . " ... This can take a while ...\n";
    $sam = Bio::DB::Sam->new( -bam => $bam );
    $sam->coverage2BedGraph($bgFH);
    $bgFH->close();
    undef($bgFH);
}

### make BedGraph in parallel
sub makeBedGraphParallel {
    my $bam = shift;

    $sam = Bio::DB::Sam->new( -bam => $bam, -split => 1 );
    $sam->max_pileup_cnt(40_000);    # libbam.a has set this to 8000 as default

    my ( $fname, $fpath, $fsuffix ) = fileparse( $bam, "_sorted.bam" );

    # parallel Jobs
    print STDERR "\nCreating BedGraph can take a while ...\n";
    print STDERR "Running $max_threads processes for bedGraph geneartion\n";
    my $pm = Parallel::ForkManager->new($max_threads);

    $pm->run_on_finish( \&bgChrDone );
    $pm->run_on_wait( \&bgTellStatus, 0.5 );

    # loop through the chromosomes
    $nrOfBGTargets = $sam->n_targets;
    for my $tid ( 0 .. $nrOfBGTargets - 1 ) {

        printf( "processing %i sequences: %i - %02.2f%% done ...\r", $nrOfBGTargets, $bgCount, $bgPCTdone );

        # run each chromosome in a separate process
        $pm->start and next;

        ### in child ###
        # clone the Bam object
        $sam->clone;

        # open chromosome wig files
        my ( $bgFile, $bgFH );
        $bgFile = nicePath($fpath . "/" . $fname . "_" . $tid . "_bg_tmp");
        $bgFH = IO::File->new( $bgFile, O_RDWR | O_CREAT | O_TRUNC ) || die( $bgFile . ": " . $! );

        # do the job
        bg_on_chromosome( $bgFH, $tid );

        # finished with this chromosome
        $bgFH->close;
        undef($bgFH);
        $pm->finish;
    }
    $pm->wait_all_children;

    # merge the files
    print STDERR "\nmerging chromosome bedGraph files\n";

    my $fnameSuffix = $fpath . "/" . $fname;
    my @bgFiles     = glob "$fnameSuffix\_*\_bg_tmp";
    die("bedGraph chunk files!\n") unless (@bgFiles);

    my $bgFile = nicePath($fpath . "/" . $fname . ".bg");

    my $bgFH = IO::File->new( $bgFile, O_RDWR | O_CREAT | O_TRUNC ) || die( $bgFile . ": " . $! );

    combineBGchunks( $bgFH, \@bgFiles );
    $bgFH->close();
    print STDERR "bedGraph created: " . $bgFile . "\n";

}

sub bg_on_chromosome {
    my $fh  = shift;
    my $tid = shift;

    my $header  = $sam->header;
    my $index   = $sam->bam_index;
    my $seqids  = $header->target_name;
    my $lengths = $header->target_len;
    my $b       = $sam->bam;

    my $convertor = getConvertor($bgScale);

    my $seqid = $seqids->[$tid];
    my $len   = $lengths->[$tid];

    my $sec_start = -1;
    my $last_val  = -1;

    for ( my $start = 0 ; $start <= $len ; $start += DUMP_INTERVAL ) {
        my $end = $start + DUMP_INTERVAL;
        $end = $len if $end > $len;
        my $coverage = $index->coverage( $b, $tid, $start, $end );
        for ( my $i = 0 ; $i < @$coverage ; $i++ ) {
            if ( $last_val == -1 ) {
                $sec_start = 0;
                $last_val  = $coverage->[$i];
            }
            if ( $last_val != $coverage->[$i] ) {
                print $fh $seqid, "\t", $sec_start, "\t", $start + $i, "\t", &$convertor($last_val), "\n"
                  unless $last_val < $minBGcov;
                $sec_start = $start + $i;
                $last_val  = $coverage->[$i];
            }
            elsif ( $start + $i == $len - 1 ) {
                print $fh $seqid, "\t", $sec_start, "\t", $start + $i, "\t", &$convertor($last_val), "\n"
                  unless $last_val < $minBGcov;
            }
        }
    }
}

sub bgTellStatus {
    printf( "processing %i sequences: %02.2f%% done ...\r", $nrOfBGTargets, $bgPCTdone );
}

sub bgChrDone {
    $bgPCTdone = ++$bgCount * 100 / $nrOfBGTargets;
    printf( "processing %i sequences: %02.2f%% done ...\r", $nrOfBGTargets, $bgPCTdone );
}

sub combineBGchunks {
    my $outFH   = shift;
    my $inFiles = shift;

    foreach my $chunkF ( @{$inFiles} ) {
        my $in = IO::File->new( $chunkF, "r" );
        while (<$in>) { $outFH->print($_) }
        $in->close();
        unlink($chunkF);
    }
}

sub getConvertor {
    my $ctype = shift;

    my $convertor;

    if ( $ctype eq 'LOG2' ) {
        $convertor = sub {
            return 0 if $_[0] == 0;
            return log( $_[0] ) / LOG2;
          }
    }
    elsif ( $ctype eq 'LOG10' ) {
        $convertor = sub {
            return 0 if $_[0] == 0;
            return log( $_[0] ) / LOG10;
          }
    }
    else {
        $convertor = sub {
            return $_[0];
          }
    }
    return $convertor;
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
      map {
        sprintf( "%.1f", ( ( $mBiasData->[$_] / ( ( $mBiasReadData->[$_] != 0 ) ? $mBiasReadData->[$_] : 1 ) ) * 100 ) )
      } 0 .. $#$mBiasReadData;

    my @mCratioHQ =
      map {
        sprintf( "%.1f",
                 ( ( $mBiasDataHQ->[$_] / ( ( $mBiasReadData->[$_] != 0 ) ? $mBiasReadData->[$_] : 1 ) ) * 100 ) )
      } 0 .. $#$mBiasReadData;

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
    my $font = [ 'LiberationSans-Regular', 'verdana', 'arial', "gdMediumBoldFont" ];

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

############################## END sssssssttttttttttttaaaaaaaaaaaaaaa


sub extractSpliceSites {
    my $genes;
    my $trans;
    
    my $gtfFH = IO::File->new( $GTF, O_RDONLY | O_EXCL ) || die( $GTF . ": " . $! );
    my $line;
    
    while ( defined( $line = $gtfFH->getline() ) ) {
        chomp($line);

        next if $line =~ /^#/;

        if ( index( $line, '#' ) != -1 ) {
            $line = ( split( '#', $line ) )[0];
            chomp($line);
            $line = trim($line);
        }
        my ( $chrom, $source, $feature, $left, $right, $score, $strand, $frame, $values ) = split( /\t/, $line );

        next if ( ( $feature ne 'exon' ) || ( $left >= $right ) );

        my $values_dict;

        foreach my $a ( split( /;\s?/, $values ) ) {
            my ( $attr, $val ) = split( /\s+/, $a, 2 );
            $val =~ tr/\"//d;
            $values_dict->{$attr} = $val;
        }

        next if ( ( !defined( $values_dict->{'gene_id'} ) ) || ( !defined( $values_dict->{'transcript_id'} ) ) );

        my $transcript_id = $values_dict->{'transcript_id'};

        if ( !defined( $trans->{$transcript_id} ) ) {
            $trans->{$transcript_id} = [ $chrom, $strand, [ [ $left, $right ] ] ];
            push( @{ $genes->{ $values_dict->{'gene_id'} } }, $transcript_id );
        }
        else {
            push( @{ $trans->{$transcript_id}->[2] }, [ $left, $right ] );
        }
    }
    
    $gtfFH->close();
    undef($gtfFH);

    foreach my $tran ( keys( %{$trans} ) ) {
        my $t       = $trans->{$tran};
        my $chrom   = $t->[0];
        my $strand  = $t->[1];
        my $exons   = $t->[2];
        my @exons_s = sort { ( $a->[0] <=> $b->[0] ) || ( $a->[1] <=> $b->[1] ) } @{$exons};

        my $tmp_exons = [ $exons_s[0] ];

        for ( my $i = 1 ; $i < scalar(@exons_s) ; $i++ ) {
            if ( ( $exons_s[$i]->[0] - $tmp_exons->[-1]->[1] ) <= 5 ) {
                $tmp_exons->[-1]->[1] = $exons_s[$i]->[1];
            }
            else {
                push( @$tmp_exons, $exons_s[$i] );
            }
        }

        $trans->{$tran} = [ $chrom, $strand, $tmp_exons ];
    }

    my @junctions;
    my $seen;

    foreach my $tran ( keys( %{$trans} ) ) {
        my $t      = $trans->{$tran};
        my $chrom  = $t->[0];
        my $strand = $t->[1];
        my $exons  = $t->[2];

        my $key;
        for ( my $i = 1 ; $i < scalar( @{$exons} ) ; $i++ ) {
            $key = join( '#', ( $chrom, $exons->[ $i - 1 ][1], $exons->[$i][0], $strand ) );
            push( @junctions, [ $chrom, $exons->[ $i - 1 ][1], $exons->[$i][0], $strand ] ) unless $seen->{$key};
            $seen->{$key} = 1;
        }
    }
    
    my @junctions_s = sort { ( $a->[0] cmp $b->[0] ) ||  ( $a->[1] <=> $b->[1] ) || ( $a->[2] <=> $b->[2] ) } @junctions;

    my $spliceSitesFH   = IO::File->new( $spliceSitesFile, O_RDWR | O_CREAT | O_TRUNC ) || die( $spliceSitesFile . ": " . $! );

    foreach my $j (@junctions_s) {
        my $chrom  = $j->[0];
        my $left   = $j->[1] - 1;
        my $right  = $j->[2] - 1;
        my $strand = $j->[3];

        $spliceSitesFH->print( join( "\t", ( $chrom, $left, $right, $strand ) ) . "\n" );
    }
    
    $spliceSitesFH->close();
    undef($spliceSitesFH);
}

sub trim { $_[0] =~ s/^\s+|\s+$//g; return $_[0] }

sub lockF {
    my ($fh) = $_[0];
    flock( $fh, LOCK_EX ) or die "Cannot lock resultfile - $!\n";
}

sub unlockF {
    my ($fh) = $_[0];
    flock( $fh, LOCK_UN ) or die "Cannot unlock resultfile - $!\n";
}

sub mergeMappedSeqHash {
    my ( $a, $b ) = @_;

    my $c;

    for my $href ( $a, $b ) {
        while ( my ( $k, $v ) = each %$href ) {
            while ( my ( $sk, $sv ) = each %{ $href->{$k} } ) {
                $c->{$k}->{$sk} = $sv;
            }
        }
    }

    return $c;
}


sub getHISAT2version {

    my $versionCheckCmd = $hisat2_cmd . " --version";
    my $pid = open( my $VERSIONCHECK, "-|", $versionCheckCmd )
      || die( "\n\b\aERROR could not determine the HISAT2 version: ", $! );

    $VERSIONCHECK->autoflush();
    my $versionStr = $VERSIONCHECK->getline();
    waitpid( $pid, 0 );
    close($VERSIONCHECK);

    if ( !defined($versionStr) ) {
        die( "Could not get HISAT2 version by running: " . $versionCheckCmd );
    }
    chomp($versionStr);

    $versionStr =~ s/(.+)(hisat2-align-s\ version\ )(.+)/$3/;
    return ($versionStr);
}

sub checkHiSat2 {
    $hisat2Version = getHISAT2version();
    if ( !exists( $supportedHiSat2Versions{$hisat2Version} ) ) {
        say STDERR "\nERROR: HISAT2 version "
          . $hisat2Version
          . " not supported please install one of the following versions:\n";
        foreach my $sv ( keys(%supportedHiSat2Versions) ) {
            say STDERR $sv;
        }
        exit(1);
    }
    return (1);
}

sub checkFastQsort {
    my $canSort = getFastqSortIDNcap();
    if ( !($canSort) ) {
        say STDERR
          "\nERROR: installed fastq-sort has not the required alphanumerically read identifier sorting capability.\n"
          . "Please install the provided fastq-sort tool or get the latest version from: https://github.com/dcjones/fastq-tools\n";
        exit(1);
    }
    return (1);
}

sub nicePath {
    my $rawPath = shift;
    
    $rawPath =~ s/\/+/\//g;
    return($rawPath);
}

sub usage {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self <runmode> [-h] [-m] [--version]

Required <runmode> any of:
    mkbsidx     :   Generate the HISAT2 BS index.
    align       :   Align bs reads to a reference genome.

Options:
    --version   :   Print the program version and exit.
    -h|help     :   Print the program help information.
    -m|man      :   Print a detailed documentation.
    
EOF
}

sub usage_mkbsidx {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self mkbsidx [-fa] [-id] [-h] [-m]

Required all of :
    -fa|fasta           : Fasta file(s) to use for BS index generation.
                          Use a comma separated file list or expression
                          (?, *, [0-9], [a-z], {a1,a2,..an}) if more than one
                          fasta file. If using an expression pattern, please put
                          single quotes arround the -fa argument, e.g:

                              -fa '/genome/chrs/chr[1-8].fa'

    -id|bsidxdir        : Directory where to store the BS index.

Options:                          
    -hisat2build|hs2b   : Path to the HiSat2 "hisat2-build" program.
                          (default: hisat2-build from the meRanTK installation
                           or systems PATH)
                          
    -t|threads          : number of CPUs/threads to run

    --version           : Print the program version and exit.
    -h|help             : Print the program help information.
    -m|man              : Print a detailed documentation.

EOF
}

sub usage_align {
    my $self = basename($0);

    print STDERR<<EOF
USAGE: $self align [-f|-r] [-id] [-h] [-m]

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

Options:
    -GTF                  : GTF gene model annotations and/or known
                            transcripts for building a transcriptome index.

    -illuminaQC|-iqc      : Filter reads that did not pass the Illumina QC.
                            Only relevant if you have Illumina 1.8+ reads.
                            (default: not set)

    -forceDir|-fDir       : Filter reads that did not pass did not pass the internal
                            directionality check:
                                FWDreads: #C > #G && #C > #T && #A > #G)
                                REVreads: #G > #C && #T > #C && #G > #A)
                            (default: not set)

    -first|-fn            : Process only this many reads per input fastq file
                            (default: process all reads)

    -outdir|-o            : Directory where results get stored
                            (default: current directory)

    -sam|-S               : Name of the SAM file for uniq and resolved alignments
                            (default: meRanGh_[timestamp].sam )

    -unalDir|-ud          : Directory where unaligned reads get stored
                            (default: outdir)

                            Note: if -hisat2un|-un is not set unaligned reads
                            will not get stored

    -threads|-t           : Use max. this many CPUs to process data
                            (default: 1)

    -hisat2cmd|-hisat2    : Path to hisat2
                            (default: use hisat2 from the meRanTK installation or in your system PATH)

    -id|-bsidxdir         : Path to bsindex directory created in 'mkbsidx' runMode.
                            This directory holds the '+' and '-' strand bs index
                            (default: use BS_HISAT2_IDX environment variable)

    -samMM|-MM            : Save multimappers? If set multimappers will be stored
                            in SAM format '\$sam_multimappers.sam'
                            (default: not set)

    -ommitBAM|-ob         : Do not create an sorted and indexed BAM file
                            (default: not set)

    -deleteSAM|-ds        : Delete the SAM files after conversion to BAM format
                            (default: not set)

    -deleteBAMus|-dbus    : Delete the unsorted BAM files after sorting BAM.
                            (default: not set)
                            
    -hisat2un|-un        : Report unaligned reads. See also -unalDir|-ud
                            (default: not set)

    -mbiasplot|-mbp       : Create an m-bias plot, that shows potentially biased
                            read positions
                            (default: not set)

    -mbiasQS|-mbQS        : Quality score for a high quality m-bias plot. This plot
                            considers only basecalls with a quality score equal or
                            higher than specified by this option.
                            (default: 30)

    -mkbg|-bg             : Generate a BEDgraph file from the aligned reads.
                            ! This can take a while !
                            (default: not set)

    -minbgCov|-mbgc       : If '-bg' is set, '-mbgc' defines the minimum coverage that we
                            should consider in the BEDgraph output?
                            (default: 1)

    -fixMateOverlap|-fmo  : The sequenced fragment and read lengths might be such that
                            alignments for the two mates from a pair overlap each other.
                            If '-fmo' is set, deduplicate alignment subregions that are
                            covered by both, forward and reverse, reads of the same
                            read pair. Only relevant for paired end reads.
                            (default: not set)

    -hardClipMO|-hcmo     : If '-fmo' is set, hardclip instead of softclip the overlaping
                            sequence parts.
                            (default: not set = softclip)

    -dovetail|-dt         : If '-dt' is set, "dovetailing" read pairs in paired end mode
                            are allowed and will not be discordant.
                            (default: not set = dovetailing reads are not aligned)

    -bgScale|-bgs         : Generate a BEDgraph in log [log2|log10] scale
                            (default: not set, no scaling)

    -max-edit-dist|-e     : Maximum edit distance to allow for a valid alignment
                            (default: 10)

    -max-mm-rate|-mmr     : Maximum mismatch ratio (mismatches over read alignement length)
                            [0 <= mmr < 1]
                            (default: 0.05)

    -hisat2_trim5
                          : see hisat2 manual for -trim5 option
                            (default: 0)    

    -hisat2_trim3
                          : see hisat2 manual for -trim3 option
                            (default: 0)    
    -hisat2_n-ceil
                          : see hisat2 manual for -n-ceil option
                            (default: L,0,0.15)    
    -hisat2_ignore-quals
                          : see hisat2 manual for -ignore-quals option
                            (default: not set)    

    -hisat2_nofw
                          : see hisat2 manual for -nofw option
                            (default: not set)    

    -hisat2_norc
                          : see hisat2 manual for -norc option
                            (default: not set)    

    -hisat2_pen-canonsplice
                          : see hisat2 manual for -pen-canonsplice option
                            (default: 0)    

    -hisat2_pen-noncanonsplice
                          : see hisat2 manual for -pen-noncanonsplice option
                            (default: 12)    

    -hisat2_pen-canintronlen
                          : see hisat2 manual for -pen-canintronlen option
                            (default: G,-8,1)    

    -hisat2_pen-noncanintronlen
                          : see hisat2 manual for -pen-noncanintronlen option
                            (default: G,-8,1)    

    -hisat2_min-intronlen
                          : see hisat2 manual for -min-intronlen option
                            (default: 20)        

    -hisat2_max-intronlen
                          : see hisat2 manual for -may-intronlen option
                            (default: 500000)        

    -hisat2_novel-splicesite-outfile
                          : see hisat2 manual for -novel-splicesite-outfile option
                            (default: not set)        

    -hisat2_novel-splicesite-infile
                          : see hisat2 manual for -novel-splicesite-infile option
                            (default: not set)        

    -hisat2_no-temp-splicesite
                          : see hisat2 manual for -no-temp-splicesite option
                            (default: not set)        

    -hisat2_no-spliced-alignment
                          : see hisat2 manual for -no-spliced-alignment option
                            (default: not set)        

    -hisat2_rna-strandness
                          : see hisat2 manual for -rna-strandness option
                            (default: FR)        

    -hisat2_transcriptome-mapping-only
                          : see hisat2 manual for -transcriptome-mapping-only option
                            (default: not set)        

    -hisat2_remove-chrname
                          : see hisat2 manual for -remove-chrname option
                            (default: not set)        

    -hisat2_add-chrname
                          : see hisat2 manual for -add-chrname option
                            (default: not set)        

    -hisat2_omit-sec-seq
                          : see hisat2 manual for -omit-sec-seq option
                            (default: not set)        

    -hisat2_offrate
                          : see hisat2 manual for -offrate option
                            (default: not set)        

    -hisat2_ma
                          : see hisat2 manual for -ma option
                            (default: not set)        
    -hisat2_mp
                          : see hisat2 manual for -mp option
                            (default: 6,2)        

    -hisat2_sp
                          : see hisat2 manual for -sp option
                            (default: 1,2)        

    -hisat2_np
                          : see hisat2 manual for -np option
                            (default: 1)        

    -hisat2_rdg
                          : see hisat2 manual for -rdg option
                            (default: 5,3)        

    -hisat2_rfg
                          : see hisat2 manual for -rfg option
                            (default: 5,3)        

    -hisat2_min-score
                          : see hisat2 manual for -min-score option
                            (default: L,0.0,-0.2)        

    -hisat2_k
                          : see hisat2 manual for -k option
                            (default: not set)        

    -hisat2_minins
                          : see hisat2 manual for -minins option
                            (default: 0)        

    -hisat2_maxins
                          : see hisat2 manual for -maxins option
                            (default: 500)        

    -hisat2_met-file
                          : see hisat2 manual for -met-file option
                            (default: not set)        

    --version             : Print the program version and exit.
    -h|help               : Print the program help information.
    -m|man                : Print a detailed documentation.

    -debug|-d             : Print some debugging information.

EOF
}

__END__

=head1 NAME

meRanGh - RNA bisulfite short read mapping to genome

=head1 SYNOPSIS

=head2 Index generation

## Generate a genome database BS index for the aligner

 meRanGh mkbsidx \
    -t 4 \
    -fa mm10.chr1.fa,mm10.chr2.fa,[...] \
    -id /export/data/mm10/BSgenomeIDX \

 Generates an index for bisulfite mapping strand specific RNA-BSseq reads to a 
 genome database provided as fasta file(s)) (e.g. from genome assembly mm10).
 The indexer will run with max. (-t) 4 threads.

 The example above assumes that the HiSat2 index builder “hisat2-build” and
 command is found in the systems path “$PATH” or “hisat2-build”
 from the meRanTK shipped third party programs is used (see Installation 2.2, 2.3).
 Alternatively, the path to “hisat2-build” can be specified using
 the command line options (-hs2b).
 
=head2 Align directed/strand specific RNA-BSseq short reads to a genome

=over 2
 
### Single End reads

 meRanGh align \
 -o ./meRanGh2Result \
 -f ./FastqDir/01.fastq,./FastqDir/02.fastq,./FastqDir/03.fastq \
 -t 16 \
 -S RNA-BSseq.sam \
 -ud ./meRanGhUnaligned \
 -un \
 -MM \
 -dovetail \
 -GTF /export/data/mm10/mm10.gtf \
 -mbp \
 -id /export/data/mm10/BSgenomeIDX \
 -bg \
 -mbgc 10
 
 The command above maps the reads from three diffent fastq files, separated by
 commas, to a genome using the index created as indicated in the previous section.
 The -GTF option indicates that the program should also search alignments spaning
 the known splicesites. meRanGh will use the index created in the "mkbsidx" run 
 mode which was stored under "/export/data/mm10/BStranscriptomeIDX" (-id) (see above).

 The mapping process will use (-t) 16 cpus to search for valid alignments, from
 which the best one will be stored in the (-S) "RNA-BSseq.sam" result file. The
 program will save (-un) the unaligned reads in (-ud) the directory "meRanGhUnaligned":
 The alignments of multi mapper reads (-MM) will additionally be stored in a separate
 SAM file.
 The example above assumes that the path to "hisat2" was found in the systems PATH.
 Alternatively, this path can be specified by using commandline the appropriate
 option (-hisat2).
 The program will also create a Bedgraph (-bg) file for the alignments. Only positions
 with minimum 10 readcounts (-mbgc 10) will be recorded in the resulting Bedgraph.
 
 
### Paired End reads

 meRanGh align \
 -o ./meRanGh2Result \
 -f ./FastqDir/fwd01-paired.fastq,./FastqDir/fwd02-paired.fastq \
 -r ./FastqDir/rev01-paired.fastq,./FastqDir/rev02-paired.fastq \
 -t 16 \
 -S RNA-BSseq.sam \
 -ud ./meRanGhUnaligned \
 -un \
 -MM \
 -dovetail \
 -GTF /export/data/mm10/mm10.gtf \
 -mbp \
 -id /export/data/mm10/BSgenomeIDX \
 -bg \
 -mbgc 10

 When using paired end reads, one can specify the forward and reverse reads using
 the commandline options "-f" and "-r" respectively. Multiple files for each read
 direction files can be specified separated by commas. Not only the sort order of
 the forward and reverse reads has to be the same within the fastq files but also
 the order in which one specifies the forward and reverse read fastq files (see
 example aobve). Note: The paired fastq files may not have unpaired reads. If
 this is the case, one can use for example the "pairfq" (S. Evan Staton) tool to
 pair and sort the mates. In this example meRanGh allows for dovetailing 
 alignment of readpairs (-dovetail).

=back

=head1 DEPENDENCIES

 The program requires additional Perl modules and depending on your perl installation
 you might need do add the following modules:

 Bio::DB::Sam
 Parallel::ForkManager
 
 These modules should be availble via CPAN or depending on your OS via the package
 manager.

 Bio::DB:Sam requires the samtools libraries (version 0.1.10 or higher, version 1.0
 is not compSatible with Bio::DB::Sam, yet) and header files in order to compile
 successfully.
 
 Optional modules for generating m-bias plots:
 GD
 GD::Text::Align
 GD::Graph::lines
 
 In addition to these Perl modules a working installation of hisat2 (>= v.2.0.13)
 and Bowtie2 (>= v.2.2.3) is required.
 A patched fastq-sort program is included in the distribution. Or available from
 https://github.com/dcjones/fastq-tools

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
 
  mkbsidx : generate the genome database BS index for the aligner
  align   : align RNA-BSseq reads to the genome database.

=back

=head1 OPTIONS

=head2 general:

=over 2

=item -version

 Print the program version and exit.

=item -h|-help

 Print the program help information.

=item -m|-man

 Print a detailed documentation.

=back

=head2 mkbsidx:

=over 2

 run "meRanGh mkbsidx -h" to get a list of options

=back

=head2 align:

=over 2

 run "meRanGh align -h" to get a list of options

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
