#!/usr/bin/env perl
#/usr/local/bioinf/perl/perl-5.18.2/bin/perl -w
#
#  mkRefSeq2GeneMap.pl
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

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use IO::File;
use File::Basename;

my $VERSION = '1.0.0';
my $version;
my $help;

my $fastaDB = "";
my $mapFile = "";

GetOptions( 'fasta|f=s' => \$fastaDB,
            'map|m=s'   => \$mapFile, );

my $refSeqID       = "";
my $geneName       = "";
my $fastaSeq       = "";
my $geneName_old   = "";
my $refSeqID_old   = "";
my $geneLength_old = 0;
my %idHash;

usage() and exit(0) if ($help);
say $VERSION and exit(0) if ($version);

if ( !$fastaDB ) {
    usage() and exit(1);
}

my $fastaFH = IO::File->new( $fastaDB, O_RDONLY | O_EXCL ) || die( $fastaDB . ": " . $! );

my $mapFH;
if ($mapFile) {
    $mapFH = IO::File->new( $mapFile, O_RDWR | O_CREAT | O_TRUNC ) || die( $mapFile . ": " . $! );
}
else {
    open( $mapFH, ">-" );
}

while ( defined( my $fastaLine = $fastaFH->getline() ) ) {

# refSeq fastaLine example:
# >gi|568905436|ref|XM_006495550.1| PREDICTED: Mus musculus X Kell blood group precursor related family member 4 (Xkr4), transcript variant X1, mRNA

    chomp($fastaLine);
    if ( $fastaLine =~ /^\>/ ) {
        my @fields = split( '\|', $fastaLine );

        if (    ( $fields[0] ne '>gi' )
             || ( $fields[1] !~ /\d+/ )
             || ( $fields[2] ne 'ref' )
             || ( $fields[3] !~ /(XM|NM|NR|XR)(_)(\d+)(\.)(\d+)/ ) )
        {
            say STDERR "Skipping: unknown refSeq fasta header: " . $fastaLine;
            next;
        }

        $refSeqID = substr( ( split( ' ', $fastaLine ) )[0], 1 );

        $geneName = 'Unknown';

        my $desc = $fields[4];
        $desc =~ /(.+\()(.+)(\).+)/;
        $geneName = $2 if ( defined($2) );

        $geneLength_old = length($fastaSeq);

        if ( $geneLength_old != 0 ) {
            $fastaSeq                          = "";
            $idHash{$refSeqID_old}{'geneName'} = $geneName_old;
            $idHash{$refSeqID_old}{'length'}   = $geneLength_old;
            $mapFH->print( $refSeqID_old . "\t" . $geneName_old . "\t" . $geneLength_old . "\n" );
        }
    }
    else {
        $refSeqID_old = $refSeqID;
        $geneName_old = $geneName;
        $fastaSeq .= $fastaLine;
    }
}
$fastaFH->close();

# Flush last sequence
$geneLength_old = length($fastaSeq);

if ( $geneLength_old != 0 ) {
    $fastaSeq                          = "";
    $idHash{$refSeqID_old}{'geneName'} = $geneName_old;
    $idHash{$refSeqID_old}{'length'}   = $geneLength_old;
    $mapFH->print( $refSeqID_old . "\t" . $geneName_old . "\t" . $geneLength_old . "\n" );
}
$mapFH->close();


undef($mapFH);
undef($fastaFH);

sub usage {
    my $self = basename($0);

    print STDERR<<EOF
mkRefSeq2GeneMap.pl creates a RefSeq-ID to gene name (offical Gene Symbol) map file with transcript length
information in the following format:

                   gi|geneID|ref|refseqID|<tab>Genesymbol<tab>sequence length


As input any mRNA refSeq fasta file provided by the NCBI can be used, e.g:

  <ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.rna.fna.gz>

USAGE: $self [options] ][-h] [--version]

Required options any of:
    -fasta|-f             : FASTA file for which to create the map file.

Options:

    -map|-m               : Map file where to store the mapping
			   (default: STDOUT)

    --version             : Print the program version and exit.
    
    -h|help               : Print the program help information.
    
EOF
}
