#!/bin/perl

use strict ;
use warnings ;

# Run TCGA datas with given row ids
die "Usage: a.pl TCGA_bam_folder < row_ids" if ( @ARGV == 0 ) ;

my $BAMPath = $ARGV[0] ;

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
}

while ( <STDIN> )
{
	my @cols = split /\s+/ ;
	#my $file = $cols[1] ;
	#my $folder = $cols[0] ;
	my $BAMId = $cols[0] ;
	if ( -e "/liulab/galib/mouse_scRNAseq_margaret/data/Margaret_CD70/${BAMId}.bam" )
	{
		mkdir $BAMId ;
		system_call( "/liulab/galib/mouse_scRNAseq_margaret/tool/TRUST4/run-trust4".
			 " -b /liulab/galib/mouse_scRNAseq_margaret/data/Margaret_CD70/${BAMId}.bam".
			 " -t 8 -o /liulab/galib/mouse_scRNAseq_margaret/data/trust4/$BAMId/$BAMId".
			 " -f /liulab/galib/mouse_scRNAseq_margaret/tool/TRUST4/mouse/GRCm38_bcrtcr.fa".
			 " --ref /liulab/galib/mouse_scRNAseq_margaret/tool/TRUST4/mouse/mouse_IMGT+C.fa".
			 " --barcode CB" ) ;
	}
	else
	{
		print STDERR "$BAMId does not exist\n" ;
	}
}
