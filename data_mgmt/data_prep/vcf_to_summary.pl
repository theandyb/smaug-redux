#!/usr/local/bin/perl

################################################################################
# Get initial summary information from input vcfs
#
# To avoid the possibility of modifying raw vcf data, includes option to copy
# raw vcfs to a working project directory
################################################################################
use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use File::Find::Rule;
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin; # data_prep location
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

#my @macs = ("common", "singletons", "doubletons", "all");
my @macs = ("common", "singletons");
#my @macs = ("all");

my $analysisdir = $config->{analysisdir};
my $inputdir = $config->{inputdir};
my $vcftoolsdir = $config->{vcftoolsdir};
my $rawvcfdir = $config->{rawvcfdir};
my $rawvcfext = $config->{rawvcfext};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

# Default summary directory
my $outdir="$analysisdir/summaries";
make_path($outdir);

# Running script with argument 'copy'
# copies VCFs from another directory to avoid overwriting raw data
my $makecopy='';
if(defined($ARGV[0])){
  $makecopy = $ARGV[0];
}

# Specify project folder for VCFs
my $vcfdir="$analysisdir/vcfs";
make_path($vcfdir);

################################################################################
# Reads raw vcfs from original location, annotates with necessary info,
# and copies to input directory. Original vcfs are preserved and not modified.
################################################################################
if ($makecopy eq "copy") {
	for my $mac (@macs){
		my @rawvcfs = File::Find::Rule->file()
	                            ->name("*.$rawvcfext")
	                            ->maxdepth(3)
	                            ->in($rawvcfdir);

		foreach my $rawvcf (@rawvcfs) {
			my $filename=fileparse($rawvcf);

	    if($filename !~ /chrX/){
	      my @parts = split(/\.vcf.gz/, $filename);
	      my $basename = $parts[0];
	      my @nameparts = split(/\./, $basename);
	      my $i = $nameparts[0];
	      $i =~ s/chr//g;

	      my $newvcf = "$analysisdir/vcfs/$basename.ma.aa.$mac.vcf.gz";
	      my $ancestral = "$analysisdir/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz";
	      my $fasta = "$analysisdir/reference_data/human_g1k_v37/chr$i.fasta.gz";

	      # first command extracts singletons with filter PASS, including any that
	      # occur in multiallelic sites
	      my $maparse = "perl $relpath/ma_parse.pl --i $rawvcf";
          my $maparse2 = "perl $relpath/ma_parse_double.pl --i $rawvcf";
          my $maparse3 = "perl $relpath/ma_parse_all.pl --i $rawvcf";

	      # second command fills ancestral allele to AA field
	      my $aaparse = "perl $vcftoolsdir/perl/fill-aa -a $ancestral";

	      # last command adds Motif and Category info fields
	      my $infoparse = "perl $relpath/fill_motif.pl -a $fasta";

	      # pipe commands and execute
	      my $pipe;
	      if($mac eq "singletons"){
	        $pipe = "$maparse | $aaparse | $infoparse | bgzip -c > $newvcf";
	      } elsif($mac eq "common"){
	        my $filter = "bcftools view -i 'AC>=10' -f PASS";
	        $pipe = "$filter $rawvcf | $infoparse | bgzip -c > $newvcf";
          } elsif($mac eq "doubletons"){
            $pipe = "$maparse2 | $aaparse | $infoparse | bgzip -c > $newvcf";
          } elsif($mac eq "all"){
            $pipe = "$maparse3 | $aaparse | $infoparse | bgzip -c > $newvcf";
          }

	      print STDERR "Input file: $rawvcf\n";
	      print STDERR "Writing to: $newvcf...\n";
	      forkExecWait($pipe);
	      print STDERR "Done\n";
	    }
		}

	  print STDERR "Operation complete\n";
	}
}

################################################################################
# Scans input directory for vcfs and outputs summary file to get per-chromosome
# summary files from bcftools. If running without the 'copy' argument, assumes
# vcfs in the input directory have already been parsed for multiallelic sites
################################################################################
my $script = 1;

if ($script==1){
	for my $mac (@macs) {
	  my @vcfs = File::Find::Rule->file()
	                            ->name("*$mac.vcf.gz")
	                            ->maxdepth(1)
	                            ->in($vcfdir);
	  my $header;
	  if($mac eq "common"){
			$header = "\"CHR\tPOS\tREF\tALT\tAN\tMotif\tCategory\tAC\"";
		} elsif ($mac eq "singletons"){
			$header = "\"CHR\tPOS\tREF\tALT\tAA\tAN\tMotif\tCategory\"";
        } elsif ($mac eq "doubletons"){
            $header = "\"CHR\tPOS\tREF\tALT\tAA\tAN\tMotif\tCategory\tAC\"";
        } elsif ($mac eq "all"){
            $header = "\"CHR\tPOS\tREF\tALT\tAA\tAN\tMotif\tCategory\tAC\"";
        }

	  my $summout = "$outdir/$mac.full.summary";
	  my $headercmd = "echo $header > $summout";
	  forkExecWait($headercmd);

		foreach my $file (@vcfs) {
	    print STDERR "Getting summary for $file...\n";
	    unless(-e "$file.tbi"){
	      print STDERR "$file not indexed--indexing now...";
	      my $tabixcmd = "tabix -p vcf $file";
	      forkExecWait($tabixcmd);
	      print STDERR "Done\n";
	    }

	    my $filename=fileparse($file);
	    my @parts = split(/\./, $filename);
	    my $chr = $parts[0];
	    $chr =~ s/chr//;
	    my $bcfquery;
	    my $outputcols;
			if($mac eq "common"){
				# $bcfquery = "bcftools query -i 'AC>=10 && FILTER=\"PASS\"' -r $chr"; # sanity check
	      $bcfquery = "bcftools query -r $chr";
	      $outputcols = "'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AN\t%INFO/Motif\t%INFO/Category\t%AC\n'";
			} elsif ($mac eq "singletons"){
				# $bcfquery = "bcftools query -i 'AC=1 && FILTER=\"PASS\"' -r $chr"; # sanity check
	      $bcfquery = "bcftools query -r $chr";
	      $outputcols = "'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/AN\t%INFO/Motif\t%INFO/Category\n'";
            } elsif ( ($mac eq "doubletons") || ($mac eq "all") ){
                $bcfquery = "bcftools query -r $chr";
                $outputcols = "'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/AN\t%INFO/Motif\t%INFO/Category\t%AC\n'";
            }

	    # my $cmd = "$bcfquery -f $outputcols $file > $outdir/chr$chr.summary";
	    my $cmd = "$bcfquery -f $outputcols $file >> $summout";
			forkExecWait($cmd);
		}
	  print STDERR "Operation complete\n";
	}
}
