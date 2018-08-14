#!/bin/bash

#############################################################################
# Script does computationally heavier formatting of the reference data
#############################################################################
verbose=1
usage()
{
    echo "usage: prep_ref_data [[-s] | [-h]]"
}
status()
{
	if [ "$verbose" = "1" ]; then
		printf "$1"
	fi
}
# Main
while [ "$1" != "" ]; do
	case $1 in
		-s | --silent )
			verbose=0
			;;
		-h | --help )
			usage
			exit
			;;
		* )
			usage
			exit 1
	esac
	shift
done
curdir=${PWD}
cd reference_data

#############################################################################
# hg19 chromosome lengths
#############################################################################

# Make fixed-width windows
status "Building fixed-width windows...\n"
status "1000kb..."
bedtools makewindows -g "hg19.genome" -w 1000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "genome.1000kb.sorted.bed"
status "Done!\n5000kb..."
bedtools makewindows -g "hg19.genome" -w 5000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "genome.5000kb.sorted.bed"
status "Done!\n100kb..."
bedtools makewindows -g "hg19.genome" -w 100000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "genome.100kb.sorted.bed"
status "Done!\n10kb..."
bedtools makewindows -g "hg19.genome" -w 10000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "genome.10kb.sorted.bed"
status "Done!\nFull genome..."
bedtools makewindows -g "hg19.genome" -w 3000000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "genome.full.sorted.bed"
status "Done!\n"

#############################################################################
# GC content in 10kb windows
#############################################################################
status "\nGC content in 10kb windows..."
sed s/chr// "genome.10kb.sorted.bed" | bedtools nuc -fi "human_g1k_v37/human_g1k_v37.fasta" -bed - > "gc10kb.bed"
status "Done!\n"

#############################################################################
# Compress and index reference genomes
#############################################################################
# v 37
status "Compress and index reference genome...\n"
for i in `seq 1 22`; do
	status "Chr ${i}..."
	samtools faidx "human_g1k_v37/human_g1k_v37.fasta" $i | bgzip -c > "human_g1k_v37/chr$i.fasta.gz"
	samtools faidx "human_g1k_v37/chr$i.fasta.gz"
	status "Done!\n"
done

# mask v 37
status "Compress and index masked genome...\n"
perl -ane 'if(/\>/){$a++;print ">$a dna:chromosome\n"}else{print;}' "human_g1k_v37_mask/human_g1k_v37.premask.fasta" > "human_g1k_v37_mask/human_g1k_v37.mask.fasta"

rm -f "human_g1k_v37_mask/human_g1k_v37.premask.fasta"

for i in `seq 1 22`; do
	status "Chr ${i}..."
	samtools faidx "human_g1k_v37_mask/human_g1k_v37.mask.fasta" $i | bgzip -c > "human_g1k_v37_mask/chr$i.fasta.gz"
	samtools faidx "human_g1k_v37_mask/chr$i.fasta.gz"
	status "Done!\n"
done

# ancestral genome
status "Compress and index ancestral genome...\n"
for i in `seq 1 22`; do
	status "human_ancestor_GRCh37_e59/human_ancestor_$i.fa..."
	cat "human_ancestor_GRCh37_e59/human_ancestor_$i.fa" | sed "s,^>.*,>$i," | \
		bgzip -c > "human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz"
	samtools faidx "human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz"
	status "Done!\n"
done
cd $curdir
