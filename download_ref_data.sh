# Function declarations
usage()
{
    echo "usage: download_ref_data [[-s] | [-h]]"
}
status()
{
	if [ "$verbose" = "1" ]; then
		printf "$1"
	fi
}
# Main
verbose=1
wgetOpt=
tarOpt="v"
while [ "$1" != "" ]; do
	case $1 in
		-s | --silent )
			verbose=0
			wgetOpt="-nv"
			tarOpt=""
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
# Create directory structure
declare -a dirs=("reference_data"
	"summaries"
	"output"
	"slurm"
	"vcfs"
	"images"
	"motif_counts"
	"singletons")

## now loop through the above array
for i in "${dirs[@]}"
do
   if [ ! -d "$i" ]; then
		 mkdir "$i"
	 fi
done

cd reference_data

#############################################################################
# hg19 chromosome lengths
#############################################################################
status "Downloading hg19 chromosome lengths..."
echo "chrom\tsize" > "hg19.genome"
curl -s "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes" >> "hg19.genome"
status "Done!\n"

##############################################################################
## 1000G strict mask
##############################################################################
status "Downloading 1000G strict mask..."
curl -s  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed" | bedtools complement -i - -g "hg19.genome" | bedtools sort | awk 'match($1, /chr[0-9]+$/) {print $0}' > "testmask2.bed"
status "Done!\n"

#############################################################################
# Reference genomes in fasta format
#############################################################################
# v37
status "\nReference genome in fasta format..."
mkdir "human_g1k_v37"
curl -s "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz" | gunzip -c > "human_g1k_v37/human_g1k_v37.fasta"
status "Done!\n"

# mask v37
status "\nMasked genome..."
mkdir "human_g1k_v37_mask"
bedtools maskfasta -fi "human_g1k_v37/human_g1k_v37.fasta" -bed "testmask2.bed" -fo "human_g1k_v37_mask/human_g1k_v37.premask.fasta"
status "Done!\n"

# Ancestral genome
status "\nAncestral genome..."
curl -s "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2" > "human_ancestor_GRCh37_e59.tar.bz2"

tar -${tarOpt}jxf "human_ancestor_GRCh37_e59.tar.bz2"
status "Done!\n"

#############################################################################
# CpG islands
#############################################################################
status "\nCpG islands..."
curl -s  "http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt" | awk 'NR>1' | sort -k1,1 -k2,2n > "cpg_islands_sorted.bed"
status "Done!\n"

#############################################################################
# Lamin-associated domains
#############################################################################
status "\nLamin-associated domains..."
curl -s  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz" | gunzip | awk '{print $2"\t"$3"\t"$4}' | bedtools sort -i - > "lamin_B1_LADS2.bed"
status "Done!\n"
#############################################################################
# DNase hypersensitive sites
#############################################################################
status "\nDNase hypersensitive sites..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz" | gunzip | cut -f1-3 | sort -k1,1V -k2,2n - > "DHS.bed"
status "Done!\n"
#############################################################################
# Replication timing
#############################################################################
status "\nReplication timing..."
wget http://mccarrolllab.org/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip; unzip Koren-et-al-Table-S2.zip; mv 'Koren et al Table S2.txt' lymph_rep_time.txt; rm Koren-et-al-Table-S2.zip
status "Done!\n"
#############################################################################
# Recombination rate
#############################################################################
status "\nRecombination rates..."
curl -s "http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw" > "SexAveraged.bw"
bigWigToWig "SexAveraged.bw" "SexAveraged.wig"
echo "CHR\tSTART\tEND\tRATE" > "recomb_rate.bed"
awk 'NR>1' "SexAveraged.wig" | cat | sort -k1,1 -V -s >> "recomb_rate.bed"
sed -i '/^#bed/d' ./recomb_rate.bed
status "Done!\n"

#############################################################################
# Histone marks
#############################################################################
status "\nHistone marks..."
wget $wgetOpt -r -nd -P . --accept-regex 'E062' https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/

gunzip *.broadPeak.gz
for f in *.broadPeak; do
	mv -- "$f" "${f%.broadPeak}.bed"
done

for i in E062*.bed; do
	bedtools sort -i $i > sort.$i
done
status "Done!\n"
#############################################################################
# de novo mutations
#############################################################################
status "\nDe novo mutations..."
mkdir DNMs
# GoNL
curl -s "https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5.2/GoNL_DNMs.txt" > "DNMs/GoNL_DNMs.txt"

# ITMI
curl -s "https://media.nature.com/original/nature-assets/ng/journal/v48/n8/extref/ng.3597-S3.xlsx" > "DNMs/goldmann_2016_dnms.xlsx"

status "Done!\n"
#############################################################################
# RefSeq v69 exons
# originally downloaded via UCSC genome browser;
# this command directly downloads the file used in analyses
#############################################################################
status "\nRefSeq v69 exons..."
curl -s "http://mutation.sph.umich.edu/hg19/GRCh37_RefSeq_sorted.bed" >  "GRCh37_RefSeq_sorted.bed"
status "Done!\n"
#############################################################################
# cytobands
#############################################################################
status "\nCytobands..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip > "cytoBand.txt"

# get mask pct per band; used in downstream analysis
bedtools coverage -a "testmask2.bed" -b "cytoBand.txt" > "testcov.bed"
status "Done!\n"
#############################################################################
# Human Accelerated Regions
#############################################################################
status "\nHuman Accelerated Regions..."
curl -s "ftp://ftp.broadinstitute.org/pub/assemblies/mammals/29mammals/2xHARs.bed" > "2xHARs.bed"

curl -s "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver" > liftOver
chmod +x liftOver

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz" > "hg18ToHg19.over.chain.gz"

"liftOver" "2xHARs.bed" "hg18ToHg19.over.chain.gz" "2xHARs.hg19.bed" "unlifted.bed"

bedtools sort -i "2xHARs.hg19.bed" > "2xHARs.hg19.sort.bed"
status "Done!\n"
#############################################################################
# Aggarwala & Voight rates
#############################################################################
status "\nAggarwala & Voight rates..."
curl -s "https://media.nature.com/original/nature-assets/ng/journal/v48/n4/extref/ng.3511-S2.xlsx" > "AV_rates.xlsx"
status "Done!\n"
cd ../
