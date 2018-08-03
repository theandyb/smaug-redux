if [ ! -d "reference_data" ]; then
  mkdir reference_data
fi

cd reference_data

#############################################################################
# hg19 chromosome lengths
#############################################################################
curl -s "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes" > "hg19.genome"

##############################################################################
## 1000G strict mask
##############################################################################
curl -s  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed" | bedtools complement -i - -g "hg19.genome" | bedtools sort | awk 'match($1, /chr[0-9]+$/) {print $0}' > "testmask2.bed"

#############################################################################
# Reference genomes in fasta format
#############################################################################
# v37
mkdir "human_g1k_v37"
curl -s "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz" | gunzip -c > "human_g1k_v37/human_g1k_v37.fasta"

# mask v37
mkdir "human_g1k_v37_mask"
bedtools maskfasta -fi "human_g1k_v37/human_g1k_v37.fasta" -bed "testmask2.bed" -fo "human_g1k_v37_mask/human_g1k_v37.premask.fasta"

# Ancestral genome
curl -s "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2" > "human_ancestor_GRCh37_e59.tar.bz2"

tar -vjxf "human_ancestor_GRCh37_e59.tar.bz2"
#############################################################################
# GC content in 10kb windows
#############################################################################
bedtools nuc -fi "human_g1k_v37/human_g1k_v37.fasta" -bed <(sed s/chr// "genome.10kb.sorted.bed") > "gc10kb.bed"

#############################################################################
# CpG islands
#############################################################################
curl -s  "http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt" | awk 'NR>1' | sort -k1,1 -k2,2n > "cpg_islands_sorted.bed"

#############################################################################
# Lamin-associated domains
#############################################################################
curl -s  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz" | gunzip | awk 'NR>1 {print $2"\t"$3"\t"$4}' | bedtools sort -i - > "lamin_B1_LADS2.bed"

#############################################################################
# DNase hypersensitive sites
#############################################################################
curl -s "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz" | gunzip | cut -f1-3 | bedtools sort -i - > "DHS.bed"

#############################################################################
# Replication timing
#############################################################################
curl -s "http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip" | gunzip > "lymph_rep_time.txt"

#############################################################################
# Recombination rate
#############################################################################
curl -s "http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw" > "SexAveraged.bw"
bigWigToWig "SexAveraged.bw" "SexAveraged.wig"
echo "CHR\tSTART\tEND\tRATE" > "recomb_rate.bed"
awk 'NR>1' "SexAveraged.wig" | cat >> "recomb_rate.bed"

#############################################################################
# Histone marks
#############################################################################
wget -r -nd -P . --accept-regex 'E062' https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/

gunzip *.broadPeak.gz
for f in *.broadPeak; do
	mv -- "$f" "${f%.broadPeak}.bed"
done

for i in E062*.bed; do
	bedtools sort -i $i > sort.$i
done

#############################################################################
# de novo mutations
#############################################################################
mkdir DNMs
# GoNL
curl -s "https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5.2/GoNL_DNMs.txt" > "DNMs/GoNL_DNMs.txt"

# ITMI
curl -s "http://www.nature.com/ng/journal/v48/n8/extref/ng.3597-S3.xlsx" > "DNMs/goldmann_2016_dnms.xlsx"

#############################################################################
# RefSeq v69 exons
# originally downloaded via UCSC genome browser;
# this command directly downloads the file used in analyses
#############################################################################
curl -s "http://mutation.sph.umich.edu/hg19/GRCh37_RefSeq_sorted.bed" >  "GRCh37_RefSeq_sorted.bed"

#############################################################################
# cytobands
#############################################################################
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip > "cytoBand.txt"

# get mask pct per band; used in downstream analysis
bedtools coverage -a "testmask2.bed" -b "cytoBand.txt" > "testcov.bed"

#############################################################################
# Human Accelerated Regions
#############################################################################
curl -s "ftp://ftp.broadinstitute.org/pub/assemblies/mammals/29mammals/2xHARs.bed" > "2xHARs.bed"

curl -s "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver" > liftOver
chmod +x liftOver

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz" > "hg18ToHg19.over.chain.gz"

"liftOver" "2xHARs.bed" "hg18ToHg19.over.chain.gz" "2xHARs.hg19.bed" "unlifted.bed"

bedtools sort -i "2xHARs.hg19.bed" > "2xHARs.hg19.sort.bed"

#############################################################################
# Aggarwala & Voight rates
#############################################################################
curl -s "https://media.nature.com/original/nature-assets/ng/journal/v48/n4/extref/ng.3511-S2.xlsx" > "AV_rates.xlsx"

cd ../
