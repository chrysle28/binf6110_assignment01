#!/bin/bash

## Script file for code used in genome assembly

# Quality control
NanoPlot --fastq SRR32410565.fastq -o nanoplot_results

# Longread assembly
flye -t 8 --nano-hq SRR32410565.fastq -o flye_results
quast assembly.fasta -r ref.fasta -o quast_results

# Alignment
minimap2 -ax asm5 -t 6 ref.fasta assembly.fasta > assembly_ref.sam # assembly to reference
minimap2 -ax lr:hq -t 6 ref.fasta SRR32410565.fastq  > raw_ref.sam # raw reads to reference

# Variant calls
samtools view -bS raw_ref.sam > raw_ref.bam
samtools sort raw_ref.bam -o raw_ref_sorted.bam
samtools faidx ref.fasta
samtools index raw_ref_sorted.bam

run_clair3.sh --bam_fn=/home/ken/binf6110/assign_01/raw_ref_sorted.bam --ref_fn=ref.fasta \
--threads=6 --platform="ont" --model_path="/home/ken/miniconda3/envs/clair3/bin/models/r1041_e82_400bps_sup_v500" \
--include_all_ctgs --output=/home/ken/binf_6110/assign_01/clair3_results

bcftools stats clair3_results/merge_output.vcf > variants.stats

# Obtaining variant stats (e.g. number of variants in gene)
awk -F'\t' 'NR<=10 {print $1 "\t" $4-1 "\t" $5}' genomic.gff > gene_coords.bed #coords were manually edited
bcftools view -R gene_coords.bed merge_output.vcf.gz | grep -v "^#" | wc -l
bcftools view -v snps -R gene_coords.bed merge_output.vcf.gz | grep -v "^#" | wc -l
bcftools view -v indels -R gene_coords.bed merge_output.vcf.gz | grep -v "^#" | wc -l

# File format conversion for visualizations
# Creating karyotype and links files for assembly-to-ref visualization
samtools faidx ref.fasta
awk '{print "chr -",$1,$1,0,$2,"grey"}' ref.fasta.fai > karyotype.ref.txt
samtools faidx assembly.fasta
awk '{print "chr -",$1,$1,0,$2,"blue"}' assembly.fasta.fai > karyotype.asm.txt
cat karyotype.ref.txt karyotype.asm.txt > karyotype.txt

awk '{print $6, $8, $9, $1, $3, $4}' assembly_ref.paf > links.txt

# Creating files for raw-to-ref visualization
# Creating karyotype file
samtools faidx ref.fasta
awk '{print "chr - "$1" "$1" 0 "$2" grey"}' ref.fasta.fai > genome.txt

# Counting variants
awk -F"\t" '$3=="gene"{print $1, $4-1, $5, $9}' OFS="\t" genomic.gff > genes.bed
grep -v '^#' merge_output.vcf | awk '{print $1, $2-1, $2}' OFS="\t" > variants.bed
bedtools intersect -a genes.bed -b variants.bed -c > genes_variant_counts.bed

# Filtering genes based on number of variants (min. 100)
awk '$5 >= 100' genes_variant_counts.bed > genes_high_variants.bed
awk '{print $1, $2, $3, "green"}' OFS="\t" genes_high_variants.bed > genes_blocks.txt
awk '{
    split($4,a,";")
    name="unknown"
    for(i in a){
        if(a[i] ~ /^Name=/){split(a[i],b,"="); name=b[2]}
    }
    print $1, $2, $3, name
}' OFS="\t" genes_high_variants.bed > genes_labels.txt

# Counting SNPs and indels
awk '{print $1" "$2" "$2" 1 "color=red"}' merge_output.vcf > snps_circos.txt
awk '{print $1" "$2" "$2" 1 "color= blue"}' merge_output.vcf > indels_circos.txt

# Variant density heatmap
bedtools makewindows -g ref.fasta.fai -w 10000 > genome_windows.bed
awk '{print $1"\t"$2-1"\t"$2}' merge_output.vcf > variants.bed
bedtools coverage -a genome_windows.bed -b variants.bed > variant_density.txt
awk '{print $1" "$2" "$3" "$5}' variant_density.txt > variant_density_heatmap.txt

# Creating coords file for dotplot
awk '{
  if ($5 == "+")
    print $8, $3, $9, $4;
  else
    print $8, $4, $9, $3;
}' assembly_ref.paf > dotplot.coords
