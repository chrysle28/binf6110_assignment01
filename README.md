# Genome assembly and comparison of Salmonella enterica
The objective of this assignment is to generate a de novo assembly of the _Salmonella enterica_ genome from long sequence reads, align the assembly to an established reference genome, and conduct and visualize variant calls.

## **1 | Introduction**
Genome assembly methods aim to generate complete chromosome-level sequences from raw sequence reads, which in turn supports key biological initiatives [(Collins, 2018)]([^1]). The presence of complete genomes allows for the linking of SNPs to phenotypic effects, the determination of organism variants, the creation of high-quality reference genomes, and myriad other studies [(Dida & Yi, 2021)]([^2]).  

A genome can be assembled from a reference, this being reference-based assembly, or from scratch, as in _de novo_ assembly. As its name implies, _de novo_ assembly eschews the use of a reference sequence, and thus the quality of the assembly depends on the continuity of the contigs [(Dida & Yi, 2021)]([^2]). Consequently, _de novo_ assembly may struggle with genomes containing short and inaccurate reads, or an abundance of repeats, leading to inaccurate and/or fragmented reconstructions [(Dida & Yi, 2021;]([^2]),[Yang et al., 2025)]([^3]). Notably, reference-based assembly was found to outperform _de novo_ assembly for plasmid reconstruction [(Li et al., 2022)]([^4]). However, these challenges mainly stem from the use of short reads in de novo assembly. 

Short-read sequencing is more cost-effective, but assembly is inherently limited by the length of the reads; as such, it is difficult to construct chromosome-level assemblies using Sanger and second-generation sequencing [(Yang et al., 2025;]([^3]),[Amarasinghe et al., 2020)]([^5]). In contrast, long reads can not only refine _de novo_ construction of genomes, but also improve the detection of variants [(Yang et al., 2025)]([^3]). Long reads based on third-generation sequencing had higher error rates and could not handle complex genomic regions (e.g. tandem repeats) very well [(Yang et al., 2025)]([^3]), but next-generation sequencing methods such as PacBio and Oxford Nanopore Technologies (ONT) have addressed these issues. These methods offer longer read lengths suitable for de novo genome assembly with little to no trade-off in per-read accuracy, as evidenced by PacBio HiFi reads boasting an error rate of less than 0.1%-1% [(Yang et al., 2025;]([^3]),[Amarasinghe et al., 2020)]([^5]). Therefore, genome assembly has pivoted into the use of long read sequences. However, it should be noted that despite the many technological and methodological advancements being made, genome assembly still faces challenges such as high costs, biological constraints clashing with current technologies (e.g. ultra-long tandem repeats and polyploid genomes), and current assemblers being too time-intensive and inefficient [(Yang et al., 2025)]([^6]). 

Highlighting the last point, it is critical that the pipelines and tools used are not only suitable for the type of assembly being performed (i.e. _de novo_ vs reference-based), but also the type of read being used (i.e. short-read vs long-read). Focusing on de novo assemblies using long-reads, assemblers such as Flye and Canu have been the most used and tested. Both use a different approach to assembly - Canu uses the Overlap Layout Consensus Method (OLC) wherein contigs are created from overlaps between the reads [(Koren et al., 2017)]([^7]), while Flye uses a hybrid approach involving repeat graphs (i.e. disjointigs) [(Kolmogorov et al., 2019)]([^8]). Many systematic reviews have found that Flye consistently outperforms Canu for assemblies using ONT reads in metrics such as computational speed, sequence identity, contiguity, misassembly count, and gene identification [(Dida & Yi, 2021;]([^2])[Cosma et al., 2023;]([^9])[Boostrom et al., 2022)]([^10]). Moreover, Canu produced the most fragmented reads out of all long-read assemblers tested, and also sometimes failed to complete with certain inputs (e.g. reads processed with Filtlong) [(Boostrom et al., 2022)]([^10]). Given Flye’s consistency in terms of producing the most accurate assemblies from long reads, it will be the main tool used for _de novo_ assembly of the _Salmonella enterica_ genome.

## **2 | Methods**
### 2.1 | Description of Data
Long read sequences of a _Salmonella enterica_ isolate (accession SRR32410565) were used for a _de novo_ assembly of a genome. Reads were generated via Oxford Nanopore, with R10 chemistry, producing an expected accuracy of Q20+ and N50 of 5-15kb. Raw read data was obtained in FASTQ format from https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR32410565. The reference _Salmonella enterica_ genome, along with its GFF3 annotation file, was obtained through the NCBI genome database (assembly ASM694v2) in FASTA format.

### 2.2 | Inspection and Quality Control
Long read sequences were inspected for quality and length using a combination of Nanoplot (v1.46.2) and Filtlong (v0.3.1). Nanoplot was used to generate sequence statistics, which informed the subsequent processing of reads with Filtlong. The initial high quality of the reads did not necessitate the filtering of any sequences.
```
NanoPlot --fastq <input.fastq>  -o <output_dir>
NanoPlot --fastq SRR32410565.fastq  -o nanoplot_results

filtlong --min_length <#> –-keep_percent <#> input.fastq > output.fastq
filtlong --min_length 1000 -q 15 –-keep_percent 90 input.fastq > q20_filt.fastq
```
> The values for the min_length and keep_percent parameters were determined through inspection of the sequence statistics generated by Nanoplot.

### 2.3 | Longread Assembly
After quality control, the genome was assembled from the long reads using Flye (v2.9.6).
```
flye -t <#> –genome-size <#> –-nano-hq input.fastq -o output_dir
flye -t 6 --nano-hq SRR32410565.fastq -o flye_results
```
> Default parameters were be used. The –nano-hq flag was employed due to the expected accuracy of the raw reads (Q20+), and the reported average accuracy of 18.9.

### 2.4 | Alignment and Verification
The de novo assembly and the raw reads were aligned with the reference genome using Minimap2 (v2.30). Two different presets were used: one for mapping of reads to the reference, and one for complete genome-to-genome alignment.
```
minimap2 -ax asm5 -t 6 GCF_000006945.2_ASM694v2_genomic.fna assembly.fasta > assembly_ref.sam
minimap2 -ax lr:hq -t 6 GCF_000006945.2_ASM694v2_genomic.fna SRR32410565.fastq  > raw_ref.sam
```
> The asm5 parameter was employed for the alignment of the assembly to the reference as it is expected that the assembly and reference will have less than 5% divergence.

The quality of each alignment was checked using QUAST (v5.3.0), and the results were visualized with the bundled Icarus viewer.

### 2.5 | Variant Calling
The SAM files generated from the alignment were converted to BAM format, sorted, and then indexed using samtools (v1.23).
```
samtools view -bS raw_ref.sam > raw_ref.bam
samtools sort raw_ref.bam -o raw_ref_sorted.bam
samtools faidx GCF_000006945.2_ASM694v2_genomic.fna
samtools index raw_ref_sorted.bam
```
Variant calling was performed using Clair3 (v1.2.0), with the input files being the sorted raw reads in BAM format, and the reference genome in FASTA format.
```
run_clair3.sh --bam_fn=/home/ken/binf6110/assign_01/raw_ref_sorted.bam --ref_fn=GCF_000006945.2_ASM694v2_genomic.fna --threads=6 --platform="ont" --model_path="/home/ken/miniconda3/envs/clair3/bin/models/r1041_e82_400bps_sup_v500" --include_all_ctgs --output=/home/ken/binf_6110/assign_01/clair3_results
```
> Default parameters for variant calling of Oxford Nanopore reads were used.

The resulting VCF file was converted to a file containing the variant stats using bcftools (v1.23).
```
bcftools stats merge_output.vcf > variants.stats
```


### 2.7 | Visualization
The _de novo_ assembly was visualized using Bandage (v0.9.0), with the input file being the assembly graphs generated by Flye.

Visualization of the assembly against the reference, as well as the variants, were done via Circos (v0.69-10). Required files (e.g. karyotype, links) were generated from the PAF file made by Minimap2, the GFF3 annotation file of the reference, and the VCF file generated by Clair3. Code used to generate the files and the Circos configuration file can be found in the circos folder.

Variants were visualized using IGV (v2.19.7). The following files were inputted: the reference genome FASTA (ASM694v2) and GFF annotation file, the assembly BAM, and the read alignment BAM.


Counting variants
```
awk -F"\t" '$3=="gene"{print $1, $4-1, $5, $9}' OFS="\t" genomic.gff > genes.bed
grep -v '^#' variants.vcf | awk '{print $1, $2-1, $2}' OFS="\t" > variants.bed
bedtools intersect -a genes.bed -b variants.bed -c > genes_variant_counts.bed
```
Filtering
```
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
```

## References
[^1]: Collins, A. (2018). The challenge of genome sequence assembly. _Open Bioinformatics Journal, 11_(1), 231-239.
[^2]: Dida, F., & Yi, G. (2021). Empirical evaluation of methods for de novo genome assembly. _PeerJ Computer Science, 7_, e636.
[^3]: Yang, Y., Du, W., Li, Y., Lei, J., & Pan, W. (2025). Recent advances and challenges in de novo genome assembly. _Genomics Communications, 2_(1)
[^4]: Li, I. C., Yu, G. Y., Huang, J. F., Chen, Z. W., & Chou, C. H. (2022). Comparison of reference-based assembly and De novo assembly for bacterial plasmid reconstruction and AMR gene localization in Salmonella enterica Serovar Schwarzengrund isolates. _Microorganisms, 10_(2), 227.
[^5]: Amarasinghe, S. L., Su, S., Dong, X., Zappia, L., Ritchie, M. E., & Gouil, Q. (2020). Opportunities and challenges in long-read sequencing data analysis. _Genome biology, 21_(1), 30.
[^6]: Yang, Y., Du, W., Li, Y., Lei, J., & Pan, W. (2025). Recent advances and challenges in de novo genome assembly. _Genomics Communications, 2_(1).
[^7]: Koren, S., Walenz, B. P., Berlin, K., Miller, J. R., Bergman, N. H., & Phillippy, A. M. (2017). Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. _Genome research, 27_(5), 722-736.
[^8]: Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. _Nature biotechnology, 37_(5), 540-546.
[^9]: Cosma, B. M., Shirali Hossein Zade, R., Jordan, E. N., van Lent, P., Peng, C., Pillay, S., & Abeel, T. (2023). Evaluating long-read de novo assembly tools for eukaryotic genomes: insights and considerations. _GigaScience, 12_, giad100.
[^10]: Boostrom, I., Portal, E. A., Spiller, O. B., Walsh, T. R., & Sands, K. (2022). Comparing long-read assemblers to explore the potential of a sustainable low-cost, low-infrastructure approach to sequence antimicrobial resistant bacteria with oxford nanopore sequencing. _Frontiers in Microbiology, 13_, 796465.

