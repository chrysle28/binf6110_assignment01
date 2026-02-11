# Genome assembly and comparison of Salmonella enterica
The objective of this project is to generate a de novo assembly of the _Salmonella enterica_ genome from long sequence reads, align the assembly to an established reference genome, and conduct and visualize variant calls.


## **1 | Introduction**
Genome assembly methods aim to generate complete chromosome-level sequences from raw sequence reads, which in turn supports key biological initiatives [(Collins, 2018)]([^1]). The presence of complete genomes allows for the linking of SNPs to phenotypic effects, the determination of organism variants, the creation of high-quality reference genomes, and myriad other studies [(Dida & Yi, 2021)]([^2]).  

A genome can be assembled from a reference, this being reference-based assembly, or from scratch, as in _de novo_ assembly. As its name implies, _de novo_ assembly eschews the use of a reference sequence, and thus the quality of the assembly depends on the continuity of the contigs [(Dida & Yi, 2021)]([^2]). Consequently, _de novo_ assembly may struggle with genomes containing short and inaccurate reads, or an abundance of repeats, leading to inaccurate and/or fragmented reconstructions [(Dida & Yi, 2021]([^2]); [Yang et al., 2025)]([^3]). Notably, reference-based assembly was found to outperform _de novo_ assembly for plasmid reconstruction [(Li et al., 2022)]([^4]). However, these challenges mainly stem from the use of short reads in de novo assembly. 

Short-read sequencing is more cost-effective, but assembly is inherently limited by the length of the reads; as such, it is difficult to construct chromosome-level assemblies using Sanger and second-generation sequencing [(Yang et al., 2025]([^3]); [Amarasinghe et al., 2020)]([^5]). In contrast, long reads can not only refine _de novo_ construction of genomes, but also improve the detection of variants [(Yang et al., 2025)]([^3]). Long reads based on third-generation sequencing had higher error rates and could not handle complex genomic regions (e.g. tandem repeats) very well [(Yang et al., 2025)]([^3]), but next-generation sequencing methods such as PacBio and Oxford Nanopore Technologies (ONT) have addressed these issues. These methods offer longer read lengths suitable for de novo genome assembly with little to no trade-off in per-read accuracy, as evidenced by PacBio HiFi reads boasting an error rate of less than 0.1%-1% [(Yang et al., 2025]([^3]); [Amarasinghe et al., 2020)]([^5]). Therefore, genome assembly has pivoted into the use of long read sequences. However, it should be noted that despite the many technological and methodological advancements being made, genome assembly still faces challenges such as high costs, biological constraints clashing with current technologies (e.g. ultra-long tandem repeats and polyploid genomes), and current assemblers being too time-intensive and inefficient [(Yang et al., 2025)]([^6]). 

Highlighting the last point, it is critical that the pipelines and tools used are not only suitable for the type of assembly being performed (i.e. _de novo_ vs reference-based), but also the type of read being used (i.e. short-read vs long-read). Focusing on de novo assemblies using long-reads, assemblers such as Flye and Canu have been the most used and tested. Both use a different approach to assembly - Canu uses the Overlap Layout Consensus Method (OLC) wherein contigs are created from overlaps between the reads [(Koren et al., 2017)]([^7]), while Flye uses a hybrid approach involving repeat graphs (i.e. disjointigs) [(Kolmogorov et al., 2019)]([^8]). Many systematic reviews have found that Flye consistently outperforms Canu for assemblies using ONT reads in metrics such as computational speed, sequence identity, contiguity, misassembly count, and gene identification [(Dida & Yi, 2021]([^2]); [Cosma et al., 2023]([^9]); [Boostrom et al., 2022)]([^10]). Moreover, Canu produced the most fragmented reads out of all long-read assemblers tested, and also sometimes failed to complete with certain inputs (e.g. reads processed with Filtlong) [(Boostrom et al., 2022)]([^10]). Given Flye’s consistency in terms of producing the most accurate assemblies from long reads, it will be the main tool used for _de novo_ assembly of the _Salmonella enterica_ genome.


## **2 | Methods**
### 2.1 | Description of Data
Long read sequences of a _Salmonella enterica_ isolate (accession SRR32410565) were used for a _de novo_ assembly of a genome. Reads were generated via Oxford Nanopore, with R10 chemistry, producing an expected accuracy of Q20+ and N50 of 5-15kb. Raw read data was obtained in FASTQ format from https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR32410565. The reference _Salmonella enterica_ genome, along with its GFF3 annotation file, was obtained through the NCBI genome database (assembly ASM694v2) in FASTA format. The reference genome file was renamed to ref.fasta for the rest of the analysis.

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

The quality of the assembly was checked using QUAST (v5.3.0), and the results were visualized with the bundled Icarus viewer. 
```
quast assembly.fasta -r ref.fasta -o quast_results
```

### 2.4 | Alignment 
The de novo assembly and the raw reads were aligned with the reference genome using Minimap2 (v2.30). Two different presets were used: one for mapping of reads to the reference, and one for complete genome-to-genome alignment.
```
minimap2 -ax asm5 -t 6 ref.fasta assembly.fasta > assembly_ref.sam
minimap2 -ax lr:hq -t 6 ref.fasta SRR32410565.fastq  > raw_ref.sam
```
> The asm5 parameter was employed for the alignment of the assembly to the reference as it was expected that the assembly and reference will have less than 5% divergence.

### 2.5 | Variant Calling
The SAM files generated from the alignment were converted to BAM format, sorted, and then indexed using samtools (v1.23). Coverage statistics were calculated using samtools (v1.23).
```
samtools view -bS raw_ref.sam > raw_ref.bam
samtools sort raw_ref.bam -o raw_ref_sorted.bam
samtools faidx ref.fasta
samtools index raw_ref_sorted.bam
samtools coverage raw_ref_sorted.bam
```
Variant calling was performed using Clair3 (v1.2.0), with the input files being the sorted raw reads in BAM format, and the reference genome in FASTA format.
```
run_clair3.sh --bam_fn=/home/ken/binf6110/assign_01/raw_ref_sorted.bam --ref_fn=ref.fasta --threads=6 --platform="ont" --model_path="/home/ken/miniconda3/envs/clair3/bin/models/r1041_e82_400bps_sup_v500" --include_all_ctgs --output=/home/ken/binf_6110/assign_01/clair3_results
```
> Default parameters for variant calling of Oxford Nanopore reads were used.

The resulting VCF file was converted to a file containing the variant stats using bcftools (v1.23).
```
bcftools stats merge_output.vcf > variants.stats
```

### 2.7 | Visualization
The _de novo_ assembly was visualized using Bandage (v0.9.0), with the input file being the assembly graphs generated by Flye. An alternative visualization of the assembly against the reference was done using gnuplot (v5.0). Visualization of the assembly against the reference, as well as the variants, were done via Circos (v0.69-10). Required files (e.g. karyotype, links) were generated from the PAF file made by Minimap2, the GFF3 annotation file of the reference, and the VCF file generated by Clair3. Code used to generate the files and the Circos configuration file can be found in the code folder. Variants were visualized using IGV (v2.19.7). The following files were inputted: the reference genome FASTA (ASM694v2) and GFF annotation file, the assembly BAM, and the read alignment BAM.


## **3 | Results**
**Quality of _de novo_ Genome Assembly and Alignment to Reference**

Assembly of the long read sequences with Flye resulted in three contigs, with lengths of 3,318,776 bp, 1,676,977 bp, and 109,059 bp (Fig. 1, Fig. 2A, Fig. 3, Fig. 4). The assembly's total length was 5,104,812 bp, compared to the reference genome's length of 4,951,383 bp. The largest continuous alignment in the assembly was 953,687 bp. The bandage plot reports an additional sequence 6,269 bp long, which represents repeated regions which were not resolved. QUAST analysis reported that the total number of aligned bases in the assembly was 4,746,095, with a duplication ratio of 1.002, while the genome fraction was 95.669%. Moreover, QUAST showed that there were 35 total misassemblies (25 relocations and 10 local), with an overall misassembled contigs length  of 4,995,753 bp. (Fig. 2B). There were 27.39 mismatches per 100 kbp, and 3.81 indels per 100 kbp.

<img width="926" height="578" alt="bandage_plot" src="https://github.com/user-attachments/assets/bbb104ba-8bbd-4036-8b74-73c1ef718a49" />

Fig. 1: Bandage plot of the _de novo_ genome assembly generated via Flye. Chromosome contigs are labeled as edge_1 and edge_2 (green and brown loops, respectively) while the plasmid is labeled as edge_4. Edge_3 represents repeat regions. 

<img width="1350" height="455" alt="quast_figs" src="https://github.com/user-attachments/assets/efa92638-a82f-4d7d-a721-d89f34c663ff" />

Fig. 2: Icarus output of the assembly (top track) aligned against the reference (bottom track). **A)** Size comparison of the assembled genome against the reference, with the chromosome contigs highlighted in red, and the plasmid contig in grey **B)** Alignment of assembly contigs against the reference, with misassemblies marked by grey spheres. The top track displays overlaps between the assembled contigs.

The dot plot of the assembly (Fig. 3) displays 3 diagonal and colinear segments, which each correspond to one contig aligning to the reference genome. None of the segments slope downwards, indicating that there were no large inversions and that the contigs are in the same orientation. The small fragments at the bottom left likely represent the repeat regions. The Circos plot (Fig. 4) displays all 3 contigs and their links to respective regions in the reference. Contigs 1 and 2  consistently map to NC_003197.2 (the reference chromosome) while contig 4 maps fairly well to NC_003277.2 (the reference plasmid).

<img width="2000" height="2000" alt="dotplot" src="https://github.com/user-attachments/assets/40166f74-44d0-46a9-b875-915bcc30ede7" />

Fig. 3: Dotplot of the alignment between the assembly and the reference genome, with contigs plotted in blue. The two longer contigs represent the chromosome, while the shortest contig represents the plasmid.

<img width="3000" height="3000" alt="circos_mm" src="https://github.com/user-attachments/assets/83b7c169-fb4d-4767-b518-35461a9dbc5c" />

Fig. 4: Circos plot visualizing links between the contigs of the assembly (blue) and the reference genome (grey).


**Assembly-to-Reference and Reads-to-Reference Alignments**

| Name  | Coverage (%) | Mean Depth | Mean base quality | Mean mapping quality |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| **Assembly to Reference** |
| NC_003197.2  | 97.4  | 0.98  | 255  | 60  |
| NC_003277.2  | 0  | 0  | 0  | 0  |
| **Reads to Reference** |
| NC_003197.2  | 97.8  | 150.89  | 41.4  | 59.5  |
| NC_003277.2  | 43.1  | 81.66  | 42.2  | 44.7  |

Table 1: Statistics of the assembly-to-reference alignment and reads-to-reference alignment generated by samtools.

In the assembly to reference alignment, the chromosome (NC_003197.2) had a coverage of 97.42%, a mean depth of 0.98, a mean base quality of 255, and a mean mapping quality of 60. On the other hand, the plasmid (NC_003277.2) had 0% coverage, and consequently, 0 in all other metrics. In the raw reads to reference alignment, the chromosme had a coverage of 97.8% and a mean depth of 150.89, while the plasmid had a coverage of 43.1% and a mean depth of 81.66. While the two differed greatly in coverage and mean depth, they had fairly close mean base qualities and mean mapping qualities (41.4, 59.5 and 42.2, 44.7, respectively).

**Visualization of Variants**

After variant calling, Clair3 reported 10,076 variants, with 8991 being single nucleotide polymorphisms (SNPs) and 1114 being indels. SNPs appeared fairly distributed across both the chromosome and the plasmid, while indels appear to cluster in specfic regions within the genome (Fig. 5). Filtering genes with more than 100 variants reveals that a majority of these genes are located within the plasmid. Only 3 genes within the chromsome (STM2628, STM1009, and STM1022) contained more than 100 variants. The greatest density of indels appears to be located within the plasmid, while SNPs occur more freqeuntly in the chromosomal region.

<img width="3000" height="3000" alt="circos_variants" src="https://github.com/user-attachments/assets/9c4486b5-7507-4e0a-a8f4-8c9f87f7a2c0" />

Fig. 5: Circos plot of variants (i.e. SNPs, insertions, and deletions) after aligning raw reads to the reference genome. Red circles represent the SNPs, while blue triangles represent indels. Genes which contain more than 100 variants in the reads are highlighted in green. The variant density heatmap is binned into windows of 10 kb, with red and blue regions representing greater variant density of SNPs and indels, respectively.

Visualizing one of the high-variant chromosomal genes (STM1009) in IGV reveals specific details about the variants within the gene (Fig. 6). STM1009 has a total of 211 variants (148 SNPs, 63 indels). Some variants lead to missense mutations and potential frameshift mutations, but the majority of SNPs appear lead to silent mutations. In contrast, another chromosomal gene (nrdD) only has 1 variant, a C to T substitution. (Fig. 7). A high-variant gene in the plasmid (repA2) has 185 variants (180 SNPs, 5 indels), with many leading to potential missense, frameshift, or nonsense mutations (Fig. 8).

<img width="1000" height="494" alt="STM1009" src="https://github.com/user-attachments/assets/09b94b56-262a-42ab-a474-6794a1f6f32b" />

Fig. 6: IGV view of STM1009, a gene in the chromosome, showing variants in the raw reads compared to the reference genome. Examples of deletions (red box), SNPs (orange box), and insertions (blue box) are outlined. Missense mutations are highlighted in green in the sequence track, while stop codons are highlighted in red. The deletion followed by a stop codon may represent a frameshift mutation.

<img width="1903" height="912" alt="nrdD" src="https://github.com/user-attachments/assets/d4e21c54-0544-4a00-b811-4550ed1161fe" />

Fig. 7: IGV view of nrdD, a gene in the chromosome, showing one variant (an SNP).

<img width="1000" height="467" alt="repA2" src="https://github.com/user-attachments/assets/d49d3cc5-c97a-4308-a042-4a5105e14311" />

Fig. 8: IGV view of repA2, a gene in the pSLT plasmid, showing variants compared to the reference genome. Indels are represented by the purple lines, while SNPs are represented by orange, green, blue, and red lines. Altered proteins and stop codons are highlighted in the sequence track in green and red, respectively.


## **4 | Discussion**
**Comparison of Assembly to Reference**

**Signifcance of Variants**


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

