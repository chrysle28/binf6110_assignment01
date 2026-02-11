# Genome assembly and comparison of Salmonella enterica
The objective of this project is to generate a de novo assembly of the _Salmonella enterica_ genome from long sequence reads, align the assembly to an established reference genome, and conduct and visualize variant calls.


## **1 | Introduction**
Genome assembly methods aim to generate complete chromosome-level sequences from raw sequence reads, which in turn supports key biological initiatives (Collins, 2018). The presence of complete genomes allows for the linking of SNPs to phenotypic effects, the determination of organism variants, the creation of high-quality reference genomes, and myriad other studies (Dida & Yi, 2021).  

A genome can be assembled from a reference, this being reference-based assembly, or from scratch, as in _de novo_ assembly. As its name implies, _de novo_ assembly eschews the use of a reference sequence, and thus the quality of the assembly depends on the continuity of the contigs (Dida & Yi, 2021). Consequently, _de novo_ assembly may struggle with genomes containing short and inaccurate reads, or an abundance of repeats, leading to inaccurate and/or fragmented reconstructions (Dida & Yi, 2021; Yang et al., 2025). Notably, reference-based assembly was found to outperform _de novo_ assembly for plasmid reconstruction (Li et al., 2022). However, these challenges mainly stem from the use of short reads in de novo assembly. 

Short-read sequencing is more cost-effective, but assembly is inherently limited by the length of the reads; as such, it is difficult to construct chromosome-level assemblies using Sanger and second-generation sequencing (Yang et al., 2025; Amarasinghe et al., 2020). In contrast, long reads can not only refine _de novo_ construction of genomes, but also improve the detection of variants (Yang et al., 2025). Long reads based on third-generation sequencing had higher error rates and could not handle complex genomic regions (e.g. tandem repeats) very well (Yang et al., 2025), but next-generation sequencing methods such as PacBio and Oxford Nanopore Technologies (ONT) have addressed these issues. These methods offer longer read lengths suitable for de novo genome assembly with little to no trade-off in per-read accuracy, as evidenced by PacBio HiFi reads boasting an error rate of less than 0.1%-1% (Yang et al., 2025; Amarasinghe et al., 2020). Therefore, genome assembly has pivoted into the use of long read sequences. However, it should be noted that despite the many technological and methodological advancements being made, genome assembly still faces challenges such as high costs, biological constraints clashing with current technologies (e.g. ultra-long tandem repeats and polyploid genomes), and current assemblers being too time-intensive and inefficient (Yang et al., 2025). 

Highlighting the last point, it is critical that the pipelines and tools used are not only suitable for the type of assembly being performed (i.e. _de novo_ vs reference-based), but also the type of read being used (i.e. short-read vs long-read). Focusing on de novo assemblies using long-reads, assemblers such as Flye and Canu have been the most used and tested. Both use a different approach to assembly - Canu uses the Overlap Layout Consensus Method (OLC) wherein contigs are created from overlaps between the reads (Koren et al., 2017), while Flye uses a hybrid approach involving repeat graphs (i.e. disjointigs) (Kolmogorov et al., 2019). Many systematic reviews have found that Flye consistently outperforms Canu for assemblies using ONT reads in metrics such as computational speed, sequence identity, contiguity, misassembly count, and gene identification (Dida & Yi, 2021; Cosma et al., 2023; Boostrom et al., 2022). Moreover, Canu produced the most fragmented reads out of all long-read assemblers tested, and also sometimes failed to complete with certain inputs (e.g. reads processed with Filtlong) (Boostrom et al., 2022). Given Flye’s consistency in terms of producing the most accurate assemblies from long reads, it will be the main tool used for _de novo_ assembly of the _Salmonella enterica_ genome.


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

Assembly of the long read sequences with Flye resulted in three contigs, with lengths of 3,318,776 bp, 1,676,977 bp, and 109,059 bp (Fig. 1, Fig. 2A, Fig. 3, Fig. 4). The assembly's total length was 5,104,812 bp, compared to the reference genome's length of 4,951,383 bp. The largest continuous alignment in the assembly was 953,687 bp. The bandage plot reports an additional sequence 6,269 bp long, which represents repeated regions which were not resolved. QUAST analysis reported that the total number of aligned bases in the assembly was 4,746,095, with a duplication ratio of 1.002, while the genome fraction was 95.669%. Moreover, QUAST showed that there were 35 total misassemblies (25 relocations and 10 local), with an overall misassembled contigs length of 4,995,753 bp. (Fig. 2B). There were 27.39 mismatches per 100 kbp, and 3.81 indels per 100 kbp. Contig 4, representing the assembled plasmid, is unaligned with the reference plasmid (Fig 2).

<img width="926" height="578" alt="bandage_plot" src="https://github.com/user-attachments/assets/bbb104ba-8bbd-4036-8b74-73c1ef718a49" />

_Fig. 1_: Bandage plot of the _de novo_ genome assembly generated via Flye. Chromosome contigs are labeled as edge_1 and edge_2 (green and brown loops, respectively) while the plasmid is labeled as edge_4. Edge_3 represents repeat regions. 

<img width="1350" height="455" alt="quast_figs" src="https://github.com/user-attachments/assets/efa92638-a82f-4d7d-a721-d89f34c663ff" />

_Fig. 2_: Icarus output of the assembly (top track) aligned against the reference (bottom track). **A)** Size comparison of the assembled genome against the reference, with the chromosome contigs highlighted in red, and the plasmid contig in grey **B)** Alignment of assembly contigs against the reference, with misassemblies marked by grey spheres. The top track displays overlaps between the assembled contigs.

The dot plot of the assembly (Fig. 3) displays 3 diagonal and colinear segments, which each correspond to one contig aligning to the reference genome. None of the segments slope downwards, indicating that there were no large inversions and that the contigs are in the same orientation. The small fragments at the bottom left likely represent the repeat regions. The Circos plot (Fig. 4) displays all 3 contigs and their links to respective regions in the reference. Contigs 1 and 2 consistently map to NC_003197.2 (the reference chromosome) while contig 4 shows only sparse links to NC_003277.2 (the reference plasmid).

<img width="2000" height="2000" alt="dotplot" src="https://github.com/user-attachments/assets/40166f74-44d0-46a9-b875-915bcc30ede7" />

_Fig. 3_: Dotplot of the alignment between the assembly and the reference genome, with contigs plotted in blue. The two longer contigs represent the chromosome, while the shortest contig represents the plasmid.

<img width="3000" height="3000" alt="circos_mm" src="https://github.com/user-attachments/assets/83b7c169-fb4d-4767-b518-35461a9dbc5c" />

_Fig. 4_: Circos plot visualizing links between the contigs of the assembly (blue) and the reference genome (grey).


**Assembly-to-Reference and Reads-to-Reference Alignments**

| Name  | Coverage (%) | Mean Depth | Mean base quality | Mean mapping quality |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| **Assembly to Reference** |
| NC_003197.2  | 97.4  | 0.98  | 255  | 60  |
| NC_003277.2  | 0  | 0  | 0  | 0  |
| **Reads to Reference** |
| NC_003197.2  | 97.8  | 150.89  | 41.4  | 59.5  |
| NC_003277.2  | 43.1  | 81.66  | 42.2  | 44.7  |

_Table 1_: Statistics of the assembly-to-reference alignment and reads-to-reference alignment generated by samtools.

In the assembly to reference alignment, the chromosome (NC_003197.2) had a coverage of 97.42%, a mean depth of 0.98, a mean base quality of 255, and a mean mapping quality of 60. On the other hand, the plasmid (NC_003277.2) had 0% coverage, and consequently, 0 in all other metrics. In the raw reads to reference alignment, the chromosme had a coverage of 97.8% and a mean depth of 150.89, while the plasmid had a coverage of 43.1% and a mean depth of 81.66. While the two differed greatly in coverage and mean depth, they had fairly close mean base qualities and mean mapping qualities (41.4, 59.5 and 42.2, 44.7, respectively).

**Visualization of Variants**

After variant calling, Clair3 reported 10,076 variants, with 8991 being single nucleotide polymorphisms (SNPs) and 1114 being indels. 2917 variants (2454 SNPs, 463 indels) are within the chromosome, while the remaining 7159 variants (6579 SNPs, 580 indels) are located within the plasmid.  SNPs appear fairly distributed across both the chromosome and the plasmid, while indels appear to cluster in specfic regions within the genome (Fig. 5). Filtering genes with more than 100 variants reveals that a majority of these genes are located within the plasmid. Only 3 genes within the chromsome (STM2628, STM1009, and STM1022) contained more than 100 variants. The greatest density of indels appears to be located within the plasmid, while SNPs occur more freqeuntly in the chromosomal region.

<img width="3000" height="3000" alt="circos_variants" src="https://github.com/user-attachments/assets/9c4486b5-7507-4e0a-a8f4-8c9f87f7a2c0" />

_Fig. 5_: Circos plot of variants (i.e. SNPs, insertions, and deletions) after aligning raw reads to the reference genome. Red circles represent the SNPs, while blue triangles represent indels. Genes which contain more than 100 variants in the reads are highlighted in green. The variant density heatmap is binned into windows of 10 kb, with red and blue regions representing greater variant density of SNPs and indels, respectively.

Visualizing one of the high-variant chromosomal genes (STM1009) in IGV reveals specific details about the variants within the gene (Fig. 6). STM1009 has a total of 211 variants (148 SNPs, 63 indels). Some variants lead to missense mutations and potential frameshift mutations, but the majority of SNPs appear lead to silent mutations. In contrast, another chromosomal gene (nrdD) only has 1 variant, a C to T substitution. (Fig. 7). A high-variant gene in the plasmid (repA2) has 185 variants (180 SNPs, 5 indels), with many leading to potential missense, frameshift, or nonsense mutations (Fig. 8).

<img width="1000" height="494" alt="STM1009" src="https://github.com/user-attachments/assets/09b94b56-262a-42ab-a474-6794a1f6f32b" />

_Fig. 6_: IGV view of STM1009, a gene in the chromosome, showing variants in the raw reads compared to the reference genome. Examples of deletions (red box), SNPs (orange box), and insertions (blue box) are outlined. Missense mutations are highlighted in green in the sequence track, while stop codons are highlighted in red. The deletion followed by a stop codon may represent a frameshift mutation.

<img width="1903" height="912" alt="nrdD" src="https://github.com/user-attachments/assets/d4e21c54-0544-4a00-b811-4550ed1161fe" />

_Fig. 7_: IGV view of nrdD, a gene in the chromosome, showing one variant (an SNP).

<img width="1000" height="467" alt="repA2" src="https://github.com/user-attachments/assets/d49d3cc5-c97a-4308-a042-4a5105e14311" />

_Fig. 8_: IGV view of repA2, a gene in the pSLT plasmid, showing variants compared to the reference genome. Indels are represented by the purple lines, while SNPs are represented by orange, green, blue, and red lines. Altered proteins and stop codons are highlighted in the sequence track in green and red, respectively.


## **4 | Discussion**

**Assembly vs Reference**

The _de novo_ genome assembly of the _Salmonella enterica_ isolate managed to be fairly successful, even without additional filtering or polishing after assembly. The relatively high genome fraction (95.669) to the reference reflects the fact that the assembly was mostly contiguous. A full resolution of the genome was not achieved due to the repeated regions which were filtered out. The longest two contigs (~3.3 Mb and ~1.7 Mb) map well to the reference chromosome with a coverage of 97.4%, and have relatively high mean base quality and mean mapping quality (as compared to the reads to reference alignment). Although QUAST reported 35 misassemblies, visualization of the assembly to reference alignment revealed that there were no translocations nor inversions, and the two contigs closely align to the chromosome.

In contrast, samtools statistics for the plasmid report 0% and 43.1% coverage (for the assembly and raw reads, respectively), which coincides with its lack of alignment to the reference as reported by QUAST. While the chromosome and the plasmid had fairly similar base and mapping qualities, the difference in their coverage may indicate that the assembled plasmid diverges from the reference. The reference plasmid, pSLT, is a virulence plasmid 94Kb long, and variants have been reported from fusion/recombination events with other plasmids (Hiley, Graham & Jennison, 2019). The assembled plasmid has a longer sequence (~109Kb), which may reflect the fact that it is a variant of pSLT that arose from gene duplication. Plasmids do not only contribute to their hosts' virulence, but also play a role in spreading antibiotic resistance and other genes (Rychlik, Gregorova & Hradecka, 2006; Kubasova et al., 2014). The mechanism for doing so has been found to depend on allelic variations within the plasmid genes (Yue & Schifferli, 2014). This aligns well with the fact that variant density is highest in the assembled plasmid (~71% of variants), with these variants potentially leading to changes in its virulence, transmission, or resistance. One of the high-variant genes in the plasmid, repA2, is involved in plasmid replication (Payne et al., 2019), lending more support to the idea that the numerous variants within the assembled plasmid has made it sufficiently different from the pSLT plasmid. These variants may then affect the overall virulence of the _Salmonella enterica_ isolate.

**Significance of Variants**

On a broader scale, the assembled genome contained 10,076 variants (8991 SNPs, 1114 indels). The general ratio between the SNPs and indels, with SNPs making up a greater proportion of variants, aligns well with literature (Chen et al., 2009), though it is lower than the average ratio (19.61 for bacterial genomes). Although the plasmid contained more high-variant genes (i.e. more than 100 variants), the high-variant genes located within the chromosome (STM1009, STM1022, and STM2626) has interesting implications for the _Salmonella enterica_ strain. These three chromosomal genes encode Gifsy-1/Gifsy-2, which are phages that contribute to the virulence of the host bacteria (Lawley et al., 2006; Mohammed & Cormican, 2015). In particular, Gifsy-1, which is encoded by STM2626, is involved in the control of bacterial proliferation in tissues (Lawley et al., 2006). Since these genes are involved in the virulence of _Salmonella enterica_, it is reasonable that they contain the most variants. Potential missense or frameshift mutations, as seen with STM1009 (Fig. 6), could lead to changes in the host bacteria, whether it be an increase or decrease in its virulence. Mutations in other genes, such as the deletion of the gidA gene (Shippy et al. 2012) or SNPs in the envZ gene (Ko & Choi, 2021) have been shown to affect _Salmonella enterica's_ virulence or morphology. In this assembly, no variants were detected in any of these genes, signifying that this specific isolate could potentially reflect a more virulent strain, due to its high variant density in its virulent plasmid and its chromosomal virulence genes. However, it is key to note that even if a gene contains many variants, mutations that lead to an altered or truncated protein may not follow. As reported by Ko & Chou (2021), although an SNP in the nrdD gene leads to an amino acid change, the function of the protein remained unaffected. In this _Salmonella enterica_ strain, an SNP was observed in the nrdD gene (Fig. 7), but it did not lead to any changes in the amino acid sequence.


## References

Amarasinghe, S. L., Su, S., Dong, X., Zappia, L., Ritchie, M. E., & Gouil, Q. (2020). Opportunities and challenges in long-read sequencing data analysis. _Genome biology, 21_(1), 30.

Boostrom, I., Portal, E. A., Spiller, O. B., Walsh, T. R., & Sands, K. (2022). Comparing long-read assemblers to explore the potential of a sustainable low-cost, low-infrastructure approach to sequence antimicrobial resistant bacteria with oxford nanopore sequencing. _Frontiers in Microbiology, 13_, 796465.

Chen, J. Q., Wu, Y., Yang, H., Bergelson, J., Kreitman, M., & Tian, D. (2009). Variation in the ratio of nucleotide substitution and indel rates across genomes in mammals and bacteria. _Molecular biology and evolution, 26(7)_, 1523-1531.

Collins, A. (2018). The challenge of genome sequence assembly. _Open Bioinformatics Journal, 11_(1), 231-239.

Cosma, B. M., Shirali Hossein Zade, R., Jordan, E. N., van Lent, P., Peng, C., Pillay, S., & Abeel, T. (2023). Evaluating long-read de novo assembly tools for eukaryotic genomes: insights and considerations. _GigaScience, 12_, giad100.

Dida, F., & Yi, G. (2021). Empirical evaluation of methods for de novo genome assembly. _PeerJ Computer Science, 7_, e636.

Hiley, L., Graham, R. M., & Jennison, A. V. (2019). Genetic characterisation of variants of the virulence plasmid, pSLT, in Salmonella enterica serovar Typhimurium provides evidence of a variety of evolutionary directions consistent with vertical rather than horizontal transmission. _Plos one, 14(4)_, e0215207.

Lawley, T. D., Chan, K., Thompson, L. J., Kim, C. C., Govoni, G. R., & Monack, D. M. (2006). Genome-wide screen for Salmonella genes required for long-term systemic infection of the mouse. _PLoS pathogens, 2(2)_, e11.

Li, I. C., Yu, G. Y., Huang, J. F., Chen, Z. W., & Chou, C. H. (2022). Comparison of reference-based assembly and De novo assembly for bacterial plasmid reconstruction and AMR gene localization in Salmonella enterica Serovar Schwarzengrund isolates. _Microorganisms, 10_(2), 227.

Ko, D., & Choi, S. H. (2021). Comparative genomics reveals an SNP potentially leading to phenotypic diversity of Salmonella enterica serovar Enteritidis. _Microbial Genomics, 7(5)_, 000572.

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. _Nature biotechnology, 37_(5), 540-546.

Koren, S., Walenz, B. P., Berlin, K., Miller, J. R., Bergman, N. H., & Phillippy, A. M. (2017). Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. _Genome research, 27_(5), 722-736.

Kubasova, T., Matiasovicova, J., Rychlik, I., & Juricova, H. (2014). Complete sequence of multidrug resistance p9134 plasmid and its variants including natural recombinant with the virulence plasmid of Salmonella serovar Typhimurium. _Plasmid, 76_, 8-14.

Payne, M., Octavia, S., Luu, L. D. W., Sotomayor-Castillo, C., Wang, Q., Tay, A. C. Y., ... & Lan, R. (2019). Enhancing genomics-based outbreak detection of endemic Salmonella enterica serovar Typhimurium using dynamic thresholds._ Microbial Genomics, 7(6)_, 000310.

Mohammed, M., & Cormican, M. (2015). Whole genome sequencing provides possible explanations for the difference in phage susceptibility among two Salmonella Typhimurium phage types (DT8 and DT30) associated with a single foodborne outbreak. _BMC research notes, 8(1)_, 728.

Rychlik, I., Gregorova, D., & Hradecka, H. (2006). Distribution and function of plasmids in Salmonella enterica. _Veterinary microbiology, 112(1)_, 1-10.

Shippy, D. C., Heintz, J. A., Albrecht, R. M., Eakley, N. M., & Fadl, A. A. (2012). Deletion of glucose-inhibited division (gidA) gene alters the morphological and replication characteristics of Salmonella enterica Serovar typhimurium. _Archives of microbiology, 194(6)_, 405-412.

Yang, Y., Du, W., Li, Y., Lei, J., & Pan, W. (2025). Recent advances and challenges in de novo genome assembly. _Genomics Communications, 2_(1).

Yue, M., & Schifferli, D. M. (2014). Allelic variation in Salmonella: an underappreciated driver of adaptation and virulence. _Frontiers in microbiology, 4_, 419.






