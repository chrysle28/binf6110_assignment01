# Genome assembly and comparison of Salmonella enterica
The objective of this assignment is to generate a de novo assembly of the _Salmonella enterica_ genome from long sequence reads, align the assembly to an established reference genome, and conduct and visualize variant calls.

## **1 | Introduction**
Genome assembly methods aim to generate complete chromosome-level sequences from raw sequence reads, which in turn supports key biological initiatives [^1]. The presence of complete genomes allows for the linking of SNPs to phenotypic effects, the determination of organism variants, the creation of high-quality reference genomes, and myriad other studies. [^2].  

A genome can be assembled from a reference, this being reference-based assembly, or from scratch, as in _de novo_ assembly. As its name implies, _de novo_ assembly eschews the use of a reference sequence, and thus the quality of the assembly depends on the continuity of the contigs [^2]. Consequently, _de novo_ assembly may struggle with genomes containing short and inaccurate reads, or an abundance of repeats, leading to inaccurate and/or fragmented reconstructions [^2,3]. Notably, reference-based assembly was found to outperform _de novo_ assembly for plasmid reconstruction [4].However, these challenges mainly stem from the use of short reads in de novo assembly. 

Short-read sequencing is more cost-effective, but assembly is inherently limited by the length of the reads; as such, it is difficult to construct chromosome-level assemblies using Sanger and second-generation sequencing [^3,5]. In contrast, long reads can not only refine _de novo_ construction of genomes, but also improve the detection of variants [^3]. Long reads based on third-generation sequencing had higher error rates and could not handle complex genomic regions (e.g. tandem repeats) very well [^3], but next-generation sequencing methods such as PacBio and Oxford Nanopore Technologies (ONT) have addressed these issues. These methods offer longer read lengths suitable for de novo genome assembly with little to no trade-off in per-read accuracy, as evidenced by PacBio HiFi reads boasting an error rate of less than 0.1%-1% [^3,5]. Therefore, genome assembly has pivoted into the use of long read sequences. However, it should be noted that despite the many technological and methodological advancements being made, genome assembly still faces challenges such as high costs, biological constraints clashing with current technologies (e.g. ultra-long tandem repeats and polyploid genomes), and current assemblers being too time-intensive and inefficient [^6]. 

Highlighting the last point, it is critical that the pipelines and tools used are not only suitable for the type of assembly being performed (i.e. _de novo_ vs reference-based), but also the type of read being used (i.e. short-read vs long-read). Focusing on de novo assemblies using long-reads, assemblers such as Flye and Canu have been the most used and tested. Both use a different approach to assembly - Canu uses the Overlap Layout Consensus Method (OLC) wherein contigs are created from overlaps between the reads [^7], while Flye uses a hybrid approach involving repeat graphs (i.e. disjointigs) [^8]. Many systematic reviews have found that Flye consistently outperforms Canu for assemblies using ONT reads in metrics such as computational speed, sequence identity, contiguity, misassembly count, and gene identification [^2,9,10]. Moreover, Canu produced the most fragmented reads out of all long-read assemblers tested, and also sometimes failed to complete with certain inputs (e.g. reads processed with Filtlong) [^10]. Given Flye’s consistency in terms of producing the most accurate assemblies from long reads, it will be the main tool used for _de novo_ assembly of the _Salmonella enterica_ genome.

## **2 | Methods**
### 2.1 | Description of Data
Long read sequences of a _Salmonella enterica_ isolate (accession SRR32410565) will be assembled into a genome. Reads were generated via Oxford Nanopore, with R10 chemistry, producing an expected accuracy of Q20+ and N50 of 5-15kb. Raw read data will be obtained in FASTQ format from https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR32410565. The reference _Salmonella enterica_ genome will be obtained through the NCBI genome database (assembly ASM694v2) in FASTA format.

### 2.2 | Inspection and Quality Control
Long read sequences will be inspected for quality and length using a combination of [Nanoplot (v1.46.2)](https://github.com/wdecoster/NanoPlot) and [Filtlong (v0.3.1)](https://github.com/rrwick/Filtlong). Nanoplot will be used to generate sequence statistics, which will inform the subsequent processing of reads with Filtlong. Sequences below a certain length and accuracy threshold will be removed.
```
NanoPlot -fastq <input.fastq>  -o <output_dir>
filtlong -min_length <#> –keep_percent <#> input.fastq > output.fastq
```
> The values for the min_length and keep_percent parameters will be determined through inspection of the sequence statistics generated by Nanoplot.

### 2.3 | Longread Assembly
After filtering and quality control, the genome will be assembled from the long reads using [Flye (v2.9.6)](https://github.com/mikolmogorov/Flye).
```
flye -t <#> –genome-size <#> –nano-hq input.fastq -o output_dir
```
> Default parameters will be used, and the value used for –genome-size will be 809.3m. The –nano-hq flag will be employed due to the raw reads expected accuracy of Q20+.

### 2.4 | Alignment
The de novo assembly will be aligned with the reference genome using [Minimap2 (v2.30)](https://github.com/lh3/minimap2). Two different presets will be used: one for mapping of reads to the reference, and one for complete genome-to-genome alignment.
```
minimap2 -ax map-ont -t <#> assembly.fasta reference.fasta > output.sam
minimap2 -ax asm5 -t <#> assembly.fasta reference.fasta > output.sam
```
> The asm5 parameter will be employed as it is expected that the assembly and reference will have less than 5% divergence.

### 2.5 | Variant Calling
The SAM files from the previous step will be converted to BAM format, sorted, and then indexed using [samtools (v1.23)](https://github.com/samtools/samtools).
```
samtools sort query.bam -o query.sorted.bam
samtools index query.sorted.bam
```
Variant calling will be performed using bcftools.
```
bcftools mpileup -Ou –max-depth <#> -f reference.fasta query.bam | bcftools call -mv -Ob -o calls.bcf
```
> Default parameters for variant calling will be used.

### 2.6 | Visualization
Visualization of the assembled genome will be done using [IGV (v2.19.7)](https://github.com/igvteam/igv). The following files will be inputted: reference genome FASTA (ASM694v2), assembly BAM, and read alignment BAM.

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

