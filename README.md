#### Overview
OWGS: An optimized whole genome sequencing pipeline for the annotation and analysis of bacterial genomes.

#### Installation
```
git clone https://github.com/lakinsm/wgs-meg.git
cd wgs-meg
./multi_sample_control_script.sh
```

#### Dependencies
* [iMetAMOS](http://metamos.readthedocs.io/en/v1.5rc3/content/imetamos.html) - Genome Assembly
* [Mauve](http://www.bioinformatics.org/wiki/Mauve) - Multiple Genome Alignment
* [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml) - Genome Annotation
* [CISA](http://sb.nhri.org.tw/CISA/en/Instruction) - Contig Integrator
* [BBMap](https://wiki.gacrc.uga.edu/wiki/BBMap) - Insert Size Estimator
* [BWA](http://bio-bwa.sourceforge.net) - Sequence Aligner
* [Freebayes](https://github.com/ekg/freebayes) - Variant Caller
* [SAMtools/BCFtools](https://samtools.github.io/bcftools/) - VCF Manipulation
* [BLASTn & Diamond BLASTx](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - Search Nucleotide Databases
* [kSNP3](https://sourceforge.net/projects/ksnp/) - KMER Selection

#### Usage
Usage: ./multi_sample_control_script.sh [options]

Pipeline Options:

    -i      DIR     Input directory of paired fastq files
    -a      BOOL    Run assembly pipeline on this sample
    -o      DIR     Directory for output of important files
    -s      STR     Output prefix for output files
    -spp    STR     Species identifier string
    -t      INT     Number of threads to use where applicable
    -td     DIR     Temp directory for intermediate dirs/files
    -T      INT     Qual threshold for N-masking in kSNP

#### Example
TODO: Put stuff here

#### Output
TODO: Put stuff here
