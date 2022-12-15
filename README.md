# gblr: genotyping by long reads

Allele calling/genotyping approaches using long reads for repetitive polymorphic loci.

## gblr.py
There are three overall models to choose from: 

#### alignment
Aligns each read to each known allele and uses edit distances to determine the most likely allele/genotype. Outputs a list of alleles/genotypes and the probabilistic scores for each, where the largest score indiciates the called allele/genotype.

#### consensus
Generates a partial-order alignment graph from all of the reads and from there determines the one or two consensus sequences that best describe the genotype. The consensus sequences are compared to the sequences for known alleles and if there is no match, a 'novel allele' flag is used to indicate the allele with the lowest edit distance to the consensus sequence. Only one genotype is output with and there is no likelihood score.

#### combined
Performs the alignment model, then the consensus model to check for novel alleles in the top genotype. Outputs a list of alleles/genotypes and the probabilistic scores for each, where the largest score indiciates the called allele/genotype, except novel alleles are given an arbitrary score of 1 to bring them to the top of the list. If the most likely genotype determined by the consensus model is homozygous but the most likely genotype determined by the alignment model is heterozygous, a 'potential false allele' flag is used.

## gap-filter.py
When dealing with noisy reads due to PCR laddering of known repeat units, indels are counted and those that occur in less than 25% of the reads are ignored. Indels of an unexpected length are also ignored. Outputs a list of reads to keep that can be fed into `samtools view -N' to filter a bam file.
