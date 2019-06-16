## Introduction

Dna-nn implements a proof-of-concept deep-learning model to learn relatively
simple features on DNA sequences. So far it has been trained to identify the
(ATTCC)n and alpha satellite repeats, which occupy over 5% of the human genome.
Taking the RepeatMasker annotations as the ground truth, dna-nn is able to
achieve high classification accuracy. It can annotate a human PacBio assembly
in less than 1.5 hours using 16 CPUs, while RepeatMasker may take several days
on 32 CPUs. Dna-nn is a practical tool to find these two types of
satellites.

Dna-nn may have potentials to learn other types of sequence features. It can
accurately identify Alu repeats as well. However, it has low sensitivity to
Beta satellites and fails to learn L1 repeats even with more hidden neurons.


## Installation

Dna-nn is implemented in C and includes the source code of the [KANN][kann]
deep-learning framework. The only external dependency is [zlib][zlib]. To
compile,
```sh
git clone https://github.com/lh3/dna-nn
cd dna-nn
make
```


## Usage

### Applying a trained model

To find (ATTCC)n and alpha satellites for long contigs,
```sh
./dna-brnn -Ai models/attcc-alpha.knm -t16 seq.fa > seq.bed
```
The output is a BED file. A label `1` on the 4th column indicates the interval
is a region of (AATTC)n ; label `2` indicates a region of alpha satellites.

### Training

Training dna-nn requires sequences in the FASTQ format, where each "base
quality" indicates the label of the corresponding base.

The following command lines shows how we generate the pre-trained model
`attcc-alpha.knm`.
```sh
# Install the k8 javascript shell (if this has not been done)
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8              # or copy it to a directory on your $PATH
# Run RepeatMasker to generate truth data
RepeatMasker -species human -pa 16 -e ncbi -xsmall -small -dir . train.fa
# The last column indicates the label of each region in the output BED
./k8 parse-rm.js train.fa.out > train.rm.bed
# Generate training data in FASTQ. Base qualities indicate labels.
./gen-fq -m2 train.fa train.rm.bed > train.lb2.fq
# Training. We trained 10 models with different random seeds
./dna-brnn -t8 -n32 -b5m -m50 -s14 -o attcc-alpha.knm train.lb2.fq
```

### Evaluation

With truth annotations in the FASTQ format, we can evaluate the accuracy of
a model with
```sh
./dna-brnn -Ei models/attcc-alpha.knm -t16 seq.fa > /dev/null
```
The stderr output gives the accuracy for each label.


## Citing dna-nn

If you use dna-nn in your work, please cite its [paper][pub]:
```txt
Li, H (2019) Identifying centromeric satellites with dna-brnn,
*Bioinformatics*, doi:10.1093/bioinformatics/btz264
```

[zlib]: https://zlib.net/
[kann]: https://github.com/attractivechaos/kann
[pub]: https://www.ncbi.nlm.nih.gov/pubmed/30989183
