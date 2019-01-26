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
./dna-brnn -Ai models/attcc-alpha.knm -t16 asm.fa > asm.bed
```
The output is a BED file. A label `1` on the 4th column indicates the interval
is a region of (AATTC)n ; label `2` indicates a region of alpha satellites.

### Training

Training and evaluation are more complex. I will document that part if there is
general interest.


## Citing dna-nn

If you use dna-nn in your work, please cite its [preprint][pub]:
```txt
Li, H (2019) Identifying centromeric satellites with dna-brnn. arXiv:1901.07327.
```

[zlib]: https://zlib.net/
[kann]: https://github.com/attractivechaos/kann
[pub]: https://arxiv.org/abs/1901.07327
