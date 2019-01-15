## ATTCC and Alpha centromeric satellites

* Purpose: to identify (ATTCC)n and alpha centromeric [satellites][sat] in
  human
* Model file: `attcc-alpha.knm`
* Lables: `1` for the (ATTCC)n microsatellite; `2` for the alpha satellite.
* Training data: chr11 in GRCh37, alpha satellites on other chromosomes and
  GRCh37 decoy. Annotations on GRCh37 was acquired from UCSC. Annotations on
  decoys were generated with [RepeatMasker][rm]-4.0.8, [rmBlast][rmb]-2.6.0+
  and [RepBase][rpb]-23.11.
* Training command line: `dna-brnn -n32 -b5m -m50 -B256 -s14 -o attcc-alpha.knm
  train.fq.gz`
* Validation data: GRCh38 decoy
* Testing data: randomly selected 326Mb sequences from CHM1\_FC\_P6

[sat]: https://en.wikipedia.org/wiki/Satellite_DNA
[rm]: http://www.repeatmasker.org/
[rmb]: http://www.repeatmasker.org/RMBlast.html
[rpb]: https://www.girinst.org/repbase/
