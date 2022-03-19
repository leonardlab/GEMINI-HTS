
NGS Analysis Code
=================

This code use used for analyzing High Throughput Sequencing (HTS) data, here referred to as Next-Gen Sequencing (NGS) data.

To run/use this code:

- python3 is required

- this requires the "Segment Processor" currently implemented as part
  of SPATS (on the `segments` branch), check it out:

```
$ git clone https://github.com/LucksLab/spats
$ git checkout segments
```

- set the `NGS_HOME` environment variable to point to the root of this repository

- set the `SPATS_HOME` environment variable to point to the root of the SPATS repository

- create a folder that holds the R1/R2 FASTQ files from the experiment

- configure your experiment: see the `ngs_experiment2_keys()` method
  in `ngs/analysis.py` for an example of how to specify the R1/R2
  files for the experiment. in the current code, it expects the files
  to be named `s1_full_r1.fastq`, `s1_full_r2.fastq`,
  `s2_full_r1.fastq`, etc.

- execute `analysis.sh` in that same folder, e.g.,
```
$ $NGS_HOME/bin/analysis.sh
```
