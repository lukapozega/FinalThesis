# Final BSc thesis (computer science - 2018/2019)

Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Installation

Tool can be installed by running bash script

```bash
bash install.sh
```

## Usage

You can run tool either by passing reads and genome to script

```bash
bash run.sh <reference genome> <reads file>
```

or by passing all needed files to tool directly

```bash
build/bin/assembly <PAF file> <reference file> <repeats file> <reads file>
```

Reads file need to be in FASTA or FASTQ format, PAF file is generated from Minimap tool and repeats file can be generated either by RED tool or by tool written by Sara Bakić: https://github.com/lbcb-edu/BSc-thesis-18-19/tree/sbakic

## Data

Tool can be tested on an Oxford Nanopore Technologies data set obtained by sequencing the Escherichia coli K-12 substr. MG1655 genome. The data set is freely available from Loman Labs [here](https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta), while the reference genome is freely available from NCBI [here](https://bit.ly/2PCYHWr).

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.
