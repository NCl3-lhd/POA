# Minipoa: A minimizer-based method for fast and memory-efficient partial order alignment

[![Downloads](https://anaconda.org/malab/minipoa/badges/downloads.svg)](https://anaconda.org/malab/minipoa)
[![License](https://anaconda.org/malab/minipoa/badges/license.svg)](https://anaconda.org/malab/minipoa)
[![Platforms](https://anaconda.org/malab/minipoa/badges/platforms.svg)](https://anaconda.org/malab/minipoa)

[Minipoa: A minimizer-based method for fast and memory-efficient partial order alignment.](https://doi.org/10.64898/2026.02.18.706716)

## Install

Install via conda
```bash
conda install -c malab minipoa
```
Or, make from source
```bash
# Git clone this repository
git clone https://github.com/NCl3-lhd/minipoa.git
cd minipoa && mkdir build && cd build
# AVX2(default)
cmake -DENABLE_AVX2=ON -DENABLE_AVX512=OFF -DENABLE_SSE2=OFF .. && make
# AVX512
cmake -DENABLE_AVX512=ON -DENABLE_AVX2=OFF -DENABLE_SSE2=OFF .. && make
# SSE2
cmake -DENABLE_SSE2=ON -DENABLE_AVX2=OFF -DENABLE_AVX512=OFF .. && make
```

## Quick start

Generate a consensus sequence (Sequencing mode):
```bash
minipoa data/mtDNA.fasta > cons.fasta
```

Perform multiple sequence alignment (MSA mode):
```bash
minipoa data/mtDNA.fasta -S -r1 -t thread > mtDNA.fasta 
```

Output the sequence graph in GFA format:
```bash
minipoa data/mtDNA.fasta -S -r2 -t thread > mtDNA.gfa
```

View the full list of parameters:
```bash
minipoa -h
```
## Input
Minipoa supports input in FASTA, FASTQ, gzipped FASTA(.fa.gz), and gzipped FASTQ(.fq.gz) formats. It incrementally construct a alignment graph by input sequences. Optionally, an existing GFA file can be provided via `-i` to initialize the alignment graph.
```bash
minipoa input.fasta -i input.gfa -S -t thread -r2 > output.gfa
```
## Output
Minipoa provides three output modes, which can be selected using the `-r` parameter. Please note that the `-r` option is independent and is not affected by any other parameters.
* `-r0` : Output the consensus sequence in FASTA format (default)
* `-r1` : Output the multiple sequence alignment:
* `-r2` : Output the sequence graph in GFA format
## Advanced alignment control

To accommodate different datasets and balance speed with accuracy, minipoa provides the following advanced command-line parameters for fine-tuning:

### Alignment Scoring Parameters
You can customize the dynamic programming scoring scheme using the following parameters:
* `-M` : Match score
* `-X` : Mismatch penalty
* `-O` : Gap open penalty
* `-E` : Gap extension penalty
```bash
minipoa input.fasta -M 2 -X -4 -O -4 -E -2 > output.fasta
```
### Banding Strategy Control

`-B` : Enable the adaptive band strategy (static banding is enabled by default).
```bash
# Enable adaptive band strategy in Sequencing mode
minipoa input.fasta -B > output.fasta
```
The band length is automatically calculated based on the specified parameters and the query sequence length using the formula: `len = b + 1 / f * query_length`.
```bash
# Narrow the bandwidth under adaptive band strategy
minipoa input.fasta -B -b 10 -f 100 > output.fasta
```
`f = 0` : If you explicitly set `f = 0`, minipoa will not apply any banding strategy and will perform full dynamic programming across the entire matrix.

```bash
minipoa input.fasta -f 0 > output.fasta
```

### Anchor Chain Optimization & Acceleration
`-S` : Enable anchor chain optimization. Please note that enabling this option will **force-enable** the adaptive band strategy to maintain alignment robustness.

`-W` : Adjust the window distance  parameter for anchors.
* **Recommended Range**: Values between `500` and `query_length / 24` are generally appropriate.
* **Default Behavior**: To prevent incorrect alignments caused by false-positive anchors, the software conservatively defaults to the maximum value in this range.
* **Performance Tuning**: If you want to further accelerate the alignment process, you can manually specify a smaller `-W` value. This will increase speed, though it may result in a very slight loss of accuracy.
```bash
minipoa input.fasta -S -W 500 -t thread > output.fasta
```

## Contact

For any questions or issues, please contact me at ncl3.lhd@gmail.com.