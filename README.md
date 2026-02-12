# Minipoa: A minimizer-based method for fast and memory-efficient partial order alignment
## How to use?

1. Install

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

2. General usage
   ```bash
   # Sequencing mode for consensus calling
   minipoa input.fasta > output.fasta
   # MSA mode for multiple sequence alignment
   minipoa input.fasta -S -r1 > output.fasta
   # Generate graph information in GFA format
   minipoa input.fasta -S -r2 > output.gfa
   # Show help message
   minipoa -h
   ```


## License

minipoa is released under the MIT license. See the file LICENSE for more details.

## Contact

If you have any problems, please email me.(3300236038@qq.com)