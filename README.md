# minimap
> A respository for minipoa project

## How to use?

1. git clone this repository
   ```bash
   git clone https://github.com/NCl3-lhd/minipoa.git
   ```

2. make
   ```bash
   cd minipoa && mkdir build && cd build
   # AVX2(default)
   cmake -DENABLE_AVX2=ON -DENABLE_AVX512=OFF -DENABLE_SSE2=OFF .. && make
   # AVX512
   cmake -DENABLE_AVX512=ON -DENABLE_AVX2=OFF -DENABLE_SSE2=OFF .. && make
   # SSE2
   cmake -DENABLE_SSE2=ON -DENABLE_AVX2=OFF -DENABLE_AVX512=OFF .. && make
   ```

3. run
   ```bash
   # consensus calling
   ./bin/minipoa ./data/mtDNA.fasta > mtDNA.cons
   # multiple sequence alignment
   ./bin/minipoa ./data/mtDNA.fasta -S -t 12 -r1 > mtDNA.msa
   ```

4. clean
   ```bash
   make clean
   ```

## License

minipoa is released under the MIT license. See the file LICENSE for more details.

## Contact

If you have any problems, please email me.(3300236038@qq.com)