# POA
> A respository for POA project

## How to use?

1. git clone this repository
   ```bash
   git clone https://github.com/NCl3-lhd/POA.git
   ```

2. make
   ```bash
   cd POA
   cmake -S . -B build
   cmake --build build
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
   cmake --build build --target clean
   ```

## License

minipoa is released under the MIT license. See the file LICENSE for more details.

## Contact

If you have any problems, please email me.(3300236038@qq.com)