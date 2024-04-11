# **DEGAP**## Dynamic Elongation of a Genome Assembly Path**DEGAP degaps gaps!****DEGAP** is a novel gap-filling software that resolves gap regions by utilizing the dual advantages of accuracy and length of high-fidelity (HiFi) reads.## System requirements**DEGAP** has been developed with the Python programming language under a Linux environment. The following packages/tools (older versions might not work properly as expected) are required:Biopython (version 1.80)Hifiasm (0.16.1-r375)Minimap2 (2.17-r941)MUMmer (4.0.0beta2)Pysam (version 0.20.0)SAMtools (version 1.6)## Tutorial**DEGAP** provides 2 modes (GapFiller, CtgLinker) to fill the gaps in difference scenarios.GapFiller is designed for one specific gap, of which the exactly left and right sequences between this gap are already know. Gapfiller is designed to fill the gap by elongating the sequence from only one direction.> python ./bin/DEGAP.py --mode gapfiller --seqleft ./path/gapLeftSequence.fasta --seqright ./path/gapRightSequence.fasta --reads ./path/HiFiReads.fasta -o ./path/Output --flag leftCtgLinker is designed for the assembly with gaps. If sequences in an assembly are unordered, Ctglinker can try to not only fill the gaps, but also link the potentially neighbored sequences by elongating the edges.> python ./bin/DEGAP.py --mode ctglinker --ctgseq ./path/contigs.fasta --reads ./path/HiFiReads.fasta --out ./path/Output --filterDepth 0.3**DEGAP** also provides parameter --filterDepth num to filter the HiFi reads and --filterDepthOnly to only filter the HiFi reads but not execute the elongation process.> python ./bin/DEGAP.py --mode ctglinker --ctgseq ./path/contigs.fasta --reads ./path/HiFiReads.fasta --out ./path/Output --filterDepthOnly 0.3> python ./bin/DEGAP.py --mode gapfiller --seqleft ./path/gapLeftSequence.fasta --seqright ./path/gapRightSequence.fasta --reads ./path/HiFiReads.fasta -o ./path/Output --filterDepthOnly 0.3More details on how to run DEGAP will be added shortly. Any comments and suggestions, please contact <jzhang@mail.hzau.edu.cn>.## How to cite?Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. Briefings in Bioinformatics, 2024, in press. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)