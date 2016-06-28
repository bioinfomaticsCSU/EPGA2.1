License
=========

Copyright (C) 2014 Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


De Novo Assembler
=================
1) Introduction

	EPGA2 updates some modules in EPGA which can improve memory efficiency in genome asssembly. 
	EPGA2 adopts new scaffolding method which can accurately determine orientations and orders of contigs. 
	The read library for EPGA2 should be paired-end or mate-paired reads (fastq). Read length shorter than 50bp and coverage larger than 100.
	
2) Before installing and running
	
	Users should install Bowtie2 and Samtools firstly and add them to your PATH. EPGA2 uses Bowtie2 for mapping read to contigs in the step of scaffolding, and uses Samtools for converting ".sam" file to ".bam" file. 
	Users can download Bowtie2 from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml 
	Samtools is available from http://samtools.sourceforge.net/index.shtml

3) Installing

	EPGA is written C++ and therefore will require a machine with GNU C++ pre-installed.
	Create a main directory (eg:EPGA2). Copy all source code to this directory.
	Type "make all".

3) Running.

	Run command line: 
	
	ulimit -n 1100 //this command is used for BCALM 
	perl EPGA.pl library.txt kmerLength threadNumber

		<library.txt>:
			Each line represents one read library.
			The first column is the first mate read file (*.fastq);
			The second column is the second mate read file (*.fastq);
			The third column is length of read;
			The fourth column is insert size of read library;
			The fifth column is standard deviation of insert size;
			The sixth column represents whether the read library is mate-paired (0 denotes paired-end reads, 1 denotes mate-paired reads);
		<kmerLength>:
			One integer (<32) which should be shorter than read length.
		<threadNumber>:
			Thread number of program.

4) Output:

	The contig set and scaffold set produced by EPGA2 are named result_ContigSet.fa and result_ScaffoldSet.fa respectively.

5) Example：

	one line in library.txt:
	
		frag_1.fastq frag_2.fastq 101 180 20 0

References
=================

To cite EPGA please use the following citation:

	[1] J. Luo, J. Wang* , Z. Zhang, F.X. Wu, M. Li, and Y. Pan. EPGA: de novo assembly using the distributions of reads and insert size. Bioinformatics, 2015, 31(6):825-833.
	[2] J. Luo, J. Wang* , Z. Zhang, F.X. Wu, M. Li, and Y. Pan. EPGA2：memory-efficient de novo assembler. Bioinformatics, 2015, 31(24):3988-3990.


