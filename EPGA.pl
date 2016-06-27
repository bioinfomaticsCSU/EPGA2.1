use strict;
	
	if(@ARGV!=3){
		print "Usage: perl EPGA.pl [library] [kmer_length] [thread_number]\n";
		die "too less or too much parameter!\n";
	}
	
	my $library_number = 0;
	my @library_name;
	my @library_insertsize;
	my @library_readLength;
	my @library_std;
	my @library_orientation;
	my @library_kmer;
	my @library_sam;
	my @library_bam;
	
	my $library_information = shift;
	my $kmer_length = shift;
	my $thread_number = shift;
	my $line;
	my $time_output = "time.txt";
	open(TIME, ">$time_output")or die $!;
	open(INPUT,$library_information)or die $!;
	print TIME "bless start--:";
	print TIME scalar localtime;
	print TIME "\n";
	while($line=<INPUT>){
	    chomp $line;
		my @infor = split /\s+/, $line;
		$library_name[$library_number*2]=$infor[0];
		$library_name[$library_number*2+1]=$infor[1];
		$library_readLength[$library_number] = $infor[2];
		$library_insertsize[$library_number]=$infor[3];
		$library_std[$library_number]=$infor[4];
		$library_orientation[$library_number]=$infor[5];
		$library_kmer[$library_number] = "library_"."$library_number";
		$library_number++;
	}
	
	my @temp;
	@temp = ("chmod +x bcalm dsk bless GetKmerHash KmerToDot parse_results SimplePathToGraph");
	`@temp`;
	
	@temp = ("ulimit -n 1100");
	`@temp`;
	
	
	for(my $k=0;$k<$library_number;$k++){
        my $output_prefix = "library_"."$k";
        @temp = ("./bless -read1 $library_name[$k*2] -read2 $library_name[$k*2+1] -prefix $output_prefix -kmerlength $kmer_length");
		`@temp`;
		unlink glob $output_prefix.".h*";
		unlink glob $output_prefix.".b*";
		unlink glob $output_prefix.".l*";
		$library_name[$k*2] = $output_prefix.".1.corrected.fastq";
		$library_name[$k*2+1] = $output_prefix.".2.corrected.fastq";
	}
	print TIME "DSK start--:";
	print TIME scalar localtime;
	print TIME "\n";

	for(my $k=0;$k<$library_number;$k++){
	    my $dsk_library_name = "dsk_library_name";
		open(DSK, ">$dsk_library_name")or die $!;
	    print DSK "$library_name[$k*2]\n";
		print DSK "$library_name[$k*2+1]\n";
        @temp = ("./dsk $dsk_library_name $kmer_length -t 2 -m 1");
		`@temp`;
		my $output_prefix = $library_kmer[$k].".dot";
		@temp = ("./parse_results $dsk_library_name".".solid_kmers_binary > $output_prefix");
		`@temp`;
		unlink "$dsk_library_name".".solid_kmers_binary";
		unlink "$dsk_library_name".".reads_binary";
		unlink $dsk_library_name;
	}
	print TIME "GetKmerHash start--:";
	print TIME scalar localtime;
	print TIME "\n";
	my $kmer_hash_count = 0;
	my $temp_file;
	for(my $k=0;$k<$library_number;$k++){
        $temp_file = $library_kmer[$k].".dot";
	    @temp = ("cat $temp_file |wc -l");
		my $temp = `@temp`;
        $kmer_hash_count = $kmer_hash_count + $temp;
	}
	$kmer_hash_count = int($kmer_hash_count*1.2);
	
	my $all_kmer = "allkmer.dot";
	my $kmer_hash = "kmerhash.b";
	my $argument = $library_kmer[0].".dot";
	for(my $k=1;$k<$library_number;$k++){
	    $argument = $argument.' '.$library_kmer[$k].".dot";
	}
	$argument = $argument.' '.$kmer_hash_count.' '.$kmer_length.' '.$thread_number.' '.$all_kmer.' '.$kmer_hash;
	@temp = ("./GetKmerHash $argument");
	`@temp`;
	
	print TIME "DSK k+1-mer start--:";
	print TIME scalar localtime;
	print TIME "\n";
	my $new_all_kmer = "all_k+1.dot";
	my $dsk_library_name = "dsk_library_name";
	open(DSK, ">$dsk_library_name")or die $!;
	for(my $k=0;$k<$library_number;$k++){
	    print DSK "$library_name[$k*2]\n";
		print DSK "$library_name[$k*2+1]\n";
	}
	$kmer_length = $kmer_length + 1;
	@temp = ("./dsk $dsk_library_name $kmer_length -t 2 -m 1");
	`@temp`;
	@temp = ("./parse_results $dsk_library_name".".solid_kmers_binary > $new_all_kmer");
	`@temp`;
	unlink "$dsk_library_name".".reads_binary";
	unlink "$dsk_library_name".".solid_kmers_binary";
	unlink $dsk_library_name;

	$kmer_length = $kmer_length - 1;
	
	my $temp_output = "temp.dot";
	print TIME "./KmerToDot $new_all_kmer $temp_output\n";
	@temp = ("./KmerToDot $new_all_kmer $temp_output");
	`@temp`;
	unlink $new_all_kmer;
	rename $temp_output, $new_all_kmer;
	
	print TIME "bcalm start--:";
	print TIME scalar localtime;
	print TIME "\n";
	
	my $debruijn_graph = "DBGraph.fa";
	my $simple_path = "path.dot";
	print TIME "./bcalm $all_kmer -s $new_all_kmer -o $simple_path\n";
    @temp = ("./bcalm $all_kmer -s $new_all_kmer -o $simple_path");
	`@temp`;
	print TIME "SimplePathToGraph start--:";
	print TIME scalar localtime;
	print TIME "\n";
	@temp = ("cat $simple_path |wc -l");
	my $temp = `@temp`;
	my $graph_node = int($temp);

	$kmer_length--;
	
	print TIME "./SimplePathToGraph $simple_path $debruijn_graph $graph_node $kmer_length $thread_number\n";
	@temp = ("./SimplePathToGraph $simple_path $debruijn_graph $graph_node $kmer_length $thread_number");
	`@temp`;
	
	print TIME "EPGA start--:";
	print TIME scalar localtime;
	print TIME "\n";
	$kmer_length++;

	$temp_file = $library_name[0].' '.$library_name[1].' '.$library_insertsize[0].' '.$library_std[0].' '.$library_orientation[0];
	for(my $k=1;$k<$library_number;$k++){
        $temp_file = $temp_file.' '.$library_name[2*$k].' '.$library_name[2*$k+1].' '.$library_insertsize[$k].' '.$library_std[$k].' '.$library_orientation[$k];
	}
	print TIME "$kmer_hash_count\n";
	print TIME "./EPGA $temp_file $debruijn_graph $kmer_hash $kmer_hash_count $kmer_length $thread_number\n";
	@temp = ("./EPGA $temp_file $debruijn_graph $kmer_hash $kmer_hash_count $kmer_length $thread_number");
	`@temp`;
	
	unlink glob "allkmer.dot";
	unlink glob "all_k+1.dot";
	unlink glob "path.dot";
	unlink glob "CompactContig.fa";
	unlink glob "contigSet.fa";
	unlink glob "nonInclude.fa";
	unlink glob "overlap11.fa";
	unlink glob "SubContig.fa";
	
	
	for(my $k=0;$k<$library_number;$k++){
		unlink glob "library_"."$k".".dot";
	}
	
	$temp_file = "result_ContigSet.fa";
	
	for(my $k=0;$k<$library_number;$k++){
		$library_sam[2*$k] = "library_$k"."_left.sam";
		$library_sam[2*$k+1] = "library_$k"."_right.sam";
		$library_bam[2*$k] = "library_$k"."_left.bam";
		$library_bam[2*$k+1] = "library_$k"."_right.bam";
		@temp = ("bowtie2-build $temp_file contigs");
		`@temp`;
		@temp = ("bowtie2 -x contigs $library_name[2*$k] -S $library_sam[2*$k]");
		`@temp`;
		@temp = ("bowtie2 -x contigs $library_name[2*$k+1] -S $library_sam[2*$k+1]");
		`@temp`;
		@temp = ("samtools view -Sb $library_sam[2*$k] > $library_bam[2*$k]");
		`@temp`;
		@temp = ("samtools view -Sb $library_sam[2*$k+1] > $library_bam[2*$k+1]");
		`@temp`;
		unlink glob $library_sam[2*$k];
		unlink glob $library_sam[2*$k+1];
	}

	unlink glob "contigs.*";
	
	my $percentage = 0.07;
	my $pairedNumber = 2;
	my $weightMin = 0.2;
	my $paired = 0;
	my $weightMethod = 0;
	for(my $k=0;$k<$library_number;$k++){
		if($library_orientation[$k] == 0){
			$paired = 1;
		}else{
			$paired = 0;
		}
		$temp_file = $temp_file.' '.$library_bam[2*$k].' '.$library_bam[2*$k+1].' '.$library_readLength[$k].' '.$library_insertsize[$k].' '.$percentage.' '.$weightMin.' '.$pairedNumber.' '.$paired.' '.$weightMethod;
	}
	
	print TIME "./boss $temp_file result\n";
	@temp = ("./boss $temp_file result");
	`@temp`;
	
	print TIME "end--:";
	print TIME scalar localtime;
	print TIME "\n";
	
	
