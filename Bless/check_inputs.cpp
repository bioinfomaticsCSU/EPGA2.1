#include "check_inputs.hpp"



//----------------------------------------------------------------------
// check_read_file
//----------------------------------------------------------------------
void C_check_read::check_read_file(const C_arg& c_inst_args, C_time& c_inst_time) {
   // start measuring run-time
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_check_read_file = asctime(localtime(&rawtime));

   f_log.open(c_inst_args.log_file_name.c_str(), std::fstream::app);

   if (f_log.is_open()) {
      std::cout << "Checking input read files" << std::endl;

      f_log << "Checking input read files" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // paired
   if (c_inst_args.paired_read == true) {
      check_read_file_fastq_paired(c_inst_args);
   }
   // single
   else {
      check_read_file_fastq_single(c_inst_args);
   }

   f_log.close();

   time(&rawtime);
   c_inst_time.end_check_read_file = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// check_read_file_fastq_single
//----------------------------------------------------------------------
void C_check_read::check_read_file_fastq_single(const C_arg& c_inst_args) {
   // open input read files
   std::ifstream f_read(c_inst_args.read_file_name.c_str());

   // initialize variables
   std::size_t num_reads(0);

   bool set_read_length(false);

   max_quality_score = -1000;
   min_quality_score = 1000;

   std::string line;

   // header
   getline(f_read, line);

   while(!f_read.eof()) {
      // increment the number of reads
      num_reads++;

      // sequence
      getline(f_read, line);
      
      // check read length
      if (set_read_length == false) {
         read_length = line.length();

         if (read_length < c_inst_args.kmer_length) {
            std::cout << std::endl << "ERROR: K-mer length(" << c_inst_args.kmer_length << ") is longer than read length(" << read_length << ")" << std::endl << std::endl;
            exit(EXIT_FAILURE);            
         }

         set_read_length = true;
      }
      else if (line.length() != read_length) {
         std::cout << std::endl << "ERROR: Read length is different" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      
      // connector
      getline(f_read, line);

      // quality score
      getline(f_read, line);
      for (std::size_t it = 0; it < line.length(); it++) {
         if ((int)line[it] > max_quality_score) {
            max_quality_score = (int)line[it];
         }

         if ((int)line[it] < min_quality_score) {
            min_quality_score = (int)line[it];
         }
      }

      // header
      getline(f_read, line);
   }


   if (min_quality_score < MIN_SCORE) {
      std::cout << "ERROR: Illegal quality score" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if ((min_quality_score >= MIN_SCORE) && (min_quality_score < MIN_64_SCORE)) {
      quality_score_offset = PHRED33;
   }
   else {
      quality_score_offset = PHRED64;
   }

   f_read.close();

   std::cout << "     Number of reads     : " << num_reads << std::endl;
   std::cout << "     Read length         : " << read_length << std::endl;
   std::cout << "     Quality score offset: " << quality_score_offset << std::endl;
   std::cout << "     Checking input read files: done" << std::endl << std::endl;

   f_log << "     Number of reads     : " << num_reads << std::endl;
   f_log << "     Read length         : " << read_length << std::endl;
   f_log << "     Quality score offset: " << quality_score_offset << std::endl;
   f_log << "     Checking input read files: done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// check_read_file_fastq_paired
//----------------------------------------------------------------------
void C_check_read::check_read_file_fastq_paired(const C_arg& c_inst_args) {
   // open input read files
   std::ifstream f_read_1(c_inst_args.read_file_name1.c_str());
   std::ifstream f_read_2(c_inst_args.read_file_name2.c_str());

   // initialize variables
   std::size_t num_reads_1(0);
   std::size_t num_reads_2(0);
   std::size_t read_length_1(0);
   std::size_t read_length_2(0);

   bool set_read_length(false);

   max_quality_score = -1000;
   min_quality_score = 1000;

   std::string line;

   // header
   getline(f_read_1, line);

   while(!f_read_1.eof()) {
      // increment the number of reads
      num_reads_1++;

      // sequence
      getline(f_read_1, line);
      
      // check read length
      if (set_read_length == false) {
         read_length_1 = line.length();

         if (read_length_1 < c_inst_args.kmer_length) {
            std::cout << std::endl << "ERROR: K-mer length(" << c_inst_args.kmer_length << ") is longer than read length(" << read_length_1 << ")" << std::endl << std::endl;
            exit(EXIT_FAILURE);            
         }

         set_read_length = true;
      }
      else if (line.length() != read_length_1) {
         std::cout << std::endl << "ERROR: Read length is different(1st file)" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      
      // connector
      getline(f_read_1, line);

      // quality score
      getline(f_read_1, line);
      for (std::size_t it = 0; it < line.length(); it++) {
         if ((int)line[it] > max_quality_score) {
            max_quality_score = (int)line[it];
         }

         if ((int)line[it] < min_quality_score) {
            min_quality_score = (int)line[it];
         }
      }

      // header
      getline(f_read_1, line);
   }

   // read the 2nd file
   read_length_2 = read_length_1;

   // header
   getline(f_read_2, line);

   while(!f_read_2.eof()) {
      // increment the number of reads
      num_reads_2++;

      // sequence
      getline(f_read_2, line);

      if (line.length() != read_length_2) {
         std::cout << std::endl << "ERROR: Read length is different(2nd file)" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // connector
      getline(f_read_2, line);

      // quality score
      getline(f_read_2, line);
      for (std::size_t it = 0; it < line.length(); it++) {
         if ((int)line[it] > max_quality_score) {
            max_quality_score = (int)line[it];
         }

         if ((int)line[it] < min_quality_score) {
            min_quality_score = (int)line[it];
         }
      }

      // header
      getline(f_read_2, line);
   }

   if (min_quality_score < MIN_SCORE) {
      std::cout << "ERROR: Illegal quality score" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if ((min_quality_score >= MIN_SCORE) && (min_quality_score < MIN_64_SCORE)) {
      quality_score_offset = PHRED33;
   }
   else {
      quality_score_offset = PHRED64;
   }

   // compair the number of reads in two input files
   read_length = read_length_1;
   if (num_reads_1 == num_reads_2) {
      num_reads = num_reads_1 + num_reads_2;
   }
   else {
      std::cout << "ERROR: Number of lines in two input files are not same" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // compair the read length in two input files
   if (read_length_1 == read_length_2) {
      read_length = read_length_1;
   }
   else {
      std::cout << "ERROR: The read length of two input files are not same " << read_length_1 << " vs. " << read_length_2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_read_1.close();
   f_read_2.close();

   std::cout << "     Number of reads     : " << num_reads << std::endl;
   std::cout << "     Read length         : " << read_length << std::endl;
   std::cout << "     Quality score offset: " << quality_score_offset << std::endl;
   std::cout << "     Checking input read files: done" << std::endl << std::endl;

   f_log << "     Number of reads     : " << num_reads << std::endl;
   f_log << "     Read length         : " << read_length << std::endl;
   f_log << "     Quality score offset: " << quality_score_offset << std::endl;
   f_log << "     Checking input read files: done" << std::endl << std::endl;
}
