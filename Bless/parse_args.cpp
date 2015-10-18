#include "parse_args.hpp"



//----------------------------------------------------------------------
// read_args
//----------------------------------------------------------------------
void C_arg::read_args() {
   if (num_args == 1) {
      print_usage();
      exit(EXIT_SUCCESS);
   }

   //--------------------------------------------------
   // parse arguments
   //--------------------------------------------------
   for (int it_arg = 1; it_arg < num_args; it_arg++) {
      if (strcmp(args[it_arg], "-help") == 0) {
         print_usage();
         exit(EXIT_SUCCESS);
      }
      else if (strcmp(args[it_arg], "-read") == 0) {
         if (it_arg <= num_args - 2) {
            read_file_name = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The read file name is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-read1") == 0) {
         if (it_arg <= num_args - 2) {
            read_file_name1 = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The first read file name is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-read2") == 0) {
         if (it_arg <= num_args - 2) {
            read_file_name2 = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The second read file name is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-prefix") == 0) {
         if (it_arg <= num_args - 2) {
            prefix = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The prefix is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-kmerlength") == 0) {
         if (it_arg <= num_args - 2) {
            kmer_length = atoi(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The k-mer length is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-fpr") == 0) {
         if (it_arg <= num_args - 2) {
            target_false_positive_prob = atof(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The false positvie probability is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-seed") == 0) {
         if (it_arg <= num_args - 2) {
            random_seed = atol(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The seed is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-count") == 0) {
         if (it_arg <= num_args - 2) {
            kmer_occurrence_threshold = atoi(args[it_arg + 1]);
            set_kmer_occurrence_threshold = true;
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The k-mer occurrence threshold is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-cluster") == 0) {
         if (it_arg <= num_args - 2) {
            num_clusters = atoi(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The number of clusters is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-extend") == 0) {
         if (it_arg <= num_args - 2) {
            extend = atoi(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The max extension is not specified" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else if (strcmp(args[it_arg], "-nowrite") == 0) {
         nowrite = true;
      }
      else if (strcmp(args[it_arg], "-debug") == 0) {
         debug = true;
      }
      else if (strcmp(args[it_arg], "-noverify") == 0) {
         verify = false;
      }
      else {
         std::cout << std::endl << "ERROR: Illegal option " << args[it_arg] << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }


   //--------------------------------------------------
   // check options
   //--------------------------------------------------
   // paried-end input
   if (read_file_name.empty()) {
      paired_read = true;

      if (read_file_name1.empty()) {
         std::cout << std::endl << "ERROR: The first read file name is not specified" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      else {
         std::ifstream f_tmp;
         f_tmp.open(read_file_name1.c_str());
         if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << read_file_name1 << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
         f_tmp.close();
      }

      if (read_file_name2.empty()) {
         std::cout << std::endl << "ERROR: The second read file name is not specified" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      else {
         std::ifstream f_tmp;
         f_tmp.open(read_file_name2.c_str());
         if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << read_file_name2 << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
         f_tmp.close();
      }
   }
   // single-end input
   else {
      paired_read = false;

      std::ifstream f_tmp;
      f_tmp.open(read_file_name.c_str());
      if (f_tmp.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << read_file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      f_tmp.close();

      if (read_file_name1.empty() == false) {
         std::cout << std::endl << "ERROR: -read1 cannot be used with -read" << std::endl;
         exit(EXIT_FAILURE);
      }
      if (read_file_name2.empty() == false) {
         std::cout << std::endl << "ERROR: -read2 cannot be used with -read" << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   if (prefix.empty()) {
      std::cout << std::endl << "ERROR: The prefix is not specified" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (kmer_length == 0) {
      std::cout << std::endl << "ERROR: k-mer length not specified" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if ((kmer_length % 2) == 0) {
      std::cout << std::endl << "ERROR: k-mer length should be an odd number" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if (kmer_length < MIN_KMER_LENGTH) {
      std::cout << std::endl << "ERROR: k-mer length should be >= " << MIN_KMER_LENGTH << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if (kmer_length > (MAX_BOOST_ARRAY_SIZE * BPS_PER_BYTE)) {
      std::cout << std::endl << "ERROR: k-mer length should be < " << MAX_BOOST_ARRAY_SIZE * BPS_PER_BYTE << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else {
      kmer_middle_index = (std::size_t)(kmer_length / 2);
   }

   if (target_false_positive_prob == 0.0) {
      std::cout << std::endl << "ERROR: Target false positive probability not specified" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if (target_false_positive_prob < 0 || target_false_positive_prob > 1) {
      std::cout << std::endl << "ERROR: False positive rate should be between 0 and 1" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (set_kmer_occurrence_threshold == true) {
      if (kmer_occurrence_threshold < 2) {
         std::cout << std::endl << "ERROR: Bit vector max count should be >= 2" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      else if (kmer_occurrence_threshold > BIT_VECTOR_MAX_COUNT) {
         std::cout << std::endl << "ERROR: Bit-vector counter width cannot exceed " << BIT_VECTOR_MAX_COUNT << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   if (num_clusters < 1) {
      std::cout << std::endl << "ERROR: The number of clusters should be >= 1" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else {
      double log_out = log10(num_clusters);
      num_clusters_digit = (std::size_t)std::floor(log_out) + 1;
   }

   if (extend < 1) {
      std::cout << std::endl << "ERROR: The max extension should be >= 1" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (debug == true) {
      nowrite = false;
   }

   residue_kmer = kmer_length % BPS_PER_BYTE;
   remainder_kmer = kmer_length / BPS_PER_BYTE;
   num_bytes_per_kmer = kmer_length / BPS_PER_BYTE;
   if (residue_kmer != 0) {
      num_bytes_per_kmer++;
   }

   //--------------------------------------------------
   // output file names
   //--------------------------------------------------
   // set a log file name
   log_file_name = prefix + ".log";

   // set original read test result file names
   bf_querying_result_org_file_name  = prefix + ".org.bf-querying-result";
   bf_querying_result_org_file_name1 = prefix + ".org.1.bf-querying-result";
   bf_querying_result_org_file_name2 = prefix + ".org.2.bf-querying-result";

   // set corrected read test result file names
   bf_querying_result_mod_file_name  = prefix + ".mod.bf-querying-result";
   bf_querying_result_mod_file_name1 = prefix + ".mod.1.bf-querying-result";
   bf_querying_result_mod_file_name2 = prefix + ".mod.2.bf-querying-result";

   // set error_correction information file names
   error_correction_info_file_name  = prefix + ".error-correction";
   error_correction_info_file_name1 = prefix + ".1.error-correction";
   error_correction_info_file_name2 = prefix + ".2.error-correction";

   verified_error_correction_info_file_name  = prefix + ".verified.error-correction";
   verified_error_correction_info_file_name1 = prefix + ".1.verified.error-correction";
   verified_error_correction_info_file_name2 = prefix + ".2.verified.error-correction";

   // set error_correction information file names (text version)
   error_correction_info_txt_file_name  = prefix + ".error-correction.txt";
   error_correction_info_txt_file_name1 = prefix + ".1.error-correction.txt";
   error_correction_info_txt_file_name2 = prefix + ".2.error-correction.txt";

   verified_error_correction_info_txt_file_name  = prefix + ".verified.error-correction.txt";
   verified_error_correction_info_txt_file_name1 = prefix + ".1.verified.error-correction.txt";
   verified_error_correction_info_txt_file_name2 = prefix + ".2.verified.error-correction.txt";

   // set error_correction_direction_file_name
   error_correction_direction_file_name  = prefix + ".error-correction-direction";
   error_correction_direction_file_name1 = prefix + ".error-correction-direction1";
   error_correction_direction_file_name2 = prefix + ".error-correction-direction2";

   // set a tef file name
   tef_file_name = prefix + ".tef";

   // set a hisgram file name
   histo_file_name = prefix + ".histo";

   // set error_corrected read file names
   corrected_read_file_name  = prefix + ".corrected.fastq";
   corrected_read_file_name1 = prefix + ".1.corrected.fastq";
   corrected_read_file_name2 = prefix + ".2.corrected.fastq";

   //--------------------------------------------------
   // print options
   //--------------------------------------------------
   std::ofstream f_log;
   f_log.open(log_file_name.c_str());

   if (f_log.is_open()) {
      std::cout << std::endl;
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << "AUTHOR : Yun Heo" << std::endl;
      std::cout << "VERSION: " VERSION << std::endl;
      std::cout << "DATE   : " DATE << std::endl;
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << std::endl;
      std::cout << "Parsing arguments is finished" << std::endl;
      if (read_file_name.empty()) {
         std::cout << "     1st Read File Name            : " << read_file_name1 << std::endl;
         std::cout << "     2nd Read File Name            : " << read_file_name2 << std::endl;
      }
      else {
         std::cout << "     Read File Name                : " << read_file_name << std::endl;
      }
      std::cout << "     Log File Name                 : " << log_file_name << std::endl;
      std::cout << "     K-mer Length                  : " << kmer_length << std::endl;
      std::cout << "     Target False Positive Prob.   : " << target_false_positive_prob << std::endl;
      std::cout << "     Random Seed                   : " << random_seed << std::endl;
      if (set_kmer_occurrence_threshold == true) {
         std::cout << "     K-mer Occurence Threshold     : " << kmer_occurrence_threshold << std::endl;
      }
      else {
         std::cout << "     K-mer Occurence Threshold     : Not specified" << std::endl;
      }
      std::cout << "     Number of Clusters            : " << num_clusters << std::endl;
      std::cout << std::endl;

      f_log << "----------------------------------------------------------------------" << std::endl;
      f_log << "AUTHOR : Yun Heo" << std::endl;
      f_log << "VERSION: " VERSION << std::endl;
      f_log << "DATE   : " DATE << std::endl;
      f_log << "----------------------------------------------------------------------" << std::endl;
      f_log << std::endl;
      f_log << "Parsing arguments is finished" << std::endl;
      if (read_file_name.empty()) {
         f_log << "     1st Read File Name            : " << read_file_name1 << std::endl;
         f_log << "     2nd Read File Name            : " << read_file_name2 << std::endl;
      }
      else {
         f_log << "     Read File Name                : " << read_file_name << std::endl;
      }
      f_log << "     Log File Name                 : " << log_file_name << std::endl;
      f_log << "     K-mer Length                  : " << kmer_length << std::endl;
      f_log << "     Target False Positive Prob.   : " << target_false_positive_prob << std::endl;
      f_log << "     Random Seed                   : " << random_seed << std::endl;
      if (set_kmer_occurrence_threshold == true) {
         f_log << "     K-mer Occurence Threshold     : " << kmer_occurrence_threshold << std::endl;
      }
      else {
         f_log << "     K-mer Occurence Threshold     : Not specified" << std::endl;
      }
      f_log << "     Number of Clusters            : " << num_clusters << std::endl;
      f_log << std::endl;

      f_log.close();
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
}

//----------------------------------------------------------------------
// print_usage
//----------------------------------------------------------------------
void C_arg::print_usage() {
   std::cout << std::endl;
   std::cout << "--------------------------------------------------------------------------------" << std::endl;
   std::cout << "AUTHOR : Yun Heo" << std::endl;
   std::cout << "VERSION: " VERSION << std::endl;
   std::cout << "DATE   : " DATE << std::endl;
   std::cout << "--------------------------------------------------------------------------------" << std::endl;
   std::cout << std::endl;
   std::cout << "USAGE: " << args[0] << " <ARGUMENTS>" << std::endl;
   std::cout << std::endl;
   std::cout << "ARGUMENT               DESCRIPTION                           MANDATORY   DEFAULT" << std::endl;
   std::cout << "--------------------------------------------------------------------------------" << std::endl;
   std::cout << "-read    <file name>   input file (single-end reads)         Y                  " << std::endl;
   std::cout << "-read1   <file name>   1st input file (paired-end reads)     Y                  " << std::endl;
   std::cout << "-read2   <file name>   2nd input file (paired-end reads)     Y                  " << std::endl;
   std::cout << "-prefix     <prefix>   outfile prefix                        Y                  " << std::endl;
   std::cout << "-kmerlength <number>   k-mer length                          Y                  " << std::endl;
   std::cout << "-fpr        <number>   target false positive probability     N           0.001  " << std::endl;
   std::cout << "-seed       <number>   random number seed                    N           0      " << std::endl;
   std::cout << "-cluster    <number>   number of clusters                    N           100    " << std::endl;
   std::cout << "-count      <number>   k-mer occurrence threshold            N                  " << std::endl;
   std::cout << "-nowrite               no output read                        N                  " << std::endl;
   std::cout << "-noverify              no false positive removal             N                  " << std::endl;
   std::cout << "--------------------------------------------------------------------------------" << std::endl;
   std::cout << std::endl;
   std::cout << "EXAMPLE (paired-end): " << args[0] << " -read1 in1.fastq -read2 in2.fastq -prefix directory/prefix -kmerlength 31" << std::endl;
   std::cout << "EXAMPLE (single-end): " << args[0] << " -read in.fastq -prefix directory/prefix -kmerlength 31" << std::endl;
   std::cout << std::endl;
}
