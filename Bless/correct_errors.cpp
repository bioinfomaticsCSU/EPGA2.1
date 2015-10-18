#include "correct_errors.hpp"



//----------------------------------------------------------------------
// generate_bloom_filter
//----------------------------------------------------------------------
void C_correct_errors::generate_bloom_filter(const C_arg& c_inst_args, C_time& c_inst_time) {
   // open a log file
   f_log.open(c_inst_args.log_file_name.c_str(), std::fstream::app);

   if (f_log.is_open()) {
      std::cout << "Generating a Bloom filters" << std::endl;

      f_log     << "Generating a Bloom filters" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // find optimal parameters of for this Bloom filter
   c_inst_bloom_filter.find_optimal_parameters(num_unique_solid_kmers, c_inst_args.target_false_positive_prob);

   // generate unique hash functions
   c_inst_bloom_filter.generate_unique_hash_func(c_inst_args.random_seed);

   // Initialize Bloom filter
   c_inst_bloom_filter.generate_bit_vector();

   bloom_type bit_vector_size_mb(c_inst_bloom_filter.bit_vector_width_byte / (1024 * 1024));
   if (bit_vector_size_mb < 1) {
      bit_vector_size_mb = 1;
   }

   std::cout << "     Bloom filter" << std::endl;
   std::cout << "     Number of keys          : " << num_unique_solid_kmers << std::endl;
   std::cout << "     Bit-vector size         : " << bit_vector_size_mb << "MB" << std::endl;
   std::cout << "     Number of hash functions: " << c_inst_bloom_filter.num_hash_func << std::endl;
   std::cout << "     Generating a Bloom filter: done" << std::endl << std::endl;

   f_log     << "     Bloom filter" << std::endl;
   f_log     << "     Number of keys          : " << num_unique_solid_kmers << std::endl;
   f_log     << "     Bit-vector size         : " << bit_vector_size_mb << "MB" << std::endl;
   f_log     << "     Number of hash functions: " << c_inst_bloom_filter.num_hash_func << std::endl;
   f_log     << "     Generating a Bloom filter: done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// program_kmers
//----------------------------------------------------------------------
void C_correct_errors::program_kmers(const C_arg& c_inst_args, C_time& c_inst_time) {
   if (f_log.is_open()) {
      std::cout << "Programming k-mers to the Bloom filter" << std::endl;

      f_log     << "Programming k-mers to the Bloom filter" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_program_kmers_into_bf = asctime(localtime(&rawtime));

   program_kmers_fastq(c_inst_args);

   time(&rawtime);
   c_inst_time.end_program_kmers_into_bf = asctime(localtime(&rawtime));

   std::cout << "     Number of programmed k-mers: " << c_inst_bloom_filter.num_elements << std::endl;
   std::cout << "     Programming k-mers to the Bloom filter: done" << std::endl << std::endl;

   f_log     << "     Number of programmed k-mers: " << c_inst_bloom_filter.num_elements << std::endl;
   f_log     << "     Programming k-mers to the Bloom filter: done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// program_kmers_fastq
//----------------------------------------------------------------------
void C_correct_errors::program_kmers_fastq(const C_arg& c_inst_args) {
   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      program_kmers_cluster(it_files, c_inst_args.num_bytes_per_kmer);
   }
}



//----------------------------------------------------------------------
// program_kmers_cluster
//----------------------------------------------------------------------
void C_correct_errors::program_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_per_kmer) {
   std::ifstream f_unique_solid_kmers;
   f_unique_solid_kmers.open(distributed_unique_solid_kmer_file_names[file_index].c_str(), std::ios::binary);
   if (f_unique_solid_kmers.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   char buffer;
   f_unique_solid_kmers.get(buffer);

   while (!f_unique_solid_kmers.eof()) {
      std::string kmer;
      kmer.resize(num_bytes_per_kmer);

      // read a unique solid k-mer
      for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
         if (!f_unique_solid_kmers.good()) {
            std::cout << std::endl << "ERROR: Number of bytes in " << distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Number of bytes in " << distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
         else {
            kmer[it_byte] = buffer;

            f_unique_solid_kmers.get(buffer);
         }
      }

      // program the unique solid k-mer into the Bloom filter
      c_inst_bloom_filter.program(kmer);
   }

   f_unique_solid_kmers.close();
}



//----------------------------------------------------------------------
// reverse_complement
//----------------------------------------------------------------------
inline void C_correct_errors::reverse_complement(const std::string& read, std::string& read_rc) {
   read_rc = read;
   std::reverse(read_rc.begin(), read_rc.end());

   for (std::string::iterator it = read_rc.begin(); it != read_rc.end(); it++) {
      switch (*it) {
         case 'A' :
            *it = 'T';
            break;
         case 'C' :
            *it = 'G';
            break;
         case 'G' :
            *it = 'C';
            break;
         case 'T' :
            *it = 'A';
            break;
         default :
            std::cout << std::endl << "ERROR: Illegal character " << *it << " (reverse_complement)" << std::endl << std::endl;
            exit(EXIT_FAILURE);
            break;
      }
   }
}



//----------------------------------------------------------------------
// remove_new_line
//----------------------------------------------------------------------
std::string C_correct_errors::remove_new_line(std::string in_string) {
   in_string.erase(std::remove(in_string.begin(), in_string.end(), '\n'), in_string.end());
   return in_string;
}



//----------------------------------------------------------------------
// summarize_outputs
//----------------------------------------------------------------------
void C_correct_errors::summarize_outputs(const C_arg& c_inst_args, C_time& c_inst_time) {
   if (f_log.is_open()) {
      std::string last_ending_time;

      //--------------------------------------------------
      // stdout
      //--------------------------------------------------
      std::cout << "Running Time" << std::endl;
      std::cout << "     Parsing arguments" << std::endl;
      std::cout << "          Start:" << remove_new_line(c_inst_time.start_parse_args) << std::endl;
      std::cout << "          End  :" << remove_new_line(c_inst_time.end_parse_args) << std::endl;

      std::cout << "     Checking reads" << std::endl;
      std::cout << "          Start:" << remove_new_line(c_inst_time.end_parse_args) << std::endl;
      std::cout << "          End  :" << remove_new_line(c_inst_time.end_check_read_file) << std::endl;

      std::cout << "     Distributing k-mers" << std::endl;
      std::cout << "          Start:" << remove_new_line(c_inst_time.end_check_read_file) << std::endl;
      std::cout << "          End  :" << remove_new_line(c_inst_time.end_distribute_kmers) << std::endl;

      std::cout << "     Counting the number of unique solid k-mers" << std::endl;
      std::cout << "          Start:" << remove_new_line(c_inst_time.end_distribute_kmers) << std::endl;
      std::cout << "          End  :" << remove_new_line(c_inst_time.end_count_kmers) << std::endl;

      std::cout << "     Programming k-mers into the Bloom filter" << std::endl;
      std::cout << "          Start:" << remove_new_line(c_inst_time.end_count_kmers) << std::endl;
      std::cout << "          End  :" << remove_new_line(c_inst_time.end_program_kmers_into_bf) << std::endl;

      std::cout << "     Correcting errors in reads" << std::endl;
      std::cout << "          Start:" << remove_new_line(c_inst_time.end_program_kmers_into_bf) << std::endl;
      std::cout << "          End  :" << remove_new_line(c_inst_time.end_correct_errors_in_reads) << std::endl;
      last_ending_time = c_inst_time.end_correct_errors_in_reads;

      if (c_inst_args.debug == true) {
         std::cout << "     Writing querying results of original reads" << std::endl;
         std::cout << "          Start:" << remove_new_line(last_ending_time) << std::endl;
         std::cout << "          End  :" << remove_new_line(c_inst_time.end_mark_unique_solid_kmers_org) << std::endl;
         last_ending_time = c_inst_time.end_mark_unique_solid_kmers_org;
      }

      if (c_inst_args.verify == true) {
         std::cout << "     Finding false positve candidate k-mers" << std::endl;
         std::cout << "          Start:" << remove_new_line(last_ending_time) << std::endl;
         std::cout << "          End  :" << remove_new_line(c_inst_time.end_find_false_positives) << std::endl;

         std::cout << "     Counting false positve candidate k-mers" << std::endl;
         std::cout << "          Start:" << remove_new_line(c_inst_time.end_find_false_positives) << std::endl;
         std::cout << "          End  :" << remove_new_line(c_inst_time.end_count_false_positives) << std::endl;

         std::cout << "     Correcting false positve k-mers" << std::endl;
         std::cout << "          Start:" << remove_new_line(c_inst_time.end_count_false_positives) << std::endl;
         std::cout << "          End  :" << remove_new_line(c_inst_time.end_correct_false_positives) << std::endl;
         last_ending_time = c_inst_time.end_correct_false_positives;
      }

      if (c_inst_args.nowrite == false) {
         std::cout << "     Writing corrected reads" << std::endl;
         std::cout << "          Start:" << remove_new_line(last_ending_time) << std::endl;
         std::cout << "          End  :" << remove_new_line(c_inst_time.end_write_corrected_reads) << std::endl;
         last_ending_time = c_inst_time.end_write_corrected_reads;
      }

      std::cout << "     Writing a TEF" << std::endl;
      std::cout << "          Start:" << remove_new_line(last_ending_time) << std::endl;
      std::cout << "          End  :" << remove_new_line(c_inst_time.end_write_tef) << std::endl;

      std::cout << std::endl;
      std::cout << "The program is successfully completed" << std::endl << std::endl;

      //--------------------------------------------------
      // log file
      //--------------------------------------------------
      f_log << "Running Time" << std::endl;
      f_log << "     Parsing arguments" << std::endl;
      f_log << "          Start:" << remove_new_line(c_inst_time.start_parse_args) << std::endl;
      f_log << "          End  :" << remove_new_line(c_inst_time.end_parse_args) << std::endl;

      f_log << "     Checking reads" << std::endl;
      f_log << "          Start:" << remove_new_line(c_inst_time.end_parse_args) << std::endl;
      f_log << "          End  :" << remove_new_line(c_inst_time.end_check_read_file) << std::endl;

      f_log << "     Distributing k-mers" << std::endl;
      f_log << "          Start:" << remove_new_line(c_inst_time.end_check_read_file) << std::endl;
      f_log << "          End  :" << remove_new_line(c_inst_time.end_distribute_kmers) << std::endl;

      f_log << "     Counting the number of unique solid k-mers" << std::endl;
      f_log << "          Start:" << remove_new_line(c_inst_time.end_distribute_kmers) << std::endl;
      f_log << "          End  :" << remove_new_line(c_inst_time.end_count_kmers) << std::endl;

      f_log << "     Programming k-mers into the Bloom filter" << std::endl;
      f_log << "          Start:" << remove_new_line(c_inst_time.end_count_kmers) << std::endl;
      f_log << "          End  :" << remove_new_line(c_inst_time.end_program_kmers_into_bf) << std::endl;

      f_log << "     Correcting errors in reads" << std::endl;
      f_log << "          Start:" << remove_new_line(c_inst_time.end_program_kmers_into_bf) << std::endl;
      f_log << "          End  :" << remove_new_line(c_inst_time.end_correct_errors_in_reads) << std::endl;
      last_ending_time = c_inst_time.end_correct_errors_in_reads;

      if (c_inst_args.debug == true) {
         f_log << "     Writing querying results of original reads" << std::endl;
         f_log << "          Start:" << remove_new_line(last_ending_time) << std::endl;
         f_log << "          End  :" << remove_new_line(c_inst_time.end_mark_unique_solid_kmers_org) << std::endl;
         last_ending_time = c_inst_time.end_mark_unique_solid_kmers_org;
      }

      if (c_inst_args.verify == true) {
         f_log << "     Finding false positve candidate k-mers" << std::endl;
         f_log << "          Start:" << remove_new_line(last_ending_time) << std::endl;
         f_log << "          End  :" << remove_new_line(c_inst_time.end_find_false_positives) << std::endl;

         f_log << "     Counting false positve candidate k-mers" << std::endl;
         f_log << "          Start:" << remove_new_line(c_inst_time.end_find_false_positives) << std::endl;
         f_log << "          End  :" << remove_new_line(c_inst_time.end_count_false_positives) << std::endl;

         f_log << "     Correcting false positve k-mers" << std::endl;
         f_log << "          Start:" << remove_new_line(c_inst_time.end_count_false_positives) << std::endl;
         f_log << "          End  :" << remove_new_line(c_inst_time.end_correct_false_positives) << std::endl;
         last_ending_time = c_inst_time.end_correct_false_positives;
      }

      if (c_inst_args.nowrite == false) {
         f_log << "     Writing corrected reads" << std::endl;
         f_log << "          Start:" << remove_new_line(last_ending_time) << std::endl;
         f_log << "          End  :" << remove_new_line(c_inst_time.end_write_corrected_reads) << std::endl;
         last_ending_time = c_inst_time.end_write_corrected_reads;
      }

      f_log << "     Writing a TEF" << std::endl;
      f_log << "          Start:" << remove_new_line(last_ending_time) << std::endl;
      f_log << "          End  :" << remove_new_line(c_inst_time.end_write_tef) << std::endl;

      f_log << std::endl;
      f_log << "The program is successfully completed" << std::endl << std::endl;
      f_log.close();
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
}



//----------------------------------------------------------------------
// correct_errors_in_reads
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_correct_errors_in_reads = asctime(localtime(&rawtime));

   if (f_log.is_open()) {
      std::cout << "Correcting errors in reads" << std::endl;

      f_log     << "Correcting errors in reads" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (c_inst_args.paired_read == true) {
      correct_errors_in_reads_paired_fastq(c_inst_args);
   }
   else {
      correct_errors_in_reads_single_fastq(c_inst_args);
   }

   std::cout << "     Number of corrected errors: " << std::setw(10) << num_corrected_errors_step1_1 + num_corrected_errors_step1_2 + num_corrected_errors_step1_3 + num_corrected_errors_step2_1 + num_corrected_errors_step2_2 << std::endl;
   if (c_inst_args.debug == true) {
      std::cout << "     Number of corrected errors (step 1-1): " << std::setw(10) << num_corrected_errors_step1_1 << std::endl;
      std::cout << "     Number of corrected errors (step 1-2): " << std::setw(10) << num_corrected_errors_step1_2 << std::endl;
      std::cout << "     Number of corrected errors (step 1-3): " << std::setw(10) << num_corrected_errors_step1_3 << std::endl;
      std::cout << "     Number of corrected errors (step 2-1): " << std::setw(10) << num_corrected_errors_step2_1 << std::endl;
      std::cout << "     Number of corrected errors (step 2-2): " << std::setw(10) << num_corrected_errors_step2_2 << std::endl;
   }
   std::cout << "     Number of corrected reads : " << std::setw(10) << num_corrected_reads << std::endl;
   std::cout << "     Correcting errors in reads: done" << std::endl << std::endl;

   f_log     << "     Number of corrected errors: " << std::setw(10) << num_corrected_errors_step1_1 + num_corrected_errors_step1_2 + num_corrected_errors_step1_3 + num_corrected_errors_step2_1 + num_corrected_errors_step2_2 << std::endl;
   if (c_inst_args.debug == true) {
      f_log     << "     Number of corrected errors (step 1-1): " << std::setw(10) << num_corrected_errors_step1_1 << std::endl;
      f_log     << "     Number of corrected errors (step 1-2): " << std::setw(10) << num_corrected_errors_step1_2 << std::endl;
      f_log     << "     Number of corrected errors (step 1-3): " << std::setw(10) << num_corrected_errors_step1_3 << std::endl;
      f_log     << "     Number of corrected errors (step 2-1): " << std::setw(10) << num_corrected_errors_step2_1 << std::endl;
      f_log     << "     Number of corrected errors (step 2-2): " << std::setw(10) << num_corrected_errors_step2_2 << std::endl;
   }
   f_log     << "     Number of corrected reads : " << std::setw(10) << num_corrected_reads << std::endl;
   f_log     << "     Correcting errors in reads: done" << std::endl << std::endl;

   time(&rawtime);
   c_inst_time.end_correct_errors_in_reads = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// correct_errors_in_reads_single_fastq
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads_single_fastq(const C_arg& c_inst_args) {
   // open error correction information files
   std::ofstream f_error_correction;
   f_error_correction.open(c_inst_args.error_correction_info_file_name.c_str(), std::ios::binary);

   // check error correction information files
   if (f_error_correction.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open error correction information files (text file)
   std::ofstream f_error_correction_txt;

   if (c_inst_args.debug == true) {
      f_error_correction_txt.open(c_inst_args.error_correction_info_txt_file_name.c_str());

      if (f_error_correction_txt.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_txt_file_name << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_txt_file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open the error correction direction file
   std::ofstream f_error_correction_direction;

   if (c_inst_args.verify == true) {
      f_error_correction_direction.open(c_inst_args.error_correction_direction_file_name.c_str());
      if (f_error_correction_direction.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open read files
   f_read.open(c_inst_args.read_file_name.c_str());

   // iterate reads
   if (f_read.is_open()) {
      std::string header;

      std::string sequence;
      std::string sequence_tmp;

      std::string connector;

      std::string quality_score;

      // header
      getline(f_read, header);

      if (c_inst_args.debug == true) {
         f_error_correction_txt << header << std::endl;
      }

      //--------------------------------------------------
      // read_length needs to be calculate all the time for reads with variable length
      //--------------------------------------------------
      // number of bytes for storing reads
      std::size_t num_byte_per_read((std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE)));

      // write remaining bps into the vectors
      std::size_t residue_read(read_length % BPS_PER_BYTE);

      // initialize variables
      num_corrected_errors_step1_1 = 0;
      num_corrected_errors_step1_2 = 0;
      num_corrected_errors_step1_3 = 0;
      num_corrected_errors_step2_1 = 0;
      num_corrected_errors_step2_2 = 0;
      num_corrected_reads = 0;

      std::size_t prev_num_corrected_errors_step1_1(0);
      std::size_t prev_num_corrected_errors_step1_2(0);
      std::size_t prev_num_corrected_errors_step1_3(0);
      std::size_t prev_num_corrected_errors_step2_1(0);
      std::size_t prev_num_corrected_errors_step2_2(0);

      std::size_t sum_errors(0);
      std::size_t prev_sum_errors(0);

      while(!f_read.eof()) {
         // DNA sequence
         getline(f_read, sequence);

         // "+"
         getline(f_read, connector);

         // quality score
         getline(f_read, quality_score);

         // change sequences to upper case
         transform(sequence.begin(), sequence.end(), sequence.begin(), toupper);

         // substitute Ns other characters
         std::replace(sequence.begin(), sequence.end(), 'N', SUBST_CHAR);

         // storage for the modification of the reads
         // # of entries = read length
         std::string sequence_modification(read_length, '0');

         std::string sequence_modification_zero(read_length, '0');

         //----------------------------------------------------------------------
         // correct errors in a read
         //----------------------------------------------------------------------
         bool too_many_errors(false);

         correct_errors_in_a_read_fastq(
                                        sequence,
                                        sequence_modification,
                                        quality_score,
                                        c_inst_args.kmer_length,
                                        c_inst_args.kmer_middle_index,
                                        c_inst_args.num_bytes_per_kmer,
                                        c_inst_args.residue_kmer,
                                        c_inst_args.extend,
                                        c_inst_args.verify
                                       );

         // update num_corrected_reads
         too_many_errors = false;
         sum_errors = num_corrected_errors_step1_1 + num_corrected_errors_step1_2 + num_corrected_errors_step1_3 + num_corrected_errors_step2_1 + num_corrected_errors_step2_2;
         if ((sum_errors - prev_sum_errors) > (sequence.length() * MAX_ERROR_RATE)) {
            num_corrected_errors_step1_1 = prev_num_corrected_errors_step1_1;
            num_corrected_errors_step1_2 = prev_num_corrected_errors_step1_2;
            num_corrected_errors_step1_3 = prev_num_corrected_errors_step1_3;
            num_corrected_errors_step2_1 = prev_num_corrected_errors_step2_1;
            num_corrected_errors_step2_2 = prev_num_corrected_errors_step2_2;

            sum_errors = prev_sum_errors;

            too_many_errors = true;
         }
         else if (sum_errors > prev_sum_errors) {
            prev_num_corrected_errors_step1_1 = num_corrected_errors_step1_1;
            prev_num_corrected_errors_step1_2 = num_corrected_errors_step1_2;
            prev_num_corrected_errors_step1_3 = num_corrected_errors_step1_3;
            prev_num_corrected_errors_step2_1 = num_corrected_errors_step2_1;
            prev_num_corrected_errors_step2_2 = num_corrected_errors_step2_2;

            prev_sum_errors = sum_errors;

            num_corrected_reads++;
         }

         //----------------------------------------------------------------------
         // write error correction information to files
         //----------------------------------------------------------------------
         std::vector<char> write_buffer;

         char buffer;

         // initialize buffers
         buffer &= ZERO;

         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            if (too_many_errors == true) {
               encode_correction_info(buffer, sequence[it_mod], '0');
            }
            else {
               if (c_inst_args.verify == true) {
                  // the modified base is in lower case
                  if (islower(sequence_modification[it_mod])) {
                     f_error_correction_direction.put('0');
                     sequence_modification[it_mod] = toupper(sequence_modification[it_mod]);
                  }
                  // the modified base is in upper case
                  else if (isupper(sequence_modification[it_mod])) {
                     f_error_correction_direction.put('1');
                  }
               }

               encode_correction_info(buffer, sequence[it_mod], sequence_modification[it_mod]);
            }

            if ((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) {
               // update write_buffer*
               write_buffer.push_back(buffer);

               // initialize buffer
               buffer &= ZERO;
            }
         }

         // read_length is a multiple of BPS_PER_BYTE
         // do nothing
         if (residue_read == 0) {
         }
         // read_length is not a multiple of BPS_PER_BYTE
         else {
            for (std::size_t it_fill = 0; it_fill < (BPS_PER_BYTE - residue_read); it_fill++) {
               buffer = buffer << 2;
            }

            // fill the last entry
            write_buffer.push_back(buffer);
         }

         // write vectors to output files
         f_error_correction.write((const char*)&write_buffer[0], num_byte_per_read);

         // initialize vectors
         write_buffer.clear();

         if (c_inst_args.debug == true) {
            if (too_many_errors == true) {
               f_error_correction_txt << sequence_modification_zero << std::endl;
            }
            else {
               f_error_correction_txt << sequence_modification << std::endl;
            }
         }

         // header
         getline(f_read, header);

         if ((c_inst_args.debug == true) && (header.length() > 0)) {
            f_error_correction_txt << header << std::endl;
         }
      }
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // close read files
   f_read.close();

   // close error correction information files
   f_error_correction.close();

   if (c_inst_args.debug == true) {
      f_error_correction_txt.close();
   }

   // close the error correction direction file
   if (c_inst_args.verify == true) {
      f_error_correction_direction.close();
   }
}



//----------------------------------------------------------------------
// correct_errors_in_reads_paired_fastq
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads_paired_fastq(const C_arg& c_inst_args) {
   // open error correction information files
   std::ofstream f_error_correction1;
   std::ofstream f_error_correction2;
   f_error_correction1.open(c_inst_args.error_correction_info_file_name1.c_str(), std::ios::binary);
   f_error_correction2.open(c_inst_args.error_correction_info_file_name2.c_str(), std::ios::binary);

   // check error correction information files
   if (f_error_correction1.is_open()) {
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_error_correction2.is_open()) {
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name2 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open error correction information files (text file)
   std::ofstream f_error_correction_txt1;
   std::ofstream f_error_correction_txt2;

   if (c_inst_args.debug == true) {
      f_error_correction_txt1.open(c_inst_args.error_correction_info_txt_file_name1.c_str());
      f_error_correction_txt2.open(c_inst_args.error_correction_info_txt_file_name2.c_str());

      if (f_error_correction_txt1.is_open()) {
      } else {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_txt_file_name1 << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_txt_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      if (f_error_correction_txt2.is_open()) {
      } else {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_txt_file_name2 << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_txt_file_name2 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open the error correction direction file
   std::ofstream f_error_correction_direction1;
   std::ofstream f_error_correction_direction2;

   if (c_inst_args.verify == true) {
      f_error_correction_direction1.open(c_inst_args.error_correction_direction_file_name1.c_str());
      if (f_error_correction_direction1.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name1 << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      f_error_correction_direction2.open(c_inst_args.error_correction_direction_file_name2.c_str());
      if (f_error_correction_direction2.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name2 << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name2 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open read files
   f_read_1.open(c_inst_args.read_file_name1.c_str());
   f_read_2.open(c_inst_args.read_file_name2.c_str());

   // iterate reads
   if (f_read_1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if (f_read_2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else {
      std::string header_1st;
      std::string header_2nd;

      std::string sequence_1st;
      std::string sequence_2nd;
      std::string sequence_1st_tmp;
      std::string sequence_2nd_tmp;

      std::string connector_1st;
      std::string connector_2nd;

      std::string quality_score_1st;
      std::string quality_score_2nd;

      // header
      getline(f_read_1, header_1st);
      getline(f_read_2, header_2nd);

      if (c_inst_args.debug == true) {
         f_error_correction_txt1 << header_1st << std::endl;
         f_error_correction_txt2 << header_2nd << std::endl;
      }

      //--------------------------------------------------
      // read_length needs to be calculate all the time for reads with variable length
      //--------------------------------------------------
      // number of bytes for storing reads
      std::size_t num_byte_per_read((std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE)));

      // write remaining bps into the vectors
      std::size_t residue_read(read_length % BPS_PER_BYTE);

      // initialize variables
      num_corrected_errors_step1_1 = 0;
      num_corrected_errors_step1_2 = 0;
      num_corrected_errors_step1_3 = 0;
      num_corrected_errors_step2_1 = 0;
      num_corrected_errors_step2_2 = 0;
      num_corrected_reads = 0;

      std::size_t prev_num_corrected_errors_step1_1(0);
      std::size_t prev_num_corrected_errors_step1_2(0);
      std::size_t prev_num_corrected_errors_step1_3(0);
      std::size_t prev_num_corrected_errors_step2_1(0);
      std::size_t prev_num_corrected_errors_step2_2(0);

      std::size_t sum_errors(0);
      std::size_t prev_sum_errors(0);

      while(!f_read_1.eof()) {
         // DNA sequence
         getline(f_read_1, sequence_1st);
         getline(f_read_2, sequence_2nd);

         // "+"
         getline(f_read_1, connector_1st);
         getline(f_read_2, connector_2nd);

         // quality score
         getline(f_read_1, quality_score_1st);
         getline(f_read_2, quality_score_2nd);

         // change sequences to upper case
         transform(sequence_1st.begin(), sequence_1st.end(), sequence_1st.begin(), toupper);
         transform(sequence_2nd.begin(), sequence_2nd.end(), sequence_2nd.begin(), toupper);

         // substitute Ns other characters
         std::replace(sequence_1st.begin(), sequence_1st.end(), 'N', SUBST_CHAR);
         std::replace(sequence_2nd.begin(), sequence_2nd.end(), 'N', SUBST_CHAR);

         // temporary storage of reads
         // sequence_1st_tmp = sequence_1st;
         // sequence_2nd_tmp = sequence_2nd;

         // storage for the modification of the reads
         // # of entries = read length
         std::string sequence_1st_modification(read_length, '0');
         std::string sequence_2nd_modification(read_length, '0');

         std::string sequence_modification_zero(read_length, '0');

         // # of entry = (# of k-mers) - (padding in both direction) = l - k + 1
         // std::vector<std::size_t> num_covered_1st(num_kmers, 0);
         // std::vector<std::size_t> num_covered_2nd(num_kmers, 0);

         // k-mer modification flag vector
         // # of entry = # of k-mers = l - k + 1
         // std::vector<bool> flag_modified_1st(num_kmers, false);
         // std::vector<bool> flag_modified_2nd(num_kmers, false);

         // bool modification_made_1st(false);
         // bool modification_made_2nd(false);

         // std::size_t accum_errors_1st(0);
         // std::size_t accum_errors_2nd(0);

         //----------------------------------------------------------------------
         // correct errors in a read
         //----------------------------------------------------------------------
         bool too_many_errors_1st(false);
         bool too_many_errors_2nd(false);

         correct_errors_in_a_read_fastq(
                                        sequence_1st,
                                        sequence_1st_modification,
                                        quality_score_1st,
                                        c_inst_args.kmer_length,
                                        c_inst_args.kmer_middle_index,
                                        c_inst_args.num_bytes_per_kmer,
                                        c_inst_args.residue_kmer,
                                        c_inst_args.extend,
                                        c_inst_args.verify
                                       );

         // update num_corrected_reads
         too_many_errors_1st = false;
         sum_errors = num_corrected_errors_step1_1 + num_corrected_errors_step1_2 + num_corrected_errors_step1_3 + num_corrected_errors_step2_1 + num_corrected_errors_step2_2;
         if ((sum_errors - prev_sum_errors) > (sequence_1st.length() * MAX_ERROR_RATE)) {
            num_corrected_errors_step1_1 = prev_num_corrected_errors_step1_1;
            num_corrected_errors_step1_2 = prev_num_corrected_errors_step1_2;
            num_corrected_errors_step1_3 = prev_num_corrected_errors_step1_3;
            num_corrected_errors_step2_1 = prev_num_corrected_errors_step2_1;
            num_corrected_errors_step2_2 = prev_num_corrected_errors_step2_2;

            sum_errors = prev_sum_errors;

            too_many_errors_1st = true;
         }
         else if (sum_errors > prev_sum_errors) {
            prev_num_corrected_errors_step1_1 = num_corrected_errors_step1_1;
            prev_num_corrected_errors_step1_2 = num_corrected_errors_step1_2;
            prev_num_corrected_errors_step1_3 = num_corrected_errors_step1_3;
            prev_num_corrected_errors_step2_1 = num_corrected_errors_step2_1;
            prev_num_corrected_errors_step2_2 = num_corrected_errors_step2_2;

            prev_sum_errors = sum_errors;

            num_corrected_reads++;
         }

         correct_errors_in_a_read_fastq(
                                        sequence_2nd,
                                        sequence_2nd_modification,
                                        quality_score_2nd,
                                        c_inst_args.kmer_length,
                                        c_inst_args.kmer_middle_index,
                                        c_inst_args.num_bytes_per_kmer,
                                        c_inst_args.residue_kmer,
                                        c_inst_args.extend,
                                        c_inst_args.verify
                                       );

         // update num_corrected_reads
         too_many_errors_2nd = false;
         sum_errors = num_corrected_errors_step1_1 + num_corrected_errors_step1_2 + num_corrected_errors_step1_3 + num_corrected_errors_step2_1 + num_corrected_errors_step2_2;
         if ((sum_errors - prev_sum_errors) > (sequence_2nd.length() * MAX_ERROR_RATE)) {
            num_corrected_errors_step1_1 = prev_num_corrected_errors_step1_1;
            num_corrected_errors_step1_2 = prev_num_corrected_errors_step1_2;
            num_corrected_errors_step1_3 = prev_num_corrected_errors_step1_3;
            num_corrected_errors_step2_1 = prev_num_corrected_errors_step2_1;
            num_corrected_errors_step2_2 = prev_num_corrected_errors_step2_2;

            sum_errors = prev_sum_errors;

            too_many_errors_2nd = true;
         }
         else if (sum_errors > prev_sum_errors) {
            prev_num_corrected_errors_step1_1 = num_corrected_errors_step1_1;
            prev_num_corrected_errors_step1_2 = num_corrected_errors_step1_2;
            prev_num_corrected_errors_step1_3 = num_corrected_errors_step1_3;
            prev_num_corrected_errors_step2_1 = num_corrected_errors_step2_1;
            prev_num_corrected_errors_step2_2 = num_corrected_errors_step2_2;

            prev_sum_errors = sum_errors;

            num_corrected_reads++;
         }

         //----------------------------------------------------------------------
         // write error correction information to files
         //----------------------------------------------------------------------
         std::vector<char> write_buffer1;
         std::vector<char> write_buffer2;

         char buffer1;
         char buffer2;

         // initialize buffers
         buffer1 &= ZERO;
         buffer2 &= ZERO;

         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            // forward read
            if (too_many_errors_1st == true) {
               encode_correction_info(buffer1, sequence_1st[it_mod], '0');
            }
            else {
               if (c_inst_args.verify == true) {
                  // the modified base is in lower case
                  if (islower(sequence_1st_modification[it_mod])) {
                     f_error_correction_direction1.put('0');
                     sequence_1st_modification[it_mod] = toupper(sequence_1st_modification[it_mod]);
                  }
                  // the modified base is in upper case
                  else if (isupper(sequence_1st_modification[it_mod])) {
                     f_error_correction_direction1.put('1');
                  }
               }

               encode_correction_info(buffer1, sequence_1st[it_mod], sequence_1st_modification[it_mod]);
            }

            // reverse read
            if (too_many_errors_2nd == true) {
               encode_correction_info(buffer2, sequence_2nd[it_mod], '0');
            }
            else {
               if (c_inst_args.verify == true) {
                  // the modified base is in lower case
                  if (islower(sequence_2nd_modification[it_mod])) {
                     f_error_correction_direction2.put('0');
                     sequence_2nd_modification[it_mod] = toupper(sequence_2nd_modification[it_mod]);
                  }
                  // the modified base is in upper case
                  else if (isupper(sequence_2nd_modification[it_mod])) {
                     f_error_correction_direction2.put('1');
                  }
               }

               encode_correction_info(buffer2, sequence_2nd[it_mod], sequence_2nd_modification[it_mod]);
            }

            if ((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) {
               // update write_buffer*
               write_buffer1.push_back(buffer1);
               write_buffer2.push_back(buffer2);

               // initialize buffer
               buffer1 &= ZERO;
               buffer2 &= ZERO;
            }
         }

         // read_length is a multiple of BPS_PER_BYTE
         // do nothing
         if (residue_read == 0) {
         }
         // read_length is not a multiple of BPS_PER_BYTE
         else {
            for (std::size_t it_fill = 0; it_fill < (BPS_PER_BYTE - residue_read); it_fill++) {
               buffer1 = buffer1 << 2;
               buffer2 = buffer2 << 2;
            }

            // fill the last entry
            write_buffer1.push_back(buffer1);
            write_buffer2.push_back(buffer2);
         }

         // write vectors to output files
         f_error_correction1.write((const char*)&write_buffer1[0], num_byte_per_read);
         f_error_correction2.write((const char*)&write_buffer2[0], num_byte_per_read);

         // initialize vectors
         write_buffer1.clear();
         write_buffer2.clear();

         if (c_inst_args.debug == true) {
            if (too_many_errors_1st == true) {
               f_error_correction_txt1 << sequence_modification_zero << std::endl;
            }
            else {
               f_error_correction_txt1 << sequence_1st_modification << std::endl;
            }

            if (too_many_errors_2nd == true) {
               f_error_correction_txt2 << sequence_modification_zero << std::endl;
            }
            else {
               f_error_correction_txt2 << sequence_2nd_modification << std::endl;
            }
         }

         // header
         getline(f_read_1, header_1st);
         getline(f_read_2, header_2nd);

         if ((c_inst_args.debug == true) && (header_1st.length() > 0)) {
            f_error_correction_txt1 << header_1st << std::endl;
            f_error_correction_txt2 << header_2nd << std::endl;
         }
      }
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // close error correction information files
   f_error_correction1.close();
   f_error_correction2.close();

   if (c_inst_args.debug == true) {
      f_error_correction_txt1.close();
      f_error_correction_txt2.close();
   }

   // close the error correction direction file
   if (c_inst_args.verify == true) {
      f_error_correction_direction1.close();
      f_error_correction_direction2.close();
   }
}



//----------------------------------------------------------------------
// correct_errors_in_a_read_fastq
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_in_a_read_fastq(const std::string& sequence, std::string& sequence_modification, const std::string& quality_score, const std::size_t kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& max_extension, const bool& verify) {
   //--------------------------------------------------
   // STEP 0-0: find solid k-mers in this read
   //--------------------------------------------------
   // variables
   std::size_t num_kmers(sequence.length() - kmer_length + 1);

   std::vector< std::pair<std::size_t, std::size_t> > solid_regions;

   bool is_solid_kmer_prev(false);

   // find solid regions
   std::pair<std::size_t, std::size_t> new_solid_region;
   for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
      std::string current_kmer(sequence.substr(it_kmer, kmer_length));

      // k-mer is solid
      if (query_text(current_kmer, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         // start point of a solid region
         if (is_solid_kmer_prev == false) {
            new_solid_region.first = it_kmer;
            is_solid_kmer_prev = true;
         }
      }
      else {
         // end point of a solid region
         if (is_solid_kmer_prev == true) {
            new_solid_region.second = it_kmer - 1;

            if (new_solid_region.second < new_solid_region.first) {
               std::cout << std::endl << "ERROR: The second index is smaller than the first" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: The second index is smaller than the first" << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }
            solid_regions.push_back(new_solid_region);

            is_solid_kmer_prev = false;
         }
      }
   }

   // last solid region
   if (is_solid_kmer_prev == true) {
      new_solid_region.second = num_kmers - 1;
      solid_regions.push_back(new_solid_region);
   }

   //--------------------------------------------------
   // STEP 0-1: remove short solid regions beside short non-solid regions
   //--------------------------------------------------
   /*
   if (solid_regions.size() > 0) {
      std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
      solid_regions_tmp.push_back(solid_regions[0]);

      if (solid_regions.size() > 1) {
         for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
            // short non-solid region to the left
            if ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < kmer_length) {
               if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) < MIN_SOLID_LENGTH) {
                  // (non-solid region length < kmer_length) &&
                  // (solid_regions[it_region] is too short) &&
                  // (solid_regions[it_region].second is not the last base of the read)
                  // -> remove solid_regions[it_region + 1]
                  // k = 10
                  // SSSSSNNNNNNNSNNSSSSS
                  //             ^ Remove it
                  // do nothing
                  if (solid_regions[it_region].second != (sequence.length() - 1)) {
                  }
                  else {
                     solid_regions_tmp.push_back(solid_regions[it_region]);
                  }
               }
               else {
                  solid_regions_tmp.push_back(solid_regions[it_region]);
               }
            }
            else {
               solid_regions_tmp.push_back(solid_regions[it_region]);
            }
         }
      }
      solid_regions = solid_regions_tmp;
   }
   */

   //--------------------------------------------------
   // STEP 0-2: remove short solid regions
   //--------------------------------------------------
   if (solid_regions.size() > 0) {
      std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;

      for (std::size_t it_region = 0; it_region < solid_regions.size(); it_region++) {
         if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) >= MIN_SOLID_LENGTH) {
            solid_regions_tmp.push_back(solid_regions[it_region]);
         }
      }

      solid_regions = solid_regions_tmp;
   }

   //--------------------------------------------------
   // STEP 0-3: remove short non-solid regions
   //--------------------------------------------------
   if (solid_regions.size() > 0) {
      std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
      solid_regions_tmp.push_back(solid_regions[0]);

      if (solid_regions.size() > 1) {
         for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
            if ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < MIN_NON_SOLID_LENGTH) {
               solid_regions_tmp[solid_regions_tmp.size() - 1].second = solid_regions[it_region].second;
            }
            else {
               solid_regions_tmp.push_back(solid_regions[it_region]);
            }
         }
      }
      solid_regions = solid_regions_tmp;
   }

   //--------------------------------------------------
   // STEP 0-4: reduce the size of solid regions
   //--------------------------------------------------
   if (solid_regions.size() > 1) {
      for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
         // (length of a non-solid region < kmer_length) && (length of a non-solid region >= kmer_length - FP_SUSPECT_LENGTH(default: 1))
         if (((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < kmer_length) &&
             ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) >= kmer_length - FP_SUSPECT_LENGTH)) {
            // length of the right solid region > FP_SUSPECT_LENGTH(default: 1)
            if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) > FP_SUSPECT_LENGTH) {
               solid_regions[it_region].first += FP_SUSPECT_LENGTH;
            }

            // length of the left solid region > FP_SUSPECT_LENGTH(default: 1)
            if ((solid_regions[it_region - 1].second - solid_regions[it_region - 1].first + 1) > FP_SUSPECT_LENGTH) {
               solid_regions[it_region - 1].second -= FP_SUSPECT_LENGTH;
            }
         }
      }
   }

   //--------------------------------------------------
   // STEP 0-5: remove a solid region that makes a non-solid reiong shorter than k
   //--------------------------------------------------
   if (solid_regions.size() == 2) {
      // the first solid region starts from the first k-mer
      if (solid_regions[0].first == 0) {
         // the distance between two regions is shorter than k
         if ((solid_regions[1].first - solid_regions[0].second) < (kmer_length + 1)) {
            // remove the second solid region
            solid_regions.erase(solid_regions.begin() + 1);
         }
      }
      // the second solid region ends in the last k-mer
      else if (solid_regions[1].second == (sequence.length() - kmer_length)) {
         // the distance between two regions is shorter than k
         if ((solid_regions[1].first - solid_regions[0].second) < (kmer_length + 1)) {
            // the length of the second solid region is >= 10% of the sequence length
            if ((solid_regions[1].second - solid_regions[1].first + 1) >= (sequence.length() * 0.1)) {
               // the length of the first solid region is < 10% of the sequence length
               if ((solid_regions[0].second - solid_regions[0].first + 1) < (sequence.length() * 0.1)) {
                  // remove the second solid region
                  solid_regions.erase(solid_regions.begin());
               }
            }
         }
      }
   }

   //--------------------------------------------------
   // STEP 0-6: check the quality scores of right side of each solid k-mer region
   //--------------------------------------------------
   // at least one solid region
   if (solid_regions.size() > 0) {
      // 1 - (n - 1) solid region
      for (std::size_t it_sr = 0; it_sr < (solid_regions.size() - 1); it_sr++) {
         // sufficient solid regions length
         if ((solid_regions[it_sr].second - solid_regions[it_sr].first) > SOLID_REGION_ADJUST_RANGE) {
            for (std::size_t it_adjust = solid_regions[it_sr].second; it_adjust > (solid_regions[it_sr].second - SOLID_REGION_ADJUST_RANGE); it_adjust--) {
               // low quality score
               if ((((std::size_t)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < QS_CUTOFF) ||
                   (((std::size_t)quality_score[it_adjust] - quality_score_offset) < QS_CUTOFF)
                  ) {
                  solid_regions[it_sr].second = it_adjust - 1;
                  break;
               }
            }
         }
      }

      // last solid region
      std::size_t index_solid_region(solid_regions.size() - 1);

      // non-solid k-mers exist at the 3-prime end
      if (solid_regions[index_solid_region].second < (sequence.length() - kmer_length)) {
         // sufficient solid regions length
         if ((solid_regions[index_solid_region].second - solid_regions[index_solid_region].first) > SOLID_REGION_ADJUST_RANGE) {
            for (std::size_t it_adjust = solid_regions[index_solid_region].second; it_adjust > (solid_regions[index_solid_region].second - SOLID_REGION_ADJUST_RANGE); it_adjust--) {
               // low quality score
               if (((std::size_t)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < QS_CUTOFF) {
                  solid_regions[index_solid_region].second = it_adjust - 1;
                  break;
               }
            }
         }
      }

      // non-solid k-mers exist at the 5-prime end
      if (solid_regions[0].first > 0) {
         // sufficient solid regions length
         if ((solid_regions[0].second - solid_regions[0].first) > SOLID_REGION_ADJUST_RANGE) {
            for (std::size_t it_adjust = solid_regions[0].first; it_adjust < (solid_regions[0].first + SOLID_REGION_ADJUST_RANGE); it_adjust++) {
               // low quality score
               if (((std::size_t)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < QS_CUTOFF) {
                  solid_regions[0].first = it_adjust + 1;
                  break;
               }
            }
         }
      }
   }

   //--------------------------------------------------
   // STEP 0-7: check whether a non-solid region < k still exists
   //--------------------------------------------------
   bool short_non_solid_region(false);
   if (solid_regions.size() > 1) {
      for (std::size_t it_sr = 1; it_sr < (solid_regions.size() - 1); it_sr++) {
         if ((solid_regions[it_sr].first - solid_regions[it_sr - 1].second) <= kmer_length) {
            short_non_solid_region = true;
            break;
         }
      }
   }

   //--------------------------------------------------
   // correct errors
   //--------------------------------------------------
   std::string sequence_modified(sequence);

   if ((solid_regions.size() > 0) && (short_non_solid_region == false)) {
      //--------------------------------------------------
      // STEP 1-1: Correct errors between solid regions
      //--------------------------------------------------
      if (solid_regions.size() > 1) {
         // for each solid region
         for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
            if ((((solid_regions[it_region].first - 1) - (solid_regions[it_region - 1].second + 1)) + 1) >= kmer_length) {
               correct_errors_between_solid_regions(
                                                    sequence_modified,
                                                    quality_score,
                                                    (solid_regions[it_region - 1].second + 1),
                                                    (solid_regions[it_region].first - 1),
                                                    kmer_length,
                                                    kmer_middle_index,
                                                    num_bytes_per_kmer,
                                                    residue_kmer,
                                                    sequence_modification
                                                   );
            }
            else {
            }
         }
      }

      //--------------------------------------------------
      // STEP 1-2: Correct errors in the 5' end
      //--------------------------------------------------
      // number of solid regions is >= 1
      if (solid_regions.size() >= 1) {
         // the first solid region does not start from the 0-th k-mer in a read
         if (solid_regions[0].first > 0) {
            correct_errors_5_prime_end(
                                       sequence_modified,
                                       quality_score,
                                       solid_regions[0].first - 1,
                                       kmer_length,
                                       kmer_middle_index,
                                       num_bytes_per_kmer,
                                       residue_kmer,
                                       sequence_modification,
                                       max_extension,
                                       verify
                                      );
         }
      }

      //--------------------------------------------------
      // STEP 1-3: Correct errors in the 3' end
      //--------------------------------------------------
      // number of solid regions is >= 1
      if (solid_regions.size() >= 1) {
         // the last solid region does not end in the last k-mer in a read
         if (solid_regions[solid_regions.size() - 1].second < (sequence.length() - kmer_length)) {
            correct_errors_3_prime_end(
                                       sequence_modified,
                                       quality_score,
                                       solid_regions[solid_regions.size() - 1].second + 1,
                                       kmer_length,
                                       kmer_middle_index,
                                       num_bytes_per_kmer,
                                       residue_kmer,
                                       sequence_modification,
                                       max_extension
                                      );
         }
      }
   }
   //--------------------------------------------------
   // no solid region or short weak regions
   //--------------------------------------------------
   else {
      //--------------------------------------------------
      // STEP 2-1: Correct errors in the first k-mer
      //--------------------------------------------------
      // find potentially wrong bases
      std::vector<C_candidate_path> candidate_path_vector_tmp;

      correct_errors_first_kmer(
                                sequence_modified,
                                quality_score,
                                kmer_length,
                                kmer_middle_index,
                                num_bytes_per_kmer,
                                residue_kmer,
                                sequence_modification,
                                candidate_path_vector_tmp
                               );

      // candidiate_path_vector_tmp: differently modified versions of the first k-mer

      // filter some candidates by extending the first k-mer to the left
      std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

      if (candidate_path_vector_tmp.size() > 0) {
         // each path
         for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size(); it_candidates++) {
            // no modified path
            if (candidate_path_vector_tmp[it_candidates].modified_bases.size() == 0) {
               candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
            }
            // check the index of the first modified base
            // extension is needed
            //else if (candidate_path_vector_tmp[it_candidates].modified_bases[0].first < (MAX_EXTENSION - 1)) {
            else if (candidate_path_vector_tmp[it_candidates].modified_bases[0].first < (kmer_length - 1)) {
               bool extension_success(false);
               solid_first_kmer(
                                candidate_path_vector_tmp[it_candidates],
                                sequence_modified,
                                kmer_length,
                                kmer_middle_index,
                                num_bytes_per_kmer,
                                residue_kmer,
                                max_extension,
                                extension_success
                               );

               if (extension_success == true) {
                  candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
               }
            }
            // extension is not needed
            else {
               candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
            }
         }
      }

      // candidiate_path_vector_tmp_tmp: solid k-mers in candidate_path_vector_tmp

      // candidates in candidiate_path_vector_tmp_tmp are moved to candidate_path_vector_tmp
      candidate_path_vector_tmp = candidate_path_vector_tmp_tmp;
      candidate_path_vector_tmp_tmp.clear();

      //--------------------------------------------------
      // STEP 2-2: extend candidate paths to the right
      //--------------------------------------------------
      if (candidate_path_vector_tmp.size() > 0) {
         // each path
         for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size(); it_candidates++) {
            bool correction_success(false);

            extend_first_kmer_to_right(
                                       sequence_modified,
                                       quality_score,
                                       candidate_path_vector_tmp[it_candidates],
                                       kmer_length,
                                       kmer_middle_index,
                                       num_bytes_per_kmer,
                                       residue_kmer,
                                       max_extension,
                                       correction_success
                                      );

            // add this path to candidate_path_vector_tmp_tmp if its correction succeeds
            if (correction_success == true) {
               candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
            }
         }
      }

      // candidiate_path_vector_tmp_tmp: successfully right extended candidates

      //--------------------------------------------------
      // STEP 2-3: choose a final one in candidate_path_vector_tmp_tmp if possible
      //--------------------------------------------------
      // compare quality scores of candidate paths
      // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
      if (candidate_path_vector_tmp_tmp.size() > 1) {
         // each path
         std::vector<C_candidate_path>::iterator it_path;
         std::vector<C_candidate_path>::iterator it_path_1st;
         std::vector<C_candidate_path>::iterator it_path_2nd;

         std::size_t qs_1st(INIT_MIN_QS);
         std::size_t qs_2nd(INIT_MIN_QS);

         // each candidate path
         for (it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); it_path++) {
            // each modification
            for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
               // add quality scores of modified bases
               if (sequence_modified[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
                  (*it_path).sum_qs += ((std::size_t)quality_score[(*it_path).modified_bases[it_mod].first] - quality_score_offset);
               }
            }

            // compare quality scores of each path
            if ((*it_path).sum_qs <= qs_1st) {
               qs_2nd = qs_1st;
               qs_1st = (*it_path).sum_qs;

               it_path_2nd = it_path_1st;
               it_path_1st = it_path;
            }
         }

         // use the 1st path
         // correction succeeds
         if (qs_1st >= qs_2nd + MIN_QS_DIFF) {
            for (std::size_t it_base = 0; it_base < (*it_path_1st).modified_bases.size(); it_base++) {
               // filter out the bases that are equal to the original ones
               if (sequence_modified[(*it_path_1st).modified_bases[it_base].first] != (*it_path_1st).modified_bases[it_base].second) {
                  sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
                  sequence_modified[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;

                  if (candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first < kmer_length) {
                     num_corrected_errors_step2_1++;
                  }
                  else {
                     num_corrected_errors_step2_2++;
                  }
               }
            }
         }
      }
      // only one path
      // correction succeeds
      else if (candidate_path_vector_tmp_tmp.size() == 1) {
         for (std::size_t it_base = 0; it_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_base++) {
            // filter out the bases that are equal to the original ones
            if (sequence_modified[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] != candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second) {
               sequence_modification[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second;
               sequence_modified[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second;

               if (candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first < kmer_length) {
                  num_corrected_errors_step2_1++;
               }
               else {
                  num_corrected_errors_step2_2++;
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_between_solid_regions
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_between_solid_regions(std::string& sequence, const std::string& quality_score, const std::size_t& index_start, const std::size_t& index_end, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification) {
   //--------------------------------------------------
   // from i-th region to (i + 1)-th region
   //--------------------------------------------------
   // i-th solid region | non-solid | (i+1)-th solid region
   // --------------------------------------- read
   //              |-----|                    (index_start)-th k-mer
   //                              |-----|    (index_end)-th k-mer: last non-solid k-mer
   //                         |-----|         (index_last_mod)-th k-mer = (index_end - kmer_length + 1)-th k-mer: last k-mer that can be modified
   //                               |----|    This regions should not be modified
   //--------------------------------------------------
   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;

   // index of the k-mer that can be modified
   // k-mers that are overlapped with a solid regioin cannot be modified
   std::size_t index_last_mod(index_end - kmer_length + 1);

   // make an initial k-mer
   std::string kmer_initial(sequence.substr(index_start, kmer_length));

   // each alternative neocletide
   for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         // generate a new path
         C_candidate_path candidate_path;

         if (sequence[index_start + kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
            std::pair<std::size_t, char> pair_tmp;
            pair_tmp.first  = index_start + kmer_length - 1;
            pair_tmp.second = NEOCLEOTIDE[it_alter];

            candidate_path.modified_bases.push_back(pair_tmp);
         }

         // if this k-mer is the last k-mer that can be modified
         // running extend_a_kmer_right is not needed any more
         if (index_start == index_last_mod) {
            candidate_path_vector_tmp.push_back(candidate_path);
         }
         else {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_a_kmer(
                          kmer_initial,
                          sequence,
                          kmer_length,
                          kmer_middle_index,
                          num_bytes_per_kmer,
                          residue_kmer,
                          index_start,
                          index_last_mod,
                          candidate_path,
                          candidate_path_vector_tmp
                         );
         }
      }
   }

   std::vector<C_candidate_path> candidate_path_vector;

   // check the solidness of k-mers between index_last_mod and index_end
   bool all_solid_wo_modification(false);

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      if ((*it_path).modified_bases.size() == 0) {
         all_solid_wo_modification = true;
         break;
      }
      else {
         // checking is needed
         std::size_t index_last_modified_base((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

         if (index_last_modified_base > index_last_mod) {
            // generate a temporary sequence
            std::string sequence_tmp(sequence);
            for (std::size_t it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
               sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
            }

            // check k-mers
            std::size_t num_success(0);
            for (std::size_t it_check = index_last_mod; it_check <= index_last_modified_base; it_check++) {
               std::string kmer_current(sequence_tmp.substr(it_check, kmer_length));

               if (query_text(kmer_current, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
                  num_success++;
               }
               else {
                  break;
               }
            }

            if (num_success == (index_last_modified_base - index_last_mod + 1)) {
               candidate_path_vector.push_back(*it_path);
            }
         }
         // checking is not needed
         else {
            candidate_path_vector.push_back(*it_path);
         }
      }
   }

   // all k-mers are solid without any modification
   // do nothing
   if (all_solid_wo_modification == true) {
   }
   // compare quality scores of candidate paths
   // if the number of paths in candidate_path_vector is larger than 1
   else if (candidate_path_vector.size() > 1) {
      // each path
      std::vector<C_candidate_path>::iterator it_path;
      std::vector<C_candidate_path>::iterator it_path_1st;
      std::vector<C_candidate_path>::iterator it_path_2nd;

      std::size_t qs_1st(INIT_MIN_QS);
      std::size_t qs_2nd(INIT_MIN_QS);

      for (it_path = candidate_path_vector.begin(); it_path != candidate_path_vector.end(); it_path++) {
         // each modification
         for (std::vector< std::pair<std::size_t, char> >::iterator it_mod = (*it_path).modified_bases.begin(); it_mod != (*it_path).modified_bases.end(); it_mod++) {
            // add quality scores of modified bases
            (*it_path).sum_qs += ((std::size_t)quality_score[(*it_mod).first] - quality_score_offset);
         }

         // compare quality scores of each path
         if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = qs_1st;
            qs_1st = (*it_path).sum_qs;

            it_path_2nd = it_path_1st;
            it_path_1st = it_path;
         }
      }

      // use the 1st path
      if (qs_1st >= qs_2nd + MIN_QS_DIFF) {
         // each modification
         for (std::vector< std::pair<std::size_t, char> >::iterator it_base = (*it_path_1st).modified_bases.begin(); it_base != (*it_path_1st).modified_bases.end(); it_base++) {
            // update sequence_modification
            sequence_modification[(*it_base).first] = (*it_base).second;
            sequence[(*it_base).first] = (*it_base).second;
            num_corrected_errors_step1_1++;
         }
      }
   }
   // only one path
   else if (candidate_path_vector.size() == 1) {
      // each modification
      for (std::size_t it_base = 0; it_base < candidate_path_vector[0].modified_bases.size(); it_base++) {
         // update sequence_modification
         sequence_modification[candidate_path_vector[0].modified_bases[it_base].first] = candidate_path_vector[0].modified_bases[it_base].second;
         sequence[candidate_path_vector[0].modified_bases[it_base].first] = candidate_path_vector[0].modified_bases[it_base].second;
         num_corrected_errors_step1_1++;
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_5_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_5_prime_end(std::string& sequence, const std::string& quality_score, const std::size_t& index_start, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification, const std::size_t& max_extension, const bool& verify) {
   // |  non-solid  | 1st solid region
   // |--------------------------------------| read
   //         |-----|                          (index_start)-th k-mer
   //--------------------------------------------------
   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;

   // make an initial k-mer
   std::string kmer_initial(sequence.substr(index_start, kmer_length));

   // each alternative neocletide
   for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[0] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         // generate a new path
         C_candidate_path candidate_path;

         if (sequence[index_start] != NEOCLEOTIDE[it_alter]) {
            std::pair<std::size_t, char> pair_tmp;
            pair_tmp.first  = index_start;
            pair_tmp.second = NEOCLEOTIDE[it_alter];

            candidate_path.modified_bases.push_back(pair_tmp);
         }

         // if this k-mer is the first k-mer in a read
         // running extend_a_kmer_5_prime_end is not needed any more
         if (index_start == 0) {
            candidate_path_vector_tmp.push_back(candidate_path);
         }
         else if (index_start > 0) {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_a_kmer_5_prime_end(
                                      kmer_initial,
                                      sequence,
                                      kmer_length,
                                      kmer_middle_index,
                                      num_bytes_per_kmer,
                                      residue_kmer,
                                      index_start,
                                      candidate_path,
                                      candidate_path_vector_tmp
                                     );
         }
      }
      else {
      }
   }

   std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

   // check the solidness of the leftmost k-mers of each modified base
   bool all_solid_wo_modification(false);

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      if ((*it_path).modified_bases.size() == 0) {
         all_solid_wo_modification = true;
         break;
      }
      else {
         std::string sequence_tmp(sequence);

         // index_smallest_modified
         std::size_t index_smallest_modified((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

         // number of bases that should be extended
         std::size_t extend_amount;

         // calculate extend_amount
         // no extension is needed
         // kmer_length = 11, max_extension = 5
         // |0|0|0|0|0|0|0|0|0|0|1|1|1|-
         // |0|1|2|3|4|5|6|7|8|9|0|1|2|-
         // |<------------------->|      k = 11
         // |--------------------------- read
         //                     |<------ index_smallest_modified >= 10
         if (index_smallest_modified >= kmer_length - 1) {
            candidate_path_vector_tmp_tmp.push_back(*it_path);
         }
         // extension is needed
         else {
            // applied the modified bases to sequence_tmp
            for (std::size_t it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
               // modify sequence_tmp
               sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
            }

            // determine the number of extensions
            // extension amount = kmer_length - index_smallest_modified - 1
            // kmer_length = 11, max_extension = 5
            // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
            // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
            //           |<------------------->|      k = 11
            //           |--------------------------- read
            //                     |<------->|        (index_smallest_modified < 10) AND (index_smallest_modified >= 5)
            //     |<------------------->|            index_smallest_modified = 7 -> extend_amount = 3
            if (index_smallest_modified >= kmer_length - max_extension - 1) {
               extend_amount = kmer_length - index_smallest_modified - 1;
            }
            // kmer_length = 11, max_extension = 5
            // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
            // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
            //           |<------------------->|      k = 11
            //           |--------------------------- read
            //           |<------->|                  index_smallest_modified < 5
            else {
               extend_amount = max_extension;
            }

            bool extension_success(false);

            // generate an initial k-mer
            std::string kmer_initial(sequence_tmp.substr(0, kmer_length - 1));
            kmer_initial = '0' + kmer_initial;

            // each alternative neocletide
            for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
               // make a change
               kmer_initial[0] = NEOCLEOTIDE[it_alter];

               // kmer_initial is solid
               if (query_text(kmer_initial, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
                  // if extend_amount == 1
                  // running extend_out_left is not needed any more
                  if (extend_amount == 1) {
                     extension_success = true;
                     break;
                  }
                  else {
                     // trace  this kmer recursively and update candidate_path_vector_tmp
                     extend_out_left(
                                     kmer_initial,
                                     kmer_length,
                                     kmer_middle_index,
                                     num_bytes_per_kmer,
                                     residue_kmer,
                                     1,
                                     extend_amount,
                                     extension_success
                                    );
                  }
               }
            }

            if (extension_success == true) {
               candidate_path_vector_tmp_tmp.push_back(*it_path);
            }
         }
      }
   }

   // all k-mers are solid without any modification
   // do nothing
   if (all_solid_wo_modification == true) {
   }
   // compare quality scores of candidate paths
   // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
   else if (candidate_path_vector_tmp_tmp.size() > 1) {
      // each path
      std::vector<C_candidate_path>::iterator it_path;
      std::vector<C_candidate_path>::iterator it_path_1st;
      std::vector<C_candidate_path>::iterator it_path_2nd;

      std::size_t qs_1st(INIT_MIN_QS);
      std::size_t qs_2nd(INIT_MIN_QS);

      // each candidate path
      for (it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); it_path++) {
         // each modification
         for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
            // add quality scores of modified bases
            if (sequence[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
               (*it_path).sum_qs += ((std::size_t)quality_score[(*it_path).modified_bases[it_mod].first] - quality_score_offset);
            }
         }

         // compare quality scores of each path
         if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = qs_1st;
            qs_1st = (*it_path).sum_qs;

            it_path_2nd = it_path_1st;
            it_path_1st = it_path;
         }
      }

      // use the 1st path
      if (qs_1st >= qs_2nd + MIN_QS_DIFF) {
         // update sequence_modification
         for (std::size_t it_base = 0; it_base < (*it_path_1st).modified_bases.size(); it_base++) {
            if (verify == true) {
               if ((*it_path_1st).modified_bases[it_base].first <= (sequence.length() - kmer_length)) {
                  sequence_modification[(*it_path_1st).modified_bases[it_base].first] = tolower((*it_path_1st).modified_bases[it_base].second);
               }
               else {
                  std::cout << std::endl << "ERROR: This algorithm cannot change " << it_base << "-th base" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: This algorithm cannot change " << it_base << "-th base" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
            else {
               sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
            }

            sequence[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
            num_corrected_errors_step1_2++;
         }
      }
   }
   // only one path
   else if (candidate_path_vector_tmp_tmp.size() == 1) {
      std::vector<C_candidate_path>::iterator it_path(candidate_path_vector_tmp_tmp.begin());

      // update sequence_modification
      for (std::size_t it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
         if (verify == true) {
            if ((*it_path).modified_bases[it_base].first <= (sequence.length() - kmer_length)) {
               sequence_modification[(*it_path).modified_bases[it_base].first] = tolower((*it_path).modified_bases[it_base].second);
            }
            else {
               std::cout << std::endl << "ERROR: This algorithm cannot change " << it_base << "-th base" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: This algorithm cannot change " << it_base << "-th base" << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }
         }
         else {
            sequence_modification[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
         }

         sequence[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
         num_corrected_errors_step1_2++;
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_3_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_3_prime_end(std::string& sequence, const std::string& quality_score, const std::size_t& index_start, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification, const std::size_t& max_extension) {
   //  last solid region | non-solid region |
   // --------------------------------------| read
   //               |-----|                   (index_start)-th k-mer
   //--------------------------------------------------
   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;

   // make an initial k-mer
   std::string kmer_initial(sequence.substr(index_start, kmer_length));

   // each alternative neocletide
   for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         // generate a new path
         C_candidate_path candidate_path;

         if (sequence[index_start + kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
            std::pair<std::size_t, char> pair_tmp;
            pair_tmp.first  = index_start + kmer_length - 1;
            pair_tmp.second = NEOCLEOTIDE[it_alter];

            candidate_path.modified_bases.push_back(pair_tmp);
         }

         // if this k-mer is the last k-mer in a read
         // running extend_a_kmer_3_prime_end is not needed any more
         if (index_start == (sequence.length() - kmer_length)) {
            candidate_path_vector_tmp.push_back(candidate_path);
         }
         else if (index_start < (sequence.length() - kmer_length)) {
         // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_a_kmer_3_prime_end(
                                      kmer_initial,
                                      sequence,
                                      kmer_length,
                                      kmer_middle_index,
                                      num_bytes_per_kmer,
                                      residue_kmer,
                                      index_start,
                                      candidate_path,
                                      candidate_path_vector_tmp
                                     );
         }
      }
   }

   std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

   // check the solidness of the rightmost k-mers of each modified base
   bool all_solid_wo_modification(false);

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      if ((*it_path).modified_bases.size() == 0) {
         all_solid_wo_modification = true;
         break;
      }
      else {
         std::string sequence_tmp(sequence);

         // index_largest_modified
         std::size_t index_largest_modified((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

         // number of bases that should be extended
         std::size_t extend_amount;

         // calculate extend_amount
         // no extension is needed
         // sequence.length() = 20, kmer_length = 11, max_extension = 5
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|
         //     |<------------------->| k = 11
         // --------------------------| read
         // ----->|                     index_largest_modified <= 9
         if (index_largest_modified <= sequence.length() - kmer_length) {
            candidate_path_vector_tmp_tmp.push_back(*it_path);
         }
         // extension is needed
         else {
            // applied the modified bases to sequence_tmp
            for (std::size_t it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
               // modify sequence_tmp
               sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
            }

            // determine the number of extensions
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
            //     |<------------------->|           k = 11
            // --------------------------|           read
            //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
            //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
            if (index_largest_modified <= sequence.length() + max_extension - kmer_length) {
               extend_amount = kmer_length - (sequence.length() - index_largest_modified);
            }
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
            //     |<------------------->|           k = 11
            // --------------------------|           read
            //                 |<------->|           index_largest_modified > 15
            else {
               extend_amount = max_extension;
            }

            bool extension_success(false);

            // generate an initial k-mer
            // sequence.length() = 20, kmer_length = 11
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|
            //       |<----------------->| kmer_length - 1 = 10
            // --------------------------| read
            //       |-|                   20 - 11 + 1 = 10
            std::string kmer_initial(sequence_tmp.substr(sequence.length() - kmer_length + 1, kmer_length - 1));
            kmer_initial = kmer_initial + '0';

            // each alternative neocletide
            for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
               // make a change
               kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

               // kmer_initial is solid
               if (query_text(kmer_initial, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
                  // if extend_amount == 1
                  // running extend_out_right is not needed any more
                  if (extend_amount == 1) {
                     extension_success = true;
                     break;
                  }
                  else {
                     // trace  this kmer recursively and update candidate_path_vector_tmp
                     extend_out_right(
                                      kmer_initial,
                                      kmer_length,
                                      kmer_middle_index,
                                      num_bytes_per_kmer,
                                      residue_kmer,
                                      1,
                                      extend_amount,
                                      extension_success
                                     );
                  }
               }
            }

            if (extension_success == true) {
               candidate_path_vector_tmp_tmp.push_back(*it_path);
            }
         }
      }
   }

   // all k-mers are solid without any modification
   // do nothing
   if (all_solid_wo_modification == true) {
   }
   // compare quality scores of candidate paths
   // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
   else if (candidate_path_vector_tmp_tmp.size() > 1) {
      // each path
      std::vector<C_candidate_path>::iterator it_path;
      std::vector<C_candidate_path>::iterator it_path_1st;
      std::vector<C_candidate_path>::iterator it_path_2nd;

      std::size_t qs_1st(INIT_MIN_QS);
      std::size_t qs_2nd(INIT_MIN_QS);

      // each candidate path
      for (it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); it_path++) {
         // each modification
         for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
            // add quality scores of modified bases
            (*it_path).sum_qs += ((std::size_t)quality_score[(*it_path).modified_bases[it_mod].first] - quality_score_offset);
         }

         // compare quality scores of each path
         if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = qs_1st;
            qs_1st = (*it_path).sum_qs;

            it_path_2nd = it_path_1st;
            it_path_1st = it_path;
         }
      }

      // use the 1st path
      if (qs_1st >= qs_2nd + MIN_QS_DIFF) {
         // update sequence_modification
         for (std::size_t it_base = 0; it_base < (*it_path_1st).modified_bases.size(); it_base++) {
            sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
            sequence[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
            num_corrected_errors_step1_3++;
         }
      }
   }
   // only one path
   else if (candidate_path_vector_tmp_tmp.size() == 1) {
      std::vector<C_candidate_path>::iterator it_path(candidate_path_vector_tmp_tmp.begin());

      // update sequence_modification
      for (std::size_t it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
         sequence_modification[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
         sequence[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
         num_corrected_errors_step1_3++;
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_first_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_first_kmer(const std::string& sequence, const std::string& quality_score, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification, std::vector<C_candidate_path>& candidate_path_vector) {
   std::string first_kmer(sequence.substr(0, kmer_length));

   std::vector<std::size_t> low_qs_indexes;

   for (std::size_t it_bases = 0; it_bases < kmer_length; it_bases++) {
      if (((std::size_t)quality_score[it_bases] - quality_score_offset) < QS_CUTOFF) {
         low_qs_indexes.push_back(it_bases);
      }
   }

   // correct errors if the number of low-quality bases is smaller than the threshold
   if ((low_qs_indexes.size() <= MAX_LOW_QS_BASES) && (low_qs_indexes.size() > 0)) {
      C_candidate_path candidate_path;

      check_first_kmer(
                       first_kmer,
                       candidate_path,
                       low_qs_indexes,
                       candidate_path_vector,
                       0,
                       kmer_length,
                       kmer_middle_index,
                       num_bytes_per_kmer,
                       residue_kmer
                      );

      // no candidate path is found
      if (candidate_path_vector.size() == 0) {
         for (std::size_t it_bases = 0; it_bases < kmer_length; it_bases++) {
            std::string kmer_tmp(first_kmer);

            // each alternative neocletide
            for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
               // not equal to the original character
               if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                  // generate a new k-mer
                  kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                  // add kmer_tmp to candidate_path_tmp if it is solid
                  if (query_text(kmer_tmp, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
                     // generate a new candidate path
                     candidate_path.clear_path();

                     std::pair<std::size_t, char> pair_tmp;
                     pair_tmp.first = it_bases;
                     pair_tmp.second = NEOCLEOTIDE[it_alter];
                     candidate_path.modified_bases.push_back(pair_tmp);

                     candidate_path_vector.push_back(candidate_path);
                  }
               }
            }
         }
      }
   }
   // no low-quality base or too many low-quality bases
   else {
      if (query_text(first_kmer, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         C_candidate_path candidate_path;

         candidate_path_vector.push_back(candidate_path);
      }
      else {
         for (std::size_t it_bases = 0; it_bases < kmer_length; it_bases++) {
            std::string kmer_tmp(first_kmer);

            // each alternative neocletide
            for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
               // not equal to the original character
               if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                  // generate a new k-mer
                  kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                  // add kmer_tmp to candidate_path_tmp if it is solid
                  if (query_text(kmer_tmp, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
                     // generate a new candidate path
                     C_candidate_path candidate_path;

                     std::pair<std::size_t, char> pair_tmp;
                     pair_tmp.first = it_bases;
                     pair_tmp.second = NEOCLEOTIDE[it_alter];
                     candidate_path.modified_bases.push_back(pair_tmp);

                     candidate_path_vector.push_back(candidate_path);
                  }
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// check_first_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::check_first_kmer(const std::string& kmer, const C_candidate_path& candidate_path_in, const std::vector<std::size_t>& low_qs_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& index, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer) {
   std::string new_kmer(kmer);

   for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
      // make a new k-mer
      new_kmer[low_qs_indexes[index]] = NEOCLEOTIDE[it_alter];

      C_candidate_path candidate_path_next(candidate_path_in);

      if (kmer[low_qs_indexes[index]] != NEOCLEOTIDE[it_alter]) {
         std::pair<std::size_t, char> pair_tmp;
         pair_tmp.first = low_qs_indexes[index];
         pair_tmp.second = NEOCLEOTIDE[it_alter];
         candidate_path_next.modified_bases.push_back(pair_tmp);
      }

      if (index == low_qs_indexes.size() - 1) {
         if (query_text(new_kmer, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
            candidate_path_vector.push_back(candidate_path_next);
         }
      }
      else {
         check_first_kmer(
                          new_kmer,
                          candidate_path_next,
                          low_qs_indexes,
                          candidate_path_vector,
                          index + 1,
                          kmer_length,
                          kmer_middle_index,
                          num_bytes_per_kmer,
                          residue_kmer
                         );
      }
   }
}



//----------------------------------------------------------------------
// solid_first_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::solid_first_kmer(const C_candidate_path& candidate_path, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& max_extension, bool& extension_success) {
   // index_smallest_modified
   std::size_t index_smallest_modified(candidate_path.modified_bases[0].first);

   // number of bases that should be extended
   std::size_t extend_amount;

   // applied the modified bases to first_kmer
   std::string first_kmer(sequence.substr(0, kmer_length));
   for (std::size_t it_base = 0; it_base < candidate_path.modified_bases.size(); it_base++) {
      first_kmer[candidate_path.modified_bases[it_base].first] = candidate_path.modified_bases[it_base].second;
   }

   // determine the number of extensions
   // extension amount = kmer_length - index_smallest_modified - 1
   // kmer_length = 11, max_extension = 5
   // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
   // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
   //           |<------------------->|      k = 11
   //           |--------------------------- read
   //                     |<------->|        (index_smallest_modified < 10) AND (index_smallest_modified >= 5)
   //     |<------------------->|            index_smallest_modified = 7 -> extend_amount = 3
   if (index_smallest_modified >= kmer_length - max_extension - 1) {
      extend_amount = kmer_length - index_smallest_modified - 1;
   }
   // kmer_length = 11, max_extension = 5
   // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
   // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
   //           |<------------------->|      k = 11
   //           |--------------------------- read
   //           |<------->|                  index_smallest_modified < 5
   else {
      extend_amount = max_extension;
   }

   // generate an initial k-mer
   std::string kmer_initial(first_kmer.substr(0, kmer_length - 1));
   kmer_initial = '0' + kmer_initial;

   // each alternative neocletide
   for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[0] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         // if extend_amount == 1
         // running extend_out_left is not needed any more
         if (extend_amount == 1) {
            extension_success = true;
            break;
         }
         else {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_out_left(
                            kmer_initial,
                            kmer_length,
                            kmer_middle_index,
                            num_bytes_per_kmer,
                            residue_kmer,
                            1,
                            extend_amount,
                            extension_success
                           );
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_first_kmer_to_right
//----------------------------------------------------------------------
inline void C_correct_errors::extend_first_kmer_to_right(const std::string& sequence, const std::string& quality_score, C_candidate_path& candidate_path_in, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& max_extension, bool& correction_success) {
   // generate the first k-mer
   std::string first_kmer(sequence.substr(0, kmer_length + 1));
   for (std::size_t it_base = 0; it_base < candidate_path_in.modified_bases.size(); it_base++) {
      first_kmer[candidate_path_in.modified_bases[it_base].first] = candidate_path_in.modified_bases[it_base].second;
   }

   // generate the second k-mer
   std::string second_kmer(first_kmer.substr(1, kmer_length));

   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;


   // second_kmer is solid
   if (query_text(second_kmer, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
      // if this k-mer is the last k-mer in a read
      // running extend_a_kmer_3_prime_end is not needed any more
      if ((sequence.length() - kmer_length) == 1) {
         candidate_path_vector_tmp.push_back(candidate_path_in);
      }
      else if ((sequence.length() - kmer_length) > 1) {
         // trace  this kmer recursively and update candidate_path_vector_tmp
         extend_a_kmer_3_prime_end(
                                   second_kmer,
                                   sequence,
                                   kmer_length,
                                   kmer_middle_index,
                                   num_bytes_per_kmer,
                                   residue_kmer,
                                   1,
                                   candidate_path_in,
                                   candidate_path_vector_tmp
                                  );
      }
   }
   // second_kmer is not solid
   else {
      // each alternative neocletide
      for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
         // not equal to the original character
         if (sequence[kmer_length] != NEOCLEOTIDE[it_alter]) {
            // make a change
            second_kmer[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // new second_kmer is solid
            if (query_text(second_kmer, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
               // generate a new path
               C_candidate_path new_candidate_path(candidate_path_in);

               std::pair<std::size_t, char> pair_tmp;
               pair_tmp.first  = kmer_length;
               pair_tmp.second = NEOCLEOTIDE[it_alter];

               new_candidate_path.modified_bases.push_back(pair_tmp);

               // if this k-mer is the last k-mer in a read
               // running extend_a_kmer_3_prime_end is not needed any more
               if ((sequence.length() - kmer_length) == 1) {
                  candidate_path_vector_tmp.push_back(new_candidate_path);
               }
               else if ((sequence.length() - kmer_length) > 1) {
                  // trace  this kmer recursively and update candidate_path_vector_tmp
                  extend_a_kmer_3_prime_end(
                                            second_kmer,
                                            sequence,
                                            kmer_length,
                                            kmer_middle_index,
                                            num_bytes_per_kmer,
                                            residue_kmer,
                                            1,
                                            new_candidate_path,
                                            candidate_path_vector_tmp
                                           );
               }
            }
         }
      }
   }

   // check the solidness of the rightmost k-mers of each modified base
   std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      std::string sequence_tmp(sequence);

      // index_largest_modified
      std::size_t index_largest_modified((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

      // number of bases that should be extended
      std::size_t extend_amount;

      // calculate extend_amount
      // no extension is needed
      // sequence.length() = 20, kmer_length = 11, max_extension = 5
      // |0|0|0|1|1|1|1|1|1|1|1|1|1|
      // |7|8|9|0|1|2|3|4|5|6|7|8|9|
      //     |<------------------->| k = 11
      // --------------------------| read
      // ----->|                     index_largest_modified <= 9
      if (index_largest_modified <= sequence.length() - kmer_length) {
         candidate_path_vector_tmp_tmp.push_back(*it_path);
      }
      // extension is needed
      else {
         // applied the modified bases to sequence_tmp
         for (std::size_t it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
            // modify sequence_tmp
            sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
         }

         // determine the number of extensions
         // sequence.length() = 20, kmer_length = 11, max_extension = 5
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
         //     |<------------------->|           k = 11
         // --------------------------|           read
         //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
         //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
         if (index_largest_modified <= sequence.length() + max_extension - kmer_length) {
            extend_amount = kmer_length - (sequence.length() - index_largest_modified);
         }
         // sequence.length() = 20, kmer_length = 11, max_extension = 5
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
         //     |<------------------->|           k = 11
         // --------------------------|           read
         //                 |<------->|           index_largest_modified > 15
         else {
            extend_amount = max_extension;
         }

         bool extension_success(false);

         // generate an initial k-mer
         // sequence.length() = 20, kmer_length = 11
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|
         //       |<----------------->| kmer_length - 1 = 10
         // --------------------------| read
         //       |-|                   20 - 11 + 1 = 10
         std::string kmer_initial(sequence_tmp.substr(sequence.length() - kmer_length + 1, kmer_length - 1));
         kmer_initial = kmer_initial + '0';

         // each alternative neocletide
         for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            if (query_text(kmer_initial, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
               // if extend_amount == 1
               // running extend_out_right is not needed any more
               if (extend_amount == 1) {
                  extension_success = true;
                  break;
               }
               else {
                  // trace  this kmer recursively and update candidate_path_vector_tmp
                  extend_out_right(
                                   kmer_initial,
                                   kmer_length,
                                   kmer_middle_index,
                                   num_bytes_per_kmer,
                                   residue_kmer,
                                   1,
                                   extend_amount,
                                   extension_success
                                  );
               }
            }
         }

         if (extension_success == true) {
            candidate_path_vector_tmp_tmp.push_back(*it_path);
         }
      }
   }

   // initialize candidate_path_vector_tmp
   candidate_path_vector_tmp.clear();

   // compare quality scores of candidate paths
   // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
   if (candidate_path_vector_tmp_tmp.size() > 1) {
      // each path
      std::vector<C_candidate_path>::iterator it_path;
      std::vector<C_candidate_path>::iterator it_path_1st;
      std::vector<C_candidate_path>::iterator it_path_2nd;

      std::size_t qs_1st(INIT_MIN_QS);
      std::size_t qs_2nd(INIT_MIN_QS);

      // each candidate path
      for (it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); it_path++) {
         // each modification
         for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
            // add quality scores of modified bases
            if (sequence[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
               (*it_path).sum_qs += ((std::size_t)quality_score[(*it_path).modified_bases[it_mod].first] - quality_score_offset);
            }
         }

         // compare quality scores of each path
         if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = qs_1st;
            qs_1st = (*it_path).sum_qs;

            it_path_2nd = it_path_1st;
            it_path_1st = it_path;
         }
      }

      // use the 1st path
      // correction succeeds
      if (qs_1st >= qs_2nd + MIN_QS_DIFF) {
         correction_success = true;
         candidate_path_in = *it_path_1st;
      }
   }
   // only one path
   // correction succeeds
   else if (candidate_path_vector_tmp_tmp.size() == 1) {
      correction_success = true;
      candidate_path_in = candidate_path_vector_tmp_tmp[0];
   }
}



//----------------------------------------------------------------------
// extend_a_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::extend_a_kmer(const std::string& kmer, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& index_kmer, const std::size_t& index_last_mod, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector) {
   // generate a new k-mer
   std::string kmer_new(kmer.substr(1, kmer_length - 1));
   kmer_new.push_back(sequence[index_kmer + kmer_length]);

   // kmer_new is a solid k-mer
   if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
      // if this k-mer is the last k-mer that can be modified
      // running extend_a_kmer_right is not needed any more
      if ((index_kmer + 1) == index_last_mod) {
         candidate_path_vector.push_back(current_path);
      }
      else {
         extend_a_kmer(
                       kmer_new,
                       sequence,
                       kmer_length,
                       kmer_middle_index,
                       num_bytes_per_kmer,
                       residue_kmer,
                       index_kmer + 1,
                       index_last_mod,
                       current_path,
                       candidate_path_vector
                      );
      }
   }
   else {
      // each alternative neocletide
      for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
         // not equal to the original character
         if (sequence[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter]) {
            // make a change
            kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_new is solid
            if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
               // generate a new path
               C_candidate_path temporary_path(current_path);

               std::pair<std::size_t, char> pair_tmp;
               pair_tmp.first  = index_kmer + kmer_length;
               pair_tmp.second = NEOCLEOTIDE[it_alter];

               temporary_path.modified_bases.push_back(pair_tmp);

               // if this k-mer is the last k-mer that can be modified
               // running extend_a_kmer_right is not needed any more
               if ((index_kmer + 1) == index_last_mod) {
                  candidate_path_vector.push_back(temporary_path);
               }
               else {
                  // trace  this kmer recursively and update candidate_path_vector
                  extend_a_kmer(
                                kmer_new,
                                sequence,
                                kmer_length,
                                kmer_middle_index,
                                num_bytes_per_kmer,
                                residue_kmer,
                                index_kmer + 1,
                                index_last_mod,
                                temporary_path,
                                candidate_path_vector
                               );
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_a_kmer_5_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::extend_a_kmer_5_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector) {
   // generate a new k-mer
   std::string kmer_new(kmer.substr(0, kmer_length - 1));
   kmer_new = sequence[index_kmer - 1] + kmer_new;

   // kmer_new is a solid k-mer
   if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
      // if this k-mer is the first k-mer in a read
      // running extend_a_kmer_5_prime_end is not needed any more
      if ((index_kmer - 1) == 0) {
         candidate_path_vector.push_back(current_path);
      }
      else if ((index_kmer - 1) > 0) {
         extend_a_kmer_5_prime_end(
                                   kmer_new,
                                   sequence,
                                   kmer_length,
                                   kmer_middle_index,
                                   num_bytes_per_kmer,
                                   residue_kmer,
                                   index_kmer - 1,
                                   current_path,
                                   candidate_path_vector
                                  );
      }
   }
   else {
      // each alternative neocletide
      for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
         // not equal to the original character
         if (sequence[index_kmer - 1] != NEOCLEOTIDE[it_alter]) {
            // make a change
            kmer_new[0] = NEOCLEOTIDE[it_alter];

            // kmer_new is solid
            if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
               // generate a new path
               C_candidate_path temporary_path(current_path);

               std::pair<std::size_t, char> pair_tmp;
               pair_tmp.first  = index_kmer - 1;
               pair_tmp.second = NEOCLEOTIDE[it_alter];

               temporary_path.modified_bases.push_back(pair_tmp);

               // if this k-mer is the first k-mer in a read
               // running extend_a_kmer_5_prime_end is not needed any more
               if ((index_kmer - 1) == 0) {
                  candidate_path_vector.push_back(temporary_path);
               }
               else if ((index_kmer - 1) > 0) {
                  // trace  this kmer recursively and update candidate_path_vector
                  extend_a_kmer_5_prime_end(
                                            kmer_new,
                                            sequence,
                                            kmer_length,
                                            kmer_middle_index,
                                            num_bytes_per_kmer,
                                            residue_kmer,
                                            index_kmer - 1,
                                            temporary_path,
                                            candidate_path_vector
                                           );
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_a_kmer_3_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::extend_a_kmer_3_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector) {
   // generate a new k-mer
   std::string kmer_new(kmer.substr(1, kmer_length - 1));
   kmer_new = kmer_new + sequence[index_kmer + kmer_length];

   // kmer_new is a solid k-mer
   if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
      // if this k-mer is the last k-mer in a read
      // running extend_a_kmer_3_prime_end is not needed any more
      if ((index_kmer + 1) == (sequence.length() - kmer_length)) {
         candidate_path_vector.push_back(current_path);
      }
      else if ((index_kmer + 1) < (sequence.length() - kmer_length)) {
         extend_a_kmer_3_prime_end(
                                   kmer_new,
                                   sequence,
                                   kmer_length,
                                   kmer_middle_index,
                                   num_bytes_per_kmer,
                                   residue_kmer,
                                   index_kmer + 1,
                                   current_path,
                                   candidate_path_vector
                                  );
      }
   }
   else {
      // each alternative neocletide
      for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
         // not equal to the original character
         if (sequence[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter]) {
            // make a change
            kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_new is solid
            if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
               // generate a new path
               C_candidate_path temporary_path(current_path);

               std::pair<std::size_t, char> pair_tmp;
               pair_tmp.first  = index_kmer + kmer_length;
               pair_tmp.second = NEOCLEOTIDE[it_alter];

               temporary_path.modified_bases.push_back(pair_tmp);

               // if this k-mer is the last k-mer in a read
               // running extend_a_kmer_3_prime_end is not needed any more
               if ((index_kmer + 1) == (sequence.length() - kmer_length)) {
                  candidate_path_vector.push_back(temporary_path);
               }
               else if ((index_kmer + 1) < (sequence.length() - kmer_length)) {
                  // trace  this kmer recursively and update candidate_path_vector
                  extend_a_kmer_3_prime_end(
                                            kmer_new,
                                            sequence,
                                            kmer_length,
                                            kmer_middle_index,
                                            num_bytes_per_kmer,
                                            residue_kmer,
                                            index_kmer + 1,
                                            temporary_path,
                                            candidate_path_vector
                                           );
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_out_left
//----------------------------------------------------------------------
inline void C_correct_errors::extend_out_left(const std::string& kmer, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success) {
   // generate a new k-mer
   std::string kmer_new(kmer.substr(0, kmer_length - 1));
   kmer_new = '0' + kmer_new;

   // each alternative neocletide
   for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
      // generate kmer_new
      kmer_new[0] = NEOCLEOTIDE[it_alter];

      // kmer_new is solid
      if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         // if current num_extend = extend_amount
         // running extend_out_left is not needed any more
         if ((num_extend + 1) == extend_amount) {
            extension_success = true;
            break;
         }
         else {
            // trace  this kmer recursively
            extend_out_left(
                            kmer_new,
                            kmer_length,
                            kmer_middle_index,
                            num_bytes_per_kmer,
                            residue_kmer,
                            num_extend + 1,
                            extend_amount,
                            extension_success
                           );
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_out_right
//----------------------------------------------------------------------
inline void C_correct_errors::extend_out_right(const std::string& kmer, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success) {
   // generate a new k-mer
   std::string kmer_new(kmer.substr(1, kmer_length - 1));
   kmer_new = kmer_new + '0';

   // each alternative neocletide
   for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
      // generate kmer_new
      kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

      // kmer_new is solid
      if (query_text(kmer_new, kmer_middle_index, num_bytes_per_kmer, kmer_length, residue_kmer) == true) {
         // if current num_extend = extend_amount
         // running extend_out_right is not needed any more
         if ((num_extend + 1) == extend_amount) {
            extension_success = true;
            break;
         }
         else {
            // trace  this kmer recursively
            extend_out_right(
                             kmer_new,
                             kmer_length,
                             kmer_middle_index,
                             num_bytes_per_kmer,
                             residue_kmer,
                             num_extend + 1,
                             extend_amount,
                             extension_success
                            );
         }
      }
   }
}



//----------------------------------------------------------------------
// encode_error_correction_info
//----------------------------------------------------------------------
inline void C_correct_errors::encode_correction_info(char& buffer, const char& char_in_read, const char& char_in_info) {
   switch (char_in_info) {
      // no modification
      case '0':
         buffer = buffer << 2;
         break;
      case 'A':
         switch (char_in_read) {
            // A -> A
            // error
            case 'A':
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
            // C -> A
            // 11
            case 'C':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            // G -> A
            // 10
            case 'G':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               // buffer |= BIT1;
               break;
            // T -> A
            // 01
            case 'T':
               buffer = buffer << 1;
               // buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            default:
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      case 'C':
         switch (char_in_read) {
            // A -> C
            // 01
            case 'A':
               buffer = buffer << 1;
               // buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            // C -> C
            // error
            case 'C':
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
            // G -> C
            // 11
            case 'G':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            // T -> C
            // 10
            case 'T':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               // buffer |= BIT1;
               break;
            default:
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      case 'G':
         switch (char_in_read) {
            // A -> G
            // 10
            case 'A':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               // buffer |= BIT1;
               break;
            // C -> G
            // 01
            case 'C':
               buffer = buffer << 1;
               // buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            // G -> G
            // error
            case 'G':
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
            // T -> G
            // 11
            case 'T':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            default:
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      case 'T':
         switch (char_in_read) {
            // A -> T
            // 11
            case 'A':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            // C -> T
            // 10
            case 'C':
               buffer = buffer << 1;
               buffer |= BIT1;
               buffer = buffer << 1;
               // buffer |= BIT1;
               break;
            // G -> T
            // 01
            case 'G':
               buffer = buffer << 1;
               // buffer |= BIT1;
               buffer = buffer << 1;
               buffer |= BIT1;
               break;
            // T -> T
            // error
            case 'T':
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
            default:
               std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " " << char_in_read << " (encode_error_info)" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      default:
         std::cout << std::endl << "ERROR: Illegal character " << char_in_info << " (encode_error_info)" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal character " << char_in_info << " (encode_error_info)" << std::endl << std::endl;
         exit(EXIT_FAILURE);
   }
}



//----------------------------------------------------------------------
// decode_correction_info
//----------------------------------------------------------------------
inline char C_correct_errors::decode_correction_info(const char& first, const char& second, const char& read) {
   switch (first) {
      case '1' :
         switch (second) {
            // 11
            case '1' :
               switch (read) {
                  case 'A' :
                     return 'T';
                     break;
                  case 'C' :
                     return 'A';
                     break;
                  case 'G' :
                     return 'C';
                     break;
                  case 'T' :
                     return 'G';
                     break;
                  default :
                     std::cout << std::endl << "ERROR: Illegal character " << read  << " (decode_correction_info)" << std::endl << std::endl;
                     f_log     << std::endl << "ERROR: Illegal character " << read  << " (decode_correction_info)" << std::endl << std::endl;
                     exit(EXIT_FAILURE);
                     break;
               }
               break;
            // 10
            case '0' :
               switch (read) {
                  case 'A' :
                     return 'G';
                     break;
                  case 'C' :
                     return 'T';
                     break;
                  case 'G' :
                     return 'A';
                     break;
                  case 'T' :
                     return 'C';
                     break;
                  default :
                     std::cout << std::endl << "ERROR: Illegal character " << read  << " (decode_correction_info)" << std::endl << std::endl;
                     f_log     << std::endl << "ERROR: Illegal character " << read  << " (decode_correction_info)" << std::endl << std::endl;
                     exit(EXIT_FAILURE);
                     break;
               }
               break;
            default :
               std::cout << std::endl << "ERROR: Illegal result in correction information file " << second << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal result in correction information file " << second << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      case '0' :
         switch (second) {
            // 01
            case '1' :
               switch (read) {
                  case 'A' :
                     return 'C';
                     break;
                  case 'C' :
                     return 'G';
                     break;
                  case 'G' :
                     return 'T';
                     break;
                  case 'T' :
                     return 'A';
                     break;
                  default :
                     std::cout << std::endl << "ERROR: Illegal character " << read  << " (decode_correction_info)" << std::endl << std::endl;
                     f_log     << std::endl << "ERROR: Illegal character " << read  << " (decode_correction_info)" << std::endl << std::endl;
                     exit(EXIT_FAILURE);
                     break;
               }
               break;
            default :
               std::cout << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      default :
         std::cout << std::endl << "ERROR: Illegal result in correction information file " << first << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal result in correction information file " << first << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }
}



//----------------------------------------------------------------------
// write_corrected_reads
//----------------------------------------------------------------------
void C_correct_errors::write_corrected_reads(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_write_corrected_reads = asctime(localtime(&rawtime));

   if (c_inst_args.paired_read == true) {
      write_corrected_reads_paired_fastq(c_inst_args);
   }
   else {
      write_corrected_reads_single_fastq(c_inst_args);
   }

   time(&rawtime);
   c_inst_time.end_write_corrected_reads = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// write_corrected_reads_single_fastq
//----------------------------------------------------------------------
void C_correct_errors::write_corrected_reads_single_fastq(const C_arg& c_inst_args) {
   // open a log file
   f_log.open(c_inst_args.log_file_name.c_str(), std::fstream::app);

   if (f_log.is_open()) {
      std::cout << "Writing corrected reads into files" << std::endl;

      f_log     << "Writing corrected reads into files" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open error correction information files
   std::ifstream f_error_correction_info;

   if (c_inst_args.verify == true) {
      f_error_correction_info.open(c_inst_args.verified_error_correction_info_file_name.c_str(), std::ios::binary);

      if (f_error_correction_info.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }
   else {
      f_error_correction_info.open(c_inst_args.error_correction_info_file_name.c_str(), std::ios::binary);

      if (f_error_correction_info.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open corrected read files
   std::ofstream f_corrected_read;
   f_corrected_read.open(c_inst_args.corrected_read_file_name.c_str());

   if (f_corrected_read.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read.open(c_inst_args.read_file_name.c_str());

   // process reads
   if (f_read.is_open()) {
      std::string line;

      std::string read;

      // header
      getline(f_read, line);

      // write the first head lines
      if (line.length() > 0) {
         f_corrected_read << line << std::endl;
      }

      // number of bytes for storing reads
      // std::size_t num_byte_per_read;
      // num_byte_per_read = (std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE));

      // std::string modification_buffer(num_byte_per_read, '0');

      while(!f_read.eof()) {
         // DNA sequence
         getline(f_read, read);

         // change sequences to upper case
         transform(read.begin(), read.end(), read.begin(), toupper);

         // substitute Ns other characters
         std::replace(read.begin(), read.end(), 'N', SUBST_CHAR);

         // modified reads
         std::string read_modified(read);

         // read modification information from files
         char buffer;

         if (!f_error_correction_info.get(buffer)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         // make new reads
         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            if ((buffer & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }

            buffer = buffer << 1;

            if ((buffer & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }

            buffer = buffer << 1;

            if ((first != '0') || (second != '0')) {
                read_modified[it_mod] = decode_correction_info(first, second, read[it_mod]);
            }

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info.get(buffer)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         // write new reads to output files
         f_corrected_read << read_modified << std::endl;

         // "+"
         getline(f_read, line);
         f_corrected_read << line << std::endl;

         // quality score
         getline(f_read, line);
         f_corrected_read << line << std::endl;

         // header
         getline(f_read, line);

         if (line.length() > 0) {
            f_corrected_read << line << std::endl;
         }
      }
   }

   // close read files
   f_read.close();

   // close error correction information files
   f_error_correction_info.close();

   // close corrected reads
   f_corrected_read.close();

   std::cout << "     Writing corrected reads into files: done" << std::endl << std::endl;

   f_log     << "     Writing corrected reads into files: done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// write_corrected_reads_paired_fastq
//----------------------------------------------------------------------
void C_correct_errors::write_corrected_reads_paired_fastq(const C_arg& c_inst_args) {
   // open a log file
   f_log.open(c_inst_args.log_file_name.c_str(), std::fstream::app);

   if (f_log.is_open()) {
      std::cout << "Writing corrected reads into files" << std::endl;

      f_log     << "Writing corrected reads into files" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open error correction information files
   std::ifstream f_error_correction_info1;
   std::ifstream f_error_correction_info2;

   if (c_inst_args.verify == true) {
      f_error_correction_info1.open(c_inst_args.verified_error_correction_info_file_name1.c_str(), std::ios::binary);
      f_error_correction_info2.open(c_inst_args.verified_error_correction_info_file_name2.c_str(), std::ios::binary);

      if (f_error_correction_info1.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      if (f_error_correction_info2.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name2 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }
   else {
      f_error_correction_info1.open(c_inst_args.error_correction_info_file_name1.c_str(), std::ios::binary);
      f_error_correction_info2.open(c_inst_args.error_correction_info_file_name2.c_str(), std::ios::binary);

      if (f_error_correction_info1.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      if (f_error_correction_info2.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name2 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open corrected read files
   std::ofstream f_corrected_read1;
   std::ofstream f_corrected_read2;
   f_corrected_read1.open(c_inst_args.corrected_read_file_name1.c_str());
   f_corrected_read2.open(c_inst_args.corrected_read_file_name2.c_str());

   if (f_corrected_read1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_corrected_read2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read_1.open(c_inst_args.read_file_name1.c_str());
   f_read_2.open(c_inst_args.read_file_name2.c_str());

   // process reads
   if (f_read_1.is_open() && f_read_2.is_open()) {
      std::string line_1;
      std::string line_2;

      std::string read_1st;
      std::string read_2nd;

      // header
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      // write the first head lines
      if (line_1.length() > 0) {
         f_corrected_read1 << line_1 << std::endl;
      }
      if (line_2.length() > 0) {
         f_corrected_read2 << line_2 << std::endl;
      }

      // number of bytes for storing reads
      // std::size_t num_byte_per_read;
      // num_byte_per_read = (std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE));

      // std::string modification_buffer1(num_byte_per_read, '0');
      // std::string modification_buffer2(num_byte_per_read, '0');

      while(!f_read_1.eof()) {
         // DNA sequence
         getline(f_read_1, read_1st);
         getline(f_read_2, read_2nd);

         // change sequences to upper case
         transform(read_1st.begin(), read_1st.end(), read_1st.begin(), toupper);
         transform(read_2nd.begin(), read_2nd.end(), read_2nd.begin(), toupper);

         // substitute Ns other characters
         std::replace(read_1st.begin(), read_1st.end(), 'N', SUBST_CHAR);
         std::replace(read_2nd.begin(), read_2nd.end(), 'N', SUBST_CHAR);

         // modified reads
         std::string read_1st_modified(read_1st);
         std::string read_2nd_modified(read_2nd);

         // read modification information from files
         char buffer1;
         char buffer2;

         if (!f_error_correction_info1.get(buffer1)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
         if (!f_error_correction_info2.get(buffer2)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         // make new reads
         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            //----------------------------------------------------------------------
            // forward read
            //----------------------------------------------------------------------
            if ((buffer1 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }

            buffer1 = buffer1 << 1;

            if ((buffer1 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }

            buffer1 = buffer1 << 1;

            if ((first != '0') || (second != '0')) {
                read_1st_modified[it_mod] = decode_correction_info(first, second, read_1st[it_mod]);
            }

            //----------------------------------------------------------------------
            // reverse read
            //----------------------------------------------------------------------
            if ((buffer2 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }

            buffer2 = buffer2 << 1;

            if ((buffer2 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }

            buffer2 = buffer2 << 1;

            if ((first != '0') || (second != '0')) {
                read_2nd_modified[it_mod] = decode_correction_info(first, second, read_2nd[it_mod]);
            }

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info1.get(buffer1)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               if (!f_error_correction_info2.get(buffer2)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         // write new reads to output files
         f_corrected_read1 << read_1st_modified << std::endl;
         f_corrected_read2 << read_2nd_modified << std::endl;

         // "+"
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);
         f_corrected_read1 << line_1 << std::endl;
         f_corrected_read2 << line_2 << std::endl;

         // quality score
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);
         f_corrected_read1 << line_1 << std::endl;
         f_corrected_read2 << line_2 << std::endl;

         // header
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         if (line_1.length() > 0) {
            f_corrected_read1 << line_1 << std::endl;
         }
         if (line_2.length() > 0) {
            f_corrected_read2 << line_2 << std::endl;
         }
      }
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // close error correction information files
   f_error_correction_info1.close();
   f_error_correction_info2.close();

   // close corrected reads
   f_corrected_read1.close();
   f_corrected_read2.close();

   std::cout << "     Writing corrected reads into files: done" << std::endl << std::endl;

   f_log     << "     Writing corrected reads into files: done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// decode_a_char
//----------------------------------------------------------------------
inline char C_correct_errors::decode_a_char(char& in) {
   //----------
   // A: 00
   // C: 01
   // G: 10
   // T: 11
   //----------

   // first character: 1
   if ((in & BIT8) == BIT8) {
      in = in << 1;

      // second character: 1
      if ((in & BIT8) == BIT8) {
         in = in << 1;
         return 'T';
      }
      // second character: 0
      else {
         in = in << 1;
         return 'G';
      }
   }
   // second character: 0
   else {
      in = in << 1;

      // second character: 1
      if ((in & BIT8) == BIT8) {
         in = in << 1;
         return 'C';
      }
      // second character: 0
      else {
         in = in << 1;
         return 'A';
      }

      in = in << 1;
   }
}



//----------------------------------------------------------------------
// query_text
//----------------------------------------------------------------------
inline bool C_correct_errors::query_text(const std::string& kmer, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& kmer_length, const std::size_t& residue_kmer) {
   //--------------------------------------------------
   // reverse complement or not
   //--------------------------------------------------
   std::string kmer_internal;

   switch (kmer[kmer_middle_index]) {
      case 'A' :
         kmer_internal = kmer;
         break;
      case 'C' :
         kmer_internal = kmer;
         break;
      case 'G' :
         reverse_complement(kmer, kmer_internal);
         break;
      case 'T' :
         reverse_complement(kmer, kmer_internal);
         break;
      default :
         std::cout << std::endl << "ERROR: Illegal character " << kmer[kmer_middle_index] << " in the input" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal character " << kmer[kmer_middle_index] << " in the input" << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }

   //--------------------------------------------------
   // encode kmer_internal
   //--------------------------------------------------
   // encode each character in kmer_internal
   std::string encoded_kmer_internal;

   char buffer;
   buffer &= ZERO;

   // encode most characters
   for (std::size_t it_char = 0; it_char < kmer_length; it_char++) {
      encode_a_char(kmer_internal[it_char], buffer);

      if ((it_char % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) {
         // update write_buffer
         encoded_kmer_internal.push_back(buffer);

         // initialize buffer
         buffer &= ZERO;
      }
   }

   // encode remaining characters if needed
   if (residue_kmer != 0) {
      encoded_kmer_internal.push_back(buffer);
   }

   return(c_inst_bloom_filter.query(encoded_kmer_internal));
}



//----------------------------------------------------------------------
// encode_a_char
//----------------------------------------------------------------------
inline void C_correct_errors::encode_a_char(const char& one_nt, char& buffer) {
   //----------
   // A: 00
   // C: 01
   // G: 10
   // T: 11
   //----------

   switch (one_nt) {
      // 00
      case 'A':
         buffer = buffer << 1;
         // buffer |= BIT1;
         buffer = buffer << 1;
         // buffer |= BIT1;
         break;
      // 01
      case 'C':
         buffer = buffer << 1;
         // buffer |= BIT1;
         buffer = buffer << 1;
         buffer |= BIT1;
         break;
      // 10
      case 'G':
         buffer = buffer << 1;
         buffer |= BIT1;
         buffer = buffer << 1;
         // buffer |= BIT1;
         break;
      // 11
      case 'T':
         buffer = buffer << 1;
         buffer |= BIT1;
         buffer = buffer << 1;
         buffer |= BIT1;
         break;
      default:
         std::cout << std::endl << "ERROR: Illegal character " << one_nt << " in read files" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal character " << one_nt << " in read files" << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }
}



//----------------------------------------------------------------------
// mark_bloom_filter_results_org
//----------------------------------------------------------------------
void C_correct_errors::mark_bloom_filter_results_org(const C_arg& c_inst_args, C_time& c_inst_time) {
   if (f_log.is_open()) {
      std::cout << "Marking Bloom filter querying results (original reads)" << std::endl;
      f_log     << "Marking Bloom filter querying results (original reads)" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_mark_unique_solid_kmers_org = asctime(localtime(&rawtime));

   if (c_inst_args.paired_read == true) {
      mark_bloom_filter_results_org_paired_fastq(c_inst_args);
   }
   else {
      mark_bloom_filter_results_org_single_fastq(c_inst_args);
   }

   std::cout <<       "     Marking Bloom filter querying results: done (original reads)" << std::endl << std::endl;
   f_log     <<       "     Marking Bloom filter querying results: done (original reads)" << std::endl << std::endl;

   time(&rawtime);
   c_inst_time.end_mark_unique_solid_kmers_org = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// mark_bloom_filter_results_org_single_fastq
//----------------------------------------------------------------------
void C_correct_errors::mark_bloom_filter_results_org_single_fastq(const C_arg& c_inst_args) {
   // open querying result files
   std::ofstream f_bf_querying_result;
   f_bf_querying_result.open(c_inst_args.bf_querying_result_org_file_name.c_str(), std::ios::binary);

   if (f_bf_querying_result.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_querying_result_org_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read.open(c_inst_args.read_file_name.c_str());

   if(f_read.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // process reads
   std::string line;

   std::string read;

   // header
   getline(f_read, line);

   f_bf_querying_result << line << std::endl;

   std::size_t num_kmers(read_length - c_inst_args.kmer_length + 1);

   while(!f_read.eof()) {
      // DNA sequence
      getline(f_read, read);

      // change sequences to upper case
      transform(read.begin(), read.end(), read.begin(), toupper);

      // substitute Ns other characters
      std::replace(read.begin(), read.end(), 'N', SUBST_CHAR);

      // strings of modification information
      std::string test_result(c_inst_args.kmer_length - 1, 'X');

      // process all the k-mers in reads
      for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
         std::string current_kmer;

         current_kmer = read.substr(it_kmer, c_inst_args.kmer_length);

         // current_kmer_tmp is not erroneous
         if (query_text(current_kmer, c_inst_args.kmer_middle_index, c_inst_args.num_bytes_per_kmer, c_inst_args.kmer_length, c_inst_args.residue_kmer) == true) {
            test_result.push_back('1');
         }
         // current_kmer is erroneous
         else {
            test_result.push_back('0');
         }
      }

      // write original reads to output files
      f_bf_querying_result.write(&test_result[0], read_length);

      f_bf_querying_result << std::endl;

      // "+"
      getline(f_read, line);

      // quality score
      getline(f_read, line);

      // header
      getline(f_read, line);

      if (line.length() > 0) {
         f_bf_querying_result << line << std::endl;
      }
   }

   // close read files
   f_read.close();

   // close error correction information files
   f_bf_querying_result.close();
}



//----------------------------------------------------------------------
// mark_bloom_filter_results_org_paired_fastq
//----------------------------------------------------------------------
void C_correct_errors::mark_bloom_filter_results_org_paired_fastq(const C_arg& c_inst_args) {
   // open querying result files
   std::ofstream f_bf_querying_result1;
   std::ofstream f_bf_querying_result2;
   f_bf_querying_result1.open(c_inst_args.bf_querying_result_org_file_name1.c_str(), std::ios::binary);
   f_bf_querying_result2.open(c_inst_args.bf_querying_result_org_file_name2.c_str(), std::ios::binary);

   if (f_bf_querying_result1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_querying_result_org_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_bf_querying_result2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_querying_result_org_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read_1.open(c_inst_args.read_file_name1.c_str());
   f_read_2.open(c_inst_args.read_file_name2.c_str());

   if(f_read_1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if(f_read_2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // process reads
   std::string line_1;
   std::string line_2;

   std::string read_1st;
   std::string read_2nd;

   // header
   getline(f_read_1, line_1);
   getline(f_read_2, line_2);

   f_bf_querying_result1 << line_1 << std::endl;
   f_bf_querying_result2 << line_2 << std::endl;

   std::size_t num_kmers(read_length - c_inst_args.kmer_length + 1);

   while(!f_read_1.eof()) {
      // DNA sequence
      getline(f_read_1, read_1st);
      getline(f_read_2, read_2nd);

      // change sequences to upper case
      transform(read_1st.begin(), read_1st.end(), read_1st.begin(), toupper);
      transform(read_2nd.begin(), read_2nd.end(), read_2nd.begin(), toupper);

      // substitute Ns other characters
      std::replace(read_1st.begin(), read_1st.end(), 'N', SUBST_CHAR);
      std::replace(read_2nd.begin(), read_2nd.end(), 'N', SUBST_CHAR);

      // strings of modification information
      std::string test_result1(c_inst_args.kmer_length - 1, 'X');
      std::string test_result2(c_inst_args.kmer_length - 1, 'X');

      // process all the k-mers in reads
      for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
         std::string current_kmer;

         //----------------------------------------------------------------------
         // forward read
         //----------------------------------------------------------------------
         current_kmer = read_1st.substr(it_kmer, c_inst_args.kmer_length);

         // current_kmer_tmp is not erroneous
         if (query_text(current_kmer, c_inst_args.kmer_middle_index, c_inst_args.num_bytes_per_kmer, c_inst_args.kmer_length, c_inst_args.residue_kmer) == true) {
            test_result1.push_back('1');
         }
         // current_kmer is erroneous
         else {
            test_result1.push_back('0');
         }

         //----------------------------------------------------------------------
         // reverse read
         //----------------------------------------------------------------------
         current_kmer = read_2nd.substr(it_kmer, c_inst_args.kmer_length);

         // current_kmer_tmp is not erroneous
         if (query_text(current_kmer, c_inst_args.kmer_middle_index, c_inst_args.num_bytes_per_kmer, c_inst_args.kmer_length, c_inst_args.residue_kmer) == true) {
            test_result2.push_back('1');
         }
         // current_kmer is erroneous
         else {
            test_result2.push_back('0');
         }
      }

      // write original reads to output files
      f_bf_querying_result1.write(&test_result1[0], read_length);
      f_bf_querying_result2.write(&test_result2[0], read_length);

      f_bf_querying_result1 << std::endl;
      f_bf_querying_result2 << std::endl;

      // "+"
      getline(f_read_1, line_1);
      getline(f_read_2, line_1);

      // quality score
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      // header
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      if (line_1.length() > 0) {
         f_bf_querying_result1 << line_1 << std::endl;
         f_bf_querying_result2 << line_2 << std::endl;
      }
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // close error correction information files
   f_bf_querying_result1.close();
   f_bf_querying_result2.close();
}



//----------------------------------------------------------------------
// mark_bloom_filter_results_mod
//----------------------------------------------------------------------
void C_correct_errors::mark_bloom_filter_results_mod(const C_arg& c_inst_args, C_time& c_inst_time) {
   if (f_log.is_open()) {
      std::cout << "Marking Bloom filter querying results (corrected reads)" << std::endl;
      f_log     << "Marking Bloom filter querying results (corrected reads)" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_mark_unique_solid_kmers_mod = asctime(localtime(&rawtime));

   if (c_inst_args.paired_read == true) {
      mark_bloom_filter_results_mod_paired_fastq(c_inst_args);
   }
   else {
      mark_bloom_filter_results_mod_single_fastq(c_inst_args);
   }

   std::cout <<       "     Marking Bloom filter querying results: done (corrected reads)" << std::endl << std::endl;
   f_log     <<       "     Marking Bloom filter querying results: done (corrected reads)" << std::endl << std::endl;

   time(&rawtime);
   c_inst_time.end_mark_unique_solid_kmers_mod = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// mark_bloom_filter_results_mod_single_fastq
//----------------------------------------------------------------------
void C_correct_errors::mark_bloom_filter_results_mod_single_fastq(const C_arg& c_inst_args) {
   // open querying result files
   std::ofstream f_bf_querying_result;
   f_bf_querying_result.open(c_inst_args.bf_querying_result_mod_file_name.c_str(), std::ios::binary);

   if (f_bf_querying_result.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_querying_result_mod_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read.open(c_inst_args.corrected_read_file_name.c_str());

   if(f_read.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // process reads
   std::string line;

   std::string read;

   // header
   getline(f_read, line);

   f_bf_querying_result << line << std::endl;

   std::size_t num_kmers(read_length - c_inst_args.kmer_length + 1);

   while(!f_read.eof()) {
      // DNA sequence
      getline(f_read, read);

      // change sequences to upper case
      transform(read.begin(), read.end(), read.begin(), toupper);

      // substitute Ns other characters
      std::replace(read.begin(), read.end(), 'N', SUBST_CHAR);

      // strings of modification information
      std::string test_result(c_inst_args.kmer_length - 1, 'X');

      // process all the k-mers in reads
      for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
         std::string current_kmer;

         current_kmer = read.substr(it_kmer, c_inst_args.kmer_length);

         // current_kmer_tmp is not erroneous
         if (query_text(current_kmer, c_inst_args.kmer_middle_index, c_inst_args.num_bytes_per_kmer, c_inst_args.kmer_length, c_inst_args.residue_kmer) == true) {
            test_result.push_back('1');
         }
         // current_kmer is erroneous
         else {
            test_result.push_back('0');
         }
      }

      // write original reads to output files
      f_bf_querying_result.write(&test_result[0], read_length);

      f_bf_querying_result << std::endl;

      // "+"
      getline(f_read, line);

      // quality score
      getline(f_read, line);

      // header
      getline(f_read, line);

      if (line.length() > 0) {
         f_bf_querying_result << line << std::endl;
      }
   }

   // close read files
   f_read.close();

   // close error correction information files
   f_bf_querying_result.close();
}



//----------------------------------------------------------------------
// mark_bloom_filter_results_mod_paired_fastq
//----------------------------------------------------------------------
void C_correct_errors::mark_bloom_filter_results_mod_paired_fastq(const C_arg& c_inst_args) {
   // open querying result files
   std::ofstream f_bf_querying_result1;
   std::ofstream f_bf_querying_result2;
   f_bf_querying_result1.open(c_inst_args.bf_querying_result_mod_file_name1.c_str(), std::ios::binary);
   f_bf_querying_result2.open(c_inst_args.bf_querying_result_mod_file_name2.c_str(), std::ios::binary);

   if (f_bf_querying_result1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_querying_result_mod_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_bf_querying_result2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_querying_result_mod_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read_1.open(c_inst_args.corrected_read_file_name1.c_str());
   f_read_2.open(c_inst_args.corrected_read_file_name2.c_str());

   if(f_read_1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if(f_read_2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // process reads
   std::string line_1;
   std::string line_2;

   std::string read_1st;
   std::string read_2nd;

   // header
   getline(f_read_1, line_1);
   getline(f_read_2, line_2);

   f_bf_querying_result1 << line_1 << std::endl;
   f_bf_querying_result2 << line_2 << std::endl;

   std::size_t num_kmers(read_length - c_inst_args.kmer_length + 1);

   while(!f_read_1.eof()) {
      // DNA sequence
      getline(f_read_1, read_1st);
      getline(f_read_2, read_2nd);

      // change sequences to upper case
      transform(read_1st.begin(), read_1st.end(), read_1st.begin(), toupper);
      transform(read_2nd.begin(), read_2nd.end(), read_2nd.begin(), toupper);

      // substitute Ns other characters
      std::replace(read_1st.begin(), read_1st.end(), 'N', SUBST_CHAR);
      std::replace(read_2nd.begin(), read_2nd.end(), 'N', SUBST_CHAR);

      // strings of modification information
      std::string test_result1(c_inst_args.kmer_length - 1, 'X');
      std::string test_result2(c_inst_args.kmer_length - 1, 'X');

      // process all the k-mers in reads
      for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
         std::string current_kmer;

         //----------------------------------------------------------------------
         // forward read
         //----------------------------------------------------------------------
         current_kmer = read_1st.substr(it_kmer, c_inst_args.kmer_length);

         // current_kmer_tmp is not erroneous
         if (query_text(current_kmer, c_inst_args.kmer_middle_index, c_inst_args.num_bytes_per_kmer, c_inst_args.kmer_length, c_inst_args.residue_kmer) == true) {
            test_result1.push_back('1');
         }
         // current_kmer is erroneous
         else {
            test_result1.push_back('0');
         }

         //----------------------------------------------------------------------
         // reverse read
         //----------------------------------------------------------------------
         current_kmer = read_2nd.substr(it_kmer, c_inst_args.kmer_length);

         // current_kmer_tmp is not erroneous
         if (query_text(current_kmer, c_inst_args.kmer_middle_index, c_inst_args.num_bytes_per_kmer, c_inst_args.kmer_length, c_inst_args.residue_kmer) == true) {
            test_result2.push_back('1');
         }
         // current_kmer is erroneous
         else {
            test_result2.push_back('0');
         }
      }

      // write original reads to output files
      f_bf_querying_result1.write(&test_result1[0], read_length);
      f_bf_querying_result2.write(&test_result2[0], read_length);

      f_bf_querying_result1 << std::endl;
      f_bf_querying_result2 << std::endl;

      // "+"
      getline(f_read_1, line_1);
      getline(f_read_2, line_1);

      // quality score
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      // header
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      if (line_1.length() > 0) {
         f_bf_querying_result1 << line_1 << std::endl;
         f_bf_querying_result2 << line_2 << std::endl;
      }
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // close error correction information files
   f_bf_querying_result1.close();
   f_bf_querying_result2.close();
}



//----------------------------------------------------------------------
// find_false_positive_candidates
//----------------------------------------------------------------------
void C_correct_errors::find_false_positive_candidates(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_find_false_positives = asctime(localtime(&rawtime));

   if (f_log.is_open()) {
      std::cout << "Finding false positive candidate k-mers" << std::endl;
      f_log     << "Finding false positive candidate k-mers" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (c_inst_args.paired_read == true) {
      find_false_positive_candidates_in_paired_fastq(c_inst_args);
   }
   else {
      find_false_positive_candidates_in_single_fastq(c_inst_args);
   }

   std::cout << "     Finding false positive candidate k-mers: done" << std::endl << std::endl;
   f_log     << "     Finding false positive candidate k-mers: done" << std::endl << std::endl;

   time(&rawtime);
   c_inst_time.end_find_false_positives = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// find_false_positive_candidates_in_single_fastq
//----------------------------------------------------------------------
void C_correct_errors::find_false_positive_candidates_in_single_fastq(const C_arg& c_inst_args) {
   // open error correction information files
   std::ifstream f_error_correction_info;
   f_error_correction_info.open(c_inst_args.error_correction_info_file_name.c_str(), std::ios::binary);

   if (f_error_correction_info.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open the error correction direction file
   std::ifstream f_error_correction_direction;

   f_error_correction_direction.open(c_inst_args.error_correction_direction_file_name.c_str());
   if (f_error_correction_direction.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // file names for false positive candidate k-mers
   distributed_false_positive_candidate_file_names.resize(c_inst_args.num_clusters);
   distributed_false_positive_candidate_file_streams.resize(c_inst_args.num_clusters);

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      std::string file_name;
      std::string str;
      std::stringstream str_tmp;

      str_tmp << std::setfill('0') << std::setw(c_inst_args.num_clusters_digit) << it_files;
      str = str_tmp.str();

      file_name = c_inst_args.prefix + "." + str + ".false-positive-candidates";
      distributed_false_positive_candidate_file_names[it_files] = file_name;

      distributed_false_positive_candidate_file_streams[it_files] = new std::ofstream;
      distributed_false_positive_candidate_file_streams[it_files]->open(file_name.c_str(), std::ios::binary);

      if (distributed_false_positive_candidate_file_streams[it_files]->is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open read files
   f_read.open(c_inst_args.read_file_name.c_str());

   // process reads
   if (f_read.is_open()) {
      std::string line;

      std::string read;

      // header
      getline(f_read, line);

      while(!f_read.eof()) {
         // DNA sequence
         getline(f_read, read);

         // change sequences to upper case
         transform(read.begin(), read.end(), read.begin(), toupper);

         // substitute Ns other characters
         std::replace(read.begin(), read.end(), 'N', SUBST_CHAR);

         // containers for corrected base pairs
         std::string corrected_bp(read_length, '0');

         // read modification information from files
         char buffer;

         if (!f_error_correction_info.get(buffer)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         // find corrected base pairs
         for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            if ((buffer & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            buffer = buffer << 1;

            if ((buffer & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            buffer = buffer << 1;

            if ((first != '0') || (second != '0')) {
                // update a corrected read
                corrected_bp[it_mod] = decode_correction_info(first, second, read[it_mod]);
            }

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info.get(buffer)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         // modified reads
         std::string corrected_read(read);

         find_false_positive_candidates_in_read(
                                                c_inst_args.kmer_middle_index,
                                                c_inst_args.kmer_length,
                                                read_length,
                                                read,
                                                corrected_read,
                                                corrected_bp,
                                                c_inst_args.residue_kmer,
                                                c_inst_args.num_clusters,
                                                c_inst_args.num_bytes_per_kmer,
                                                f_error_correction_direction
                                               );

         // "+"
         getline(f_read, line);

         // quality score
         getline(f_read, line);

         // header
         getline(f_read, line);
      }
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // close read files
   f_read.close();

   // close false positive candidate files
   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      distributed_false_positive_candidate_file_streams[it_files]->close();

      delete distributed_false_positive_candidate_file_streams[it_files];
   }

   // close error correction information files
   f_error_correction_info.close();

   // close the error correction direction file
   f_error_correction_direction.close();
}



//----------------------------------------------------------------------
// find_false_positive_candidates_in_paired_fastq
//----------------------------------------------------------------------
void C_correct_errors::find_false_positive_candidates_in_paired_fastq(const C_arg& c_inst_args) {
   // open error correction information files
   std::ifstream f_error_correction_info1;
   std::ifstream f_error_correction_info2;
   f_error_correction_info1.open(c_inst_args.error_correction_info_file_name1.c_str(), std::ios::binary);
   f_error_correction_info2.open(c_inst_args.error_correction_info_file_name2.c_str(), std::ios::binary);

   if (f_error_correction_info1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_error_correction_info2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name2 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open the error correction direction file
   std::ifstream f_error_correction_direction1;
   std::ifstream f_error_correction_direction2;

   f_error_correction_direction1.open(c_inst_args.error_correction_direction_file_name1.c_str());
   if (f_error_correction_direction1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name1 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_error_correction_direction2.open(c_inst_args.error_correction_direction_file_name2.c_str());
   if (f_error_correction_direction2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name2 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // file names for false positive candidate k-mers
   distributed_false_positive_candidate_file_names.resize(c_inst_args.num_clusters);
   distributed_false_positive_candidate_file_streams.resize(c_inst_args.num_clusters);

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      std::string file_name;
      std::string str;
      std::stringstream str_tmp;

      str_tmp << std::setfill('0') << std::setw(c_inst_args.num_clusters_digit) << it_files;
      str = str_tmp.str();

      file_name = c_inst_args.prefix + "." + str + ".false-positive-candidates";
      distributed_false_positive_candidate_file_names[it_files] = file_name;

      distributed_false_positive_candidate_file_streams[it_files] = new std::ofstream;
      distributed_false_positive_candidate_file_streams[it_files]->open(file_name.c_str(), std::ios::binary);

      if (distributed_false_positive_candidate_file_streams[it_files]->is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open read files
   f_read_1.open(c_inst_args.read_file_name1.c_str());
   f_read_2.open(c_inst_args.read_file_name2.c_str());

   // process reads
   if (f_read_1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else if (f_read_2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   else {
      std::string line_1;
      std::string line_2;

      std::string read_1st;
      std::string read_2nd;

      // header
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      while(!f_read_1.eof()) {
         // DNA sequence
         getline(f_read_1, read_1st);
         getline(f_read_2, read_2nd);

         // change sequences to upper case
         transform(read_1st.begin(), read_1st.end(), read_1st.begin(), toupper);
         transform(read_2nd.begin(), read_2nd.end(), read_2nd.begin(), toupper);

         // substitute Ns other characters
         std::replace(read_1st.begin(), read_1st.end(), 'N', SUBST_CHAR);
         std::replace(read_2nd.begin(), read_2nd.end(), 'N', SUBST_CHAR);

         // containers for corrected base pairs
         std::string corrected_bp_1st(read_length, '0');
         std::string corrected_bp_2nd(read_length, '0');

         // read modification information from files
         char buffer1;
         char buffer2;

         if (!f_error_correction_info1.get(buffer1)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
         if (!f_error_correction_info2.get(buffer2)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         // find corrected base pairs
         for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            //----------------------------------------------------------------------
            // forward read
            //----------------------------------------------------------------------
            if ((buffer1 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            buffer1 = buffer1 << 1;

            if ((buffer1 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            buffer1 = buffer1 << 1;

            if ((first != '0') || (second != '0')) {
                // update a corrected read
                corrected_bp_1st[it_mod] = decode_correction_info(first, second, read_1st[it_mod]);
            }

            //----------------------------------------------------------------------
            // reverse read
            //----------------------------------------------------------------------
            if ((buffer2 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            buffer2 = buffer2 << 1;

            if ((buffer2 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            buffer2 = buffer2 << 1;

            if ((first != '0') || (second != '0')) {
                // update a corrected read
                corrected_bp_2nd[it_mod] = decode_correction_info(first, second, read_2nd[it_mod]);
            }

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info1.get(buffer1)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               if (!f_error_correction_info2.get(buffer2)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         // modified reads
         std::string corrected_read_1st(read_1st);
         std::string corrected_read_2nd(read_2nd);

         find_false_positive_candidates_in_read(
                                                c_inst_args.kmer_middle_index,
                                                c_inst_args.kmer_length,
                                                read_length,
                                                read_1st,
                                                corrected_read_1st,
                                                corrected_bp_1st,
                                                c_inst_args.residue_kmer,
                                                c_inst_args.num_clusters,
                                                c_inst_args.num_bytes_per_kmer,
                                                f_error_correction_direction1
                                               );

         find_false_positive_candidates_in_read(
                                                c_inst_args.kmer_middle_index,
                                                c_inst_args.kmer_length,
                                                read_length,
                                                read_2nd,
                                                corrected_read_2nd,
                                                corrected_bp_2nd,
                                                c_inst_args.residue_kmer,
                                                c_inst_args.num_clusters,
                                                c_inst_args.num_bytes_per_kmer,
                                                f_error_correction_direction2
                                               );

         // "+"
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         // quality score
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         // header
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);
      }
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // close false positive candidate files
   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      distributed_false_positive_candidate_file_streams[it_files]->close();

      delete distributed_false_positive_candidate_file_streams[it_files];
   }

   // close error correction information files
   f_error_correction_info1.close();
   f_error_correction_info2.close();

   // close the error correction direction file
   f_error_correction_direction1.close();
   f_error_correction_direction2.close();
}



//----------------------------------------------------------------------
// find_false_positive_candidates_in_read
//----------------------------------------------------------------------
inline void C_correct_errors::find_false_positive_candidates_in_read(const std::size_t& kmer_middle_index, const std::size_t& kmer_length, const std::size_t& read_length, const std::string& read, std::string& corrected_read, const std::string& corrected_bp, const std::size_t& residue_kmer, const std::size_t& num_clusters, const std::size_t& num_bytes_per_kmer, std::ifstream& f_direction) {
   char first_kmer_direction('N');
   char correction_direction;

   std::string kmer;

   // generate a corrected read
   for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
      if (corrected_bp[it_mod] != '0') {
         corrected_read[it_mod] = corrected_bp[it_mod];
      }
   }

   // find candidates
   for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
      // only modified bases
      if (corrected_bp[it_mod] != '0') {
         // get the correction direction information
         correction_direction = f_direction.get();
         if (f_direction.good() != true) {
            std::cout << std::endl << "ERROR: Fail to get a character from the error correction direction file" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Fail to get a character from the error correction direction file" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         // 0 <= it_mod < kmer_length
         if (it_mod < kmer_length) {
            if (correction_direction == '0') {
               if (first_kmer_direction == 'N') {
                  kmer = corrected_read.substr(it_mod, kmer_length);
                  find_false_positive_candidates_in_kmer(kmer, kmer_middle_index, num_clusters, kmer_length, residue_kmer, num_bytes_per_kmer);
                  first_kmer_direction = '0';
               }
               else if (first_kmer_direction == '0') {
                  kmer = corrected_read.substr(it_mod, kmer_length);
                  find_false_positive_candidates_in_kmer(kmer, kmer_middle_index, num_clusters, kmer_length, residue_kmer, num_bytes_per_kmer);
               }
               else if (first_kmer_direction == '1') {
                  std::cout << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               else {
                  std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
            else if (correction_direction == '1') {
               if (first_kmer_direction == 'N') {
                  kmer = corrected_read.substr(0, kmer_length);
                  find_false_positive_candidates_in_kmer(kmer, kmer_middle_index, num_clusters, kmer_length, residue_kmer, num_bytes_per_kmer);
                  first_kmer_direction = '1';
               }
               else if (first_kmer_direction == '0') {
                  std::cout << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               else if (first_kmer_direction == '1') {
                  // do nothing
               }
               else {
                  std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
            else {
            }
         }
         // kmer_length <= it_mod < (read_length - kmer_length)
         else if (it_mod < (read_length - kmer_length)) {
            if (correction_direction == '0') {
               kmer = corrected_read.substr(it_mod, kmer_length);
            }
            else if (correction_direction == '1') {
               kmer = corrected_read.substr((it_mod - kmer_length + 1), kmer_length);
            }
            else {
               std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }

            find_false_positive_candidates_in_kmer(kmer, kmer_middle_index, num_clusters, kmer_length, residue_kmer, num_bytes_per_kmer);
         }
         // (read_length - kmer_length) <= it_mod < read_length
         else {
            if (correction_direction == '0') {
               std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }
            else if (correction_direction == '1') {
               kmer = corrected_read.substr((it_mod - kmer_length + 1), kmer_length);
               find_false_positive_candidates_in_kmer(kmer, kmer_middle_index, num_clusters, kmer_length, residue_kmer, num_bytes_per_kmer);
            }
            else {
               std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// find_false_positive_candidates_in_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::find_false_positive_candidates_in_kmer(const std::string& kmer, const std::size_t& kmer_middle_index, const std::size_t& num_clusters, const std::size_t& kmer_length, const std::size_t& residue_kmer, const std::size_t& num_bytes_per_kmer) {
   std::string kmer_internal;

   // handle reverse complement
   switch (kmer[kmer_middle_index]) {
      case 'A' :
         kmer_internal = kmer;
         break;
      case 'C' :
         kmer_internal = kmer;
         break;
      case 'G' :
         reverse_complement(kmer, kmer_internal);
         break;
      case 'T' :
         reverse_complement(kmer, kmer_internal);
         break;
      default :
         std::cout << std::endl << "ERROR: Illegal character " << kmer[kmer_middle_index] << " in read files" << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }

   // calculate a hash value
   unsigned int original_index(c_inst_hash_functions.djb_hash(kmer_internal));
   std::size_t file_index = original_index % num_clusters;

   //--------------------------------------------------
   // write the k-mer
   //--------------------------------------------------
   // encode each character in a k-mer
   std::vector<char> encoded_kmer;
   char buffer;
   buffer &= ZERO;

   // encode most characters
   for (std::size_t it_char = 0; it_char < kmer_length; it_char++) {
      encode_a_char(kmer_internal[it_char], buffer);

      if ((it_char % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) {
         // update write_buffer
         encoded_kmer.push_back(buffer);

         // initialize buffer
         buffer &= ZERO;
      }
   }

   // encode remaining characters if needed
   if (residue_kmer != 0) {
      encoded_kmer.push_back(buffer);
   }

   // binary writing
   distributed_false_positive_candidate_file_streams[file_index]->write((const char*)&encoded_kmer[0], num_bytes_per_kmer);
}



//----------------------------------------------------------------------
// verify_false_positive_candidates
//----------------------------------------------------------------------
void C_correct_errors::verify_false_positive_candidates(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_count_false_positives = asctime(localtime(&rawtime));

   if (f_log.is_open()) {
      std::cout << "Counting false positive candidate k-mers" << std::endl;

      f_log     << "Counting false positive candidate k-mers" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   verify_false_positive_candidates_fastq(c_inst_args);

   std::cout << "     Counting false positive candidate k-mers: done" << std::endl << std::endl;
   f_log     << "     Counting false positive candidate k-mers: done" << std::endl << std::endl;

   time(&rawtime);
   c_inst_time.end_count_false_positives = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// verify_false_positive_candidates_fastq
//----------------------------------------------------------------------
void C_correct_errors::verify_false_positive_candidates_fastq(const C_arg& c_inst_args) {
   // file names for false positive candidate k-mers
   distributed_false_positive_result_file_names.resize(c_inst_args.num_clusters);

   // initialize the number of corrected base pairs and reads
   num_wrongly_corrected_errors = 0;

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      std::string file_name;
      std::string str;
      std::stringstream str_tmp;

      str_tmp << std::setfill('0') << std::setw(c_inst_args.num_clusters_digit) << it_files;
      str = str_tmp.str();

      file_name = c_inst_args.prefix + "." + str + ".false-positive-resutls";
      distributed_false_positive_result_file_names[it_files] = file_name;
   }

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      verify_false_positive_candidates_cluster(it_files, c_inst_args.num_bytes_per_kmer);
   }
}



//----------------------------------------------------------------------
// verify_false_positive_candidates_cluster
//----------------------------------------------------------------------
void C_correct_errors::verify_false_positive_candidates_cluster(const std::size_t& file_index, const std::size_t& num_bytes_per_kmer) {
   //--------------------------------------------------
   // add all the unique solid k-mers in the hash table
   //--------------------------------------------------
   std::ifstream f_unique_solid_kmers;

   f_unique_solid_kmers.open(distributed_unique_solid_kmer_file_names[file_index].c_str(), std::ios::binary);
   if (f_unique_solid_kmers.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   char buffer;
   KEY_TYPE kmer;

   google::sparse_hash_map<boost::array<char, MAX_BOOST_ARRAY_SIZE>, unsigned char> gs_unique_solid_kmers;

   std::size_t num_kmers(0);

   f_unique_solid_kmers.get(buffer);

   while (!f_unique_solid_kmers.eof()) {
      // read a k-mer
      kmer[0] = buffer;
      for (std::size_t it_byte = 1; it_byte < num_bytes_per_kmer; it_byte++) {
         f_unique_solid_kmers.get(buffer);

         if (f_unique_solid_kmers.good()) {
            kmer[it_byte] = buffer;
         }
         else {
            std::cout << std::endl << "ERROR: Number of bytes in " << distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Number of bytes in " << distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

      }

      // initialize padding area
      for (std::size_t it_byte = num_bytes_per_kmer; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
         buffer &= ZERO;
         kmer[it_byte] = buffer;
      }

      f_unique_solid_kmers.get(buffer);

      // add kmer to the google sparse hash
      google::sparse_hash_map<KEY_TYPE, unsigned char>::iterator it_hash;
      it_hash = gs_unique_solid_kmers.find(kmer);

      // kmer is not in the hash table: add it
      if (it_hash == gs_unique_solid_kmers.end()) {
         gs_unique_solid_kmers[kmer] = 1;
         num_kmers++;
      }
   }

   // check the number of programmed k-mers
   if (num_kmers != distributed_num_unique_solid_kmers[file_index]) {
      std::cout << std::endl << "ERROR: Number of unique solid k-mers in not matched" << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Number of unique solid k-mers in not matched" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_unique_solid_kmers.close();

   //--------------------------------------------------
   // check all the false positive candidate k-mers
   //--------------------------------------------------
   // false positive candidate k-mers
   std::ifstream f_false_positive_candidates;

   f_false_positive_candidates.open(distributed_false_positive_candidate_file_names[file_index].c_str(), std::ios::binary);
   if (f_false_positive_candidates.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << distributed_false_positive_candidate_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // querying result file
   std::ofstream f_false_positive_results;

   f_false_positive_results.open(distributed_false_positive_result_file_names[file_index].c_str());
   if (f_false_positive_results.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << distributed_false_positive_result_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // test each false positive candidate k-mer
   f_false_positive_candidates.get(buffer);

   while (!f_false_positive_candidates.eof()) {
      // read a k-mer
      kmer[0] = buffer;
      for (std::size_t it_byte = 1; it_byte < num_bytes_per_kmer; it_byte++) {
         f_false_positive_candidates.get(buffer);

         if (f_false_positive_candidates.good()) {
            kmer[it_byte] = buffer;
         }
         else {
            std::cout << std::endl << "ERROR: Number of bytes in " << distributed_false_positive_candidate_file_names[file_index] << " is illegal" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Number of bytes in " << distributed_false_positive_candidate_file_names[file_index] << " is illegal" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

      }

      // initialize padding area
      for (std::size_t it_byte = num_bytes_per_kmer; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
         buffer &= ZERO;
         kmer[it_byte] = buffer;
      }

      f_false_positive_candidates.get(buffer);

      // query kmer to the google sparse hash
      google::sparse_hash_map<KEY_TYPE, unsigned char>::iterator it_hash;
      it_hash = gs_unique_solid_kmers.find(kmer);
      // kmer is not in the hash table: false positive
      if (it_hash == gs_unique_solid_kmers.end()) {
         f_false_positive_results << 'F';
         num_wrongly_corrected_errors++;
      }
      else {
         f_false_positive_results << 'T';
      }
   }

   f_false_positive_candidates.close();
   f_false_positive_results.close();

   gs_unique_solid_kmers.clear();
}



//----------------------------------------------------------------------
// correct_false_positive_candidates
//----------------------------------------------------------------------
void C_correct_errors::correct_false_positives(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_correct_false_positives = asctime(localtime(&rawtime));

   if (f_log.is_open()) {
      std::cout << "Correcting false positive k-mers" << std::endl;

      f_log     << "Correcting false positive k-mers" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (c_inst_args.paired_read == true) {
      correct_false_positives_paired_fastq(c_inst_args);
   }
   else {
      correct_false_positives_single_fastq(c_inst_args);
   }

   std::cout << "     Number of wrong corrections: " << num_wrongly_corrected_errors << std::endl;
   std::cout << "     Correcting false positive candidate k-mers: done" << std::endl << std::endl;

   f_log     << "     Number of wrong corrections: " << num_wrongly_corrected_errors << std::endl;
   f_log     << "     Correcting false positive candidate k-mers: done" << std::endl << std::endl;
   f_log.close();

   time(&rawtime);
   c_inst_time.end_correct_false_positives = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// correct_false_positives_single_fastq
//----------------------------------------------------------------------
void C_correct_errors::correct_false_positives_single_fastq(const C_arg& c_inst_args) {
   // open error correction information files (read)
   std::ifstream f_error_correction_info;
   f_error_correction_info.open(c_inst_args.error_correction_info_file_name.c_str(), std::ios::binary);

   if (f_error_correction_info.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open verified error correction information files (write)
   std::ofstream f_verified_error_correction;
   f_verified_error_correction.open(c_inst_args.verified_error_correction_info_file_name.c_str(), std::ios::binary);

   if (f_verified_error_correction.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // initialize the number of corrected base pairs
   num_wrongly_corrected_errors_check = 0;

   // open error correction information files (text file)
   std::ofstream f_verified_error_correction_txt;

   // open all the false positive test result files
   distributed_false_positive_result_file_streams.resize(c_inst_args.num_clusters);

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {

      distributed_false_positive_result_file_streams[it_files] = new std::ifstream;
      distributed_false_positive_result_file_streams[it_files]->open(distributed_false_positive_result_file_names[it_files].c_str(), std::ios::binary);

      if (distributed_false_positive_result_file_streams[it_files]->is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << distributed_false_positive_result_file_names[it_files] << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << distributed_false_positive_result_file_names[it_files] << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open the error correction direction file
   std::ifstream f_error_correction_direction;
   f_error_correction_direction.open(c_inst_args.error_correction_direction_file_name.c_str());

   if (f_error_correction_direction.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (c_inst_args.debug == true) {
      f_verified_error_correction_txt.open(c_inst_args.verified_error_correction_info_txt_file_name.c_str());

      if (f_verified_error_correction_txt.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_txt_file_name << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_txt_file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open read files
   f_read.open(c_inst_args.read_file_name.c_str());

   // process reads
   if (f_read.is_open()) {
      std::string line;

      std::string read;
      std::string read_tmp;

      // header
      getline(f_read, line);

      if (c_inst_args.debug == true) {
         f_verified_error_correction_txt << line << std::endl;
      }

      // number of bytes for storing reads
      std::size_t num_byte_per_read;
      num_byte_per_read = (std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE));

      // calculate residues
      std::size_t residue_read(read_length % BPS_PER_BYTE);

      while(!f_read.eof()) {
         // DNA sequence
         getline(f_read, read);

         // change sequences to upper case
         transform(read.begin(), read.end(), read.begin(), toupper);

         // substitute Ns other characters
         std::replace(read.begin(), read.end(), 'N', SUBST_CHAR);

         // containers for corrected base pairs
         std::string corrected_bp(read_length, '0');

         // read modification information from files
         char correction_info_buffer;

         if (!f_error_correction_info.get(correction_info_buffer)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         //----------------------------------------------------------------------
         // find corrected base pairs
         //----------------------------------------------------------------------
         for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            if ((correction_info_buffer & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            correction_info_buffer = correction_info_buffer << 1;

            if ((correction_info_buffer & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            correction_info_buffer = correction_info_buffer << 1;

            if ((first != '0') || (second != '0')) {
                corrected_bp[it_mod] = decode_correction_info(first, second, read[it_mod]);
            }

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info.get(correction_info_buffer)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         //----------------------------------------------------------------------
         // check corrected reads
         //----------------------------------------------------------------------
         std::string corrected_read(read);
         std::string corrected_bp_verified(read_length, '0');

         correct_false_positives_read(
                                      c_inst_args.kmer_middle_index,
                                      c_inst_args.kmer_length,
                                      read,
                                      corrected_read,
                                      corrected_bp,
                                      corrected_bp_verified,
                                      c_inst_args.debug,
                                      c_inst_args.num_clusters,
                                      f_error_correction_direction
                                     );

         //----------------------------------------------------------------------
         // write error correction information to files
         //----------------------------------------------------------------------
         std::vector<char> write_buffer;

         char buffer;

         // initialize buffers
         buffer &= ZERO;

         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            encode_correction_info(buffer, read[it_mod], corrected_bp_verified[it_mod]);

            if ((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) {
               // update write_buffer*
               write_buffer.push_back(buffer);

               // initialize buffer
               buffer &= ZERO;
            }
         }

         // read_length is a multiple of BPS_PER_BYTE
         // do nothing
         if (residue_read == 0) {
         }
         // read_length is not a multiple of BPS_PER_BYTE
         else {
            for (std::size_t it_fill = 0; it_fill < (BPS_PER_BYTE - residue_read); it_fill++) {
               buffer = buffer << 1;
               buffer = buffer << 1;
            }

            // fill the last entry
            write_buffer.push_back(buffer);
         }

         // write vectors to output files
         f_verified_error_correction.write((const char*)&write_buffer[0], num_byte_per_read);

         // initialize vectors
         write_buffer.clear();

         if (c_inst_args.debug == true) {
            f_verified_error_correction_txt << corrected_bp_verified << std::endl;
         }

         // "+"
         getline(f_read, line);

         // quality score
         getline(f_read, line);

         // header
         getline(f_read, line);

         if ((c_inst_args.debug == true) && (line.length() > 0)) {
            f_verified_error_correction_txt << line << std::endl;
         }
      }
   }

   // check the number of corrected false positives
   if (num_wrongly_corrected_errors != num_wrongly_corrected_errors_check) {
      std::cout << std::endl << "ERROR: The number of detected false positives is not matched: " << num_wrongly_corrected_errors << " vs " << num_wrongly_corrected_errors_check << std::endl << std::endl;
      f_log     << std::endl << "ERROR: The number of detected false positives is not matched: " << num_wrongly_corrected_errors << " vs " << num_wrongly_corrected_errors_check << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // close read files
   f_read.close();

   // close error correction information files
   f_error_correction_info.close();

   // close error correction information files
   f_verified_error_correction.close();

   if (c_inst_args.debug == true) {
      f_verified_error_correction_txt.close();
   }

   // close the error correction direction file
   f_error_correction_direction.close();

   // destroy the allocated memory for the false positive check result file names
   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      delete distributed_false_positive_result_file_streams[it_files];
   }
}



//----------------------------------------------------------------------
// correct_false_positives_paired_fastq
//----------------------------------------------------------------------
void C_correct_errors::correct_false_positives_paired_fastq(const C_arg& c_inst_args) {
   // open error correction information files (read)
   std::ifstream f_error_correction_info1;
   std::ifstream f_error_correction_info2;
   f_error_correction_info1.open(c_inst_args.error_correction_info_file_name1.c_str(), std::ios::binary);
   f_error_correction_info2.open(c_inst_args.error_correction_info_file_name2.c_str(), std::ios::binary);

   if (f_error_correction_info1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_error_correction_info2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open verified error correction information files (write)
   std::ofstream f_verified_error_correction1;
   std::ofstream f_verified_error_correction2;
   f_verified_error_correction1.open(c_inst_args.verified_error_correction_info_file_name1.c_str(), std::ios::binary);
   f_verified_error_correction2.open(c_inst_args.verified_error_correction_info_file_name2.c_str(), std::ios::binary);

   if (f_verified_error_correction1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_verified_error_correction2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // initialize the number of corrected base pairs
   num_wrongly_corrected_errors_check = 0;

   // open error correction information files (text file)
   std::ofstream f_verified_error_correction_txt1;
   std::ofstream f_verified_error_correction_txt2;

   // open all the false positive test result files
   distributed_false_positive_result_file_streams.resize(c_inst_args.num_clusters);

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {

      distributed_false_positive_result_file_streams[it_files] = new std::ifstream;
      distributed_false_positive_result_file_streams[it_files]->open(distributed_false_positive_result_file_names[it_files].c_str(), std::ios::binary);

      if (distributed_false_positive_result_file_streams[it_files]->is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << distributed_false_positive_result_file_names[it_files] << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << distributed_false_positive_result_file_names[it_files] << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open the error correction direction file
   std::ifstream f_error_correction_direction1;
   std::ifstream f_error_correction_direction2;
   f_error_correction_direction1.open(c_inst_args.error_correction_direction_file_name1.c_str());
   f_error_correction_direction2.open(c_inst_args.error_correction_direction_file_name2.c_str());

   if (f_error_correction_direction1.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name1 << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
   if (f_error_correction_direction2.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name2 << std::endl << std::endl;
      f_log << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_direction_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (c_inst_args.debug == true) {
      f_verified_error_correction_txt1.open(c_inst_args.verified_error_correction_info_txt_file_name1.c_str());
      f_verified_error_correction_txt2.open(c_inst_args.verified_error_correction_info_txt_file_name2.c_str());

      if (f_verified_error_correction_txt1.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_txt_file_name1 << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_txt_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      if (f_verified_error_correction_txt2.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_txt_file_name2 << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_txt_file_name2 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open read files
   f_read_1.open(c_inst_args.read_file_name1.c_str());
   f_read_2.open(c_inst_args.read_file_name2.c_str());

   // process reads
   if (f_read_1.is_open() && f_read_2.is_open()) {
      std::string line_1;
      std::string line_2;

      std::string read_1st;
      std::string read_2nd;
      std::string read_1st_tmp;
      std::string read_2nd_tmp;

      // header
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      if (c_inst_args.debug == true) {
         f_verified_error_correction_txt1 << line_1 << std::endl;
         f_verified_error_correction_txt2 << line_2 << std::endl;
      }

      // number of bytes for storing reads
      std::size_t num_byte_per_read;
      num_byte_per_read = (std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE));

      // calculate residues
      std::size_t residue_read(read_length % BPS_PER_BYTE);

      while(!f_read_1.eof()) {
         // DNA sequence
         getline(f_read_1, read_1st);
         getline(f_read_2, read_2nd);

         // change sequences to upper case
         transform(read_1st.begin(), read_1st.end(), read_1st.begin(), toupper);
         transform(read_2nd.begin(), read_2nd.end(), read_2nd.begin(), toupper);

         // substitute Ns other characters
         std::replace(read_1st.begin(), read_1st.end(), 'N', SUBST_CHAR);
         std::replace(read_2nd.begin(), read_2nd.end(), 'N', SUBST_CHAR);

         // containers for corrected base pairs
         std::string corrected_bp_1st(read_length, '0');
         std::string corrected_bp_2nd(read_length, '0');

         // read modification information from files
         char correction_info_buffer1;
         char correction_info_buffer2;

         if (!f_error_correction_info1.get(correction_info_buffer1)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
         if (!f_error_correction_info2.get(correction_info_buffer2)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         //----------------------------------------------------------------------
         // find corrected base pairs
         //----------------------------------------------------------------------
         for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            //----------------------------------------------------------------------
            // forward read
            //----------------------------------------------------------------------
            if ((correction_info_buffer1 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            correction_info_buffer1 = correction_info_buffer1 << 1;

            if ((correction_info_buffer1 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            correction_info_buffer1 = correction_info_buffer1 << 1;

            if ((first != '0') || (second != '0')) {
                corrected_bp_1st[it_mod] = decode_correction_info(first, second, read_1st[it_mod]);
            }

            //----------------------------------------------------------------------
            // reverse read
            //----------------------------------------------------------------------
            if ((correction_info_buffer2 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            correction_info_buffer2 = correction_info_buffer2 << 1;

            if ((correction_info_buffer2 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            correction_info_buffer2 = correction_info_buffer2 << 1;

            if ((first != '0') || (second != '0')) {
                corrected_bp_2nd[it_mod] = decode_correction_info(first, second, read_2nd[it_mod]);
            }
            //----------------------------------------------------------------------

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info1.get(correction_info_buffer1)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               if (!f_error_correction_info2.get(correction_info_buffer2)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         //----------------------------------------------------------------------
         // check corrected reads
         //----------------------------------------------------------------------
         std::string corrected_read_1st(read_1st);
         std::string corrected_read_2nd(read_2nd);
         std::string corrected_bp_1st_verified(read_length, '0');
         std::string corrected_bp_2nd_verified(read_length, '0');

         correct_false_positives_read(
                                      c_inst_args.kmer_middle_index,
                                      c_inst_args.kmer_length,
                                      read_1st,
                                      corrected_read_1st,
                                      corrected_bp_1st,
                                      corrected_bp_1st_verified,
                                      c_inst_args.debug,
                                      c_inst_args.num_clusters,
                                      f_error_correction_direction1
                                     );


         correct_false_positives_read(
                                      c_inst_args.kmer_middle_index,
                                      c_inst_args.kmer_length,
                                      read_2nd,
                                      corrected_read_2nd,
                                      corrected_bp_2nd,
                                      corrected_bp_2nd_verified,
                                      c_inst_args.debug,
                                      c_inst_args.num_clusters,
                                      f_error_correction_direction2
                                     );

         //----------------------------------------------------------------------
         // write error correction information to files
         //----------------------------------------------------------------------
         std::vector<char> write_buffer1;
         std::vector<char> write_buffer2;

         char buffer1;
         char buffer2;

         // initialize buffers
         buffer1 &= ZERO;
         buffer2 &= ZERO;

         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            encode_correction_info(buffer1, read_1st[it_mod], corrected_bp_1st_verified[it_mod]);
            encode_correction_info(buffer2, read_2nd[it_mod], corrected_bp_2nd_verified[it_mod]);

            if ((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) {
               // update write_buffer*
               write_buffer1.push_back(buffer1);
               write_buffer2.push_back(buffer2);

               // initialize buffer
               buffer1 &= ZERO;
               buffer2 &= ZERO;
            }
         }

         // read_length is a multiple of BPS_PER_BYTE
         // do nothing
         if (residue_read == 0) {
         }
         // read_length is not a multiple of BPS_PER_BYTE
         else {
            for (std::size_t it_fill = 0; it_fill < (BPS_PER_BYTE - residue_read); it_fill++) {
               buffer1 = buffer1 << 1;
               buffer1 = buffer1 << 1;

               buffer2 = buffer2 << 1;
               buffer2 = buffer2 << 1;
            }

            // fill the last entry
            write_buffer1.push_back(buffer1);
            write_buffer2.push_back(buffer2);
         }

         // write vectors to output files
         f_verified_error_correction1.write((const char*)&write_buffer1[0], num_byte_per_read);
         f_verified_error_correction2.write((const char*)&write_buffer2[0], num_byte_per_read);

         // initialize vectors
         write_buffer1.clear();
         write_buffer2.clear();

         if (c_inst_args.debug == true) {
            f_verified_error_correction_txt1 << corrected_bp_1st_verified << std::endl;
            f_verified_error_correction_txt2 << corrected_bp_2nd_verified << std::endl;
         }

         // "+"
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         // quality score
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         // header
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         if ((c_inst_args.debug == true) && (line_1.length() > 0)) {
            f_verified_error_correction_txt1 << line_1 << std::endl;
            f_verified_error_correction_txt2 << line_2 << std::endl;
         }
      }
   }

   // check the number of corrected false positives
   if (num_wrongly_corrected_errors != num_wrongly_corrected_errors_check) {
      std::cout << std::endl << "ERROR: The number of detected false positives is not matched: " << num_wrongly_corrected_errors << " vs " << num_wrongly_corrected_errors_check << std::endl << std::endl;
      f_log     << std::endl << "ERROR: The number of detected false positives is not matched: " << num_wrongly_corrected_errors << " vs " << num_wrongly_corrected_errors_check << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // close error correction information files
   f_error_correction_info1.close();
   f_error_correction_info2.close();

   // close error correction information files
   f_verified_error_correction1.close();
   f_verified_error_correction2.close();

   if (c_inst_args.debug == true) {
      f_verified_error_correction_txt1.close();
      f_verified_error_correction_txt2.close();
   }

   // close the error correction direction file
   f_error_correction_direction1.close();
   f_error_correction_direction2.close();

   // destroy the allocated memory for the false positive check result file names
   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      delete distributed_false_positive_result_file_streams[it_files];
   }
}



//----------------------------------------------------------------------
// correct_false_positives_read
//----------------------------------------------------------------------
inline void C_correct_errors::correct_false_positives_read(const std::size_t& kmer_middle_index, const std::size_t& kmer_length, const std::string& read, std::string& corrected_read, const std::string& corrected_bp, std::string& corrected_bp_verified, const bool& debug, const std::size_t& num_clusters, std::ifstream& f_direction) {
   std::string kmer;

   char correction_direction;
   char first_kmer_direction('N');

   // generate a corrected read
   for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
      if (corrected_bp[it_mod] != '0') {
         corrected_read[it_mod] = corrected_bp[it_mod];
      }
   }

   // check and correct candidates
   for (std::size_t it_mod = 0; it_mod < read_length; it_mod++) {
      // only modified bases
      if (corrected_bp[it_mod] != '0') {
         // get the correction direction information
         correction_direction = f_direction.get();
         if (f_direction.good() != true) {
            std::cout << std::endl << "ERROR: Fail to get a character from the error correction direction file" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Fail to get a character from the error correction direction file" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         // 0 <= it_mod < kmer_length
         if (it_mod < kmer_length) {
            if (correction_direction == '0') {
               if (first_kmer_direction == 'N') {
                  kmer = corrected_read.substr(it_mod, kmer_length);
                  // true positive
                  if (not_false_positive(kmer, kmer_middle_index, num_clusters) == true) {
                     corrected_bp_verified[it_mod] = corrected_bp[it_mod];
                  }
                  // false positive
                  else {
                     num_wrongly_corrected_errors_check++;

                     if (debug == true) {
                        std::cout << "          False positive k-mer: " << kmer << std::endl;
                        f_log     << "          False positive k-mer: " << kmer << std::endl;
                     }
                  }
                  first_kmer_direction = '0';
               }
               else if (first_kmer_direction == '0') {
                  kmer = corrected_read.substr(it_mod, kmer_length);
                  // true positive
                  if (not_false_positive(kmer, kmer_middle_index, num_clusters) == true) {
                     corrected_bp_verified[it_mod] = corrected_bp[it_mod];
                  }
                  // false positive
                  else {
                     num_wrongly_corrected_errors_check++;

                     if (debug == true) {
                        std::cout << "          False positive k-mer: " << kmer << std::endl;
                        f_log     << "          False positive k-mer: " << kmer << std::endl;
                     }
                  }
               }
               else if (first_kmer_direction == '1') {
                  std::cout << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               else {
                  std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
            else if (correction_direction == '1') {
               if (first_kmer_direction == 'N') {
                  kmer = corrected_read.substr(0, kmer_length);
                  // true positive
                  if (not_false_positive(kmer, kmer_middle_index, num_clusters) == true) {
                     corrected_bp_verified[it_mod] = corrected_bp[it_mod];
                  }
                  // false positive
                  else {
                     num_wrongly_corrected_errors_check++;

                     if (debug == true) {
                        std::cout << "          False positive k-mer: " << kmer << std::endl;
                        f_log     << "          False positive k-mer: " << kmer << std::endl;
                     }
                  }
                  first_kmer_direction = '1';
               }
               else if (first_kmer_direction == '0') {
                  std::cout << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Correction direction mismatch" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               else if (first_kmer_direction == '1') {
                  corrected_bp_verified[it_mod] = corrected_bp[it_mod];
               }
               else {
                  std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
            else {
            }
         }
         // kmer_length <= it_mod < (read_length - kmer_length)
         else if (it_mod < (read_length - kmer_length)) {
            if (correction_direction == '0') {
               kmer = corrected_read.substr(it_mod, kmer_length);
            }
            else if (correction_direction == '1') {
               kmer = corrected_read.substr((it_mod - kmer_length + 1), kmer_length);
            }
            else {
               std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }

            // true positive
            if (not_false_positive(kmer, kmer_middle_index, num_clusters) == true) {
               corrected_bp_verified[it_mod] = corrected_bp[it_mod];
            }
            // false positive
            else {
               num_wrongly_corrected_errors_check++;

               if (debug == true) {
                  std::cout << "          False positive k-mer: " << kmer << std::endl;
                  f_log     << "          False positive k-mer: " << kmer << std::endl;
               }
            }
         }
         // (read_length - kmer_length) <= it_mod < read_length
         else {
            if (correction_direction == '0') {
               std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }
            else if (correction_direction == '1') {
               kmer = corrected_read.substr((it_mod - kmer_length + 1), kmer_length);
               // true positive
               if (not_false_positive(kmer, kmer_middle_index, num_clusters) == true) {
                  corrected_bp_verified[it_mod] = corrected_bp[it_mod];
               }
               // false positive
               else {
                  num_wrongly_corrected_errors_check++;

                  if (debug == true) {
                     std::cout << "          False positive k-mer: " << kmer << std::endl;
                     f_log     << "          False positive k-mer: " << kmer << std::endl;
                  }
               }
            }
            else {
               std::cout << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal correction direction " << correction_direction << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// not_false_positive
//----------------------------------------------------------------------
inline bool C_correct_errors::not_false_positive(const std::string& kmer, const std::size_t& kmer_middle_index, const std::size_t& num_clusters) {
   std::string kmer_internal;

   // handle reverse complement
   switch (kmer[kmer_middle_index]) {
      case 'A' :
         kmer_internal = kmer;
         break;
      case 'C' :
         kmer_internal = kmer;
         break;
      case 'G' :
         reverse_complement(kmer, kmer_internal);
         break;
      case 'T' :
         reverse_complement(kmer, kmer_internal);
         break;
      default :
         std::cout << std::endl << "ERROR: Illegal character " << kmer[kmer_middle_index] << " in read files" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal character " << kmer[kmer_middle_index] << " in read files" << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }

   // calculate a hash value
   unsigned int original_index(c_inst_hash_functions.djb_hash(kmer_internal));
   std::size_t file_index = original_index % num_clusters;

   char buffer;
   if (distributed_false_positive_result_file_streams[file_index]->get(buffer)) {
      if (buffer == 'T') {
         return true;
      }
      else if (buffer == 'F') {
         return false;
      }
      else {
         std::cout << std::endl << "ERROR: Illegal character " << buffer << " in " << distributed_false_positive_result_file_names[file_index] << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal character " << buffer << " in " << distributed_false_positive_result_file_names[file_index] << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }
   else {
      std::cout << std::endl << "ERROR: Fail to read a character from " << distributed_false_positive_result_file_names[file_index] << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Fail to read a character from " << distributed_false_positive_result_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }
}



//----------------------------------------------------------------------
// remove_kmer_files
//----------------------------------------------------------------------
void C_correct_errors::remove_kmer_files(const std::size_t& num_clusters, const bool& verify) {
   for (std::size_t it_files = 0; it_files < num_clusters; it_files++) {
      remove(distributed_unique_solid_kmer_file_names[it_files].c_str());

      if (verify == true) {
         remove(distributed_false_positive_candidate_file_names[it_files].c_str());
         remove(distributed_false_positive_result_file_names[it_files].c_str());
      }
   }
}



//----------------------------------------------------------------------
// remove_error_correction_info_files
//----------------------------------------------------------------------
void C_correct_errors::remove_error_correction_info_files(const C_arg& c_inst_args, const bool& verify) {
   if (c_inst_args.paired_read == true) {
      remove(c_inst_args.error_correction_info_file_name1.c_str());
      remove(c_inst_args.error_correction_info_file_name2.c_str());
   }
   else {
      remove(c_inst_args.error_correction_info_file_name.c_str());
   }

   if (verify == true) {
      if (c_inst_args.paired_read == true) {
         remove(c_inst_args.verified_error_correction_info_file_name1.c_str());
         remove(c_inst_args.verified_error_correction_info_file_name2.c_str());

         remove(c_inst_args.error_correction_direction_file_name1.c_str());
         remove(c_inst_args.error_correction_direction_file_name2.c_str());
      }
      else {
         remove(c_inst_args.verified_error_correction_info_file_name.c_str());

         remove(c_inst_args.error_correction_direction_file_name.c_str());
      }
   }
}



//----------------------------------------------------------------------
// write_tef
//----------------------------------------------------------------------
void C_correct_errors::write_tef(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_write_tef = asctime(localtime(&rawtime));

   if (c_inst_args.paired_read == true) {
      write_tef_paired_fastq(c_inst_args);
   }
   else {
      write_tef_single_fastq(c_inst_args);
   }

   time(&rawtime);
   c_inst_time.end_write_tef = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// write_tef_single_fastq
//----------------------------------------------------------------------
void C_correct_errors::write_tef_single_fastq(const C_arg& c_inst_args) {
   // open a log file
   std::ofstream f_log;
   f_log.open(c_inst_args.log_file_name.c_str(), std::fstream::app);

   if (f_log.is_open()) {
      std::cout << "Writing a TEF file" << std::endl;

      f_log     << "Writing a TEF file" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open error correction information files
   std::ifstream f_error_correction_info;

   if (c_inst_args.verify == true) {
      f_error_correction_info.open(c_inst_args.verified_error_correction_info_file_name.c_str(), std::ios::binary);

      if (f_error_correction_info.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }
   else {
      f_error_correction_info.open(c_inst_args.error_correction_info_file_name.c_str(), std::ios::binary);

      if (f_error_correction_info.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open a TEF file
   std::ofstream f_tef;
   f_tef.open(c_inst_args.tef_file_name.c_str());

   if (f_tef.is_open()) {
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.tef_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read.open(c_inst_args.read_file_name.c_str());

   // process reads
   if (f_read.is_open()) {
      std::string line;

      std::string read;

      // header
      getline(f_read, line);

      // take the first word in the header line
      // an initial value of index is 1 to remove "@"
      std::size_t index(1);
      std::string read_id;
      while ((index < line.length()) && (line[index] != ' ') && (line[index] != '\t') && (line[index] != '\n')) {
         read_id.push_back(line[index]);
         index++;
      }

      while(!f_read.eof()) {
         // DNA sequence
         getline(f_read, read);

         // change sequences to upper case
         transform(read.begin(), read.end(), read.begin(), toupper);

         // substitute Ns other characters
         std::replace(read.begin(), read.end(), 'N', SUBST_CHAR);

         // read modification information from files
         char buffer;

         if (!f_error_correction_info.get(buffer)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         std::vector<std::size_t> position;
         std::vector<char> true_bp;
         std::vector<char> wrong_bp;

         // collect correction information
         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            if ((buffer & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            buffer = buffer << 1;

            if ((buffer & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            buffer = buffer << 1;

            if ((first != '0') || (second != '0')) {
               position.push_back(it_mod);
               wrong_bp.push_back(nt_to_num(read[it_mod]));
               true_bp.push_back(decode_correction_info_num(first, second, read[it_mod]));
            }

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info.get(buffer)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         // write a tef information
         if (position.size() != 0) {
            // read identifier
            f_tef << read_id << " " << position.size();

            // errors
            for (std::size_t it_position = 0; it_position < position.size(); it_position++) {
               f_tef << " " << position[it_position] << " " << true_bp[it_position] << " " << wrong_bp[it_position] << " 0";
            }

            f_tef << std::endl;
         }

         // "+"
         getline(f_read, line);

         // quality score
         getline(f_read, line);

         // header
         getline(f_read, line);

         // an initial value of index is 1 to remove "@"
         index = 1;
         read_id.clear();
         while ((index < line.length()) && (line[index] != ' ') && (line[index] != '\t') && (line[index] != '\n')) {
            read_id.push_back(line[index]);
            index++;
         }
      }
   }

   // close read files
   f_read.close();

   // close error correction information files
   f_error_correction_info.close();

   // close corrected reads
   f_tef.close();

   std::cout << "     Writing a TEF file: done" << std::endl << std::endl;

   f_log     << "     Writing a TEF file: done" << std::endl << std::endl;
   f_log.close();
}



//----------------------------------------------------------------------
// write_tef_paired_fastq
//----------------------------------------------------------------------
void C_correct_errors::write_tef_paired_fastq(const C_arg& c_inst_args) {
   // open a log file
   std::ofstream f_log;
   f_log.open(c_inst_args.log_file_name.c_str(), std::fstream::app);

   if (f_log.is_open()) {
      std::cout << "Writing a TEF file" << std::endl;

      f_log     << "Writing a TEF file" << std::endl;
   } else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open error correction information files
   std::ifstream f_error_correction_info1;
   std::ifstream f_error_correction_info2;

   if (c_inst_args.verify == true) {
      f_error_correction_info1.open(c_inst_args.verified_error_correction_info_file_name1.c_str(), std::ios::binary);
      f_error_correction_info2.open(c_inst_args.verified_error_correction_info_file_name2.c_str(), std::ios::binary);

      if (f_error_correction_info1.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      if (f_error_correction_info2.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.verified_error_correction_info_file_name2 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }
   else {
      f_error_correction_info1.open(c_inst_args.error_correction_info_file_name1.c_str(), std::ios::binary);
      f_error_correction_info2.open(c_inst_args.error_correction_info_file_name2.c_str(), std::ios::binary);

      if (f_error_correction_info1.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name1 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
      if (f_error_correction_info2.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.error_correction_info_file_name2 << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // open a TEF file
   std::ofstream f_tef;
   f_tef.open(c_inst_args.tef_file_name.c_str());

   if (f_tef.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.tef_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // open read files
   f_read_1.open(c_inst_args.read_file_name1.c_str());
   f_read_2.open(c_inst_args.read_file_name2.c_str());

   // process reads
   if (f_read_1.is_open() && f_read_2.is_open()) {
      std::string line_1;
      std::string line_2;

      std::string read_1st;
      std::string read_2nd;

      // header
      getline(f_read_1, line_1);
      getline(f_read_2, line_2);

      // take the first word in the header line
      // an initial value of index is 1 to remove "@"
      std::size_t index1(1);
      std::string read_id1;
      while ((index1 < line_1.length()) && (line_1[index1] != ' ') && (line_1[index1] != '\t') && (line_1[index1] != '\n')) {
         read_id1.push_back(line_1[index1]);
         index1++;
      }

      // an initial value of index is 1 to remove "@"
      std::size_t index2(1);
      std::string read_id2;
      while ((index2 < line_2.length()) && (line_2[index2] != ' ') && (line_2[index2] != '\t') && (line_2[index2] != '\n')) {
         read_id2.push_back(line_2[index2]);
         index2++;
      }

      while(!f_read_1.eof()) {
         // DNA sequence
         getline(f_read_1, read_1st);
         getline(f_read_2, read_2nd);

         // change sequences to upper case
         transform(read_1st.begin(), read_1st.end(), read_1st.begin(), toupper);
         transform(read_2nd.begin(), read_2nd.end(), read_2nd.begin(), toupper);

         // substitute Ns other characters
         std::replace(read_1st.begin(), read_1st.end(), 'N', SUBST_CHAR);
         std::replace(read_2nd.begin(), read_2nd.end(), 'N', SUBST_CHAR);

         // read modification information from files
         char buffer1;
         char buffer2;

         if (!f_error_correction_info1.get(buffer1)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
         if (!f_error_correction_info2.get(buffer2)) {
            std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         std::vector<std::size_t> position_1st;
         std::vector<std::size_t> position_2nd;
         std::vector<char> true_bp_1st;
         std::vector<char> true_bp_2nd;
         std::vector<char> wrong_bp_1st;
         std::vector<char> wrong_bp_2nd;

         // collect correction information
         std::size_t it_mod;
         for (it_mod = 0; it_mod < read_length; it_mod++) {
            char first;
            char second;

            //----------------------------------------------------------------------
            // forward read
            //----------------------------------------------------------------------
            if ((buffer1 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }
            buffer1 = buffer1 << 1;

            if ((buffer1 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }
            buffer1 = buffer1 << 1;

            if ((first != '0') || (second != '0')) {
               position_1st.push_back(it_mod);
               wrong_bp_1st.push_back(nt_to_num(read_1st[it_mod]));
               true_bp_1st.push_back(decode_correction_info_num(first, second, read_1st[it_mod]));
            }

            //----------------------------------------------------------------------
            // reverse read
            //----------------------------------------------------------------------
            if ((buffer2 & BIT8) == BIT8) {
               first = '1';
            }
            else {
               first = '0';
            }

            buffer2 = buffer2 << 1;

            if ((buffer2 & BIT8) == BIT8) {
               second = '1';
            }
            else {
               second = '0';
            }

            buffer2 = buffer2 << 1;

            if ((first != '0') || (second != '0')) {
               position_2nd.push_back(it_mod);
               wrong_bp_2nd.push_back(nt_to_num(read_2nd[it_mod]));
               true_bp_2nd.push_back(decode_correction_info_num(first, second, read_2nd[it_mod]));
            }

            // increment indexes
            if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (read_length - 1))) {
               if (!f_error_correction_info1.get(buffer1)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name1 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               if (!f_error_correction_info2.get(buffer2)) {
                  std::cout << std::endl << "ERROR: The size of " << c_inst_args.error_correction_info_file_name2 << " is wrong" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }

         // write a tef information
         if (position_1st.size() != 0) {
            // read identifier
            f_tef << read_id1 << " " << position_1st.size();

            // errors
            for (std::size_t it_position = 0; it_position < position_1st.size(); it_position++) {
               f_tef << " " << position_1st[it_position] << " " << true_bp_1st[it_position] << " " << wrong_bp_1st[it_position] << " 0";
            }

            f_tef << std::endl;
         }

         if (position_2nd.size() != 0) {
            // read identifier
            f_tef << read_id2 << " " << position_2nd.size();

            // errors
            for (std::size_t it_position = 0; it_position < position_2nd.size(); it_position++) {
               f_tef << " " << position_2nd[it_position] << " " << true_bp_2nd[it_position] << " " << wrong_bp_2nd[it_position] << " 0";
            }

            f_tef << std::endl;
         }

         // "+"
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         // quality score
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         // header
         getline(f_read_1, line_1);
         getline(f_read_2, line_2);

         // an initial value of index is 1 to remove "@"
         index1 = 1;
         read_id1.clear();
         while ((index1 < line_1.length()) && (line_1[index1] != ' ') && (line_1[index1] != '\t') && (line_1[index1] != '\n')) {
            read_id1.push_back(line_1[index1]);
            index1++;
         }

         // an initial value of index is 1 to remove "@"
         index2 = 1;
         read_id2.clear();
         while ((index2 < line_2.length()) && (line_2[index2] != ' ') && (line_2[index2] != '\t') && (line_2[index2] != '\n')) {
            read_id2.push_back(line_2[index2]);
            index2++;
         }
      }
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // close error correction information files
   f_error_correction_info1.close();
   f_error_correction_info2.close();

   // close corrected reads
   f_tef.close();

   std::cout << "     Writing a TEF file: done" << std::endl << std::endl;

   f_log     << "     Writing a TEF file: done" << std::endl << std::endl;
   f_log.close();
}



//----------------------------------------------------------------------
// decode_correction_info_num
//----------------------------------------------------------------------
inline char C_correct_errors::decode_correction_info_num(const char& first, const char& second, const char& read) {
   switch (first) {
      case '1' :
         switch (second) {
            // 11
            case '1' :
               switch (read) {
                  case 'A' :
                     return '3';
                     break;
                  case 'C' :
                     return '0';
                     break;
                  case 'G' :
                     return '1';
                     break;
                  case 'T' :
                     return '2';
                     break;
                  default :
                     std::cout << std::endl << "ERROR: Illegal character " << read  << " in the read file" << std::endl << std::endl;
                     f_log     << std::endl << "ERROR: Illegal character " << read  << " in the read file" << std::endl << std::endl;
                     exit(EXIT_FAILURE);
                     break;
               }
               break;
            // 10
            case '0' :
               switch (read) {
                  case 'A' :
                     return '2';
                     break;
                  case 'C' :
                     return '3';
                     break;
                  case 'G' :
                     return '0';
                     break;
                  case 'T' :
                     return '1';
                     break;
                  default :
                     std::cout << std::endl << "ERROR: Illegal character " << read  << " in the read file" << std::endl << std::endl;
                     f_log     << std::endl << "ERROR: Illegal character " << read  << " in the read file" << std::endl << std::endl;
                     exit(EXIT_FAILURE);
                     break;
               }
               break;
            default :
               std::cout << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      case '0' :
         switch (second) {
            // 01
            case '1' :
               switch (read) {
                  case 'A' :
                     return '1';
                     break;
                  case 'C' :
                     return '2';
                     break;
                  case 'G' :
                     return '3';
                     break;
                  case 'T' :
                     return '0';
                     break;
                  default :
                     std::cout << std::endl << "ERROR: Illegal character " << read  << " in the read file" << std::endl << std::endl;
                     f_log     << std::endl << "ERROR: Illegal character " << read  << " in the read file" << std::endl << std::endl;
                     exit(EXIT_FAILURE);
                     break;
               }
               break;
            default :
               std::cout << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
               exit(EXIT_FAILURE);
               break;
         }
         break;
      default :
         std::cout << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal result in correction information file" << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }
}



//----------------------------------------------------------------------
// nt_to_num
//----------------------------------------------------------------------
inline char C_correct_errors::nt_to_num(const char& in_char) {
   switch (in_char) {
      case 'A' :
         return '0';
         break;
      case 'C' :
         return '1';
         break;
      case 'G' :
         return '2';
         break;
      case 'T' :
         return '3';
         break;
      default :
         std::cout << std::endl << "ERROR: Illegal character " << in_char << " (nt_to_num)" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal character " << in_char << " (nt_to_num)" << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }
}



//----------------------------------------------------------------------
// clear_path
//----------------------------------------------------------------------
void C_candidate_path::clear_path() {
   modified_bases.clear();
   sum_qs = 0;
}
