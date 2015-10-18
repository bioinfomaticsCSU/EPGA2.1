#include "count_solid_kmers.hpp"



//----------------------------------------------------------------------
// distribute_kmers
//----------------------------------------------------------------------
void C_count_solid_kmers::distribute_kmers(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_distribute_kmers = asctime(localtime(&rawtime));

   // open a log file
   f_log.open(c_inst_args.log_file_name.c_str(), std::fstream::app);

   if (f_log.is_open()) {
      std::cout << "Distributing k-mers" << std::endl;
      f_log     << "Distributing k-mers" << std::endl;
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.log_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   if (c_inst_args.paired_read == true) {
	   distribute_kmers_paired_fastq(c_inst_args, c_inst_check_reads);
   }
   else {
	   distribute_kmers_single_fastq(c_inst_args, c_inst_check_reads);
   }

   std::cout << "     Distributing k-mers: done" << std::endl << std::endl;
   f_log     << "     Distributing k-mers: done" << std::endl << std::endl;

   time(&rawtime);
   c_inst_time.end_distribute_kmers = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// distribute_kmers_single_fastq
//----------------------------------------------------------------------
void C_count_solid_kmers::distribute_kmers_single_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads) {
   // open distributed k-mer files
   distributed_kmer_file_streams.resize(c_inst_args.num_clusters);
   distributed_num_kmers.resize(c_inst_args.num_clusters);
   distributed_num_unique_kmers.resize(c_inst_args.num_clusters);
   distributed_num_unique_solid_kmers.resize(c_inst_args.num_clusters);

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      // distributed_kmer_file_streams
      std::string file_name;
      std::string str;
      std::stringstream str_tmp;

      str_tmp << std::setfill('0') << std::setw(c_inst_args.num_clusters_digit) << it_files;
      str = str_tmp.str();

      // distributed_kmer_file_names
      file_name = c_inst_args.prefix + "." + str + ".dist-kmer";
      distributed_kmer_file_names.push_back(file_name);

      distributed_kmer_file_streams[it_files] = new std::ofstream;
      distributed_kmer_file_streams[it_files]->open(file_name.c_str(), std::ios::binary);

      if (distributed_kmer_file_streams[it_files]->is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // distributed_unique_kmer_file_names
      file_name = c_inst_args.prefix + "." + str + ".dist-unique-solid-kmer";
      distributed_unique_solid_kmer_file_names.push_back(file_name);

      // distributed_num_kmers
      distributed_num_kmers[it_files] = 0;

      // distributed_num_unique_kmers
      distributed_num_unique_kmers[it_files] = 0;

      // distributed_num_unique_solid_kmers
      distributed_num_unique_solid_kmers[it_files] = 0;
   }

   // initialize the histgram of k-mer occurrences
   num_occurrences_histogram.resize(HISTOGRAM_SIZE);
   for (std::size_t it_histo = 0; it_histo < HISTOGRAM_SIZE; it_histo++) {
      num_occurrences_histogram[it_histo] = 0;
   }

   // open read files
   f_read.open(c_inst_args.read_file_name.c_str());

   // process reads
   if (f_read.is_open()) {
      std::string line;

      std::string read;
      std::string read_rc;

      // header
      getline(f_read, line);
      std::size_t num_kmers(c_inst_check_reads.read_length - c_inst_args.kmer_length + 1);

      // iterate all reads
      while(!f_read.eof()) {
         // DNA sequence
         getline(f_read, read);

         // change sequences to upper case
         transform(read.begin(), read.end(), read.begin(), toupper);

         // substitute Ns other characters
         std::replace(read.begin(), read.end(), 'N', SUBST_CHAR);

         // make reverse complement reads
         reverse_complement(read, read_rc);

         // process (l - k + 1) k-mers in reads
         for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
            std::string current_kmer;
            std::string current_kmer_rc;

            current_kmer = read.substr(it_kmer, c_inst_args.kmer_length);
            current_kmer_rc = read_rc.substr((c_inst_check_reads.read_length - c_inst_args.kmer_length - it_kmer), c_inst_args.kmer_length);

            write_a_kmer(current_kmer, current_kmer_rc, c_inst_args.kmer_middle_index, c_inst_args.kmer_length, c_inst_args.residue_kmer, c_inst_args.num_bytes_per_kmer, c_inst_args.num_clusters);
         }

         // "+"
         getline(f_read, line);

         // quality score
         getline(f_read, line);

         // header
         getline(f_read, line);
      }
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open input read files" << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open input read files" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // close read files
   f_read.close();

   // report number of k-mers in each file
   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      distributed_kmer_file_streams[it_files]->close();

      delete distributed_kmer_file_streams[it_files];

      if (c_inst_args.debug == true) {
         std::cout << "     Cluster " << std::setw(c_inst_args.num_clusters_digit) << it_files << ": " << std::setw(10) << distributed_num_kmers[it_files] << std::endl;
         f_log     << "     Cluster " << std::setw(c_inst_args.num_clusters_digit) << it_files << ": " << std::setw(10) << distributed_num_kmers[it_files] << std::endl;
      }
   }
}



//----------------------------------------------------------------------
// distribute_kmers_paired_fastq
//----------------------------------------------------------------------
void C_count_solid_kmers::distribute_kmers_paired_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads) {
   // open distributed k-mer files
   distributed_kmer_file_streams.resize(c_inst_args.num_clusters);
   distributed_num_kmers.resize(c_inst_args.num_clusters);
   distributed_num_unique_kmers.resize(c_inst_args.num_clusters);
   distributed_num_unique_solid_kmers.resize(c_inst_args.num_clusters);

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      // distributed_kmer_file_streams
      std::string file_name;
      std::string str;
      std::stringstream str_tmp;

      str_tmp << std::setfill('0') << std::setw(c_inst_args.num_clusters_digit) << it_files;
      str = str_tmp.str();

      // distributed_kmer_file_names
      file_name = c_inst_args.prefix + "." + str + ".dist-kmer";
      distributed_kmer_file_names.push_back(file_name);

      distributed_kmer_file_streams[it_files] = new std::ofstream;
      distributed_kmer_file_streams[it_files]->open(file_name.c_str(), std::ios::binary);

      if (distributed_kmer_file_streams[it_files]->is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << file_name << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // distributed_unique_kmer_file_names
      file_name = c_inst_args.prefix + "." + str + ".dist-unique-solid-kmer";
      distributed_unique_solid_kmer_file_names.push_back(file_name);

      // distributed_num_kmers
      distributed_num_kmers[it_files] = 0;

      // distributed_num_unique_kmers
      distributed_num_unique_kmers[it_files] = 0;

      // distributed_num_unique_solid_kmers
      distributed_num_unique_solid_kmers[it_files] = 0;
   }

   // initialize the histgram of k-mer occurrences
   num_occurrences_histogram.resize(HISTOGRAM_SIZE);
   for (std::size_t it_histo = 0; it_histo < HISTOGRAM_SIZE; it_histo++) {
      num_occurrences_histogram[it_histo] = 0;
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
      std::string read_1st_rc;
      std::string read_2nd_rc;

      // header
      getline(f_read_1, line_1);
      getline(f_read_2, line_1);

      std::size_t num_kmers(c_inst_check_reads.read_length - c_inst_args.kmer_length + 1);

      // iterate all reads
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

         // make reverse complement reads
         reverse_complement(read_1st, read_1st_rc);
         reverse_complement(read_2nd, read_2nd_rc);

         // process (l - k + 1) k-mers in reads
         for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
            std::string current_kmer;
            std::string current_kmer_rc;

            //----------------------------------------------------------------------
            // forward read
            //----------------------------------------------------------------------
            current_kmer = read_1st.substr(it_kmer, c_inst_args.kmer_length);
            current_kmer_rc = read_1st_rc.substr((c_inst_check_reads.read_length - c_inst_args.kmer_length - it_kmer), c_inst_args.kmer_length);

            write_a_kmer(current_kmer, current_kmer_rc, c_inst_args.kmer_middle_index, c_inst_args.kmer_length, c_inst_args.residue_kmer, c_inst_args.num_bytes_per_kmer, c_inst_args.num_clusters);

            //----------------------------------------------------------------------
            // reverse read
            //----------------------------------------------------------------------
            current_kmer = read_2nd.substr(it_kmer, c_inst_args.kmer_length);
            current_kmer_rc = read_2nd_rc.substr((c_inst_check_reads.read_length - c_inst_args.kmer_length - it_kmer), c_inst_args.kmer_length);

            write_a_kmer(current_kmer, current_kmer_rc, c_inst_args.kmer_middle_index, c_inst_args.kmer_length, c_inst_args.residue_kmer, c_inst_args.num_bytes_per_kmer, c_inst_args.num_clusters);
         }

         // "+"
         getline(f_read_1, line_1);
         getline(f_read_2, line_1);

         // quality score
         getline(f_read_1, line_1);
         getline(f_read_2, line_1);

         // header
         getline(f_read_1, line_1);
         getline(f_read_2, line_1);
      }
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open input read files" << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open input read files" << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // close read files
   f_read_1.close();
   f_read_2.close();

   // report number of k-mers in each file
   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      distributed_kmer_file_streams[it_files]->close();

      delete distributed_kmer_file_streams[it_files];

      if (c_inst_args.debug == true) {
         std::cout << "     Cluster " << std::setw(c_inst_args.num_clusters_digit) << it_files << ": " << std::setw(10) << distributed_num_kmers[it_files] << std::endl;
         f_log     << "     Cluster " << std::setw(c_inst_args.num_clusters_digit) << it_files << ": " << std::setw(10) << distributed_num_kmers[it_files] << std::endl;
      }
   }
}



//----------------------------------------------------------------------
// count_kmers
//----------------------------------------------------------------------
void C_count_solid_kmers::count_kmers(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_count_kmers = asctime(localtime(&rawtime));

   // open a log file
   std::cout << "Counting the number of k-mers" << std::endl;
   f_log     << "Counting the number of k-mers" << std::endl;

   // open read files
   if (kmer_occurrence_threshold == 0) {
      count_unique_kmers_fastq(c_inst_args, c_inst_check_reads);
      write_histogram(c_inst_args.histo_file_name, c_inst_args.set_kmer_occurrence_threshold);
      count_unique_solid_kmers_fastq(c_inst_args, c_inst_check_reads);
   }
   else {
      count_both_unique_and_solid_kmers_fastq(c_inst_args, c_inst_check_reads);
      write_histogram(c_inst_args.histo_file_name, c_inst_args.set_kmer_occurrence_threshold);
   }

   remove_kmer_files(c_inst_args.num_clusters);

   time(&rawtime);
   c_inst_time.end_count_kmers = asctime(localtime(&rawtime));

   std::cout << "     Number of unique k-mers      : " << std::setw(10) << num_unique_kmers << std::endl;
   std::cout << "     Number of unique solid k-mers: " << std::setw(10) << num_unique_solid_kmers << std::endl;
   std::cout << "     Counting the number of k-mers: done" << std::endl << std::endl;

   f_log     << "     Number of unique k-mers      : " << std::setw(10) << num_unique_kmers << std::endl;
   f_log     << "     Number of unique solid k-mers: " << std::setw(10) << num_unique_solid_kmers << std::endl;
   f_log     << "     Counting the number of k-mers: done" << std::endl << std::endl;
   f_log.close();
}



//----------------------------------------------------------------------
// count_unique_kmers_fastq
//----------------------------------------------------------------------
void C_count_solid_kmers::count_unique_kmers_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads) {
   std::size_t num_empty_characters(BPS_PER_BYTE - c_inst_args.residue_kmer);

   num_unique_kmers = 0;

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      // count k-mers in each cluster
      count_unique_kmers_cluster(it_files, c_inst_args.num_bytes_per_kmer, num_empty_characters, distributed_unique_solid_kmer_file_names[it_files]);

      // count total number of unique (solid) k-mers
      num_unique_kmers += distributed_num_unique_kmers[it_files];
   }
}



//----------------------------------------------------------------------
// count_unique_kmers_cluster
//----------------------------------------------------------------------
void C_count_solid_kmers::count_unique_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_per_kmer, const std::size_t num_empty_characters, const std::string& distributed_unique_solid_kmer_file_names) {
   std::ifstream f_in;
   f_in.open(distributed_kmer_file_names[file_index].c_str(), std::ios::binary);

   std::size_t num_kmers(0);

   if (f_in.is_open()) {
      char buffer;
      KEY_TYPE kmer;

      google::sparse_hash_map<boost::array<char, MAX_BOOST_ARRAY_SIZE>, unsigned short int> gs_count_kmers;

      // program k-mers into the google sparse hash
      f_in.get(buffer);

      while (!f_in.eof()) {
         // k-mer
         kmer[0] = buffer;
         for (std::size_t it_byte = 1; it_byte < num_bytes_per_kmer; it_byte++) {
            f_in.get(buffer);

            if (f_in.good()) {
               kmer[it_byte] = buffer;
            }
            else {
               std::cout << std::endl << "ERROR: Number of bytes in " << distributed_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Number of bytes in " << distributed_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }

         }

         // initialize padding area
         for (std::size_t it_byte = num_bytes_per_kmer; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
            buffer &= ZERO;
            kmer[it_byte] = buffer;
         }

         f_in.get(buffer);

#ifdef DEBUG12042012
         //--------------------------------------------------
         // DEBUG
         // read k-mers from files
         //--------------------------------------------------
         if (file_index == 0) {
            std::cout << "D1: ";
            for (std::size_t it_byte = 0; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
               char char_buf(kmer[it_byte]);

               if (it_byte == (num_bytes_per_kmer - 1)) {
                  std::cout << decode_a_byte(char_buf, num_empty_characters);
               }
               else {
                  std::cout << decode_a_byte(char_buf, 0);
               }
            }
            std::cout << std::endl;
         }
         //--------------------------------------------------
#endif

         // add kmer to the google sparse hash
         google::sparse_hash_map<KEY_TYPE, unsigned short int>::iterator it_hash;
         it_hash = gs_count_kmers.find(kmer);
         // kmer is not in the hash table: add it
         if (it_hash == gs_count_kmers.end()) {
            gs_count_kmers[kmer] = 1;
         }
         // kmer is in the hash table: incrment it
         else {
            gs_count_kmers[kmer]++;
         }

         num_kmers++;
      }

      if (num_kmers != distributed_num_kmers[file_index]) {
         std::cout << std::endl << "ERROR: Number of k-mers in " << distributed_kmer_file_names[file_index] << " is not matched with the original number" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Number of k-mers in " << distributed_kmer_file_names[file_index] << " is not matched with the original number" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // count unique solid k-mers
      distributed_num_unique_kmers[file_index] = gs_count_kmers.size();

      // iterate all the keys in the hash table
      google::sparse_hash_map<KEY_TYPE, unsigned short int>::iterator it_hash;
      for (it_hash = gs_count_kmers.begin(); it_hash != gs_count_kmers.end(); it_hash++) {
#ifdef DEBUG12042012
         //--------------------------------------------------
         // DEBUG
         // read the number of occurrences of each k-mer from the hash table
         //--------------------------------------------------
         if (file_index == 0) {
            std::cout << "D2: ";
            for (std::size_t it_byte = 0; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
               char char_buf(it_hash->first[it_byte]);

               // last byte
               if (it_byte == (num_bytes_per_kmer - 1)) {
                  std::cout << decode_a_byte(char_buf, num_empty_characters);
               }
               // other bytes
               else {
                  std::cout << decode_a_byte(char_buf, 0);
               }
            }

            std::cout << std::setw(10) << it_hash->second;
            std::cout << std::endl;
         }
         //--------------------------------------------------
#endif

         // update the histogram
         std::size_t num_occurrences_tmp;
         if (it_hash->second > (HISTOGRAM_SIZE - 1)) {
            num_occurrences_tmp = HISTOGRAM_SIZE - 1;
         }
         else {
            num_occurrences_tmp = it_hash->second;
         }

         num_occurrences_histogram[num_occurrences_tmp]++;
      }

      // purge the hash table
      gs_count_kmers.clear();

#ifdef DEBUG12042012
      //--------------------------------------------------
      // DEBUG
      // check unique solid k-mer files
      //--------------------------------------------------
      if (file_index == 0) {
         std::ifstream f_tmp;
         f_tmp.open(c_inst_args.distributed_unique_solid_kmer_file_names[file_index].c_str(), std::ios::binary);
         if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         char buffer;
         f_tmp.get(buffer);

         while (!f_tmp.eof()) {
            std::cout << "D4: ";
            for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
               if (!f_tmp.good()) {
                  std::cout << std::endl << "ERROR: Number of bytes in " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Number of bytes in " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               else {
                  char char_buf(buffer);
                  // last byte
                  if (it_byte == (num_bytes_per_kmer - 1)) {
                     std::cout << decode_a_byte(char_buf, num_empty_characters);
                  }
                  // other bytes
                  else {
                     std::cout << decode_a_byte(char_buf, 0);
                  }

                  f_tmp.get(buffer);
               }
            }

            std::cout << std::endl;
         }

         f_tmp.close();
      }
      //--------------------------------------------------
#endif
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << distributed_kmer_file_names[file_index] << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << distributed_kmer_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_in.close();
}



//----------------------------------------------------------------------
// count_unique_solid_kmers_fastq
//----------------------------------------------------------------------
void C_count_solid_kmers::count_unique_solid_kmers_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads) {
   std::size_t num_empty_characters(BPS_PER_BYTE - c_inst_args.residue_kmer);

   num_unique_kmers = 0;
   num_unique_solid_kmers = 0;

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      // count k-mers in each cluster
      count_unique_solid_kmers_cluster(it_files, c_inst_args.num_bytes_per_kmer, num_empty_characters, distributed_unique_solid_kmer_file_names[it_files]);

      // count total number of unique (solid) k-mers
      num_unique_kmers += distributed_num_unique_kmers[it_files];
      num_unique_solid_kmers += distributed_num_unique_solid_kmers[it_files];
   }
}



//----------------------------------------------------------------------
// count_unique_solid_kmers_cluster
//----------------------------------------------------------------------
void C_count_solid_kmers::count_unique_solid_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_per_kmer, const std::size_t num_empty_characters, const std::string& distributed_unique_solid_kmer_file_names) {
   std::ifstream f_in;
   f_in.open(distributed_kmer_file_names[file_index].c_str(), std::ios::binary);

   std::size_t num_kmers(0);

   if (f_in.is_open()) {
      char buffer;
      KEY_TYPE kmer;

      google::sparse_hash_map<boost::array<char, MAX_BOOST_ARRAY_SIZE>, unsigned short int> gs_count_kmers;

      // program k-mers into the google sparse hash
      f_in.get(buffer);

      while (!f_in.eof()) {
         // k-mer
         kmer[0] = buffer;
         for (std::size_t it_byte = 1; it_byte < num_bytes_per_kmer; it_byte++) {
            f_in.get(buffer);

            if (f_in.good()) {
               kmer[it_byte] = buffer;
            }
            else {
               std::cout << std::endl << "ERROR: Number of bytes in " << distributed_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Number of bytes in " << distributed_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }

         }

         // initialize padding area
         for (std::size_t it_byte = num_bytes_per_kmer; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
            buffer &= ZERO;
            kmer[it_byte] = buffer;
         }

         f_in.get(buffer);

#ifdef DEBUG12042012
         //--------------------------------------------------
         // DEBUG
         // read k-mers from files
         //--------------------------------------------------
         /*
         if (file_index == 0) {
            std::cout << "D1: ";
            for (std::size_t it_byte = 0; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
               char char_buf(kmer[it_byte]);

               // last byte
               if (it_byte == (num_bytes_per_kmer - 1)) {
                  std::cout << decode_a_byte(char_buf, num_empty_characters);
               }
               // other bytes
               else {
                  std::cout << decode_a_byte(char_buf, 0);
               }
            }
            std::cout << std::endl;
         }
         */
         //--------------------------------------------------
#endif

         // add kmer to the google sparse hash
         google::sparse_hash_map<KEY_TYPE, unsigned short int>::iterator it_hash;
         it_hash = gs_count_kmers.find(kmer);
         // kmer is not in the hash table: add it
         if (it_hash == gs_count_kmers.end()) {
            gs_count_kmers[kmer] = 1;
         }
         // kmer is in the hash table: incrment it
         else {
            gs_count_kmers[kmer]++;
         }

         num_kmers++;
      }

      if (num_kmers != distributed_num_kmers[file_index]) {
         std::cout << std::endl << "ERROR: Number of k-mers in " << distributed_kmer_file_names[file_index] << " is not matched with the original number" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Number of k-mers in " << distributed_kmer_file_names[file_index] << " is not matched with the original number" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // open a unique solid k-mer file
      std::ofstream f_unique_solid_kmers;
      f_unique_solid_kmers.open(distributed_unique_solid_kmer_file_names.c_str(), std::ios::binary);

      if (f_unique_solid_kmers.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << distributed_unique_solid_kmer_file_names << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << distributed_unique_solid_kmer_file_names << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // iterate all the keys in the hash table
      google::sparse_hash_map<KEY_TYPE, unsigned short int>::iterator it_hash;
      for (it_hash = gs_count_kmers.begin(); it_hash != gs_count_kmers.end(); it_hash++) {
#ifdef DEBUG12042012
         //--------------------------------------------------
         // DEBUG
         // read the number of occurrences of each k-mer from the hash table
         //--------------------------------------------------
         if (file_index == 0) {
            std::cout << "D2: ";
            for (std::size_t it_byte = 0; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
               char char_buf(it_hash->first[it_byte]);

               // last byte
               if (it_byte == (num_bytes_per_kmer - 1)) {
                  std::cout << decode_a_byte(char_buf, num_empty_characters);
               }
               // other bytes
               else {
                  std::cout << decode_a_byte(char_buf, 0);
               }
            }

            std::cout << std::setw(10) << it_hash->second;
            std::cout << std::endl;
         }
         //--------------------------------------------------
#endif

         // count the number of unique solid k-mers
         if (it_hash->second >= kmer_occurrence_threshold) {
            distributed_num_unique_solid_kmers[file_index]++;

            // copy the key to the buffer
            std::vector<char> write_buffer;
            write_buffer.resize(num_bytes_per_kmer);
            for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
               write_buffer[it_byte] = it_hash->first[it_byte];
            }

#ifdef DEBUG12052012
            //--------------------------------------------------
            // DEBUG
            // print unique solid k-mers
            //--------------------------------------------------
            // if (file_index == 0) {
               std::string string_tmp("");
               std::string string_tmp_rc;

               for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
                  char char_buf(it_hash->first[it_byte]);

                  // last byte
                  if (it_byte == (num_bytes_per_kmer - 1)) {
                     string_tmp = string_tmp + decode_a_byte(char_buf, num_empty_characters);
                     // std::cout << decode_a_byte(char_buf, num_empty_characters);
                  }
                  // other bytes
                  else {
                     string_tmp = string_tmp + decode_a_byte(char_buf, 0);
                     // std::cout << decode_a_byte(char_buf, 0);
                  }
               }
               reverse_complement(string_tmp, string_tmp_rc);
               std::cout << "D3: " << string_tmp << " " << string_tmp_rc << std::endl;
            // }
            //--------------------------------------------------
#endif

            // write an unique solid k-mer
            f_unique_solid_kmers.write((const char*)&write_buffer[0], num_bytes_per_kmer);
         }

         // update the histogram
         std::size_t num_occurrences_tmp;
         if (it_hash->second > (HISTOGRAM_SIZE - 1)) {
            num_occurrences_tmp = HISTOGRAM_SIZE - 1;
         }
         else {
            num_occurrences_tmp = it_hash->second;
         }

         num_occurrences_histogram[num_occurrences_tmp]++;
      }

      // purge the hash table
      gs_count_kmers.clear();

      f_unique_solid_kmers.close();

#ifdef DEBUG12042012
      //--------------------------------------------------
      // DEBUG
      // check unique solid k-mer files
      //--------------------------------------------------
      if (file_index == 0) {
         std::ifstream f_tmp;
         f_tmp.open(c_inst_args.distributed_unique_solid_kmer_file_names[file_index].c_str(), std::ios::binary);
         if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         char buffer;
         f_tmp.get(buffer);

         while (!f_tmp.eof()) {
            std::cout << "D4: ";
            for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
               if (!f_tmp.good()) {
                  std::cout << std::endl << "ERROR: Number of bytes in " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Number of bytes in " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               else {
                  char char_buf(buffer);
                  // last byte
                  if (it_byte == (num_bytes_per_kmer - 1)) {
                     std::cout << decode_a_byte(char_buf, num_empty_characters);
                  }
                  // other bytes
                  else {
                     std::cout << decode_a_byte(char_buf, 0);
                  }

                  f_tmp.get(buffer);
               }
            }

            std::cout << std::endl;
         }

         f_tmp.close();
      }
      //--------------------------------------------------
#endif
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << distributed_kmer_file_names[file_index] << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << distributed_kmer_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_in.close();
}



//----------------------------------------------------------------------
// count_both_unique_and_solid_kmers_fastq
//----------------------------------------------------------------------
void C_count_solid_kmers::count_both_unique_and_solid_kmers_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads) {
   std::size_t num_empty_characters(BPS_PER_BYTE - c_inst_args.residue_kmer);

   num_unique_kmers = 0;
   num_unique_solid_kmers = 0;

   for (std::size_t it_files = 0; it_files < c_inst_args.num_clusters; it_files++) {
      // count k-mers in each cluster
      count_both_unique_and_solid_kmers_cluster(it_files, c_inst_args.num_bytes_per_kmer, num_empty_characters, distributed_unique_solid_kmer_file_names[it_files]);

      // count total number of unique (solid) k-mers
      num_unique_kmers += distributed_num_unique_kmers[it_files];
      num_unique_solid_kmers += distributed_num_unique_solid_kmers[it_files];
   }
}



//----------------------------------------------------------------------
// count_both_unique_and_solid_kmers_cluster
//----------------------------------------------------------------------
void C_count_solid_kmers::count_both_unique_and_solid_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_per_kmer, const std::size_t num_empty_characters, const std::string& distributed_unique_solid_kmer_file_names) {
   std::ifstream f_in;
   f_in.open(distributed_kmer_file_names[file_index].c_str(), std::ios::binary);

   std::size_t num_kmers(0);

   if (f_in.is_open()) {
      char buffer;
      KEY_TYPE kmer;

      google::sparse_hash_map<boost::array<char, MAX_BOOST_ARRAY_SIZE>, unsigned short int> gs_count_kmers;

      // program k-mers into the google sparse hash
      f_in.get(buffer);

      while (!f_in.eof()) {
         // k-mer
         kmer[0] = buffer;
         for (std::size_t it_byte = 1; it_byte < num_bytes_per_kmer; it_byte++) {
            f_in.get(buffer);

            if (f_in.good()) {
               kmer[it_byte] = buffer;
            }
            else {
               std::cout << std::endl << "ERROR: Number of bytes in " << distributed_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
               f_log     << std::endl << "ERROR: Number of bytes in " << distributed_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
               exit(EXIT_FAILURE);
            }

         }

         // initialize padding area
         for (std::size_t it_byte = num_bytes_per_kmer; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
            buffer &= ZERO;
            kmer[it_byte] = buffer;
         }

         f_in.get(buffer);

#ifdef DEBUG12042012
         //--------------------------------------------------
         // DEBUG
         // read k-mers from files
         //--------------------------------------------------
         /*
         if (file_index == 0) {
            std::cout << "D1: ";
            for (std::size_t it_byte = 0; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
               char char_buf(kmer[it_byte]);

               // last byte
               if (it_byte == (num_bytes_per_kmer - 1)) {
                  std::cout << decode_a_byte(char_buf, num_empty_characters);
               }
               // other bytes
               else {
                  std::cout << decode_a_byte(char_buf, 0);
               }
            }
            std::cout << std::endl;
         }
         */
         //--------------------------------------------------
#endif

         // add kmer to the google sparse hash
         google::sparse_hash_map<KEY_TYPE, unsigned short int>::iterator it_hash;
         it_hash = gs_count_kmers.find(kmer);
         // kmer is not in the hash table: add it
         if (it_hash == gs_count_kmers.end()) {
            gs_count_kmers[kmer] = 1;
         }
         // kmer is in the hash table: incrment it
         else {
            gs_count_kmers[kmer]++;
         }

         num_kmers++;
      }

      if (num_kmers != distributed_num_kmers[file_index]) {
         std::cout << std::endl << "ERROR: Number of k-mers in " << distributed_kmer_file_names[file_index] << " is not matched with the original number" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Number of k-mers in " << distributed_kmer_file_names[file_index] << " is not matched with the original number" << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // count unique solid k-mers
      distributed_num_unique_kmers[file_index] = gs_count_kmers.size();

      // open a unique solid k-mer file
      std::ofstream f_unique_solid_kmers;
      f_unique_solid_kmers.open(distributed_unique_solid_kmer_file_names.c_str(), std::ios::binary);

      if (f_unique_solid_kmers.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << distributed_unique_solid_kmer_file_names << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Cannot open " << distributed_unique_solid_kmer_file_names << std::endl << std::endl;
         exit(EXIT_FAILURE);
      }

      // iterate all the keys in the hash table
      google::sparse_hash_map<KEY_TYPE, unsigned short int>::iterator it_hash;
      for (it_hash = gs_count_kmers.begin(); it_hash != gs_count_kmers.end(); it_hash++) {
#ifdef DEBUG12042012
         //--------------------------------------------------
         // DEBUG
         // read the number of occurrences of each k-mer from the hash table
         //--------------------------------------------------
         if (file_index == 0) {
            std::cout << "D2: ";
            for (std::size_t it_byte = 0; it_byte < MAX_BOOST_ARRAY_SIZE; it_byte++) {
               char char_buf(it_hash->first[it_byte]);

               // last byte
               if (it_byte == (num_bytes_per_kmer - 1)) {
                  std::cout << decode_a_byte(char_buf, num_empty_characters);
               }
               // other bytes
               else {
                  std::cout << decode_a_byte(char_buf, 0);
               }
            }

            std::cout << std::setw(10) << it_hash->second;
            std::cout << std::endl;
         }
         //--------------------------------------------------
#endif

         // count the number of unique solid k-mers
         if (it_hash->second >= kmer_occurrence_threshold) {
            distributed_num_unique_solid_kmers[file_index]++;

            // copy the key to the buffer
            std::vector<char> write_buffer;
            write_buffer.resize(num_bytes_per_kmer);
            for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
               write_buffer[it_byte] = it_hash->first[it_byte];
            }

#ifdef DEBUG12052012
            //--------------------------------------------------
            // DEBUG
            // print unique solid k-mers
            //--------------------------------------------------
            // if (file_index == 0) {
               std::string string_tmp("");
               std::string string_tmp_rc;

               for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
                  char char_buf(it_hash->first[it_byte]);

                  // last byte
                  if (it_byte == (num_bytes_per_kmer - 1)) {
                     string_tmp = string_tmp + decode_a_byte(char_buf, num_empty_characters);
                     // std::cout << decode_a_byte(char_buf, num_empty_characters);
                  }
                  // other bytes
                  else {
                     string_tmp = string_tmp + decode_a_byte(char_buf, 0);
                     // std::cout << decode_a_byte(char_buf, 0);
                  }
               }
               reverse_complement(string_tmp, string_tmp_rc);
               std::cout << "D3: " << string_tmp << " " << string_tmp_rc << std::endl;
            // }
            //--------------------------------------------------
#endif

            // write an unique solid k-mer
            f_unique_solid_kmers.write((const char*)&write_buffer[0], num_bytes_per_kmer);
         }

         // update the histogram
         std::size_t num_occurrences_tmp;
         if (it_hash->second > (HISTOGRAM_SIZE - 1)) {
            num_occurrences_tmp = HISTOGRAM_SIZE - 1;
         }
         else {
            num_occurrences_tmp = it_hash->second;
         }

         num_occurrences_histogram[num_occurrences_tmp]++;
      }

      // purge the hash table
      gs_count_kmers.clear();

      f_unique_solid_kmers.close();

#ifdef DEBUG12042012
      //--------------------------------------------------
      // DEBUG
      // check unique solid k-mer files
      //--------------------------------------------------
      if (file_index == 0) {
         std::ifstream f_tmp;
         f_tmp.open(c_inst_args.distributed_unique_solid_kmer_file_names[file_index].c_str(), std::ios::binary);
         if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Cannot open " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }

         char buffer;
         f_tmp.get(buffer);

         while (!f_tmp.eof()) {
            std::cout << "D4: ";
            for (std::size_t it_byte = 0; it_byte < num_bytes_per_kmer; it_byte++) {
               if (!f_tmp.good()) {
                  std::cout << std::endl << "ERROR: Number of bytes in " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
                  f_log     << std::endl << "ERROR: Number of bytes in " << c_inst_args.distributed_unique_solid_kmer_file_names[file_index] << " is illegal" << std::endl << std::endl;
                  exit(EXIT_FAILURE);
               }
               else {
                  char char_buf(buffer);
                  // last byte
                  if (it_byte == (num_bytes_per_kmer - 1)) {
                     std::cout << decode_a_byte(char_buf, num_empty_characters);
                  }
                  // other bytes
                  else {
                     std::cout << decode_a_byte(char_buf, 0);
                  }

                  f_tmp.get(buffer);
               }
            }

            std::cout << std::endl;
         }

         f_tmp.close();
      }
      //--------------------------------------------------
#endif
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << distributed_kmer_file_names[file_index] << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << distributed_kmer_file_names[file_index] << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_in.close();
}



//----------------------------------------------------------------------
// write_a_kmer
//----------------------------------------------------------------------
inline void C_count_solid_kmers::write_a_kmer(const std::string kmer, const std::string kmer_rc, const std::size_t kmer_middle_index, const std::size_t kmer_length, const std::size_t& residue, const std::size_t& num_bytes, const std::size_t& num_clusters) {
   //--------------------------------------------------
   // determine the output k-mer file
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
         kmer_internal = kmer_rc;
         break;
      case 'T' :
         kmer_internal = kmer_rc;
         break;
      default :
         std::cout << std::endl << "ERROR: Illegal character " << kmer[kmer_middle_index] << " in read files" << std::endl << std::endl;
         f_log     << std::endl << "ERROR: Illegal character " << kmer[kmer_middle_index] << " in read files" << std::endl << std::endl;
         exit(EXIT_FAILURE);
         break;
   }

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
         // char_index = it_char / BPS_PER_BYTE;
         encoded_kmer.push_back(buffer);

         // initialize buffer
         buffer &= ZERO;
      }
   }

   // encode remaining characters if needed
   if (residue != 0) {
      encoded_kmer.push_back(buffer);
   }

   // binary writing
   distributed_kmer_file_streams[file_index]->write((const char*)&encoded_kmer[0], num_bytes);
   distributed_num_kmers[file_index]++;

#ifdef DEBUG12042012
   //--------------------------------------------------
   // DEBUG
   // k-mer -> encoded k-mer
   //--------------------------------------------------
   if (file_index == 0) {
      std::size_t num_empty_characters(BPS_PER_BYTE - residue);
      std::cout << "D0-0: " << kmer_internal << std::endl;

      std::cout << "D0-1: ";
      for (std::size_t it_byte = 0; it_byte < num_bytes; it_byte++) {
         char char_buf(encoded_kmer[it_byte]);

         // last byte
         if (it_byte == (num_bytes - 1)) {
            std::cout << decode_a_byte(char_buf, num_empty_characters);
         }
         // other bytes
         else {
            std::cout << decode_a_byte(char_buf, 0);
         }
      }
      std::cout << std::endl;
   }
   //--------------------------------------------------
#endif
}



//----------------------------------------------------------------------
// encode_a_char
//----------------------------------------------------------------------
inline void C_count_solid_kmers::encode_a_char(const char& one_nt, char& buffer) {
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
// decode_a_byte
//----------------------------------------------------------------------
inline std::string C_count_solid_kmers::decode_a_byte(char& in, const std::size_t& num_empty_characters) {
   std::string out;

   if (num_empty_characters >= BPS_PER_BYTE) {
      std::cout << std::endl << "ERROR: The number of empty characters should be smaller than " << BPS_PER_BYTE << std::endl << std::endl;
      f_log     << std::endl << "ERROR: The number of empty characters should be smaller than " << BPS_PER_BYTE << std::endl << std::endl;
   }
   else {
      for (std::size_t it_empty = 0; it_empty < num_empty_characters; it_empty++) {
         in = in << 2;
      }

      for (std::size_t it_char = 0; it_char < (BPS_PER_BYTE - num_empty_characters); it_char++) {
         out.push_back(decode_a_char(in));
      }
   }

   return out;
}



//----------------------------------------------------------------------
// decode_a_char
//----------------------------------------------------------------------
inline char C_count_solid_kmers::decode_a_char(char& in) {
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
// reverse_complement
//----------------------------------------------------------------------
inline void C_count_solid_kmers::reverse_complement(const std::string& read, std::string& read_rc) {
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
            std::cout << std::endl << "ERROR: Illegal character " << *it << " in read files" << std::endl << std::endl;
            f_log     << std::endl << "ERROR: Illegal character " << *it << " in read files" << std::endl << std::endl;
            exit(EXIT_FAILURE);
            break;
      }
   }
}



//----------------------------------------------------------------------
// write_histogram
//----------------------------------------------------------------------
void C_count_solid_kmers::write_histogram(const std::string& histo_file_name, const bool& set_kmer_occurrence_threshold) {
   std::ofstream f_histo;
   f_histo.open(histo_file_name.c_str());

   bloom_type prev_value(0);

   if (f_histo.is_open()) {
      for (std::size_t it_histo = 1; it_histo < HISTOGRAM_SIZE - 1; it_histo++) {
         // find the valey point
         if ((num_occurrences_histogram[it_histo] > prev_value) && (it_histo > 1) && (get_valey_point == false)) {
            valey_point = it_histo - 1;
            get_valey_point = true;
         }

         // write histogram values
         f_histo << std::setw(7) << it_histo << ": " << std::setw(10) << num_occurrences_histogram[it_histo] << std::endl;

         prev_value = num_occurrences_histogram[it_histo];
      }

      f_histo << "<=" << std::setw(5) << HISTOGRAM_SIZE - 1<< ": " << std::setw(10) << num_occurrences_histogram[HISTOGRAM_SIZE - 1] << std::endl;

      if (set_kmer_occurrence_threshold == false) {
         if (get_valey_point == true) {
            std::cout << "     k-mer occurrence threshold   : " << valey_point << std::endl;
            f_log     << "     k-mer occurrence threshold   : " << valey_point << std::endl;

            kmer_occurrence_threshold = valey_point;
         }
         else {
            std::cout << "     No valey point exists in the histogram and no k-mer occurrence threshold is given" << valey_point << std::endl;
            f_log     << "     No valey point exists in the histogram and no k-mer occurrence threshold is given" << valey_point << std::endl;
            exit(EXIT_FAILURE);
         }
      }
   }
   else {
      std::cout << std::endl << "ERROR: Cannot open " << histo_file_name << std::endl << std::endl;
      f_log     << std::endl << "ERROR: Cannot open " << histo_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   f_histo.close();
}



//----------------------------------------------------------------------
// remove_kmer_files
//----------------------------------------------------------------------
void C_count_solid_kmers::remove_kmer_files(const std::size_t& num_clusters) {
   for (std::size_t it_files = 0; it_files < num_clusters; it_files++) {
      remove(distributed_kmer_file_names[it_files].c_str());
   }
}
