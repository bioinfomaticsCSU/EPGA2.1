#ifndef _COUNT_SOLID_KMERS_H
#define _COUNT_SOLID_KMERS_H



#include "parse_args.hpp"
#include "bloom_filter.hpp"
#include "check_inputs.hpp"
#include "hash_functions.hpp"



//----------------------------------------------------------------------
// C_count_solid_kmers
//----------------------------------------------------------------------
class C_count_solid_kmers {
public:
   // variables
   bloom_type num_unique_kmers;
   bloom_type num_unique_solid_kmers;
   bloom_type hash_value;

   std::vector<std::ofstream *> distributed_kmer_file_streams;

   std::vector<std::string> distributed_kmer_file_names;
   std::vector<std::string> distributed_unique_solid_kmer_file_names;

   std::vector<std::size_t> distributed_num_kmers;
   std::vector<std::size_t> distributed_num_unique_kmers;
   std::vector<std::size_t> distributed_num_unique_solid_kmers;

   std::vector<bloom_type> num_occurrences_histogram;

   std::size_t kmer_occurrence_threshold;
   std::size_t valey_point;

   bool get_valey_point;

   // constructors
   C_count_solid_kmers() :
                          num_unique_kmers(0),
                          num_unique_solid_kmers(0),
                          hash_value(0),
                          kmer_occurrence_threshold(0),
                          valey_point(0),
                          get_valey_point(false)
                         {};

   // functions
   void distribute_kmers(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads, C_time& c_inst_time);
   void count_kmers(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads, C_time& c_inst_time);

private:
   C_hash_functions c_inst_hash_functions;

   std::ifstream f_read;
   std::ifstream f_read_1;
   std::ifstream f_read_2;

   std::ofstream f_log;

   //functions
   void reverse_complement(const std::string& read, std::string& read_rc);
   void write_a_kmer(const std::string kmer, const std::string kmer_rc, const std::size_t kmer_middle_index, const std::size_t kmer_length, const std::size_t& residue, const std::size_t& num_bytes, const std::size_t& num_clusters);
   void encode_a_char(const char& one_nt, char& buffer);
   void distribute_kmers_single_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads);
   void distribute_kmers_paired_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads);
   void count_unique_kmers_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads);
   void count_unique_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_for_kmer, const std::size_t num_empty_characters, const std::string& distributed_unique_solid_kmer_file_names);
   void count_unique_solid_kmers_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads);
   void count_unique_solid_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_for_kmer, const std::size_t num_empty_characters, const std::string& distributed_unique_solid_kmer_file_names);
   void count_both_unique_and_solid_kmers_fastq(const C_arg& c_inst_args, const C_check_read& c_inst_check_reads);
   void count_both_unique_and_solid_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_for_kmer, const std::size_t num_empty_characters, const std::string& distributed_unique_solid_kmer_file_names);
   void write_histogram(const std::string& histo_file_name, const bool& set_kmer_occurrence_threshold);
   void remove_kmer_files(const std::size_t& num_clusters);

   std::string decode_a_byte(char& in, const std::size_t& num_empty_characters);

   char decode_a_char(char& in);
};



#endif
