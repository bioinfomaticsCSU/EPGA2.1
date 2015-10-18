#ifndef _CORRECT_H
#define _CORRECT_H



#include "parse_args.hpp"
#include "bloom_filter.hpp"
#include "check_inputs.hpp"
#include "hash_functions.hpp"



//----------------------------------------------------------------------
// C_candidate_path
//----------------------------------------------------------------------
class C_candidate_path {
public:
   // std::string last_kmer;
   // modified positions and modification results
   std::vector< std::pair<std::size_t, char> > modified_bases;

   // sum of quality scores of modified bases
   std::size_t sum_qs;

   // constructors
   C_candidate_path() : sum_qs(0) {};

   void clear_path();
};



//----------------------------------------------------------------------
// C_correct_errors
//----------------------------------------------------------------------
class C_correct_errors {
public:
   // variables
   bloom_type num_unique_solid_kmers;

   std::size_t read_length;
   std::size_t num_false_positive_candidate_kmers;
   std::size_t num_false_corrections;
   std::size_t quality_score_offset;

   std::vector<std::size_t> distributed_num_unique_solid_kmers;

   std::vector<std::string> distributed_unique_solid_kmer_file_names;

   C_bloom_filter_normal c_inst_bloom_filter;

   // constructors
   C_correct_errors() :
                       num_unique_solid_kmers(0),
                       read_length(0),
                       num_false_positive_candidate_kmers(0),
                       num_false_corrections(0),
                       quality_score_offset(0),
                       num_corrected_errors_step1_1(0),
                       num_corrected_errors_step1_2(0),
                       num_corrected_errors_step1_3(0),
                       num_corrected_errors_step2_1(0),
                       num_corrected_errors_step2_2(0),
                       num_corrected_reads(0),
                       num_wrongly_corrected_errors(0),
                       num_wrongly_corrected_errors_check(0)
                      {};

   // functions
   void generate_bloom_filter(const C_arg& c_inst_args, C_time& c_inst_time);
   void program_kmers(const C_arg& c_inst_args, C_time& c_inst_time);
   void summarize_outputs(const C_arg& c_inst_args, C_time& c_inst_time);
   void correct_errors_in_reads(const C_arg& c_inst_args, C_time& c_inst_time);
   void find_false_positive_candidates(const C_arg& c_inst_args, C_time& c_inst_time);
   void verify_false_positive_candidates(const C_arg& c_inst_args, C_time& c_inst_time);
   void correct_false_positives(const C_arg& c_inst_args, C_time& c_inst_time);
   void write_corrected_reads(const C_arg& c_inst_args, C_time& c_inst_time);
   void write_tef(const C_arg& c_inst_args, C_time& c_inst_time);
   void mark_bloom_filter_results_org(const C_arg& c_inst_args, C_time& c_inst_time);
   void mark_bloom_filter_results_mod(const C_arg& c_inst_args, C_time& c_inst_time);
   void remove_kmer_files(const std::size_t& num_clusters, const bool& verify);
   void remove_error_correction_info_files(const C_arg& c_inst_args, const bool& verify);

private:
   // variables
   typedef boost::unordered_map <std::string, std::string> map_str_str;
   map_str_str multiple_incoming;

   std::vector<std::string> distributed_false_positive_candidate_file_names;
   std::vector<std::string> distributed_false_positive_result_file_names;

   std::vector<std::ofstream *> distributed_false_positive_candidate_file_streams;

   std::vector<std::ifstream *> distributed_false_positive_result_file_streams;

   C_hash_functions c_inst_hash_functions;
 
   std::ifstream f_read;
   std::ifstream f_read_1;
   std::ifstream f_read_2;

   std::ofstream f_log;

   std::size_t num_corrected_errors_step1_1;
   std::size_t num_corrected_errors_step1_2;
   std::size_t num_corrected_errors_step1_3;
   std::size_t num_corrected_errors_step2_1;
   std::size_t num_corrected_errors_step2_2;
   std::size_t num_corrected_reads;

   std::size_t num_wrongly_corrected_errors;
   std::size_t num_wrongly_corrected_errors_check;

   //functions
   void program_kmers_fastq(const C_arg& c_inst_args);
   void program_kmers_cluster(const std::size_t& file_index, const std::size_t& num_bytes_per_kmer);
   void reverse_complement(const std::string& read, std::string& read_rc);
   void encode_correction_info(char& buffer, const char& char_in_read, const char& char_in_info);
   void encode_a_char(const char& one_nt, char& buffer);
   void correct_errors_in_reads_single_fastq(const C_arg& c_inst_args);
   void correct_errors_in_reads_paired_fastq(const C_arg& c_inst_args);
   void correct_errors_in_a_read_fastq(const std::string& sequence, std::string& sequence_modification, const std::string& quality_score, const std::size_t kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& max_extension, const bool& verify);
   void correct_errors_between_solid_regions(std::string& sequence, const std::string& quality_score, const std::size_t& index_start, const std::size_t& index_end, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification);
   void correct_errors_5_prime_end(std::string& sequence, const std::string& quality_score, const std::size_t& index_start, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification, const std::size_t& max_extension, const bool& verify);
   void correct_errors_3_prime_end(std::string& sequence, const std::string& quality_score, const std::size_t& index_start, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification, const std::size_t& max_extension);
   void correct_errors_first_kmer(const std::string& sequence, const std::string& quality_score, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, std::string& sequence_modification, std::vector<C_candidate_path>& candidate_path_vector);
   void extend_a_kmer(const std::string& kmer, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& index_kmer, const std::size_t& index_last_mod, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector);
   void extend_a_kmer_5_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector);
   void extend_a_kmer_3_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector);
   void extend_out_left(const std::string& kmer, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success);
   void extend_out_right(const std::string& kmer, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success);
   void find_false_positive_candidates_in_single_fastq(const C_arg& c_inst_args);
   void find_false_positive_candidates_in_paired_fastq(const C_arg& c_inst_args);
   void find_false_positive_candidates_in_read(const std::size_t& kmer_middle_index, const std::size_t& kmer_length, const std::size_t& read_length, const std::string& read, std::string& corrected_read, const std::string& corrected_bp, const std::size_t& residue_kmer, const std::size_t& num_clusters, const std::size_t& num_bytes_per_kmer, std::ifstream& f_direction);
   void find_false_positive_candidates_in_kmer(const std::string& kmer, const std::size_t& kmer_middle_index, const std::size_t& num_clusters, const std::size_t& kmer_length, const std::size_t& residue_kmer, const std::size_t& num_btes_per_kmer);
   void correct_false_positives_single_fastq(const C_arg& c_inst_args);
   void correct_false_positives_paired_fastq(const C_arg& c_inst_args);
   void correct_false_positives_read(const std::size_t& kmer_middle_index, const std::size_t& kmer_length, const std::string& read, std::string& corrected_read, const std::string& corrected_bp, std::string& corrected_bp_verified, const bool& debug, const std::size_t& num_clusters, std::ifstream& f_direction);
   void verify_false_positive_candidates_fastq(const C_arg& c_inst_args);
   void verify_false_positive_candidates_cluster(const std::size_t& file_index, const std::size_t& num_bytes_per_kmer);
   void write_corrected_reads_single_fastq(const C_arg& c_inst_args);
   void write_corrected_reads_paired_fastq(const C_arg& c_inst_args);
   void write_tef_single_fastq(const C_arg& c_inst_args);
   void write_tef_paired_fastq(const C_arg& c_inst_args);
   void mark_bloom_filter_results_org_single_fastq(const C_arg& c_inst_args);
   void mark_bloom_filter_results_org_paired_fastq(const C_arg& c_inst_args);
   void mark_bloom_filter_results_mod_single_fastq(const C_arg& c_inst_args);
   void mark_bloom_filter_results_mod_paired_fastq(const C_arg& c_inst_args);
   void check_first_kmer(const std::string& kmer, const C_candidate_path& candidate_path_in, const std::vector<std::size_t>& low_qs_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& index, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer);
   void solid_first_kmer(const C_candidate_path& candidate_path, const std::string& sequence, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& max_extension, bool& extension_success);
   void extend_first_kmer_to_right(const std::string& sequence, const std::string& quality_score, C_candidate_path& candidate_path_in, const std::size_t& kmer_length, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& residue_kmer, const std::size_t& max_extension, bool& correction_success);

   bool query_text(const std::string& kmer, const std::size_t& kmer_middle_index, const std::size_t& num_bytes_per_kmer, const std::size_t& kmer_length, const std::size_t& residue_kmer);
   bool not_false_positive(const std::string& kmer, const std::size_t& kmer_middle_index, const std::size_t& num_clusters);

   char decode_correction_info(const char& first, const char& second, const char& read);
   char decode_correction_info_num(const char& first, const char& second, const char& read);
   char nt_to_num(const char& in_char);
   char decode_a_char(char& in);

   std::string remove_new_line(std::string in_string);
   std::string decode_a_byte(char& in, const std::size_t& num_empty_characters);
};



#endif
