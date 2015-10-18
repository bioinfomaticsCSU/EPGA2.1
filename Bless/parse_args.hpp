#ifndef _PARSE_ARGS_H
#define _PARSE_ARGS_H



#include "define.hpp"
#include "time.hpp"



//----------------------------------------------------------------------
// C_arg
//----------------------------------------------------------------------
class C_arg {
public:
   // variables
   std::string read_file_name;
   std::string read_file_name1;
   std::string read_file_name2;
   std::string log_file_name;
   std::string prefix;
   std::string error_correction_info_file_name;
   std::string error_correction_info_file_name1;
   std::string error_correction_info_file_name2;
   std::string verified_error_correction_info_file_name;
   std::string verified_error_correction_info_file_name1;
   std::string verified_error_correction_info_file_name2;
   std::string corrected_read_file_name;
   std::string corrected_read_file_name1;
   std::string corrected_read_file_name2;
   std::string bf_querying_result_org_file_name;
   std::string bf_querying_result_org_file_name1;
   std::string bf_querying_result_org_file_name2;
   std::string bf_querying_result_mod_file_name;
   std::string bf_querying_result_mod_file_name1;
   std::string bf_querying_result_mod_file_name2;
   std::string error_correction_info_txt_file_name;
   std::string error_correction_info_txt_file_name1;
   std::string error_correction_info_txt_file_name2;
   std::string verified_error_correction_info_txt_file_name;
   std::string verified_error_correction_info_txt_file_name1;
   std::string verified_error_correction_info_txt_file_name2;
   std::string tef_file_name;
   std::string histo_file_name;
   std::string error_correction_direction_file_name;
   std::string error_correction_direction_file_name1;
   std::string error_correction_direction_file_name2;

   bloom_type random_seed;

   double target_false_positive_prob;

   std::size_t kmer_length;
   std::size_t kmer_middle_index;
   std::size_t kmer_occurrence_threshold;
   std::size_t num_clusters;
   std::size_t num_clusters_digit;
   std::size_t residue_kmer;
   std::size_t remainder_kmer;
   std::size_t num_bytes_per_kmer;
   std::size_t extend;

   bool paired_read;
   bool nowrite;
   bool debug;
   bool verify;
   bool set_kmer_occurrence_threshold;

   // constructors
   explicit C_arg(int argc, char** argv, C_time& c_inst_time) :
                                                               random_seed(0),
                                                               target_false_positive_prob(DEFAULT_FPR),
                                                               kmer_length(0),
                                                               kmer_middle_index(0),
                                                               kmer_occurrence_threshold(0),
                                                               num_clusters(NUM_PARTITIONS_FOR_COUNT),
                                                               num_clusters_digit(0),
                                                               residue_kmer(0),
                                                               remainder_kmer(0),
                                                               num_bytes_per_kmer(0),
                                                               extend(MAX_EXTENSION),
                                                               paired_read(true),
                                                               nowrite(false),
                                                               debug(false),
                                                               verify(true),
                                                               set_kmer_occurrence_threshold(false),
                                                               num_args(argc),
                                                               args(argv)
                                                              {
      time_t rawtime;
      time(&rawtime);
      c_inst_time.start_parse_args = asctime(localtime(&rawtime));

      read_args();

      time(&rawtime);
      c_inst_time.end_parse_args = asctime(localtime(&rawtime));
   };

private:
   // variables
   int num_args;

   char** args;

   std::ofstream f_log;

   // functions
   void read_args();
   void print_usage();
};



#endif
