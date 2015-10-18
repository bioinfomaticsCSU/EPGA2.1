#ifndef _TIME_LOCAL_H
#define _TIME_LOCAL_H



#include "define.hpp"



//----------------------------------------------------------------------
// C_time
//----------------------------------------------------------------------
class C_time {
public:
   // variables
   // parse arguments
   std::string start_parse_args;
   std::string end_parse_args;

   // check read files
   std::string start_check_read_file;
   std::string end_check_read_file;

   // distribute k-mers
   std::string start_distribute_kmers;
   std::string end_distribute_kmers;

   // count k-mers
   std::string start_count_kmers;
   std::string end_count_kmers;

   // program k-mers into the Bloom filter
   std::string start_program_kmers_into_bf;
   std::string end_program_kmers_into_bf;

   // find false positive k-mers
   std::string start_find_false_positives;
   std::string end_find_false_positives;

   // count false positive k-mers
   std::string start_count_false_positives;
   std::string end_count_false_positives;

   // correct false positive k-mers
   std::string start_correct_false_positives;
   std::string end_correct_false_positives;

   // correct errors in reads
   std::string start_correct_errors_in_reads;
   std::string end_correct_errors_in_reads;

   // write corrected reads
   std::string start_write_corrected_reads;
   std::string end_write_corrected_reads;

   // write a tef file
   std::string start_write_tef;
   std::string end_write_tef;

   // write querying results of original reads
   std::string start_mark_unique_solid_kmers_org;
   std::string end_mark_unique_solid_kmers_org;

   // write querying results of corrected reads
   std::string start_mark_unique_solid_kmers_mod;
   std::string end_mark_unique_solid_kmers_mod;

   // constructors
   C_time() {};
};



#endif
