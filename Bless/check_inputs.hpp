#ifndef _PARSE_READ_H
#define _PARSE_READ_H



#include "parse_args.hpp"
#include "time.hpp"



//----------------------------------------------------------------------
// C_check_read
//----------------------------------------------------------------------
class C_check_read {
public:
   // variables
   std::size_t num_reads;
   std::size_t read_length;
   std::size_t init_num_elements;
   std::size_t interval;
   std::size_t min_num_processed_reads;
   std::size_t quality_score_offset;

   int max_quality_score;
   int min_quality_score;

   // constructor
   C_check_read() :
                   num_reads(0),
                   read_length(0),
                   init_num_elements(0),
                   interval(0),
                   min_num_processed_reads(0),
                   quality_score_offset(0),
                   max_quality_score(-1000),
                   min_quality_score(1000)
                  {};

   // functions
   void check_read_file(const C_arg& c_inst_args, C_time& c_inst_time);

private:
   // variables
   std::ofstream f_log;

   // functions
   void check_read_file_fastq_single(const C_arg& c_inst_args);
   void check_read_file_fastq_paired(const C_arg& c_inst_args);
};



#endif
