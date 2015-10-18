#ifndef _BLOOM_FILTER_H
#define _BLOOM_FILTER_H



#include "define.hpp"



//----------------------------------------------------------------------
// C_bloom_filter_normal
//----------------------------------------------------------------------
class C_bloom_filter_normal {
public:
   // variables
   unsigned char* bit_vector;

   unsigned int num_hash_func;

   bloom_type num_elements;
   bloom_type bit_vector_width_byte;

   // constructors
   C_bloom_filter_normal() :
                            bit_vector(NULL),
                            num_hash_func(0),
                            num_elements(0),
                            bit_vector_width_byte(0),
                            bit_vector_width(0)
                           {};

   // functions
   void find_optimal_parameters(const bloom_type num_elements, const double target_false_positive_prob);
   void generate_unique_hash_func(const bloom_type random_seed);
   void generate_bit_vector();
   void min_size_bit_vector();
   void clear_bit_vector();
   void program(const std::string& element);
   bool query(const std::string& element);

private:
   // variables
   bloom_type bit_vector_width;

   std::vector<bloom_type> hash_func;

   // functions
   bloom_type hash_ap(const unsigned char* begin, std::size_t remaining_length, bloom_type hash) const;

   void reverse_complement(const std::string& read, std::string& read_rc);
};



#endif
