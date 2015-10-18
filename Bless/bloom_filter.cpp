#include "bloom_filter.hpp"



//----------------------------------------------------------------------
// C_bloom_filter_normal::find_optimal_parameters
//----------------------------------------------------------------------
void C_bloom_filter_normal::find_optimal_parameters(const bloom_type num_elements, const double target_false_positive_prob) {
   // initialize variables
   double min_bit_vector_width_element = std::numeric_limits<double>::infinity();
   double min_num_hash_func = 1.0;
   double current_bit_vector_width_element = 0.0;

   // find the number of hash functions for minimum bit-vector size
   for (double j = min_num_hash_func; j < 1000.0; ++j) {
      if ((current_bit_vector_width_element = ((- j * num_elements) / std::log(1.0 - std::pow(target_false_positive_prob, 1.0 / j)))) < min_bit_vector_width_element) {
         min_bit_vector_width_element = current_bit_vector_width_element;
         min_num_hash_func = j;
      }
   }

   num_hash_func         = static_cast<unsigned int>(min_num_hash_func);
   bit_vector_width      = static_cast<bloom_type>(min_bit_vector_width_element);
   bit_vector_width      += (((bit_vector_width % BITS_PER_CHAR) != 0) ? (BITS_PER_CHAR - (bit_vector_width % BITS_PER_CHAR)) : 0);
   bit_vector_width_byte = bit_vector_width / BITS_PER_CHAR;
}



//----------------------------------------------------------------------
// C_bloom_filter_normal::generate_unique_hash_func
//----------------------------------------------------------------------
void C_bloom_filter_normal::generate_unique_hash_func(const bloom_type random_seed) {
   // initialize hash function vector
   hash_func.clear();

   // generate hash functions
   if (num_hash_func <= PREDEF_NUM_HASH_FUNC) {
      std::copy(PREDEF_HASH_FUNC, PREDEF_HASH_FUNC + num_hash_func, std::back_inserter(hash_func));

      for (unsigned int i = 0; i < hash_func.size(); ++i) {
         hash_func[i] = hash_func[i] * hash_func[(i + 3) % hash_func.size()] + random_seed;
      }
   }
   else {
      std::copy(PREDEF_HASH_FUNC, PREDEF_HASH_FUNC + PREDEF_NUM_HASH_FUNC, std::back_inserter(hash_func));

      srand(static_cast<unsigned int>(random_seed));
      while(hash_func.size() < num_hash_func) {
         bloom_type current_hash_func = static_cast<bloom_type>(rand()) * static_cast<bloom_type>(rand());
         if (0 == current_hash_func) continue;
         if (hash_func.end() == std::find(hash_func.begin(), hash_func.end(), current_hash_func)) {
            hash_func.push_back(current_hash_func);
         }
      }
   }
}



//----------------------------------------------------------------------
// C_bloom_filter_normal::generate_bit_vector
//----------------------------------------------------------------------
void C_bloom_filter_normal::generate_bit_vector() {
   delete[] bit_vector;
   bit_vector = NULL;

   // allocate memory for bit-vector and initialize it
   bit_vector = new unsigned char[static_cast<std::size_t>(bit_vector_width_byte)];
   std::fill_n(bit_vector, bit_vector_width_byte, 0x00);

   // initialize the number of elements
   num_elements = 0;
}



//----------------------------------------------------------------------
// C_bloom_filter_normal::min_size_bit_vector
//----------------------------------------------------------------------
void C_bloom_filter_normal::min_size_bit_vector() {
   delete[] bit_vector;
//   bit_vector = NULL;

   // allocate memory for bit-vector and initialize it
//   bit_vector = new unsigned char[1];
}



//----------------------------------------------------------------------
// C_bloom_filter_normal::clear_bit_vector
//----------------------------------------------------------------------
void C_bloom_filter_normal::clear_bit_vector() {
   std::fill_n(bit_vector, bit_vector_width_byte, 0x00);
   num_elements = 0;
}



//----------------------------------------------------------------------
// C_bloom_filter_normal::hash_ap
//----------------------------------------------------------------------
inline bloom_type C_bloom_filter_normal::hash_ap(const unsigned char* begin, std::size_t remaining_length, bloom_type hash) const {
   const unsigned char* itr = begin;
   unsigned int loop = 0;

   while (remaining_length >= 8) {
      const unsigned int& i1 = *(reinterpret_cast<const unsigned int*>(itr));
      itr += sizeof(unsigned int);
      const unsigned int& i2 = *(reinterpret_cast<const unsigned int*>(itr));
      itr += sizeof(unsigned int);
      hash ^= (hash <<  7) ^ i1 * (hash >> 3) ^ (~((hash << 11) + (i2 ^ (hash >> 5))));
      remaining_length -= 8;
   }

   while (remaining_length >= 4) {
      const unsigned int& i = *(reinterpret_cast<const unsigned int*>(itr));
      if (loop & 0x01) {
         hash ^= (hash <<  7) ^ i * (hash >> 3);
      }
      else {
         hash ^= (~((hash << 11) + (i ^ (hash >> 5))));
      }

      ++loop;
      remaining_length -= 4;
      itr += sizeof(unsigned int);
   }

   while (remaining_length >= 2) {
      const unsigned short& i = *(reinterpret_cast<const unsigned short*>(itr));
      if (loop & 0x01) {
         hash ^=    (hash <<  7) ^  i * (hash >> 3);
      }
      else {
         hash ^= (~((hash << 11) + (i ^ (hash >> 5))));
      }

      ++loop;
      remaining_length -= 2;
      itr += sizeof(unsigned short);
   }

   if (remaining_length) {
      hash += ((*itr) ^ (hash * 0xA5A5A5A5A5A5A5A5ULL)) + loop;
   }

   return hash;
}



//----------------------------------------------------------------------
// reverse_complement
//----------------------------------------------------------------------
inline void C_bloom_filter_normal::reverse_complement(const std::string& read, std::string& read_rc) {
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
            exit(EXIT_FAILURE);
            break;
      }
   }
}



//----------------------------------------------------------------------
// C_bloom_filter_normal::program
//----------------------------------------------------------------------
void C_bloom_filter_normal::program(const std::string& element) {
   bloom_type original_index(0);
   bloom_type bit_index(0);

   for (std::vector<bloom_type>::iterator it_hash_func = hash_func.begin(); it_hash_func != hash_func.end(); it_hash_func++) {
      original_index = (hash_ap(reinterpret_cast<const unsigned char*>(element.c_str()), element.size(), (*it_hash_func)) % bit_vector_width);
      bit_index = original_index % BITS_PER_CHAR;

      bit_vector[original_index / BITS_PER_CHAR] |= BIT_MASK[bit_index];
   }
   num_elements++;
}



//----------------------------------------------------------------------
// C_bloom_filter_normal::query
//----------------------------------------------------------------------
bool C_bloom_filter_normal::query(const std::string& element) {
   bloom_type original_index(0);
   bloom_type bit_index(0);

   for (std::vector<bloom_type>::iterator it_hash_func = hash_func.begin(); it_hash_func != hash_func.end(); it_hash_func++) {
      original_index = (hash_ap(reinterpret_cast<const unsigned char*>(element.c_str()), element.size(), (*it_hash_func)) % bit_vector_width);
      bit_index = original_index % BITS_PER_CHAR;

      if ((bit_vector[original_index / BITS_PER_CHAR] & BIT_MASK[bit_index]) != BIT_MASK[bit_index]) {
         return false;
      }
   }

   return true;
}
