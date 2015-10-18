#include "count_solid_kmers.hpp"



//----------------------------------------------------------------------
// rs_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::rs_hash(const std::string& str) {
   unsigned int b    = 378551;
   unsigned int a    = 63689;
   unsigned int hash = 0;

   for(std::size_t i = 0; i < str.length(); i++) {
      hash = hash * a + str[i];
      a    = a * b;
   }

   return hash;
}



//----------------------------------------------------------------------
// js_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::js_hash(const std::string& str) {
   unsigned int hash = 1315423911;

   for(std::size_t i = 0; i < str.length(); i++) {
      hash ^= ((hash << 5) + str[i] + (hash >> 2));
   }

   return hash;
}



//----------------------------------------------------------------------
// pjw_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::pjw_hash(const std::string& str) {
   unsigned int BitsInUnsignedInt = (unsigned int)(sizeof(unsigned int) * 8);
   unsigned int ThreeQuarters     = (unsigned int)((BitsInUnsignedInt  * 3) / 4);
   unsigned int OneEighth         = (unsigned int)(BitsInUnsignedInt / 8);
   unsigned int HighBits          = (unsigned int)(0xFFFFFFFF) << (BitsInUnsignedInt - OneEighth);
   unsigned int hash              = 0;
   unsigned int test              = 0;

   for(std::size_t i = 0; i < str.length(); i++) {
      hash = (hash << OneEighth) + str[i];

      if((test = hash & HighBits)  != 0) {
         hash = (( hash ^ (test >> ThreeQuarters)) & (~HighBits));
      }
   }

   return hash;
}



//----------------------------------------------------------------------
// elf_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::elf_hash(const std::string& str) {
   unsigned int hash = 0;
   unsigned int x    = 0;

   for(std::size_t i = 0; i < str.length(); i++) {
      hash = (hash << 4) + str[i];
      if((x = hash & 0xF0000000L) != 0)
      {
         hash ^= (x >> 24);
      }
      hash &= ~x;
   }

   return hash;
}



//----------------------------------------------------------------------
// bkdr_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::bkdr_hash(const std::string& str) {
   unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
   unsigned int hash = 0;

   for(std::size_t i = 0; i < str.length(); i++) {
      hash = (hash * seed) + str[i];
   }

   return hash;
}



//----------------------------------------------------------------------
// sdbm_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::sdbm_hash(const std::string& str) {
   unsigned int hash = 0;

   for(std::size_t i = 0; i < str.length(); i++) {
      hash = str[i] + (hash << 6) + (hash << 16) - hash;
   }

   return hash;
}



//----------------------------------------------------------------------
// djb_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::djb_hash(const std::string& str) {
   unsigned int hash = 5381;

   for(std::size_t i = 0; i < str.length(); i++)
   {
      hash = ((hash << 5) + hash) + str[i];
   }

   return hash;
}



//----------------------------------------------------------------------
// dek_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::dek_hash(const std::string& str) {
   unsigned int hash = static_cast<unsigned int>(str.length());

   for(std::size_t i = 0; i < str.length(); i++) {
      hash = ((hash << 5) ^ (hash >> 27)) ^ str[i];
   }

   return hash;
}



//----------------------------------------------------------------------
// bp_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::bp_hash(const std::string& str) {
   unsigned int hash = 0;
   for(std::size_t i = 0; i < str.length(); i++) {
      hash = hash << 7 ^ str[i];
   }

   return hash;
}



//----------------------------------------------------------------------
// fnv_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::fnv_hash(const std::string& str) {
   const unsigned int fnv_prime = 0x811C9DC5;
   unsigned int hash = 0;
   for(std::size_t i = 0; i < str.length(); i++) {
      hash *= fnv_prime;
      hash ^= str[i];
   }

   return hash;
}



//----------------------------------------------------------------------
// ap_hash
//----------------------------------------------------------------------
unsigned int C_hash_functions::ap_hash(const std::string& str) {
   unsigned int hash = 0xAAAAAAAA;

   for(std::size_t i = 0; i < str.length(); i++) {
      hash ^= ((i & 1) == 0) ? (  (hash <<  7) ^ str[i] * (hash >> 3)) :
                               (~((hash << 11) + (str[i] ^ (hash >> 5))));
   }

   return hash;
}
