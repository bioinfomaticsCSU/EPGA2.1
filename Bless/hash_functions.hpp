#ifndef _HASH_FUNCTIONS_H
#define _HASH_FUNCTIONS_H



#include "define.hpp"



//----------------------------------------------------------------------
// C_hash_functions
//----------------------------------------------------------------------
class C_hash_functions {
public:
   // constructors
   C_hash_functions() {};

   // functions
	unsigned int rs_hash(const std::string& str);
	unsigned int js_hash(const std::string& str);
	unsigned int pjw_hash(const std::string& str);
	unsigned int elf_hash(const std::string& str);
	unsigned int bkdr_hash(const std::string& str);
	unsigned int sdbm_hash(const std::string& str);
	unsigned int djb_hash(const std::string& str);
	unsigned int dek_hash(const std::string& str);
	unsigned int bp_hash(const std::string& str);
	unsigned int fnv_hash(const std::string& str);
	unsigned int ap_hash(const std::string& str);

private:
};



#endif
