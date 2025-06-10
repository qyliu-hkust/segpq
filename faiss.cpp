#include <iostream>
#include <fstream>
#include <faiss/IndexPQ.h>
#include <cstring> 
#include <algorithm>
#include <vector>
#include "la_vector.hpp"


void load_data(char* filename, float*& data, unsigned& num, unsigned& dim) { 
  std::ifstream in(filename, std::ios::binary);	
  if (!in.is_open()) {
    std::cout << "open file error" << std::endl;
    exit(-1);
  }
  in.read((char*)&dim, 4);	
  in.seekg(0, std::ios::end);	
  std::ios::pos_type ss = in.tellg();	
  size_t fsize = (size_t)ss;
  num = (unsigned)(fsize / (dim + 1) / 4);	
  data = new float[(size_t)num * (size_t)dim];

  in.seekg(0, std::ios::beg);	
  for (size_t i = 0; i < num; i++) {
    in.seekg(4, std::ios::cur);	
    in.read((char*)(data + i * dim), dim * 4);	
  }
  in.close();
}

void pack(const uint8_t* uint8Array, uint64_t* uint64Array, size_t uint8Count) {  
    size_t uint64Count = uint8Count / 8;  
    for (size_t i = 0; i < uint64Count; ++i) {  
        uint64_t value = 0;  
        memcpy(&value, &uint8Array[i * 8], 8);  
        uint64Array[i] = value;  
    }  
}


int main(int argc, char** argv) {
  float* data_load = NULL;
  unsigned points_num, dim;
  load_data(argv[1], data_load, points_num, dim);
  std::cout << "points_num:"<< points_num << std::endl << "data dimension:" << dim << std::endl;

  faiss::IndexPQ index(128, 8, 8);

  index.train(25000, data_load);
  
  uint8_t* codes = new uint8_t[25000*8];
  index.pq.compute_codes(data_load, codes, 25000);

  uint64_t* code_ints = new uint64_t[25000];
  pack(codes, code_ints, 25000*8);

//   for (int i=0; i<10; ++i) {
//     std::cout << "============== i=" << i << " ==============" << std::endl;
//     for (int j=i*8; j<i*8+8; ++j) {
//         std::cout << std::hex << static_cast<int>(codes[j]) << " ";
//     }
//     std::cout << std::endl;
//     std::cout << std::hex << code_ints[i] << std::endl;
//   }
  std::sort(code_ints, code_ints+25000);
  std::vector<uint64_t> code_ints_vec(code_ints, code_ints+10000);


  la_vector<uint64_t, 16> v1(code_ints_vec);
  std::cout << "segs: " << v1.segments_count() << std::endl;
  std::cout << "bytes:" << v1.size_in_bytes() << std::endl; 


  delete[] codes;
  delete[] data_load;
  return 0;
}
