#include <vector>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include <unistd.h> // Linux
#include "pgm_index_build.hpp"

std::string log_filename = "";
std::string decode_type = "";
bool read_only = false;

template <uint64_t epsilon=64>
void create_collection_pgm(const string input_basename, const string output_filename) {

    typedef pgm_sequence::pgm_builder<uint32_t, epsilon> PGM_INDEX_BUILDER;
    PGM_INDEX_BUILDER index;
    if (!read_only){
        index.build_model(input_basename + ".docs");
        index.statistic(log_filename);
        index.save_model(output_filename);
        // if(!log_filename.empty()) // save covered and residual
            // index.save_covered(log_filename);
            // index.save_residual(log_filename);
    } else {
        index.load_model(output_filename);
        // index.normal_init();
        // index.data_test(input_basename + ".docs");
        // index.statistic(log_filename);
        index.huffman_encode();
    }
}

int main(int argc, const char** argv)
{
    ios::sync_with_stdio(0);
    int mandatory = 6;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <index_type> <collection_basename> <output_filename> <epsilon> <read_only> <decode_type> <test_log>" << std::endl;
        return 1;
    }

    const std::string index_type = argv[1];
    const std::string input_basename = argv[2];
    const std::string output_filename = argv[3];
    const std::string epsilonstr = argv[4];
    const std::string read_only_str = argv[5];
    decode_type = argv[6];
    log_filename = argv[7];
    const uint64_t epsilon = static_cast<uint64_t>(std::stoi(epsilonstr));
    read_only = (read_only_str == "t");

    if (index_type == "pgm")  //compress docs
        switch (epsilon)
        {
            case 4: create_collection_pgm<4>(input_basename, output_filename); break;
            case 8: create_collection_pgm<8>(input_basename, output_filename); break;
            case 16: create_collection_pgm<16>(input_basename, output_filename); break;
            case 32: create_collection_pgm<32>(input_basename, output_filename); break;
            case 62: create_collection_pgm<62>(input_basename, output_filename); break;
            case 63: create_collection_pgm<63>(input_basename, output_filename); break;
            case 64: create_collection_pgm<64>(input_basename, output_filename); break;
            case 96: create_collection_pgm<96>(input_basename, output_filename); break;
            case 126: create_collection_pgm<126>(input_basename, output_filename); break;
            case 127: create_collection_pgm<127>(input_basename, output_filename); break;
            case 128: create_collection_pgm<128>(input_basename, output_filename); break;
            case 192: create_collection_pgm<192>(input_basename, output_filename); break;
            case 254: create_collection_pgm<254>(input_basename, output_filename); break;
            case 255: create_collection_pgm<255>(input_basename, output_filename); break;
            case 256: create_collection_pgm<256>(input_basename, output_filename); break;
            case 512: create_collection_pgm<512>(input_basename, output_filename); break;
            case 1024: create_collection_pgm<1024>(input_basename, output_filename); break;
            case 2048: create_collection_pgm<2048>(input_basename, output_filename); break;
            case 4096: create_collection_pgm<4096>(input_basename, output_filename); break;
            case 8192: create_collection_pgm<8192>(input_basename, output_filename); break;
            case 16384: create_collection_pgm<16384>(input_basename, output_filename); break;
            case 32768: create_collection_pgm<32768>(input_basename, output_filename); break;
            case 65536: create_collection_pgm<65536>(input_basename, output_filename); break;
            case 131072: create_collection_pgm<131072>(input_basename, output_filename); break;
            case 262144: create_collection_pgm<262144>(input_basename, output_filename); break;
            case 524288: create_collection_pgm<524288>(input_basename, output_filename); break;
            default: logStream << "Unsupported Epsilon Value: " << epsilon << endl;break;
        }
    else
        std::cerr << "ERROR: only support pgm index " << std::endl;

    return 0;
}