#ifdef _MSC_VER
#define NOMINMAX
#endif

#include "kmer.h"
#include "argparse/argparse.hpp"
#include <filesystem>

int MAX_MER;
int MIN_MER;
int TABLE_MAX_MER;
int NUM_THREAD;
double BASELINE;

int main(int argc, char** argv) {
    argparse::ArgumentParser program("tris", "0.1.0");

    program.add_argument("MIN_MER")
            .help("minimum length of sequence to find telomere [MIN_MER >= " + std::to_string(ABS_MIN_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    program.add_argument("MAX_MER")
            .help("maximum length of sequence to find telomere [MAX_MER <= " + std::to_string(ABS_MAX_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    program.add_argument("FASTQ_LOC")
            .help("locations of FASTQ file")
            .nargs(argparse::nargs_pattern::at_least_one);

    program.add_argument("-t", "--thread")
            .help("number of threads")
            .default_value(1)
            .scan<'d', int>()
            .metavar("THREAD");

    program.add_argument("-m", "--table_max_mer")
            .help("maximum length of sequence to use table (reduce this option if memory usage is high) [TABLE_MAX_MER <= " + std::to_string(ABS_TABLE_MAX_MER) + "]")
            .default_value(12)
            .scan<'d', int>()
            .metavar("TABLE_MAX_MER");

    program.add_argument("-b", "--baseline")
            .help("BASELINE for repeat read [0.5 <= BASELINE <= 1]")
            .default_value(0.8)
            .scan<'g', double>()
            .metavar("BASELINE");

    try {
        program.parse_args(argc, argv);

        MIN_MER = program.get<int>("MIN_MER");
        MAX_MER = program.get<int>("MAX_MER");
        NUM_THREAD = program.get<int>("--thread");
        TABLE_MAX_MER = program.get<int>("--table_max_mer");
        BASELINE = program.get<double>("--baseline");

        // argument check
        if (MIN_MER > MAX_MER) {
            fprintf(stderr, "MIN_MER must not be greater than MAX_MER.\n");
            throw std::exception();
        }

        if (MIN_MER < ABS_MIN_MER) {
            fprintf(stderr, "MIN_MER must be greater than or equal to %d.\n", ABS_MIN_MER);
            throw std::exception();
        }

        if (MAX_MER > ABS_MAX_MER) {
            fprintf(stderr, "MAX_MER must be less than or equal to %d.\n", ABS_MAX_MER);
            throw std::exception();
        }

        if (TABLE_MAX_MER > ABS_TABLE_MAX_MER) {
            fprintf(stderr, "TABLE_MAX_MER must be less than or equal to %d.\n", ABS_TABLE_MAX_MER);
            throw std::exception();
        }

        if (0.5 > BASELINE || BASELINE > 1) {
            fprintf(stderr, "BASELINE must be at least 0.5 and no more than 1.\n");
            throw std::exception();
        }

        if (TABLE_MAX_MER <= 0) {
            fprintf(stderr, "TABLE_MAX_MER must be positive.\n");
            throw std::exception();
        }

        if (NUM_THREAD <= 0) {
            fprintf(stderr, "number of threads must be positive.\n");
            throw std::exception();
        }
    }
    catch (...) {
        std::cerr << program;
        return 1;
    }

    std::vector<std::string> fastq_loc_list = program.get<std::vector<std::string>>("FASTQ_LOC");
    std::vector<std::filesystem::path> fastq_path_list {};

    for (const auto& fastq_loc : fastq_loc_list) {
        std::filesystem::path fastq_path {fastq_loc};
        if (! std::filesystem::is_regular_file(fastq_loc)) {
            fprintf(stderr, "%s : file not found\n", fastq_loc.c_str());
            return 1;
        }
        fastq_path_list.push_back(fastq_path);
    }

    uint8_t **repeat_check_table = nullptr;
    uint32_t **rot_table = nullptr;
    if (MIN_MER <= TABLE_MAX_MER) {
        repeat_check_table = set_repeat_check_table();
        rot_table = set_rotation_table(repeat_check_table);
    }

    uint64_t *extract_k_mer = nullptr;
    uint128_t *extract_k_mer_128 = nullptr;
    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        extract_k_mer = set_extract_k_mer();
    }
    else {
        extract_k_mer_128 = set_extract_k_mer_128();
    }

    ThreadData* thread_data_list = new ThreadData[NUM_THREAD];

    std::vector<std::string> gz_extension_list = std::vector<std::string> {".gz", ".bgz"};
    for (const auto& fastq_path :fastq_path_list) {
        fprintf(stdout, ">%s\n", std::filesystem::canonical(fastq_path).string().c_str());

        bool is_gz = false;
        std::string fastq_ext = fastq_path.extension().string();
        for (const auto& ext : gz_extension_list) {
            if (ext == fastq_ext) {
                is_gz = true;
                break;
            }
        }

        process_kmer(fastq_path.string().c_str(), repeat_check_table, rot_table,
                     extract_k_mer, extract_k_mer_128,
                     thread_data_list, is_gz);
    }

    delete[] thread_data_list;
    return 0;
}