#ifndef TRIS_KMER_H
#define TRIS_KMER_H

#ifdef _MSC_VER
#define NOMINMAX
#endif

#define LENGTH (1 << 22)
#define MAX_SEQ 1000

#define ABS_MIN_MER 3
#define ABS_TABLE_MAX_MER 15
#define ABS_UINT64_MAX_MER 32
#define ABS_MAX_MER 64

// k_mer_data macro
#define K_MER_DATA_COUNT 0
#define K_MER_DATA_MAX 1
#define K_MER_DATA_MAX_SEQ 2

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#include <cstdint>
#include <tuple>

#include <boost/unordered_map.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <tbb/concurrent_queue.h>
#include <tbb/tbb.h>

#include <zlib.h>

extern int MIN_MER;
extern int MAX_MER;
extern int TABLE_MAX_MER;
extern int NUM_THREAD;
extern int QUEUE_SIZE;
extern double BASELINE;

typedef boost::multiprecision::uint128_t uint128_t;
typedef std::tuple<uint16_t**, uint64_t**, uint128_t**, uint32_t**> ThreadDataTuple;

typedef std::pair<int, uint128_t> KmerSeq;
typedef boost::unordered_map<KmerSeq, uint32_t> ResultMap;
typedef std::vector<std::pair<int, int>> LocationVector;

typedef boost::unordered_map<uint64_t, uint16_t> CounterMap;
typedef boost::unordered_map<uint128_t, uint16_t> CounterMap_128;

struct QueueData {
    char* buffer;
    LocationVector* loc_vector;
};

typedef tbb::concurrent_bounded_queue<QueueData> TBBQueue;
typedef tbb::task_group TaskGroup;

uint8_t** set_repeat_check_table();
uint32_t** set_rotation_table(uint8_t** repeat_check_table);

uint64_t* set_extract_k_mer();
uint128_t* set_extract_k_mer_128();

void fill_rotation_table(uint32_t* rot, int k, uint32_t seq, uint8_t** repeat_check_table);

uint16_t** set_k_mer_counter();

uint64_t** set_k_mer_data();
uint128_t** set_k_mer_data_128();

uint32_t** set_k_mer_counter_list();

class ThreadData {
private:
    uint16_t **k_mer_counter = nullptr;
    uint64_t **k_mer_data = nullptr;
    uint128_t **k_mer_data_128 = nullptr;
    uint32_t **k_mer_counter_list = nullptr;
public:
    ThreadDataTuple init_check() {
        if (k_mer_counter == nullptr) {
            if (MAX_MER <= ABS_UINT64_MAX_MER) {
                k_mer_data = set_k_mer_data();
            } else {
                k_mer_data_128 = set_k_mer_data_128();
            }

            if (MIN_MER <= TABLE_MAX_MER) {
                k_mer_counter = set_k_mer_counter();
                k_mer_counter_list = set_k_mer_counter_list();
            }
        }
        return ThreadDataTuple {k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list};
    }
};

ResultMap* buffer_task(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table);

ResultMap* read_fastq(FILE* fp, ThreadData* thread_data, uint32_t** rot_table, uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint8_t** repeat_check_table);
void read_fastq_thread(FILE* fp, TBBQueue* buffer_task_queue);

ResultMap* read_fastq_gz(gzFile fp, ThreadData* thread_data, uint32_t** rot_table, uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint8_t** repeat_check_table);
void read_fastq_gz_thread(gzFile fp, TBBQueue* buffer_task_queue);

void int_to_four(char* buffer, uint128_t seq, int n);

void process_kmer(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                  uint64_t *extract_k_mer, uint128_t *extract_k_mer_128,
                  ThreadData* thread_data_list, bool is_gz);


void k_mer(const char* seq, int st, int nd, uint32_t** rot_table, const uint64_t *extract_k_mer,
           uint16_t** k_mer_counter, CounterMap* k_mer_counter_map,
           uint64_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
           ResultMap* result,
           int16_t* k_mer_total_cnt);

void k_mer_128(const char* seq, int st, int nd, uint32_t** rot_table, const uint128_t *extract_k_mer,
               uint16_t** k_mer_counter, CounterMap_128* k_mer_counter_map,
               uint128_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
               ResultMap* result,
               int16_t* k_mer_total_cnt);

uint64_t get_rot_seq(uint64_t seq, int k);
uint128_t get_rot_seq_128(const uint128_t& seq, int k);

uint8_t get_repeat_check(uint64_t seq, int k);
uint8_t get_repeat_check(uint128_t seq, int k);

uint32_t reverse_complement(uint32_t x);
uint64_t reverse_complement(uint64_t x);
uint128_t reverse_complement_128(uint128_t x);
#endif //TRIS_KMER_H
