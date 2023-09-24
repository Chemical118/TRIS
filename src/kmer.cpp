#include "kmer.h"

#include <cstdio>
#include <cstdlib>
#include <vector>

char trans_arr[4] = {'T', 'G', 'C', 'A'};
#define O (-1)
#define A 3
#define C 2
#define G 1
#define T 0

static const int codes[256] = {
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, A, O, C, O, O, O, G, O, O, O, O, O, O, O, O,
        O, O, O, O, T, O, O, O, O, O, O, O, O, O, O, O,
        O, A, O, C, O, O, O, G, O, O, O, O, O, O, O, O,
        O, O, O, O, T, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O
};

#undef O
#undef A
#undef C
#undef G
#undef T


uint32_t reverse_complement(uint32_t x) {
    x = (x >> 16) | (x << 16);
    x = ((x >> 8) & 0x00ff00ff) | ((x & 0x00ff00ff) << 8);
    x = ((x >> 4) & 0x0f0f0f0f) | ((x & 0x0f0f0f0f) << 4);
    x = ((x >> 2) & 0x33333333) | ((x & 0x33333333) << 2);
    return ~x;
}

uint64_t reverse_complement(uint64_t x) {
    x = (x >> 32) | (x << 32);
    x = ((x >> 16) & 0x0000ffff0000ffffULL) | ((x & 0x0000ffff0000ffffULL) << 16);
    x = ((x >> 8)  & 0x00ff00ff00ff00ffULL) | ((x & 0x00ff00ff00ff00ffULL) <<  8);
    x = ((x >> 4)  & 0x0f0f0f0f0f0f0f0fULL) | ((x & 0x0f0f0f0f0f0f0f0fULL) <<  4);
    x = ((x >> 2)  & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) <<  2);
    return ~x;
}

static const uint128_t rev_cmp[5] = {uint128_t("0x00000000ffffffff00000000ffffffff"),
                                     uint128_t("0x0000ffff0000ffff0000ffff0000ffff"),
                                     uint128_t("0x00ff00ff00ff00ff00ff00ff00ff00ff"),
                                     uint128_t("0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f"),
                                     uint128_t("0x33333333333333333333333333333333")};

uint128_t reverse_complement_128(uint128_t x) {
    x = (x >> 64) | (x << 64);
    x = ((x >> 32) & rev_cmp[0]) | ((x & rev_cmp[0]) << 32);
    x = ((x >> 16) & rev_cmp[1]) | ((x & rev_cmp[1]) << 16);
    x = ((x >> 8)  & rev_cmp[2]) | ((x & rev_cmp[2]) <<  8);
    x = ((x >> 4)  & rev_cmp[3]) | ((x & rev_cmp[3]) <<  4);
    x = ((x >> 2)  & rev_cmp[4]) | ((x & rev_cmp[4]) <<  2);
    return ~x;
}

ResultMap* buffer_task(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMap* result = new ResultMap {};
    QueueData temp_task {};

    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        CounterMap* k_mer_counter_map = nullptr;
        if (TABLE_MAX_MER < MAX_MER) {
            k_mer_counter_map = new CounterMap[MAX_MER - TABLE_MAX_MER];
        }

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector == nullptr) {
                break;
            }
            for (auto& [st, nd]: *(temp_task.loc_vector)) {
                k_mer(temp_task.buffer, st, nd, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                      k_mer_data, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);
            }

            delete temp_task.loc_vector;
            free(temp_task.buffer);
        }

        delete[] k_mer_counter_map;
    } else {
        CounterMap_128* k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector == nullptr) {
                break;
            }
            for (auto& [st, nd]: *(temp_task.loc_vector)) {
                k_mer_128(temp_task.buffer, st, nd, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                          k_mer_data_128, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);
            }

            delete temp_task.loc_vector;
            free(temp_task.buffer);
        }

        delete[] k_mer_counter_map;
    }


    free(k_mer_total_cnt);
    return result;
}

ResultMap* read_fastq(FILE* fp, ThreadData* thread_data, uint32_t** rot_table, uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMap* result = new ResultMap {};

    int num = 0;
    int shift = 0;
    int bytes_read;
    int idx;

    char* buffer = (char*)malloc(sizeof(char) * LENGTH);

    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        CounterMap* k_mer_counter_map = nullptr;
        if (TABLE_MAX_MER < MAX_MER) {
            k_mer_counter_map = new CounterMap[MAX_MER - TABLE_MAX_MER];
        }

        while (1) {
            bytes_read = (int) fread(buffer + shift, 1, LENGTH - 1 - shift, fp);
            buffer[bytes_read + shift] = '\0';

            for (int i = 0; i < bytes_read + shift; i++) {
                if (buffer[i] == '\n') {
                    num += 1;
                    if ((num & 3) == 2) {
                        if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                            fprintf(stderr, "This program is designed for short-read sequencing\n");
                            exit(EXIT_FAILURE);
                        }

                        k_mer(buffer, idx + 1, i - 1, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                              k_mer_data, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);
                    }
                    idx = i;
                }
            }

            if (bytes_read <= 0) {
                if (feof(fp)) {
                    break;
                } else {
                    fprintf(stderr, "File-IO Error: %s.\n", strerror(errno));
                    exit(EXIT_FAILURE);
                }
            }

            if ((num & 3) == 1) {
                strcpy(buffer, buffer + idx + 1);
                shift = bytes_read + shift - idx - 1;
            } else {
                shift = 0;
            }
        }
        delete[] k_mer_counter_map;
    } else {
        CounterMap_128* k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];

        while (1) {
            bytes_read = (int) fread(buffer + shift, 1, LENGTH - 1 - shift, fp);
            buffer[bytes_read + shift] = '\0';

            for (int i = 0; i < bytes_read + shift; i++) {
                if (buffer[i] == '\n') {
                    num += 1;
                    if ((num & 3) == 2) {
                        if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                            fprintf(stderr, "This program is designed for short-read sequencing\n");
                            exit(EXIT_FAILURE);
                        }

                        k_mer_128(buffer, idx + 1, i - 1, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                  k_mer_data_128, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);
                    }
                    idx = i;
                }
            }

            if (bytes_read <= 0) {
                if (feof(fp)) {
                    break;
                } else {
                    fprintf(stderr, "File-IO Error: %s.\n", strerror(errno));
                    exit(EXIT_FAILURE);
                }
            }

            if ((num & 3) == 1) {
                strcpy(buffer, buffer + idx + 1);
                shift = bytes_read + shift - idx - 1;
            } else {
                shift = 0;
            }
        }
        delete[] k_mer_counter_map;
    }

    free(k_mer_total_cnt);
    free(buffer);
    return result;
}

void read_fastq_thread(FILE* fp, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bytes_read;
    int idx;
    char* buffer;
    char* buffer_new;

    buffer = (char*)malloc(sizeof(char) * LENGTH);
    while (1) {
        LocationVector* loc_vector = new LocationVector{};

        bytes_read = (int)fread(buffer + shift, 1, LENGTH - 1 - shift, fp);
        buffer[bytes_read + shift] = '\0';

        for(int i = 0; i < bytes_read + shift; i++) {
            if(buffer[i] == '\n') {
                num += 1;
                if((num & 3) == 2){
                    if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                        fprintf(stderr, "This program is designed for short-read sequencing\n");
                        exit (EXIT_FAILURE);
                    }

                    loc_vector -> emplace_back(idx + 1, i - 1);
                }
                idx = i;
            }
        }

        buffer_task_queue -> push(QueueData{buffer, loc_vector});
        if (bytes_read <= 0) {
            if (feof (fp)){
                break;
            }
            else {
                fprintf(stderr, "File-IO Error: %s.\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
        }

        buffer_new = (char*)malloc(sizeof(char) * LENGTH);
        if ((num & 3) == 1) {
            strcpy(buffer_new, buffer + idx + 1);
            shift = bytes_read + shift - idx - 1;
        }
        else {
            shift = 0;
        }

        buffer = buffer_new;
    }
}

ResultMap* read_fastq_gz(gzFile fp, ThreadData* thread_data, uint32_t** rot_table, uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr,  "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMap* result = new ResultMap {};

    int num = 0;
    int shift = 0;
    int bytes_read;
    int idx;
    char* buffer = (char*)malloc(sizeof(char) * LENGTH);

    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        CounterMap* k_mer_counter_map = nullptr;
        if (TABLE_MAX_MER < MAX_MER) {
            k_mer_counter_map = new CounterMap[MAX_MER - TABLE_MAX_MER];
        }

        while (1) {
            bytes_read = gzread(fp, buffer + shift, LENGTH - 1 - shift);
            buffer[bytes_read + shift] = '\0';

            for(int i = 0; i < bytes_read + shift; i++) {
                if(buffer[i] == '\n') {
                    num += 1;
                    if((num & 3) == 2){
                        if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                            fprintf(stderr, "This program is designed for short-read sequencing\n");
                            exit (EXIT_FAILURE);
                        }

                        k_mer(buffer, idx + 1, i - 1, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                              k_mer_data, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);
                    }
                    idx = i;
                }
            }

            if (bytes_read < LENGTH - 1 - shift) {
                if (gzeof(fp)){
                    break;
                }
                else {
                    fprintf(stderr, "File-IO Error: %s.\n", strerror(errno));
                    exit(EXIT_FAILURE);
                }
            }

            if ((num & 3) == 1) {
                strcpy(buffer, buffer + idx + 1);
                shift = bytes_read + shift - idx - 1;
            }
            else {
                shift = 0;
            }

        }
        delete[] k_mer_counter_map;
    } else {
        CounterMap_128* k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];

        while (1) {
            bytes_read = gzread(fp, buffer + shift, LENGTH - 1 - shift);
            buffer[bytes_read + shift] = '\0';

            for (int i = 0; i < bytes_read + shift; i++) {
                if (buffer[i] == '\n') {
                    num += 1;
                    if ((num & 3) == 2) {
                        if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                            fprintf(stderr, "This program is designed for short-read sequencing\n");
                            exit(EXIT_FAILURE);
                        }

                        k_mer_128(buffer, idx + 1, i - 1, rot_table, extract_k_mer_128, k_mer_counter,
                                  k_mer_counter_map,
                                  k_mer_data_128, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);
                    }
                    idx = i;
                }
            }

            if (bytes_read < LENGTH - 1 - shift) {
                if (gzeof(fp)) {
                    break;
                } else {
                    fprintf(stderr, "File-IO Error: %s.\n", strerror(errno));
                    exit(EXIT_FAILURE);
                }
            }

            if ((num & 3) == 1) {
                strcpy(buffer, buffer + idx + 1);
                shift = bytes_read + shift - idx - 1;
            } else {
                shift = 0;
            }

        }
        delete[] k_mer_counter_map;
    }

    free(k_mer_total_cnt);
    free(buffer);
    return result;
}

void read_fastq_gz_thread(gzFile fp, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bytes_read;
    int idx;
    char* buffer;
    char* buffer_new;

    buffer = (char*)malloc(sizeof(char) * LENGTH);
    while (1) {
        LocationVector* loc_vector = new LocationVector{};

        bytes_read = gzread(fp, buffer + shift, LENGTH - 1 - shift);
        buffer[bytes_read + shift] = '\0';

        for(int i = 0; i < bytes_read + shift; i++) {
            if(buffer[i] == '\n') {
                num += 1;
                if((num & 3) == 2){
                    if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                        fprintf(stderr, "This program is designed for short-read sequencing\n");
                        exit (EXIT_FAILURE);
                    }

                    loc_vector -> emplace_back(idx + 1, i - 1);
                }
                idx = i;
            }
        }

        buffer_task_queue -> push(QueueData{buffer, loc_vector});
        if (bytes_read < LENGTH - 1 - shift) {
            if (gzeof(fp)){
                break;
            }
            else {
                fprintf(stderr, "File-IO Error: %s.\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
        }

        buffer_new = (char*)malloc(sizeof(char) * LENGTH);
        if ((num & 3) == 1) {
            strcpy(buffer_new, buffer + idx + 1);
            shift = bytes_read + shift - idx - 1;
        }
        else {
            shift = 0;
        }

        buffer = buffer_new;
    }
}

void process_kmer(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                  uint64_t *extract_k_mer, uint128_t *extract_k_mer_128,
                  ThreadData* thread_data_list, bool is_gz) {
    ResultMap** result_list = (ResultMap**) malloc(sizeof(ResultMap*) * NUM_THREAD);

    TBBQueue buffer_task_queue {};
    TaskGroup tasks;

    if (NUM_THREAD > 1) {
        for (int i = 0; i < NUM_THREAD - 1; i++) {
            tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
                result_list[i] = buffer_task(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
            });
        }
    }

    if (is_gz) {
        gzFile fp = gzopen(file_name, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        if (NUM_THREAD > 1) {
            read_fastq_gz_thread(fp, &buffer_task_queue);
        } else {
            result_list[0] = read_fastq_gz(fp, thread_data_list, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        }

        gzclose(fp);
    } else {
        FILE* fp = fopen(file_name, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        if (NUM_THREAD > 1) {
            read_fastq_thread(fp, &buffer_task_queue);
        } else {
            result_list[0] = read_fastq(fp, thread_data_list, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        }

        fclose(fp);
    }

    if (NUM_THREAD > 1) {
        for (int i = 1; i < NUM_THREAD; i ++) {
            buffer_task_queue.push(QueueData{nullptr, nullptr});
        }

        if (!buffer_task_queue.empty()) {
            buffer_task_queue.push(QueueData{nullptr, nullptr});

            int i = NUM_THREAD - 1;
            tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
                result_list[i] = buffer_task(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
            });
        } else {
            result_list[NUM_THREAD - 1] = new ResultMap();
        }

        tasks.wait();
    }

    ResultMap* result = result_list[0];
    for (int i = 1; i < NUM_THREAD; i++) {
        for (auto& [k, v] : *(result_list[i])) {
            (*result)[k] += v;
        }
        delete result_list[i];
    }

    std::vector<std::pair<KmerSeq, uint32_t>> result_vector(result -> begin(), result -> end());
    delete result;

    std::sort(result_vector.begin(), result_vector.end(), [](auto &a, auto &b) {
        return a.second > b.second;
    });

    char buffer[ABS_MAX_MER + 1];
    for (auto& [k, v] : result_vector) {
        int_to_four(buffer, k.second, k.first);
        fprintf(stdout, "%d,%s,%u\n", k.first, buffer, v);
    }
}

uint8_t** set_repeat_check_table() {
    uint8_t** repeat_check_table = (uint8_t**) malloc(sizeof(uint8_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1));
    if (repeat_check_table == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    uint32_t max_ind;
    for (int i = 0; i <= MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER; i++) {
        max_ind = 1 << (2 * (i + MIN_MER));

        repeat_check_table[i] = (uint8_t*) calloc(max_ind, sizeof(uint16_t));
        if (repeat_check_table[i] == nullptr){
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }

    return repeat_check_table;
}


uint64_t* set_extract_k_mer() {
    uint64_t* extract_k_mer = (uint64_t*) malloc(sizeof(uint64_t) * (MAX_MER - MIN_MER + 1));
    for (int i = 0; i <= MAX_MER - MIN_MER; i ++) {
        extract_k_mer[i] = MIN_MER + i < 32 ? (1ULL << (2 * (i + MIN_MER))) - 1 : -1ULL;
    }
    return extract_k_mer;
}

uint128_t* set_extract_k_mer_128() {
    uint128_t* extract_k_mer = (uint128_t*) malloc(sizeof(uint128_t) * (MAX_MER - MIN_MER + 1));
    for (int i = 0; i <= MAX_MER - MIN_MER; i ++) {
        extract_k_mer[i] = MIN_MER + i < 64 ? (uint128_t(1) << (2 * (i + MIN_MER))) - 1 : uint128_t(-1);
    }
    return extract_k_mer;
}

uint32_t** set_rotation_table(uint8_t** repeat_check_table) {
    uint32_t** rot_table = (uint32_t**) malloc(sizeof(uint32_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1));
    if (rot_table == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    uint32_t max_ind;
    for (int i = 0; i <= MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER; i++) {
        max_ind = 1 << (2 * (i + MIN_MER));

        rot_table[i] = (uint32_t*) calloc(max_ind, sizeof(uint32_t));
        if (rot_table[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }

        for (uint32_t j = 0; j < max_ind; j++) {
            if (rot_table[i][j] == 0) {
                fill_rotation_table(rot_table[i], i + MIN_MER, j, repeat_check_table);
            }
        }
    }

    return rot_table;
}

uint16_t** set_k_mer_counter() {
    uint16_t** cnt = (uint16_t**) malloc(sizeof(uint16_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1));
    if (cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    uint32_t max_ind;
    for (int i = 0; i <= MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER; i++) {
        max_ind =  1 << (2 * (i + MIN_MER));

        cnt[i] = (uint16_t*) calloc(max_ind, sizeof(uint16_t));
        if (cnt[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }

    return cnt;
}

uint64_t** set_k_mer_data() {
    uint64_t** k_mer_data = (uint64_t**) malloc(sizeof(uint64_t*) * (MAX_MER - MIN_MER + 1));
    if (k_mer_data == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < MAX_MER - MIN_MER + 1; i++) {
        k_mer_data[i] = (uint64_t*) calloc(3, sizeof(uint64_t));
        if (k_mer_data[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }

    return k_mer_data;
}

uint128_t** set_k_mer_data_128() {
    uint128_t** k_mer_data = (uint128_t**) malloc(sizeof(uint128_t*) * (MAX_MER - MIN_MER + 1));
    if (k_mer_data == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < MAX_MER - MIN_MER + 1; i++) {
        k_mer_data[i] = (uint128_t*) calloc(3, sizeof(uint128_t));
        if (k_mer_data[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }

    return k_mer_data;
}

uint32_t** set_k_mer_counter_list() {
    uint32_t** k_mer_counter_list = (uint32_t**) malloc(sizeof(uint32_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1));
    if (k_mer_counter_list == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1; i++) {
        k_mer_counter_list[i] = (uint32_t*) calloc(MAX_SEQ, sizeof(uint32_t));
        if (k_mer_counter_list[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }
    return k_mer_counter_list;
}

void fill_rotation_table(uint32_t* rot, int k, uint32_t seq, uint8_t** repeat_check_table) {
    uint32_t code;
    uint32_t tmp = seq;
    uint32_t cnt = 0;
    uint16_t nut_cnt[4]={0, 0, 0, 0};
    rot[tmp] = seq;
    for (int i = 0; i < k - 1; i++) {
        code = tmp & 0x3;
        nut_cnt[code]++;

        tmp = ((tmp & 0x3) << (2 * (k - 1))) + (tmp >> 2);
        rot[tmp] = seq;
    }
    nut_cnt[tmp & 0x3]++;

    for (uint16_t i : nut_cnt) {
        if (i > 0) {
            cnt++;
        }
    }

    repeat_check_table[k - MIN_MER][seq] = cnt <= 2 ? 1 : 0;
}

uint64_t get_rot_seq(uint64_t seq, int k) {
    uint64_t tmp = seq;
    uint64_t ans = seq;
    for (int i = 0; i < k - 1; i++) {
        tmp = ((tmp & 0x3) << (2 * (k - 1))) + (tmp >> 2);
        ans = MIN(tmp, ans);
    }
    return ans;
}

uint128_t get_rot_seq_128(const uint128_t& seq, int k) {
    uint128_t tmp = seq;
    uint128_t ans = seq;
    for (int i = 0; i < k - 1; i++) {
        tmp = ((tmp & 0x3) << (2 * (k - 1))) + (tmp >> 2);
        ans = MIN(tmp, ans);
    }
    return ans;
}

uint8_t get_repeat_check(uint64_t seq, int k) {
    uint32_t cnt = 0;
    uint16_t nut_cnt[4]={0, 0, 0, 0};

    nut_cnt[seq & 0x3]++;
    for (int i = 0; i < k - 1; i++, seq >>= 2) {
        nut_cnt[seq & 0x3]++;
    }

    for (uint16_t i : nut_cnt) {
        if (i > 0) {
            cnt++;
        }
    }

    return cnt <= 2 ? 1 : 0;
}

uint8_t get_repeat_check(uint128_t seq, int k) {
    uint32_t cnt = 0;
    uint16_t nut_cnt[4]={0, 0, 0, 0};

    nut_cnt[(uint8_t) (seq & 0x3)]++;
    for (int i = 0; i < k - 1; i++, seq >>= 2) {
        nut_cnt[(uint8_t) (seq & 0x3)]++;
    }

    for (uint16_t i : nut_cnt) {
        if (i > 0) {
            cnt++;
        }
    }

    return cnt <= 2 ? 1 : 0;
}

void int_to_four(char* buffer, uint128_t seq, int n) {
    for(int i = 0; i < n; i++) {
        buffer[n - 1 - i] = trans_arr[(uint8_t) (seq & 0x3)];
        seq >>= 2;
    }
    buffer[n] = '\0';
}

void k_mer(const char* seq, int st, int nd, uint32_t** rot_table, const uint64_t *extract_k_mer,
           uint16_t** k_mer_counter, CounterMap* k_mer_counter_map,
           uint64_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
           ResultMap* result,
           int16_t* k_mer_total_cnt) {

    memset(k_mer_total_cnt, 0, sizeof(int16_t) * (MAX_MER - MIN_MER + 2));

    int len;
    int chr_code;
    int err_loc;
    uint16_t cur_cnt;
    uint32_t estimated_count;
    uint64_t tmp;
    uint64_t val;

    double target_frequency = 0;
    int target_k;
    double frequency;

    err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            len = MIN(i - err_loc, MAX_MER);
            if(len >= MIN_MER){
                k_mer_total_cnt[0]++;
                k_mer_total_cnt[len - MIN_MER + 1]--;
            }
        } else {
            err_loc = i;
        }
    }

    for(int i = 1; i <= MAX_MER - MIN_MER; i++) {
        k_mer_total_cnt[i] = (int16_t)(k_mer_total_cnt[i] + k_mer_total_cnt[i - 1]);
    }

    for (int k = MIN_MER; k <= MAX_MER; k++) {
        tmp = 0; err_loc = st - 1;
        for (int i = st ; i <= nd; i++) {
            tmp <<= 2;
            chr_code = codes[seq[i]];
            if (chr_code >= 0) {
                tmp += chr_code;
                if(i - err_loc >= k) {
                    if (k <= TABLE_MAX_MER) {
                        val = rot_table[k - MIN_MER][tmp & extract_k_mer[k - MIN_MER]];
                        cur_cnt = ++k_mer_counter[k - MIN_MER][val];

                        k_mer_counter_list[k - MIN_MER][k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++] = val;
                    } else {
                        val = get_rot_seq(tmp & extract_k_mer[k - MIN_MER], k);
                        cur_cnt = ++k_mer_counter_map[k - TABLE_MAX_MER - 1][val];
                        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++;
                    }

                    if (k_mer_data[k - MIN_MER][K_MER_DATA_MAX] < cur_cnt) {
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = cur_cnt;
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ] = val;
                    }

                    estimated_count = k_mer_data[k - MIN_MER][K_MER_DATA_MAX] + k_mer_total_cnt[k - MIN_MER] - k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
                    if(estimated_count < k_mer_total_cnt[k - MIN_MER] * BASELINE){
                        break;
                    }
                }
            } else {
                err_loc = i;
            }
        }
    }


    for (int k = MIN_MER; k <= MAX_MER; k++) {
        frequency = (double)k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double)k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
        if(frequency > target_frequency && !(k <= TABLE_MAX_MER ? repeat_check_table[k - MIN_MER][k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k))){
            target_frequency = frequency;
            target_k = k;
        }
    }

    if(target_frequency >= BASELINE){
        for (int k = MIN_MER; k <= MAX_MER; k++) {
            if (target_k == k) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result)[KmerSeq {k, MIN(rot_table[k - MIN_MER][reverse_complement(k_mer_counter_list[k - MIN_MER][i]) >> (2 * (16 - k))], k_mer_counter_list[k - MIN_MER][i])}] += val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result)[KmerSeq {k, MIN(key, get_rot_seq(reverse_complement(key >> (2 * (32 - k))), k))}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                    }
                }
                else {
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
    }
    else {
        for (int k = MIN_MER; k <= MAX_MER; k++) {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                }
            }
            else {
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
    }
}

void k_mer_128(const char* seq, int st, int nd, uint32_t** rot_table, const uint128_t *extract_k_mer,
               uint16_t** k_mer_counter, CounterMap_128* k_mer_counter_map,
               uint128_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
               ResultMap* result,
               int16_t* k_mer_total_cnt) {

    memset(k_mer_total_cnt, 0, sizeof(int16_t) * (MAX_MER - MIN_MER + 2));

    int len;
    int chr_code;
    int err_loc;
    uint16_t cur_cnt;
    uint32_t estimated_count;
    uint128_t tmp;
    uint128_t val;

    double target_frequency = 0;
    int target_k;
    double frequency;

    err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            len = MIN(i - err_loc, MAX_MER);
            if(len >= MIN_MER){
                k_mer_total_cnt[0]++;
                k_mer_total_cnt[len - MIN_MER + 1]--;
            }
        } else {
            err_loc = i;
        }
    }

    for(int i = 1; i <= MAX_MER - MIN_MER; i++) {
        k_mer_total_cnt[i] = (int16_t)(k_mer_total_cnt[i] + k_mer_total_cnt[i - 1]);
    }

    for (int k = MIN_MER; k <= MAX_MER; k++) {
        tmp = 0; err_loc = st - 1;
        for (int i = st ; i <= nd; i++) {
            tmp <<= 2;
            chr_code = codes[seq[i]];
            if (chr_code >= 0) {
                tmp += chr_code;
                if(i - err_loc >= k) {
                    if (k <= TABLE_MAX_MER) {
                        val = (uint64_t) rot_table[k - MIN_MER][(uint32_t) (tmp & extract_k_mer[k - MIN_MER])];
                        cur_cnt = ++k_mer_counter[k - MIN_MER][(uint32_t) val];

                        k_mer_counter_list[k - MIN_MER][(uint32_t) (k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++)] = (uint32_t) val;
                    } else {
                        val = get_rot_seq_128(tmp & extract_k_mer[k - MIN_MER], k);
                        cur_cnt = ++k_mer_counter_map[k - TABLE_MAX_MER - 1][val];
                        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++;
                    }

                    if (k_mer_data[k - MIN_MER][K_MER_DATA_MAX] < cur_cnt) {
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = cur_cnt;
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ] = val;
                    }

                    estimated_count = (uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] + k_mer_total_cnt[k - MIN_MER] - (uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
                    if(estimated_count < k_mer_total_cnt[k - MIN_MER] * BASELINE){
                        break;
                    }
                }
            } else {
                err_loc = i;
            }
        }
    }

    for (int k = MIN_MER; k <= MAX_MER; k++) {
        frequency = (double) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
        if (frequency > target_frequency && !(k <= TABLE_MAX_MER ? repeat_check_table[k - MIN_MER][(uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k))){
            target_frequency = frequency;
            target_k = k;
        }
    }

    if(target_frequency >= BASELINE){
        for (int k = MIN_MER; k <= MAX_MER; k++) {
            if (target_k == k) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result)[KmerSeq {k, MIN(rot_table[k - MIN_MER][reverse_complement(k_mer_counter_list[k - MIN_MER][i]) >> (2 * (16 - k))], k_mer_counter_list[k - MIN_MER][i])}] += (uint32_t) val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result)[KmerSeq {k, MIN(key, get_rot_seq_128(reverse_complement_128(key) >> (2 * (64 - k)), k))}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                    }
                }
                else {
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
    }
    else {
        for (int k = MIN_MER; k <= MAX_MER; k++) {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                }
            }
            else {
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
    }
}