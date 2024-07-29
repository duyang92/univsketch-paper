#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "../algorithms/m2d/m2d.h"
#include "../algorithms/rskt/rskt.h"
#include "../algorithms/spreadSketch/spreadsketch.h"
#include "../algorithms/univ/univsketch.h"
#include "../algorithms/vbitmap/vbitmap.h"

using namespace std;

typedef std::pair<uint32_t, uint32_t> pii;

vector<vector<pii>> data_set;
unordered_set<uint32_t> flow_set;
ifstream inf;
ofstream ouf;

void data_preprocess();
void read_packet_data();
void generate_hash_seeds(int len);
void free();
void process_args(int argc, char** argv);

int data_size;
int opt = 0;
int T_MEM = 2048;  // total size: 2048KB

uint32_t* hash_seeds;
pii* data_arr;

void test_vBitmap_ss() {
    printf("===============> Test vBitmap + ss\n");

    // vBitmap
    int m = int(1.0 * 8 * 1024 * T_MEM / 2);
    int s = 5000;
    vBitmap vb(m, s, hash_seeds);

    // ss
    unsigned long long buf_size = 500000000;
    int lgn = 32, cmdepth = 4, b = 79, c = 3, memory = 438;
    int cmwidth = 1.0 * T_MEM * 1024 * 8 / 2 / cmdepth / (memory + lgn + 8);
    DetectorSS ss(cmdepth, cmwidth, lgn, b, c, memory);

    clock_t start, end;
    double time, throughput;
    // insert
    start = clock();
    for (int i = 0; i < data_size; i++) {
        vb.insert(data_arr[i].first, data_arr[i].second);
        ss.Update(data_arr[i].first, data_arr[i].second, 1);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = data_size * 1.0 / 1000000 / time;
    printf("Insert throughput: %.6f Mpps\n", throughput);

    // query per-flow
    unordered_map<uint32_t, double> vb_result;
    start = clock();
    double sum_bits = vb.sum_bits();
    for (const uint32_t& flow : flow_set) {
        vb_result[flow] = vb.query_per_flow(flow, sum_bits);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = time * 1000 / flow_set.size();
    printf("Per-flow query speed: %.6f ms\n", throughput);
}

void test_vBitmap_ss_rskt() {
    printf("===============> Test vBitmap + ss + rskt\n");

    // vBitmap
    int m = int(1.0 * 8 * 1024 * T_MEM / 3);
    int s = 5000;
    vBitmap vb(m, s, hash_seeds);

    // ss
    unsigned long long buf_size = 500000000;
    int lgn = 32, cmdepth = 4, b = 79, c = 3, memory = 438;
    int cmwidth = 1.0 * 1024 * 8 * T_MEM / 3 / cmdepth / (memory + lgn + 8);
    DetectorSS ss(cmdepth, cmwidth, lgn, b, c, memory);

    // rskt
    m = 128;
    int HLL_size = 5;
    int w = int(1.0 * 1024 * 8 * T_MEM / 3 / HLL_size / m / 2);
    RSKT rskt(w, m, hash_seeds);

    clock_t start, end;
    double time, throughput;
    // insert
    start = clock();
    for (int i = 0; i < data_size; i++) {
        vb.insert(data_arr[i].first, data_arr[i].second);
        ss.Update(data_arr[i].first, data_arr[i].second, 1);
        rskt.insert(data_arr[i].first, data_arr[i].second);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = data_size * 1.0 / 1000000 / time;
    printf("Insert throughput: %.6f Mpps\n", throughput);

    // query per-flow
    unordered_map<uint32_t, double> rskt_result;
    start = clock();
    double sum_bits = vb.sum_bits();
    for (const uint32_t& flow : flow_set) {
        rskt_result[flow] = rskt.query_per_flow(flow);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = time * 1000 / flow_set.size();
    printf("Per-flow query speed: %.6f ms\n", throughput);
}

void test_m2d() {
    printf("===============> Test M2D\n");

    // m2d
    int hashTableMemory = 132, l = 3, heapSize = 100, d = 128;
    double remian = T_MEM * 8 * 1024 - hashTableMemory * 8 * 1024 - heapSize * 64 * (l + 1);
    double ps = 0.25;
    int w = int(remian / (d * 5 + 32 * 32) / (1 + ps + ps * ps));
    pii* dw = new pii[3]{pii(d, max(2, w)), pii(d, max(2, int(w * ps))), pii(d, max(2, int(w * ps * ps)))};
    M2D m2(l, dw, 3, heapSize, ps, hash_seeds);

    clock_t start, end;
    double time, throughput;
    // insert
    start = clock();
    for (int i = 0; i < data_size; i++) {
        m2.insert(data_arr[i].first, data_arr[i].second);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = data_size * 1.0 / 1000000 / time;
    printf("Insert throughput: %.6f Mpps\n", throughput);

    // query per-flow
    unordered_map<uint32_t, double> m2_result;
    start = clock();
    for (const uint32_t& flow : flow_set) {
        m2_result[flow] = m2.query_per_flow(flow);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = time * 1000 / flow_set.size();
    printf("Per-flow query speed: %.6f ms\n", throughput);
}

void test_univsketch() {
    printf("===============> Test univSketch\n");

    // univsketch
    double register_ratio = 7, bucket_ratio = 3;
    int l = 16, bucket_size = 4, level = 4, s = 128;
    int m = int(1.0 * T_MEM * 1024 * 8 * (register_ratio / (register_ratio + bucket_ratio)) / l);
    int m_ = int(1.0 * T_MEM * 1024 * 8 * (bucket_ratio / (register_ratio + bucket_ratio)) / bucket_size / (32 + 32));
    univSketch univ(m, l, s, m_, level, bucket_size, hash_seeds);

    clock_t start, end;
    double time, throughput;
    // insert
    start = clock();
    for (int i = 0; i < data_size; i++) {
        univ.insert(data_arr[i].first, data_arr[i].second);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = data_size * 1.0 / 1000000 / time;
    printf("Insert throughput: %.6f Mpps\n", throughput);

    // query per-flow
    unordered_map<uint32_t, double> univ_result;
    start = clock();
    for (const uint32_t& flow : flow_set) {
        univ_result[flow] = univ.query_per_flow(flow);
    }
    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    throughput = 1000 * time / flow_set.size();
    printf("Per-flow query speed: %.6f ns\n", throughput);
}

int main(int argc, char** argv) {
    process_args(argc, argv);

    generate_hash_seeds(4);

    read_packet_data();

    // vBitmap + ss
    test_vBitmap_ss();
    // vBitmap + ss + rskt
    test_vBitmap_ss_rskt();
    // m2d
    test_m2d();
    // univsketch
    test_univsketch();

    free();

    return 0;
}

void read_packet_data() {
    string path = "../data/";
    data_set.clear();

    pii temp;
    int total = 0;

    string fs[] = {"00.txt"};
    for (string f : fs) {
        vector<pii> data;
        string c_path = path + f;
        inf.open(c_path, ios::in);
        cout << "Reading " << c_path << endl;
        int cc = 0;
        while (inf >> temp.second && inf >> temp.first) {
            data.push_back({temp.first, temp.second});
            flow_set.insert(temp.first);
            cc += 1;
        }
        cout << "Read " << cc << " lines" << endl;
        data_set.push_back(data);
        total += cc;
        inf.close();
        // break;
    }

    data_size = total;

    cout << "Total data size: " << total << endl;

    data_preprocess();
}

void data_preprocess() {
    data_arr = new pii[data_size];
    memset(data_arr, 0, data_size * sizeof(pii));
    int index = 0;
    for (const auto& data : data_set) {
        for (const auto& d : data) {
            data_arr[index++] = d;
        }
    }
}

void generate_hash_seeds(int len = 8) {
    hash_seeds = new uint32_t[len];

    int count = 0;
    std::unordered_set<int> diff_ele;
    while (count < len) {
        int num = rand();
        if (diff_ele.find(num) == diff_ele.end()) {
            diff_ele.insert(num);
            hash_seeds[count++] = num;
        }
    }
}

void process_args(int argc, char** argv) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << "<Memory>" << endl;
        exit(1);
    }

    T_MEM = atoi(argv[1]);

    printf("Memory: %d\n", T_MEM);
}

void free() {
    delete[] data_arr;
    delete[] hash_seeds;
}