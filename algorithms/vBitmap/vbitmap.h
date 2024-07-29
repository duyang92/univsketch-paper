#ifndef VBITMAP_H
#define VBITMAP_H

#include <string.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../utils/MurmurHash3.h"

class vBitmap {
    typedef std::pair<uint32_t, uint32_t> pii;

   public:
    int m;
    int s;
    uint32_t* hash_seeds = nullptr;
    uint32_t* M = nullptr;

   public:
    uint32_t H(uint32_t t, uint32_t s = 31) {
        uint32_t hash_val = 0;
        char hash_input_str[5] = {0};
        memcpy(hash_input_str, &t, sizeof(uint32_t));
        MurmurHash3_x86_32(hash_input_str, 4, s, &hash_val);
        return hash_val;
    }

    uint32_t H64(uint32_t t1, uint32_t t2, uint32_t s = 31) {
        uint32_t hash_val = 0;
        uint64_t t = ((uint64_t)t1) << 32 | t2;
        char hash_input_str[9] = {0};
        memcpy(hash_input_str, &t, sizeof(uint64_t));
        MurmurHash3_x86_32(hash_input_str, 8, s, &hash_val);
        return hash_val;
    }

    vBitmap(int m_val, int s_val, uint32_t* hash_seeds_)
        : m(m_val), s(s_val) {
        hash_seeds = hash_seeds_;
        M = new uint32_t[m];
        memset(M, 0, m * sizeof(uint32_t));
    }

    ~vBitmap() {
        if (M != nullptr) {
            delete[] M;
        }
    }

    void insert(uint32_t flow, uint32_t element) {
        uint32_t i = H(element, hash_seeds[0]) % s;
        size_t m_value = H64(flow, i, hash_seeds[1]) % m;

        M[m_value] = 1;
    }

    double sum_bits(){
        double sum_bits = 0;
        for (size_t i = 0; i < m; ++i) {
            sum_bits += M[i];
        }
        return sum_bits;
    }

    double query_per_flow(const uint32_t& flow, double sum_bits) {

        double k1 = -s * log(std::max(1.0, (double)(m - sum_bits)) / m);

        int U0_0 = 0;
        for (size_t i = 0; i < s; ++i) {
            size_t m_value = H64(flow, i, hash_seeds[1]) % m;
            U0_0 += M[m_value];
        }

        if (U0_0 == s) {
            U0_0 = 1;
        } else {
            U0_0 = s - U0_0;
        }

        double k2 = -s * log((double)U0_0 / s);
        double spread = -k1 + k2;
        spread = std::max(1.0, spread);

        return spread;
    }
};

#endif