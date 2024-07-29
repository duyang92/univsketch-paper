#ifndef RSKT_H
#define RSKT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <unordered_set>
#include <string.h>

#include "../../utils/MurmurHash3.h"

class RSKT {

    typedef std::pair<uint32_t, uint32_t> pii;

   public:
    int w;
    int m;

    char seeds_type[4] = {'g', 'f', 'e', 'p'};
    std::unordered_map<char, uint32_t> hash_seeds;

    uint32_t** C = nullptr;
    uint32_t** _C = nullptr;

   public:

    uint32_t H(uint32_t t, uint32_t s = 31) {
        uint32_t hash_val = 0;
        char hash_input_str[5] = {0};
        memcpy(hash_input_str, &t, sizeof(uint32_t));
        MurmurHash3_x86_32(hash_input_str, 4, s, &hash_val);
        return hash_val;
    }

    RSKT(int _w, int _m, uint32_t* hash_seeds_)
        : w(_w), m(_m) {

        int dx = 0;
        for (char c : seeds_type) {
            hash_seeds[c] = hash_seeds_[dx++];
        }

        C = new uint32_t*[w];
        for (int i = 0; i < w; ++i) {
            C[i] = new uint32_t[m];
            memset(C[i], 0, m * sizeof(uint32_t));
        }

        _C = new uint32_t*[w];
        for (int i = 0; i < w; ++i) {
            _C[i] = new uint32_t[m];
            memset(_C[i], 0, m * sizeof(uint32_t));
        }
    }

    ~RSKT() {
        if (C != nullptr) {
            for (int i = 0; i < w; ++i) {
                delete[] C[i];
            }
            delete[] C;
        }
        if (_C != nullptr) {
            for (int i = 0; i < w; ++i) {
                delete[] _C[i];
            }
            delete[] _C;
        }
    }

    void insert(uint32_t& flow, uint32_t& element) {
        int hashbits = 32;
        int maxvalue = 31;
        int b = std::ceil(std::log2(m));
        uint32_t right_move = hashbits - b;
        uint32_t left_move = (1u << right_move) - 1;
        
        uint32_t h_e_value = H(element, hash_seeds['e']);

        uint32_t h_f = H(flow, hash_seeds['f']) % w;
        uint32_t h_e = (h_e_value >> right_move) % m;

        uint32_t q = h_e_value & left_move;
        uint32_t num = 0;
        while (q) {
            num += 1;
            q >>= 1;
        }
        num = right_move - num + 1;

        uint32_t g_f_i = H(flow ^ h_e, hash_seeds['g']) % 2;
        if (g_f_i == 0) {
            C[h_f][h_e] = std::max(num, C[h_f][h_e]);
        } else {
            _C[h_f][h_e] = std::max(num, _C[h_f][h_e]);
        }
    }

    double query_per_flow(const uint32_t& flow) {
        double alpha_m;
        if (m == 16){
            alpha_m = 0.673;
        } else if (m == 32) {
            alpha_m = 0.697;
        } else if (m == 64) {
            alpha_m = 0.709;
        } else {
            alpha_m = 0.7213 / (1 + 1.079 / m);
        }

        double small_size = 5.0 / 2 * m;
        double large_size = 1.0 / 30 * std::pow(2, 32);
        double large_temp = std::pow(2, 32);

        auto query = [&](const std::vector<int>& Lf, const std::vector<int>& _Lf) -> std::pair<double, double> {
            double sum_Lf = 0;
            double _sum_Lf = 0;
            for (int i = 0; i < m; ++i) {
                sum_Lf += std::pow(2, -Lf[i]);
                _sum_Lf += std::pow(2, -_Lf[i]);
            }
            double Z = 1.0 / sum_Lf;
            double n = alpha_m * m * m * Z;

            double _Z = 1.0 / _sum_Lf;
            double _n = alpha_m * m * m * _Z;

            if (n <= small_size) {
                int V = 0;
                for (int M : Lf) {
                    if (M == 0) {
                        V += 1;
                    }
                }
                if (V) {
                    n = m * std::log((double)m / V);
                }
            } else if (n > large_size) {
                n = -large_temp * std::log(1.0 - n / large_temp);
            }

            if (_n <= small_size) {
                int _V = 0;
                for (int _M : _Lf) {
                    if (_M == 0) {
                        _V += 1;
                    }
                }
                if (_V) {
                    _n = m * std::log((double)m / _V);
                }
            } else if (_n > large_size) {
                _n = -large_temp * std::log(1.0 - _n / large_temp);
            }

            return {n, _n};
        };

        uint32_t h_f = H(flow, hash_seeds['f']) % w;

        std::vector<int> Lf(m, 0);
        std::vector<int> _Lf(m, 0);
        for (int i = 0; i < m; ++i) {
            uint32_t g_f_i = H(flow ^ i, hash_seeds['g']) % 2;
            if (g_f_i == 0) {
                Lf[i] = C[h_f][i];
                _Lf[i] = _C[h_f][i];
            } else {
                Lf[i] = _C[h_f][i];
                _Lf[i] = C[h_f][i];
            }
        }

        auto [n, _n] = query(Lf, _Lf);
        
        return std::max(int(std::round(n - _n)), 1);
    }
};

#endif