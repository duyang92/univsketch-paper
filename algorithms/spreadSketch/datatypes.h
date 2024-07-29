#ifndef DATATYPE_H
#define DATATYPE_H

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../../utils/MurmurHash3.h"
#include "string.h"

typedef uint32_t key_tp;

typedef uint32_t val_tp;

typedef struct edge_t_s {
    uint32_t src_ip;
    uint32_t dst_ip;
} edge_tp;


/**
 * Object for hash
 */
typedef struct {
    /// overloaded operation
    long operator() (const edge_tp &k) const { 
        uint32_t hash_val = 0;
        uint64_t t = ((uint64_t)k.src_ip) << 32 | k.dst_ip;
        char hash_input_str[9] = {0};
        memcpy(hash_input_str, &t, sizeof(uint64_t));
        MurmurHash3_x86_32(hash_input_str, 8, 388650253, &hash_val);
        return hash_val;
    }
} edge_tp_hash;


/**
 * Object for equality
 */
typedef struct {
    /// overloaded operation
    bool operator() (const edge_tp &x, const edge_tp &y) const {
        return memcmp(&x, &y, 8)==0;
    }
} edge_tp_eq;

typedef std::unordered_set<key_tp> myset;

typedef std::unordered_map<key_tp, myset> mymap;

typedef std::unordered_map<key_tp, val_tp> nodemap;

typedef std::vector<std::pair<key_tp, val_tp> > myvector;

typedef std::unordered_set<edge_tp, edge_tp_hash, edge_tp_eq> edgeset;

#endif
