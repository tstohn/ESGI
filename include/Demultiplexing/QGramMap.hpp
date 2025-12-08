#pragma once

#include <iostream>
#include <string>
#include <tuple>
#include <cmath>
#include <vector>
#include <memory>

#include "helper.hpp"
#include "ankerl/unordered_dense.h"


class QGramMap
{

    public:

        QGramMap() = default;          // default constructor does very little
        bool useQgramMap = true;

        // Constructor:
        //  - barcodes: all possible barcodes
        // we calcualte the length of the qgrams given mismatches (chosen so that given MM at least 2 qgrams must match)
        inline void init(int inMismatches,
                         std::vector<std::shared_ptr<std::string>>& inBarcodes)
        {
            barcodes = inBarcodes;
            mismatches = inMismatches;

            //set length if qgrams and minimum overlap
            qgramLength = std::floor( (barcodes.at(0)->size() - mismatches) / (mismatches + 1));
            minQgramOverlap = ( barcodes.at(0)->size() - qgramLength + 1) - (mismatches * qgramLength);
            qgramPositions = barcodes.at(0)->size() - qgramLength + 1;

            POSITION_BITS = compute_position_bits(qgramPositions);
            QMER_BITS = 2 * qgramLength; // 2 bits per base (ACGT)

            POSITION_MASK = (1ULL << POSITION_BITS) - 1ULL;

            useQgramMap = use_qgram_map();

            // Build the index from all barcodes
            if(useQgramMap)
            {
                build_index();
                std::cout << "\t\tUsing qgrams of length " << std::to_string(qgramLength) << \
                " to filter barcodes that share at least " << std::to_string(minQgramOverlap) << \
                " qgrams with the potential barcode-sequence\n";
            }
        }

    // query: calcualte all qmers for a sequence seq, 
    //get all barcodes that share a qmer and filter the barcode that have a minimalQmer overlap
std::vector<std::shared_ptr<std::string>> query(const std::string& seq) const
{

    //the final result of all potential barcodes that share at least <minQgramOverlap> qgrams
    std::vector<std::shared_ptr<std::string>> potentialBarcodes;
    //temporary structure to store the potential barcodes with the number of qmers (counts) that they share with seq
    //key is an index of the barcode in barcodes
    ankerl::unordered_dense::map<int, int> counts;

    //iterate through the input sequence encoding the positional qmer in key
    for (int pos = 0; pos < static_cast<int>(qgramPositions); ++pos) 
    {
        //iterate over all possible positions due to insertions deletions
        for(int alteredPos = (pos - mismatches); alteredPos <= (pos + mismatches); ++alteredPos)
        {          

            //skip negative positions
            if(alteredPos < 0){continue;}

            //1.) calcualte the key for positional qmer if possible
            uint64_t key;
            if (!encode_and_pack(seq, static_cast<uint32_t>(pos), static_cast<uint32_t>(alteredPos), key)) 
            {
                //invalid base, like N in the sequence. It is essentailly treated as a mismatch,
                //maybe we still have enough matching qmers even without considering this qmer
                continue;
            }

            //2.) map look-up: does this positional qmer exist
            auto it = qgramMap.find(key);
            if (it == qgramMap.end()) 
            {
                continue;
            }
            const std::vector<int>& barcodeList = it->second;
            for (int barcodeIdx : barcodeList) 
            {
                ++counts[barcodeIdx];
                if (counts.at(barcodeIdx) == minQgramOverlap) 
                {
                    potentialBarcodes.push_back(barcodes.at(barcodeIdx));
                }

            }
        }
    }

    return potentialBarcodes;
}

private:

    int qgramLength = 0;
    int minQgramOverlap = 0;
    uint32_t qgramPositions = 0;
    int mismatches;

    uint32_t POSITION_BITS = 0;
    int QMER_BITS = 0;
    uint64_t POSITION_MASK = 0;

    std::vector<std::shared_ptr<std::string>> barcodes;
    //map of key(pos+qgram) to the index of barcode in the barcodes-list above
    ankerl::unordered_dense::map<uint64_t, std::vector<int>> qgramMap;

    //quick check if qgrams make any sense, in case the user supplied too many mismatches 
    //and we do not expect any overlapping qgrams skip this step
    inline bool use_qgram_map() const
    {
        return( (qgramLength > 1) &&
            (minQgramOverlap > 1) &&
            (QMER_BITS + (int)POSITION_BITS <= 64));
    }

    // Encode q-gram at position `pos` of string s, then pack qmer+pos into key.
    //e.g.: string AAAAGGGG with pos=4 creates a key for qgram=GGGG at pos 4 (positions are 0-indexed)
    //pos is the position in the string where we start extracting the qmer
    //alteredPos is the position where we could see this qgram (given certain number of MM)
    inline bool encode_and_pack(const std::string& s,
                                uint32_t pos,
                                uint32_t alteredPos,
                                uint64_t& key) const
    {

        //if we can not create the whole qgram, e.g., at the end of a fastq-file when barcodes are cut-off
        if(s.size() < static_cast<size_t>(qgramLength)){return false;}

        uint64_t code = 0;
        for (int i = 0; i < qgramLength; ++i) 
        {
            char c = s[pos + i];
            uint64_t v;
            switch (c) {
                case 'A': case 'a': v = 0; break;
                case 'C': case 'c': v = 1; break;
                case 'G': case 'g': v = 2; break;
                case 'T': case 't': v = 3; break;
                default:
                    return false;  // invalid base // e.g an N
            }
            code = (code << 2) | v;
        }

        // Pack (qmer, position) into a single uint64 key
        key = (code << POSITION_BITS) | (alteredPos & POSITION_MASK);
        return true;
    }

    // Compute minimal number of bits to store values 0..maxPos
    inline static uint32_t compute_position_bits(uint32_t maxPos)
    {
        if (maxPos <= 0) 
        {
            return 1;
        }
        uint32_t bits = 0;
        while ((1u << bits) <= maxPos)
        {
            ++bits;
        }
        return bits;
    }

    // Build the inverted qgram map: QGRAM -> barcodes with this QGRAM
    inline void build_index()
    {
        qgramMap.clear();

        for (size_t idx = 0; idx < barcodes.size(); ++idx) 
        {
            if (barcodes.at(idx)->size() < static_cast<size_t>(qgramLength)) 
            {
                std::cerr << "Barcode " << *(barcodes.at(idx)) << " is not of the same length as other barcodes: " << std::to_string(qgramLength) << "\n";
            }

            for (uint32_t pos = 0; pos < qgramPositions; ++pos) 
            {
                uint64_t key;
                //pos and alteredPos are the same: for the look-up map we expect qgrams at their exact positions
                if(!encode_and_pack(*(barcodes.at(idx)), static_cast<uint32_t>(pos),static_cast<uint32_t>(pos), key))
                {
                    //while building the qgram map barcodes SHOULD ONLY HAVE VALID BASES
                    std::cerr << " None base-character (A,C,G,T) in barcode " << *(barcodes.at(idx)) << "\n";
                    exit(EXIT_FAILURE);
                }
                qgramMap[key].push_back(idx);
            }
        }

    }

    inline void print_qgram_map()
    {

        std::cout << "==== QGRAM INDEX ====\n";
        auto decode_qmer = [&](uint64_t qmer) {
            std::string s(qgramLength, 'A');
            for (int i = qgramLength - 1; i >= 0; --i) {
                uint64_t v = qmer & 3u;   // last 2 bits
                qmer >>= 2;

                switch (v) {
                    case 0: s[i] = 'A'; break;
                    case 1: s[i] = 'C'; break;
                    case 2: s[i] = 'G'; break;
                    case 3: s[i] = 'T'; break;
                }
            }
            return s;
        };
        for (const auto& kv : qgramMap) {
            uint64_t key = kv.first;
            const auto& postings = kv.second;

            // Decode key → (qmer, pos)
            uint32_t pos  = key & POSITION_MASK;
            uint64_t qmer = key >> POSITION_BITS;

            // Decode qmer to DNA sequence
            std::string qmerSeq = decode_qmer(qmer);

            std::cout << "Key raw=0x" << std::hex << key << std::dec
                    << "  pos=" << pos
                    << "  qmer=" << qmerSeq
                    << "\n";

            std::cout << "  Barcodes:\n";
            for (int idx : postings) {
                std::cout << "    [" << idx << "] " << *(barcodes[idx]) << "\n";
            }
            std::cout << "\n";
        }
        std::cout << "=====================\n";

    }
};