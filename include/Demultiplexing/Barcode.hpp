#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <thread>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <functional>
#include <array>
#include <cmath>

#include "helper.hpp"
#include "ankerl/unordered_dense.h"
#include "QGramMap.hpp"

class Barcode;
typedef std::shared_ptr<Barcode> BarcodePtr;
typedef std::vector<BarcodePtr> BarcodeVector; 
typedef std::shared_ptr<BarcodeVector> BarcodeVectorPtr; 

constexpr bool FORWARD = false;
constexpr bool REVERSE = true;

constexpr int KMER = 2;
using KmerArray = std::array<uint8_t, 16>;

enum class PatternType 
{
    Forward,
    Reverse
};

struct mappingSolution{             
    int seq_start;
    int seq_end;
    int score = INT_MAX;
    std::string realBarcode;
    int differenceInBarcodeLength;
}; 

struct baseNum
{
    int A=0;
    int T=0;
    int C=0;
    int G=0;

    inline bool operator==(const baseNum& other) const 
    {
        return A == other.A && T == other.T &&
               C == other.C && G == other.G;
    }
};

namespace std {
template <>
struct hash<baseNum> {
    std::size_t operator()(baseNum const& b) const noexcept {
        // Pack the four counts into a 64-bit integer.
        uint64_t x = 0;
        x |= (uint64_t)(uint16_t)b.A;
        x |= (uint64_t)(uint16_t)b.T << 16;
        x |= (uint64_t)(uint16_t)b.C << 32;
        x |= (uint64_t)(uint16_t)b.G << 48;

        // Mix function (similar to splitmix64 / Murmur-style finalizer).
        x ^= x >> 33;
        x *= 0xff51afd7ed558ccdULL;
        x ^= x >> 33;
        x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= x >> 33;

        return static_cast<std::size_t>(x);
    }
};
}

// Custom hash function for std::array<int, 16>
struct ArrayHash {
    inline std::size_t operator()(const std::array<int, 16>& arr) const {
        std::size_t h = 0;
        for (int val : arr) {
            h ^= std::hash<int>{}(val) + 0x9e3779b9 + (h << 6) + (h >> 2);  // boost-style hash combine
        }
        return h;
    }
};

//new datatypes
class Barcode
{
    public:
    //per default the length of a barcode is set to 0 (unknown), only for UMIs it must be set
    Barcode(std::string name, int inMismatches, int inLength = 0) : name(name), mismatches(inMismatches), length(inLength) {}
    //virtual destructor, needed to be called when destructing classes that inherit from Barcode
    virtual ~Barcode(){};
    virtual std::shared_ptr<Barcode> clone() const = 0;

    std::string name;
    int mismatches;
    unsigned int length; //length of zero means the length of this barcode is unknown

    //reverse complement is a Barcode function that should be available globally
    static std::string generate_reverse_complement(std::string seq)
    {
        auto lambda = [](const char c) 
        {
            switch (c) 
            {
                case 'A':
                    return 'T';
                case 'G':
                    return 'C';
                case 'C':
                    return 'G';
                case 'T':
                    return 'A';
                case 'N':
                    return 'N';
                default:
                    throw std::domain_error("Invalid nucleotide.");
            }
        };
        std::string newSeq = seq;
        std::transform(seq.crbegin(), seq.crend(), newSeq.begin(), lambda);
        return newSeq;
    }
    //overwritten function to match sequence pattern(s)
    virtual bool align(std::string& matchedBarcode, const std::string& fastqLine, const unsigned int targetOffset,
                       int& targetEnd, int& delNum, int& insNum, int& substNum,
                       bool reverse = false) = 0;
    virtual std::vector<std::shared_ptr<std::string>> get_patterns() = 0;
    virtual bool is_wildcard() = 0;
    virtual bool is_constant() = 0;
    virtual bool is_stop() = 0;
    virtual bool is_dna() = 0;
    virtual bool is_read_end() = 0;

};

class ConstantBarcode final : public Barcode
{

    public:
    ConstantBarcode(std::string inPattern, int inMismatches) : Barcode(inPattern, inMismatches), pattern(inPattern)
    {
        revCompPattern = generate_reverse_complement(pattern);
        // Configure Edlib
        config = edlibNewAlignConfig(
            inMismatches,        // Maximum allowed edit distance
            EDLIB_MODE_SHW,     // Semi-global alignment (deletions in the target at the end are not penalized), 
                                // the developers refer to it also as 'Prefix method'
            EDLIB_TASK_PATH,    // Request full alignment path (M, I, D, S)
            NULL, 0             // No custom alphabet
        );
    }
    std::shared_ptr<Barcode> clone() const override {
        return std::make_shared<ConstantBarcode>(*this);
    }

    inline bool align(std::string& matchedBarcode, const std::string& fastqLine,const unsigned int targetOffset,
               int& targetEnd, int& delNum, int& insNum, int& substNum,
               bool reverse = false)
    {
        bool foundAlignment = false;

        //set the pattern to use for reverse or forward mapping
        std::string usedPattern = pattern;
        if(reverse){usedPattern = revCompPattern;}

        //if we allow for zero mismathces check if we can map it immediately
        std::string exactTarget = fastqLine.substr(targetOffset, usedPattern.length());
        if(exactTarget == usedPattern)
        {
            targetEnd = usedPattern.length();
            matchedBarcode = pattern;
            return true;
        }
        else if(mismatches == 0) //return false if we do not allowe for any mismatches
        {
            return false;
        }

        //get length of substring, length does not depend on reverse/ forward pattern
        std::string target;
        int substringLength = pattern.length()+mismatches;
        if(targetOffset + substringLength > fastqLine.size()){substringLength = fastqLine.size()-targetOffset;};
        //std::cout << "\t LENGTH: " <<substringLength << " seq: " << fastqLine.size()<< "\n";
        target = fastqLine.substr(targetOffset, substringLength);
        
        //map the pattern to the target sequence
        foundAlignment = run_alignment(usedPattern, target, targetEnd, config, delNum,  insNum, substNum);
        matchedBarcode = pattern;

        return foundAlignment;
    }

    inline std::vector<std::shared_ptr<std::string>> get_patterns()
    {
        std::vector<std::shared_ptr<std::string>> patterns = {std::make_shared<std::string>(pattern)};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return true;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

    private: 
        std::string pattern;
        std::string revCompPattern;
        EdlibAlignConfig config;
};

class VariableBarcode final : public Barcode 
{

    public:
    VariableBarcode(std::vector<std::string>& inPatterns, std::string& name, int inMismatches, bool hammingIn = false) : Barcode(name, inMismatches)
    {
        for(std::string pattern : inPatterns)
        {
            patterns.push_back(std::make_shared<std::string>(pattern));
        }

        config = edlibNewAlignConfig(
            inMismatches,        // Maximum allowed edit distance
            EDLIB_MODE_SHW,     // Semi-global alignment (deletions in the target at the end are not penalized), 
                                // the developers refer to it also as 'Prefix method'
            EDLIB_TASK_PATH,    // Request full alignment path (M, I, D, S)
            NULL, 0             // No custom alphabet
        );

        //instant look-up table for perfect barcodes
        equalLengthBarcodes = true;
        length = 0; //set length to 0, then set it for the first found abrcode and compare to others
        for (const std::shared_ptr<std::string>& pattern : patterns) 
        {
            if(length == 0){length = pattern->size();}
            else
            {
                if (pattern->size() != length) 
                {
                    equalLengthBarcodes = false;
                    length = 0; //re-set to zero (length is unknown)
                    break;
                }
            }
        }
        if(equalLengthBarcodes)
        {
            for (const std::shared_ptr<std::string>& ptr : patterns)
            {
                barcodeSet.insert(*ptr);  // Dereference the shared_ptr
            }
        }
        
        //calculate minimum conversion rates of barcodes
        //for now this is only calculated if the number of barcodes is less than 1000, 
        // to avoid billions of comaprisons for 10X data
        hamming = hammingIn;
        calculate_hamming_map();
        if(!hamming)
        {
            //calculate_baseNum_hash();
            //calculate_kmer_hash();

            //barcodes need to be same length: otherwise different minimum qgram overlaps are required
            //additionally call use_qgram_map to check if qgram map could be build and makes sense 
            // //(e.g., for many mismatches we might expect few qgrams to overlap and make it useless)
            if(equalLengthBarcodes)
            {
                std::cout << "\t\tCreate qgram-index for fast barcode filtering\n";
                qgramMap.init(mismatches, patterns);
            }
            if(!qgramMap.useQgramMap)
            {
                std::cout << "\t\tSkipping barcode-space reduction by qgram-mapping for barcode " << name <<  "\n";
            }

            //if we do not have too many patterns calculate the conversion rates
            //we can then stop mapping early if we find a barcode that is closer to the sequence than the conversion of this barcode
            //to its closest barcode-neighbor
            if(patterns.size() < 100000)
            {
                calculateConversionRate = true;
                calculate_barcode_conversionRates();
            }
        }
    }
    std::shared_ptr<Barcode> clone() const override 
    {
        return std::make_shared<VariableBarcode>(*this);
    }

    inline void calculate_hamming_map()
    {
        std::cout << "    pre-calculate 1-hamming-distace barcodes\n";
        const std::string bases = "ATGC";
        //go through all forward patterns
        for(std::shared_ptr<std::string> fwPattern : patterns)
        {
            for (size_t i = 0; i < fwPattern->size(); ++i) {
                for (char b : bases) {
                    if (b != fwPattern->at(i)) 
                    {
                        std::string mutated = *fwPattern;
                        mutated.at(i) = b;
                        //if the barcode with mismatch already exists, and we have the same for the current barcode we need to take it out,
                        //apparently we can reach it from two different barcodes and can therefore not truly know where this barcode comes from
                        auto [it, inserted] = hamming_map.emplace(mutated, fwPattern);
                        if (!inserted) 
                        {
                            it->second = nullptr;
                        }
                    }
                }
            }
        }
    }

    inline void calculate_barcode_conversionRates()
    {
        EdlibAlignConfig barcodeConversionConfig = edlibNewAlignConfig(
            -1,                 // no limit for edit distancve
            EDLIB_MODE_SHW,     // Semi-global alignment (deletions in the target at the end are not penalized), 
                                // the developers refer to it also as 'Prefix method', e.g. one barcode is ATC and another ATCATC
                                //this distance should still be zero...
            EDLIB_TASK_PATH,    // Request full alignment path (M, I, D, S)
            NULL, 0);             // No custom alphabet

        for (const std::shared_ptr<std::string>& patternPtrA : patterns) 
        {
            const std::string a = *(patternPtrA);
            int min_rate = std::numeric_limits<int>::max();
            std::string minElement = "";

            for (const std::shared_ptr<std::string>& patternPtrB : patterns) 
            {
                const std::string b = *(patternPtrB);

                if(a == b){continue;}

                int del = 0;
                int ins = 0;
                int subst = 0;
                int targetEnd;

                //for stagger barcode the first barcode has to be pattern, the secone the target, bcs. we do not cound deletions on the target
                // e.g.: barcode A and AGT should have a conversion rate of 0! So for staggered barcodes we need to test ALL barcodes and then take
                //the best match with the longest matching sequence...
                if (a.length() > b.length()) 
                {
                    run_alignment(b, a, targetEnd, barcodeConversionConfig, del, ins, subst);
                } else 
                {
                    run_alignment(a, b, targetEnd, barcodeConversionConfig, del, ins, subst);
                }

                int rate = del + ins + subst;
                if (rate < min_rate) 
                {
                    min_rate = rate;
                    minElement = b;
                }
            }

            //warning if barcodes are the same/ or one is suffix of the other
            if(min_rate == 0)
            {
                std::cout << "WARNING: The data contains barcodes with a mismatch distance of 0!!! Barcodes: " << a << ", "<< minElement << "\n" <<
                "In case these are barcodes of variable length, we map the longest barcode mapping with no errors!!\n";
            }

            //warning if minimum conversion is lower than allowed mismatches
            if(min_rate <= mismatches)
            {
                std::cout << "WARNING: \nFor barcodes in " << name << ": " << mismatches << " mismatches are allowed, but with " << min_rate <<
                " mismatches we can already convert " << a  << " into "<< minElement  <<  " (semi-global alignment with un-punished deletions in target)!!! We can still find a best-fitting barcode, but you might reconsider the choice of allowed mistmaches!\n";
            }

            pattern_conversionrates[a] = min_rate;

        }
    }

    inline bool align(std::string& matchedBarcode, const std::string& fastqLine, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)
    {

        //check if we can instantly match pattern
        if(equalLengthBarcodes && (fastqLine.size() >= (targetOffset + length)))
        {

            //1.) check for simple perfect matches of a barcode
            std::string exactTarget = fastqLine.substr(targetOffset, length);

            if(reverse){exactTarget = generate_reverse_complement(exactTarget);}
            if (barcodeSet.find(exactTarget) != barcodeSet.end()) 
            {
                targetEnd = length;
                matchedBarcode = exactTarget;

                return true;
            } 
            else if(mismatches == 0)
            {
                return false;
            }

            //2.) check match with 1hamming distance (precomputed hash-map of 1 MM barcodes)
            if(hamming_map.find(exactTarget) != hamming_map.end())
            {
                if(hamming_map.at(exactTarget) == nullptr){return false;} //if this observed barcode with MM can be reached by several barcodes
                targetEnd = length;
                matchedBarcode = *(hamming_map.at(exactTarget));
                substNum = 1;

                return true;
            }
            else if(hamming){return false;} //if we only do hamming distance mapping, stop here if we found nothing
        }

        std::vector<std::shared_ptr<std::string>> patternsToMap = patterns;

        std::string bestFoundPattern;
        bool bestFoundAlignment = false;
        int bestTargetEnd = -1;
        int bestEditDist = mismatches+1;
        bool severalMatches = false; //if there are several best-fitting solutions (only the case when the number of allowed mismatches
        //is bigger than possible barcode-conversion numbers) we discard the solution

//REPLACCE THIS WITH Q-GRAMS
        //3.) reduce possible patterns but comparing base number/ kmers
        //get the basenum hashes: based on counts of bases (considering allowed MM) - how many barcodes could fit
        
        /*
        if(equalLengthBarcodes)
        {

            std::string targetSequence = fastqLine.substr(targetOffset, length);

            //we count bases in target
            baseNum baseCountTmp;
            count_bases(targetSequence, baseCountTmp);

            //the baseNum maps map all possible baseNum keys (number for {A,C,G,T}) to the barcodes that could be described by this
            //including barcodes with MM
            std::vector<std::shared_ptr<std::string>> possiblePatterns;
            if(reverse)
            {
                auto it = rvEditBaseNumMap.find(baseCountTmp);
                if (it != rvEditBaseNumMap.end()) 
                {
                    possiblePatterns = it->second;
                } else 
                {
                    possiblePatterns.clear(); //no possible patterns found
                }
            }
            else
            {
                auto it = fwEditBaseNumMap.find(baseCountTmp);
                if (it != fwEditBaseNumMap.end()) 
                {
                    possiblePatterns = it->second;
                } else 
                {
                    possiblePatterns.clear(); //no possible patterns found
                }
            }
            //us epossible pattern only if there is NO N in the seqeunce and the lenght of the target is equal to the ones of barcodes
            //(e.g., at the end of the fastq this could not be the case) 
            if(targetSequence.find('N') == std::string::npos && !(targetSequence.size()< length) )
            {
                patternsToMap = possiblePatterns;
            }

            //reduce possible barcodes by kmer filter
            // every string is encoded in a kmer-array (e.g. [1001001110110201]), then calcualte similariy of these arrays 
            // looking for kmer-array within a minimum distance, e.g., for 2MM SET_KMERS- 2(every EDIT changes 2 kmers max)*2MM kmers must be same at least
            std::vector<std::shared_ptr<std::string>> prunedPatternsToMap; 

            prunedPatternsToMap = kmer_align_patterns(patternsToMap, targetSequence, reverse);
            patternsToMap = prunedPatternsToMap;

        }*/

        //NEW Q-GRAM APPROACH
        if(equalLengthBarcodes && qgramMap.useQgramMap)
        {
            std::string targetSequence = fastqLine.substr(targetOffset, length);
            if(reverse){targetSequence = generate_reverse_complement(targetSequence);}

            patternsToMap = qgramMap.query(targetSequence);

        }

        //do not compute super heavy alignments: 10x data showed that most barcodes have at max 1000 possibilities
        if(patternsToMap.size() > 2000)
        {
            clear_progress_bar();
            std::cerr << "Skipping fastq-line with more than 500 possibilities after qgram-mapping: " << patternsToMap.size() << std::endl;
            return false;
        }
        for(const std::shared_ptr<std::string>& patternPtr : patternsToMap)
        {

            bool foundAlignment = false;

            int delNumTmp;
            int insNumTmp;
            int substNumTmp;
            delNumTmp=insNumTmp=substNumTmp=targetEnd=0;
            //get pattern, its length can vary
            std::string usedPattern = *(patternPtr);
            if(reverse){usedPattern = generate_reverse_complement(usedPattern);}

            //define target sequence (can differ for every barcode due to its length)
            std::string target;
            int substringLength = usedPattern.length()+mismatches;
            //std::cout << "\t LENGTH: " <<substringLength << " seq: " << fastqLine.size()<< "\n";
            if(targetOffset + substringLength > fastqLine.size()){substringLength = fastqLine.size()-targetOffset;};
            target = fastqLine.substr(targetOffset, substringLength);
            
            //map the pattern to the target sequence
            foundAlignment = run_alignment(usedPattern, target, targetEnd, config, delNumTmp,  insNumTmp, substNumTmp);
            
            if(foundAlignment && (delNumTmp+insNumTmp+substNumTmp)<=bestEditDist)
            {

                //if we have variable length barcodes, and two barcodes map equally well, 
                //we store the barcode of a longer sequence
                //example: barcodes = ["ATC", "ATCATC"], if we find the barcode ATCATC we store that one and not ATC...
                //this if clause checks if the new match is shorter than odl one, if so do not store it
                if(!equalLengthBarcodes && targetEnd < bestTargetEnd)
                {
                    continue;
                }

                //check if we have several best solutions
                if((delNumTmp+insNumTmp+substNumTmp)==bestEditDist)
                {
                    severalMatches = true;
                }
                else
                {
                    severalMatches = false;
                }

                bestFoundAlignment = foundAlignment;
                //bestFoundPattern = patternsToStore.at(patternIdx); //the barcode is the TRUE forward barcode, not the reverse complement
                if(reverse)
                {bestFoundPattern = generate_reverse_complement(usedPattern);}
                else {bestFoundPattern =usedPattern;}
                bestTargetEnd = targetEnd;
                bestEditDist = (delNumTmp+insNumTmp+substNumTmp);
                delNum = delNumTmp;
                insNum = insNumTmp;
                substNum = substNumTmp;

                //if we found a new best match, check if this is already the best match we can ever get (minimal conversion dist between barcodes)
                if(calculateConversionRate)
                {
                    int minConversion = pattern_conversionrates.at(bestFoundPattern); //patternsToStore.at(patternIdx)

                    //we can do this since levenshtein distance fullfills the triangle inequality is a distance metric
                    //imagine there is a second barcode that could fit better: this second barcode must have a shorter distance to target sequence
                    //than our pattern. Now there r two options 1.) while converting pattern to target we would 'go through' the second barcode. In this
                    //case minConversion/2 is always bigger than the editDist and we don t break
                    //2.) the minconversion is the maximum conversion from pattern to second barcode and target is inbetween converting these two
                    //now if the current barcode is however closer to target (minConversion/2), then this is the best match we cna ever find...
                    if(bestEditDist < (minConversion/2))
                    {
                        break;
                    }
                }
            }

        }

        //if we have several best matches we have to return false
        if(severalMatches)
        {
            return false;
        }

        targetEnd = bestTargetEnd;
        matchedBarcode = bestFoundPattern;

        return bestFoundAlignment;
    }

    inline std::vector<std::shared_ptr<std::string>> get_patterns()
    {
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

    private:
        std::vector<std::shared_ptr<std::string>> patterns;
        std::unordered_set<std::string> barcodeSet;
        bool equalLengthBarcodes;
        bool hamming;

        //map that maps a nucleotide sequence with haming distances of 1 to its real barcodePtr
        //if only hamming of 1 is set this runs super fast, otherwise after this barcodes are aligned
        ankerl::unordered_dense::map<std::string, std::shared_ptr<std::string>> hamming_map;

        //inverse map, that maps qgrams to its barcode containing it, e.g., AGAG -> AGAGTA,TAGAGC
        QGramMap qgramMap;

        ankerl::unordered_dense::map<std::string, int> pattern_conversionrates;
        bool calculateConversionRate = false;

        EdlibAlignConfig config;
};

class WildcardBarcode final : public Barcode
{
    // TODO: move public parameter mismatches to a private and derived parameter
    //wildcardBarcode doe snot make use of mismatches yet, since anyways we do not know the sequence,
    //therefore its an unused parameter, just set for completeness as these classes derive from Barcode (initialized with mismatches, see up...)
    public:
    WildcardBarcode(int inMismatches, std::string name, int inLength) : Barcode(name, inMismatches, inLength) {}
    std::shared_ptr<Barcode> clone() const override {
        return std::make_shared<WildcardBarcode>(*this);
    }

    inline bool align(std::string& matchedBarcode, const std::string& target, const unsigned int positionInFastqLine,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)
    {
        (void)delNum; // silence unused parameter warning
        (void)insNum; // silence unused parameter warning
        (void)substNum; // silence unused parameter warning
        (void)reverse; // silence unused parameter warning

        matchedBarcode = target.substr(positionInFastqLine, length);
        targetEnd = (target.length() < length) ? target.length() : length;
        // e.g.: [AGTAGT]cccc: start=0 end=6 end is first not included idx
        return true;
    }

    inline std::vector<std::shared_ptr<std::string>> get_patterns()
    {
        std::vector<std::shared_ptr<std::string>> patterns = {std::make_shared<std::string>(std::string(length, 'X'))};
        return patterns;
    }
    bool is_wildcard(){return true;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

};

//A stop-barcode: mapping on both sides is only done up to here
//the name ('*') is not writen into the output file
class StopBarcode final : public Barcode
{
    public:
    StopBarcode(std::string inPattern, int inMismatches) : Barcode("*", inMismatches),pattern(inPattern) {}
    std::shared_ptr<Barcode> clone() const override {
        return std::make_shared<StopBarcode>(*this);
    }
    inline bool align(std::string& matchedBarcode, const std::string& target, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            (void)matchedBarcode; // silence unused parameter warning
            (void)target; // silence unused parameter warning
            (void)targetOffset; // silence unused parameter warning
            (void)targetEnd; // silence unused parameter warning
            (void)delNum; // silence unused parameter warning
            (void)insNum; // silence unused parameter warning
            (void)substNum; // silence unused parameter warning
            (void)reverse; // silence unused parameter warning

            return false;
        }
    inline std::vector<std::shared_ptr<std::string>> get_patterns()
    {
        std::vector<std::shared_ptr<std::string>> patterns = {std::make_shared<std::string>(pattern)};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return true;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

    private:
    std::string pattern; //just a string of "XXXXX"
};

//A stop-barcode: mapping on both sides is only done up to here
//the name ('*') is not writen into the output file
class ReadSeperatorBarcode final : public Barcode
{
    public:
    ReadSeperatorBarcode(std::string inPattern, int inMismatches) : Barcode("-", inMismatches),pattern(inPattern) {}
    std::shared_ptr<Barcode> clone() const override {
        return std::make_shared<ReadSeperatorBarcode>(*this);
    }
    inline bool align(std::string& matchedBarcode, const std::string& target, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            (void)matchedBarcode; // silence unused parameter warning
            (void)target; // silence unused parameter warning
            (void)targetOffset; // silence unused parameter warning
            (void)targetEnd; // silence unused parameter warning
            (void)delNum; // silence unused parameter warning
            (void)insNum; // silence unused parameter warning
            (void)substNum; // silence unused parameter warning
            (void)reverse; // silence unused parameter warning

            return false;
        }
    inline std::vector<std::shared_ptr<std::string>> get_patterns()
    {
        std::vector<std::shared_ptr<std::string>> patterns = {std::make_shared<std::string>(pattern)};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return true;}

    private:
    std::string pattern; //just a string of "XXXXX"
};

//barcode representing a genomic region. We can set the number of mismatches or when set to -1 simply run it with default
//star settings
//pattern is simply "DNA"
class DNABarcode final : public Barcode
{
    public:
    DNABarcode(int inMismatches = -1) : Barcode("DNA", inMismatches) {pattern = "DNA";}
    std::shared_ptr<Barcode> clone() const override {
        return std::make_shared<DNABarcode>(*this);
    }
    inline bool align(std::string& matchedBarcode, const std::string& target, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            (void)matchedBarcode; // silence unused parameter warning
            (void)target; // silence unused parameter warning
            (void)targetOffset; // silence unused parameter warning
            (void)targetEnd; // silence unused parameter warning
            (void)delNum; // silence unused parameter warning
            (void)insNum; // silence unused parameter warning
            (void)substNum; // silence unused parameter warning
            (void)reverse; // silence unused parameter warning

            return false;
        }
    inline std::vector<std::shared_ptr<std::string>> get_patterns()
    {
        std::vector<std::shared_ptr<std::string>> patterns = {std::make_shared<std::string>(pattern)};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return true;}
    bool is_read_end(){return false;}

    private:
    std::string pattern; //just a string of "DNA"
};

//class to handle a barcode pattern 
//can be used to iterate through the barcode, and stores additional information like:
//it stores if the pattern contains DNA barcodes which require different handling
class BarcodePattern
{
    public:
        //constructor
        BarcodePattern(bool dna, std::string name, BarcodeVectorPtr pattern) :  
        containsDNA(dna), patternName(name), barcodePattern(pattern) {}

        // Deep copy constructor for barcode pattern
        BarcodePattern(const BarcodePattern& other)
            : containsDNA(other.containsDNA),
            patternName(other.patternName),
            barcodePattern(copy_barcode_vector(other.barcodePattern)),
            independentReversePattern(copy_barcode_vector(other.independentReversePattern)) {}

        //class variables
        bool containsDNA;
        std::string patternName; //this is also the file this pattern will be written to
        BarcodeVectorPtr barcodePattern;
        BarcodeVectorPtr independentReversePattern; //in case we do not have one long pattern with a fw&rv read
        //but more two independent read that should be mapped seperately

        //class functions
        //write multiplexed lines to file (must be specific for DNA, AB-barcodes, etc.)
        //void write_demultiplexed_line(const std::vector<std::string> barcodeList, std::string dna = "");

        //barcodeVector/ iterator functions
        // Add a barcode to the barcodePattern
        inline void add_barcode(const BarcodePtr& barcode, PatternType type = PatternType::Forward) 
        {
            get_pattern(type)->push_back(barcode);
        }
        // Get the size of the barcodePattern
        inline std::size_t size(PatternType type = PatternType::Forward) const 
        {
            return get_pattern(type)->size();
        }

        // Iterator types
        using iterator = typename std::vector<BarcodePtr>::iterator;
        using const_iterator = typename std::vector<BarcodePtr>::const_iterator;
        using reverse_iterator = typename std::vector<BarcodePtr>::reverse_iterator;
        using const_reverse_iterator = typename std::vector<BarcodePtr>::const_reverse_iterator;
        // Begin and end iterators
        iterator begin(PatternType type = PatternType::Forward) {
            return get_pattern(type)->begin();
        }
        const_iterator begin(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->begin();
        }
        iterator end(PatternType type = PatternType::Forward) {
            return get_pattern(type)->end();
        }
        const_iterator end(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->end();
        }
        // Reverse iterators
        reverse_iterator rbegin(PatternType type = PatternType::Forward) {
            return get_pattern(type)->rbegin();
        }
        const_reverse_iterator rbegin(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->rbegin();
        }
        reverse_iterator rend(PatternType type = PatternType::Forward) {
            return get_pattern(type)->rend();
        }
        const_reverse_iterator rend(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->rend();
        }

    private:
        BarcodeVectorPtr get_pattern(PatternType type) const 
        {
            if(type == PatternType::Reverse){return independentReversePattern;}
            return barcodePattern; //if not reverse return forward pattern
        }

        //copy barcodeVectors of the barcodePattern
        static BarcodeVectorPtr copy_barcode_vector(const BarcodeVectorPtr& original) 
        {
            if (!original) return nullptr;
            std::shared_ptr<BarcodeVector> copy = std::make_shared<BarcodeVector>();
            copy->reserve(original->size());
            for (const BarcodePtr& barcode : *original) 
            {
                if (barcode) 
                {
                    // clone() returns a std::shared_ptr<Barcode>
                    std::shared_ptr<Barcode> clonedBarcode = barcode->clone();
                    copy->push_back(clonedBarcode);
                } else 
                {
                    // still push a nullptr if the original was nullptr
                    copy->push_back(nullptr);
                }
            }
            return copy;
        }
};

typedef std::shared_ptr<BarcodePattern> BarcodePatternPtr; 
typedef std::shared_ptr<std::vector<BarcodePatternPtr>> MultipleBarcodePatternVectorPtr; 