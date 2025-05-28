#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cassert>
#include <string_view>
#include <cstring>
#include <atomic>
#include <mutex>
#include <algorithm>

#include "edlib/edlib/include/edlib.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#define MAX(X,Y) (X<Y ? Y : X)
#define MIN(X,Y) (X<Y ? X : Y)

inline bool endWith(std::string const &fullString, std::string const &ending) 
{
    if (fullString.length() >= ending.length()) 
    {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } 
    else 
    {
        return false;
    }
}

inline unsigned long long totalNumberOfLines(std::string fileName)
{
    unsigned long long totalReads = 0;
    unsigned char buffer[1000];
    gzFile fp = gzopen(fileName.c_str(),"r");
    if(NULL == fp){
        fprintf(stderr,"Fail to open file: %s\n", fileName.c_str());
    }
    while(!gzeof(fp))
    {
        gzread(fp, buffer, 999);
        for (const unsigned char &c : buffer) 
        {
            if ( c == '\n' )
            {
                if(totalReads == ULLONG_MAX)
                {
                    std::cout << "WARNING: Analysing more than " << std::to_string(ULLONG_MAX) << " reads. There will be no status update\n";
                    gzrewind(fp);
                    return(ULLONG_MAX);
                }
                ++totalReads;
            }
        }
    }
    gzrewind(fp);

    return totalReads;
}

inline unsigned long long numberOfReads(std::string fileName)
{
    unsigned long long totalReads = 0;
    
    if(endWith(fileName, "fastq") || endWith(fileName, "fastq.gz"))
    {
        unsigned long long totalLines = totalNumberOfLines(fileName);
        if(totalLines==ULLONG_MAX){return(ULLONG_MAX);}
        totalReads = (totalLines/4);
    }
    else if(endWith(fileName, "txt"))
    {
        std::ifstream fileStream;
        fileStream.open(fileName);
        totalReads = std::count(std::istreambuf_iterator<char>(fileStream), std::istreambuf_iterator<char>(), '\n');
        fileStream.clear();
        fileStream.seekg(0);
        fileStream.close();
    }

    return totalReads;
}

inline void printProgress(double percentage) 
{
    int val = (int) (percentage*100);
    int loadLength = (int) (percentage * PBWIDTH);
    int emptyLength = PBWIDTH - loadLength;
    std::cout << "\t\r[" << std::string(loadLength, '|') << std::string(emptyLength, ' ') << "] " << val << "%" << std::flush;
}

//stores all the input parameters for the mapping tools
struct input{
    std::string inFile;

    //this is the path to the output files
    std::string outPath;
    std::string prefix;

    std::string reverseFile = "";

    std::string barcodeFile; //file of all barcode-vectors, each line sequentially representing a barcode 
    std::string mismatchFile; //file withg several lines with coma seperated list of mismathces per barcode
    std::string patternLine; //list of patterns in abstract form

    std::string barcodePatternsFile;

    //additional informations
    bool writeStats = false; 
    bool writeFailedLines = false;
    bool writeFilesOnTheFly = false;
    
    long long int fastqReadBucketSize = 10000000;
    int threads = 5;
};

struct levenshtein_value{
        unsigned int val = 0;
        int i = 0;
        int j = 0;
        levenshtein_value(int a,int b, int c):val(a), i(b), j(c){}
        levenshtein_value():val(UINT_MAX-1), i(0), j(0){} // -1 to make val+1 not start from 0
};
inline levenshtein_value min(levenshtein_value a, levenshtein_value b)
{
    if(a.val <= b.val)
    {
        return a;
    }
    return b;
}

struct frontMatrix
{
    std::vector<unsigned int> previous;
    std::vector<unsigned int> current;
    unsigned int offset;
    unsigned int d;

    frontMatrix(unsigned int m, unsigned int n): previous(m+n+3, UINT_MAX), current(m+n+3, UINT_MAX){}
};
//calculate longest common prefix of two sequences
inline int lcp(const std::string& a, const std::string& b)
{
    unsigned int lcp = 0;
    unsigned int a_len = a.length();
    unsigned int b_len = b.length();
    unsigned int len = (a_len > b_len) ? b_len : a_len;

    while( (a[lcp] == b[lcp]) && (lcp < len) )
    {
        ++lcp;
    }
    return lcp;
}
//calculates the next front
//idea: how far along all the diagonals that r within x-mismatches can i go in my edit-matrix with x-mismatches
inline void front(const std::string& a, const std::string& b, frontMatrix& f)
{
    int min = std::min(static_cast<int>(a.length()), static_cast<int>(f.d));
    int max = std::min(static_cast<int>(b.length()), static_cast<int>(f.d));

    f.previous = f.current;

    unsigned int l = 0;
    for(unsigned int i = (f.offset-min); i <= (f.offset + max); ++i)
    {
        unsigned int aVal = 0, bVal = 0, cVal = 0;
        if(f.previous.at(i-1) != UINT_MAX){ aVal = f.previous.at(i-1);}
        if(f.previous.at(i+1) != UINT_MAX){ bVal = f.previous.at(i+1)+1;}
        if(f.previous.at(i) != UINT_MAX){ cVal = f.previous.at(i)+1;}

        l = std::max(std::max(aVal, bVal), cVal);
            
            if(l >= a.length())
            {
                f.current.at(i) = a.length();
            }
            else if( (i - f.offset + l) >= b.length())
            {
                //in this case the row must be i, however it already is at this point...
                //to set explicitely uncomment
                //f.current.at(i) = i;
            }
            else
            {
                f.current.at(i) = l + lcp( a.substr(l), b.substr(i - f.offset + l) );
            }    
    }

}

//output sensitive -> O() depends on the mismatches that we allow, since we allow mostly just for a few it runs super fast...
//no backtracking implemented for now, only used to align UMIs where we do not care about alingment start, end
inline bool outputSense(const std::string& sequence, const std::string& pattern, const int& mismatches, int& score)
{
    const int m = sequence.length();
    const int n = pattern.length();
    
    frontMatrix f(m,n);
    f.offset = m + 1;
    f.d = 0;

    f.current.at(f.offset) = lcp(sequence, pattern);
    if(f.current.at(n-m+f.offset) == m)
    {
        score = f.d;
        return true;
        if(score <= mismatches){return true;}
        else
        {
            score = mismatches + 1;
            return false;
        }
    }

    ++f.d;
    while(f.d <= std::min(std::max(m,n), mismatches))
    {
        front(sequence, pattern, f);
        if(f.current.at(n-m+f.offset) == m )
        {
            score = f.d;
            if(score <= mismatches)
            {
                return true;
            }
            else
            {
                score = mismatches + 1;
                return false;
            }
        }
        ++f.d;
    }

    //no match found within limits: assign a score > mismatches
    score = mismatches + 1;
    return false;
}

inline void free_levenshtein(levenshtein_value** dist, int ls)
{
    for (int i = 0; i <= ls; ++i) 
    {
        delete[] dist[i];
    }
    delete[] dist;
}

inline std::string stripQuotes(const std::string& input) 
{
    if (input.size() >= 2 && input.front() == '\'' && input.back() == '\'') {
        return input.substr(1, input.size() - 2); // cut first and last character
    } else {
        return input;
    }
}

//running alingment with the myers bit-parallel algorithm implemented in the edlib
//https://github.com/Martinsos/edlib published here: Bioinformatics 2017 btw753. doi: 10.1093/bioinformatics/btw753

//we run the alignment in following mode:
    //we map the pattern to a longer target (pattern + maxEditDist) in case we have deletions at the beginning
    //and then map the pattern perfectly.
    //we align with a semi-global strategy, allowing deletions at the end, in case we map the pattern perfectly at the 
    //beginning

    //TODO: add an alternative strategy without saving ins, del, subst for cases where the quality is not needed
    //and running time is more important
inline bool run_alignment(std::string pattern, const std::string target, 
                          int& targetEnd,
                          EdlibAlignConfig config,
                          int& delNum, int& insNum, int& substNum)
{
    // bool to check if we found a result
    bool alignmentFound = false;

    // Run the alignment
    EdlibAlignResult result = edlibAlign(
        pattern.c_str(), pattern.length(),
        target.c_str(), target.length(),
        config
    );

    std::vector<char> alignmentPath;

    if (result.status == EDLIB_STATUS_OK && result.editDistance != -1) {
        // Translate Edlib codes (0 = M, 1 = I, 2 = D, 3 = X) into characters
        for (int i = 0; i < result.alignmentLength; ++i) 
        {
            switch (result.alignment[i]) {
                case 0: alignmentPath.push_back('M'); break; // Match
                case 1: alignmentPath.push_back('D'); break; // Insertion to target = deletion in pattern
                case 2: alignmentPath.push_back('I'); break; // Deletion in target = insertion in pattern
                case 3: alignmentPath.push_back('S'); break; // Mismatch (substitution)
            }
        }

        alignmentFound = true;

        // replace edits at end with SUBS, we prefer SUBS instead of DEL,INS at the end
        //edlib prefers DEL (or INS in target) at the end (e.g., ACGT maps to ACGAT with M-M-M-D instead of a S or I-M at the end)
        //to cut the target seqeunce however, we need to also consider the last insertion and cut [ACGA]-T and continue from T for
        //the next pattern alignment
        // example:
        // pattern:ACGT target:CACGATGTAA is I M M M D

        // COUNT EDITS (and replace last edits with S)
        bool replacinglastEdits = true; //this variables stores if we still replace last edits in alignment with substitutions
        //as soon as we encouter a first match from back, we set it to false and count the true edits
        for (int i = alignmentPath.size() - 1; i >= 0; --i) 
        {
            //stop replacing edits at end
            if (replacinglastEdits && alignmentPath[i] == 'M') 
            {
                replacinglastEdits = false; // stop if we find 'M'
            }

            //repalce edits at end with substitutions (if the last edits are not a subst assume they are)
            //e.g., the EDLIB finds deletions in the end, to find the shortest alignment, however, we do not know if those
            //last edits were deletions/ subst, just store them as subst
            if(replacinglastEdits){++substNum; ++targetEnd;}
            //count edit types
            else if (!replacinglastEdits && alignmentPath[i] == 'S') {++substNum; ++targetEnd;}
            else if (!replacinglastEdits && alignmentPath[i] == 'I') {++insNum; ++targetEnd;}
            //we do not count deletions in target end!!! As these are edits that are not present in the read
            else if (!replacinglastEdits && alignmentPath[i] == 'D') {++delNum;}
            else {++targetEnd;}
        }
    } 
    else 
    {
        alignmentFound = false;
    }

    edlibFreeAlignResult(result);

   // std::cout << "operaitons: ";
   // for(auto el : alignmentPath)
   // {
   //     std::cout << el << " ";
   // }
   // std::cout << "Target end in alignment funct: " << targetEnd << "\n";

    return alignmentFound;
}

/*
//levenshtein distance, implemented with backtracking to get start and end of alingment, however slower than output sensitive algorithm:
//used for parser so far: it has an additional flavor of unpunished deletions at the start and end of the alignment
//start is 0 indexed, end are the first indices that arre not part of the match
inline bool levenshtein(const std::string sequence, std::string pattern, const int& mismatches, int& match_start, int& match_end, int& score,
                        int& endInPattern, int& startInPattern, bool upperBoundCheck = false)
{
    int i,j,ls,la,substitutionValue, deletionValue;
    //stores the lenght of strings s1 and s2
    ls = sequence.length();
    la = pattern.length();
    //levenshtein_value dist[ls+1][la+1];
    //dynamic allcation of edit-matrix
    levenshtein_value** dist = new levenshtein_value*[ls + 1];
    for (int i = 0; i <= ls; ++i) 
    {
        dist[i] = new levenshtein_value[la + 1];
    }

    //allow unlimited deletions in beginning; start is zero
    for(i=0;i<=ls;i++) {
        levenshtein_value val(0,i-1,0);
        dist[i][0] = val;
    }
    //deletions in pattern is punished; start is zero
    for(j=0;j<=la;j++) {
        levenshtein_value val(j,0,j-1);
        dist[0][j] = val;
    }

    levenshtein_value val(0,-1,-1);
    dist[0][0] = val;

    int upperBoundCol = mismatches;
    for (i=1;i<=ls;i++) 
    {
        for(j=1;j<=la;j++) 
        {

            //Punishement for substitution
            if(sequence[i-1] == pattern[j-1]) 
            {
                substitutionValue= 0;
            } 
            else 
            {
                substitutionValue = 1;
            }
            //punishement for deletion (allow deletions on beginnign and end)
            if(j==la )
            {
                deletionValue= 0;
            } 
            else 
            {
                deletionValue = 1;
            }

            //determine new value of cell
            levenshtein_value seq_del = dist[i-1][j]; 
            seq_del.val += deletionValue;
            seq_del.i = i-1;
            seq_del.j = j;
            levenshtein_value seq_ins = dist[i][j-1];
            seq_ins.val += 1;
            seq_ins.i = i;
            seq_ins.j = j-1;
            levenshtein_value subst = dist[i-1][j-1];
            subst.val += substitutionValue;
            subst.i = i-1;
            subst.j = j-1;

            //if equal score, prefer: del  > subst > ins (deletions are prefereed for end, at the end we want only deletions intead of first deletions and then
            //subst which might go into the next barcode)
            levenshtein_value tmp1 = min(subst, seq_ins);
            levenshtein_value tmp2 = min(seq_del, tmp1);

            if(upperBoundCheck && (tmp2.val > mismatches) && (j > upperBoundCol) )
            {
                    upperBoundCol = j-1;
                    break;
            }
            
            dist[i][j] = tmp2;
            //uncomment to show the edit-matrix
            //std::cout << tmp2.val << "(" << sequence[i-1] <<  ","<< pattern[j-1] << ")" << " ";
        }
        //std::cout << "\n";
    }
    if(upperBoundCheck && (dist[ls][la].val > mismatches))
    {
        free_levenshtein(dist, ls);
        return false;
    }
    //backtracking to find match start and end
    // start and end are defined as first and last match of bases 
    //(deletion, insertion and substitution are not considered, since they could also be part of the adjacent sequences)
    int start;
    int end;
    i = ls;
    j =la;
    bool noEnd = true;
    while(j!=0 && i!=0) //only go until i==2, in this care the new value would be one. For the final result we must substract one to get to zero
    {
        int iNew = dist[i][j].i;
        int jNew = dist[i][j].j;

        if(dist[i][j].val == dist[iNew][jNew].val && i!=iNew && j!=jNew)
        {
            start = i;
            startInPattern = j;
        }

        if( (jNew<la) && noEnd && dist[i][j].val == dist[iNew][jNew].val)
        {
            noEnd = false;
            end = i;
            endInPattern = j;
        }

        i = iNew;
        j = jNew;
    }
    //in case number of mismatches is less or equal to threshold
    if((dist[ls][la]).val <= mismatches)
    {
        // asser that those values can never get zero
        assert(start!=0);
        assert(startInPattern!=0);
        // end is the first index out of sequence, since our strings in matrix are 1-indexed end is fine
        // start has to be start -= 1, to be the start in a 0-indexed string
        score = (dist[ls][la]).val;
        match_start = start-1;
        startInPattern -= 1;
        match_end = end;
        free_levenshtein(dist, ls);
        return true;
    }
    else
    {
        score = (dist[ls][la]).val;
    }

    free_levenshtein(dist, ls);
    return false;
}
*/

inline int backBarcodeMappingExtension(const std::string& sequence, const std::string& pattern, int seq_end, const int& patternEnd)
{
    int elongation = 0;
    //check if the end of sequebnces still maps for deletions
    for(int i =0; i < (pattern.length() - patternEnd); ++i)
    {
        if( (seq_end + i) >= sequence.length()){return elongation;} //special case, the sequence ends before finishing of the whole pattern
        if( sequence.at(seq_end+i) == pattern.at(patternEnd + i) )
        {
            ++elongation;
        }
        else
        {
            return elongation;
        }
    }
    return elongation;
}

inline int frontBarcodeMappingExtension(const std::string& sequence, const std::string& pattern, const int& seq_start, const int& patternStart)
{
    int elongation = 0;
    //check if the end of sequebnces still maps for deletions
    //start in pattern is 1-based: patternstart of 2 euquals one base deletion at start
    for(int i =0; i < (patternStart); ++i)
    {
        if( sequence.at(patternStart - 1 - i) == pattern.at(patternStart - 1 - i) )
        {
            ++elongation;
        }
        else
        {
            return elongation;
        }
    }
    return elongation;
}

inline void ullong_save_add(unsigned long long& a, const unsigned long long& b)
{
   if( (b>0) && ( (ULLONG_MAX - b) < a ) )
   {
       std::cerr << "ERROR ULONG OVERFLOW\n";
       a = ULLONG_MAX;
   }
   else
   {
       a += b;
   }
}

inline std::vector<std::string> splitByDelimiter(std::string line, const std::string& del)
{
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    while ((pos = line.find(del)) != std::string::npos) {
        token = line.substr(0, pos);
        tokens.push_back(token);
        line.erase(0, pos + del.length());
    }
    tokens.push_back(line);

    return tokens;
}
