#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cassert>
#include <string_view>
#include <cstring>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#define MAX(X,Y) (X<Y ? Y : X)
#define MIN(X,Y) (X<Y ? X : Y)

inline int totalNumberOfLines(std::string fileName)
{
    int totalReads = 0;
    unsigned char buffer[1000];
    gzFile fp = gzopen(fileName.c_str(),"r");
    if(NULL == fp){
        fprintf(stderr,"Fail to open file: %s\n", fileName.c_str());
    }
    while(!gzeof(fp))
    {
        gzread(fp, buffer, 999);
        for (const char &c : buffer) 
        {
            if ( c == '\n' )
            {
                ++totalReads;
            }
        }
    }
    gzrewind(fp);

    return totalReads;
}

inline bool endWith(std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

inline void printProgress(double percentage) {
    int val = (int) (percentage*100);
    int loadLength = (int) (percentage * PBWIDTH);
    int emptyLength = PBWIDTH - loadLength;
    std::cout << "\t\r[" << std::string(loadLength, '|') << std::string(emptyLength, ' ') << "] " << val << "%" << std::flush;
}

//stores all the input parameters for the parser tool
struct input{
    std::string inFile;
    std::string outFile;

    std::string barcodeFile;
    std::string mismatchLine;
    std::string patternLine;

    bool storeRealSequences = false;
    int fastqReadBucketSize = 10000000;
    int threads = 5;
};

struct fastqStats{
    //parameters that are evaluated over the whole fastq line
    //e.g. perfect match occurs only if ALL barcodes match perfectly in a fastq line
    int perfectMatches = 0;
    int noMatches = 0;
    int moderateMatches = 0;
    //parameter stating how often a barcode sequence could be matched to several sequences, can occure more than once per line
    //can only happen for vairable sequences
    int multiBarcodeMatch =0;
    //a dictionary of the number of mismatches in a barcode, in the case of a match
    std::map<std::string, std::vector<int> > mapping_dict;
};

struct levenshtein_value{
        unsigned int val = 0;
        int i = 0;
        int j = 0;
        levenshtein_value(int a,int b, int c):val(a), i(b), j(c){}
        levenshtein_value():val(UINT_MAX-1), i(0), j(0){} // -1 to make val+1 not start from 0
};
inline levenshtein_value min(levenshtein_value a, levenshtein_value b, bool& first)
{
    if(a.val <= b.val)
    {
        first = true;
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
    unsigned int min = MIN(a.length(), f.d);
    unsigned int max = MIN(b.length(), f.d);

    f.previous = f.current;

    unsigned int l = 0;
    for(unsigned int i = (f.offset-min); i <= (f.offset + max); ++i)
    {
        unsigned int aVal = 0, bVal = 0, cVal = 0;
        if(f.previous.at(i-1) != UINT_MAX){ aVal = f.previous.at(i-1);}
        if(f.previous.at(i+1) != UINT_MAX){ bVal = f.previous.at(i+1)+1;}
        if(f.previous.at(i) != UINT_MAX){ cVal = f.previous.at(i)+1;}

        l = MAX(MAX(aVal, bVal), cVal);
            
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
    const unsigned int m = sequence.length();
    const unsigned int n = pattern.length();
    
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
    while(f.d <= MIN(MAX(m,n), mismatches))
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

//levenshtein distance, implemented with backtracking to get start and end of alingment, however slower than output sensitive algorithm:
//used for parser so far: it has an additional flavor of unpunished deletions at the start and end of the alignment
inline bool levenshtein(const std::string sequence, std::string pattern, const int& mismatches, int& match_start, int& match_end, int& score,
                        bool upperBoundCheck = false)
{
    int i,j,ls,la,substitutionValue, deletionValue;
    //stores the lenght of strings s1 and s2
    ls = sequence.length();
    la = pattern.length();
    levenshtein_value dist[ls+1][la+1];

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

            //if equal score, prefer: subst > ins > del
            bool firstValueIsMin = false;
            levenshtein_value tmp1 = min(seq_ins, seq_del, firstValueIsMin);
            firstValueIsMin = false;
            levenshtein_value tmp2 = min(subst, tmp1, firstValueIsMin);

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
    if(upperBoundCheck && (dist[ls][la].val > mismatches)){return false;}
    //backtracking to find match start and end
    // start and end are defined as first and last match of bases 
    //(deletion, insertion and substitution are not considered, since they could also be part of the adjacent sequences)
    int start;
    int end;
    i = ls;
    j =la;
    bool noEnd = true;
    while(j != 1)
    {

        int iNew = dist[i][j].i;
        int jNew = dist[i][j].j;

        if(dist[i][j].val == dist[iNew][jNew].val && i!=iNew && j!=jNew)
        {
            start = iNew;
        }

        if( (jNew<la) && noEnd)
        {
            noEnd = false;
            end = i;
        }

        i = iNew;
        j = jNew;
        assert(j!=0);
    }
    if(start==0){start = i+1;}

    //in case number of mismatches is less or equal to threshold
    if((dist[ls][la]).val <= mismatches)
    {
        // end is the first index out of sequence, since our strings in matrix are 1-indexed end is fine
        // start has to be start -= 1, to be the start in a 0-indexed string
        score = (dist[ls][la]).val;
        match_start = start-1;
        match_end = end;
        return true;
    }
    else
    {
        score = (dist[ls][la]).val;
    }
    return false;
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
