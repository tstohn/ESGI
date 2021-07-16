#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cassert>

//stores all the input parameters
struct input{
    std::string inFile;
    std::string outFile;

    std::string barcodeFile;
    std::string mismatchLine;
    std::string patternLine;

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
        int val = 0;
        int i = 0;
        int j = 0;
        levenshtein_value(int a,int b, int c):val(a), i(b), j(c){}
        levenshtein_value():val(0), i(0), j(0){}
};
levenshtein_value min(levenshtein_value a, levenshtein_value b, bool& first)
{
    if(a.val <= b.val)
    {
        first = true;
        return a;
    }
    return b;
}
bool levenshtein(const std::string sequence, std::string pattern, int mismatches, int& match_start, int& match_end, int& score)
{
    int i,j,ls,la,t,substitutionValue, deletionValue;
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
            //if( j==1 | j==la )
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
            
            dist[i][j] = tmp2;
            //std::cout << tmp2.val << "(" << sequence[i-1] <<  ","<< pattern[j-1] << ")" << " ";
        }
        //std::cout << "\n";
    }

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

        if(dist[i][j].val == dist[iNew][jNew].val & i!=iNew & j!=jNew)
        {
            start = iNew;
        }

        if( (jNew<la) & noEnd)
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
    return false;
}
