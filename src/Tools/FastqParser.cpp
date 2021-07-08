#include <iostream>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

//#include <stdlib.h>
//#include <limits.h>

#include "seqtk/kseq.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

KSEQ_INIT(gzFile, gzread)

struct barcode{
    std::string antibody; //BARCODE 1
    std::string sample; //BARCODE 2
    std::string umi;
};
typedef std::vector<barcode> barcodeVector;
using namespace boost::program_options;

struct input{
    std::string inFile;
    std::string outFile;

    std::string wellBarcodes;
    std::string abBarcodes;

    std::string anchor;
    int mb = 1;
    int ma = 2;
    int threads = 5;
};

struct fastqStats{
    int perfectMatches = 0;
    int noMatches = 0;
    int moderateMatches = 0;
};

bool parse_arguments(char** argv, int argc, input& input)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&(input.inFile))->required(), "directory of files or single file in fastq(.gz) format")
            ("output,o", value<std::string>(&(input.outFile))->required(), "output file with all split barcodes")

            ("anchor,a", value<std::string>(&(input.anchor))->required(), "anchor sequence")
            ("well,w", value<std::string>(&(input.wellBarcodes))->required(), "well barcodes")
            ("antibody,ab", value<std::string>(&(input.abBarcodes))->required(), "antibody barcodes")

            ("mb", value<int>(&(input.mb))->default_value(1), "mismatches allowed in barcodes")
            ("ma", value<int>(&(input.ma))->default_value(2), "mismatches allowed in anchor sequence")
            ("t", value<int>(&(input.threads))->default_value(5), "number of threads")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }

        notify(vm);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    return true;

}
 
const char *bitap_fuzzy_bitwise_search(const char *text, const char *pattern, int k)
{
    const char *result = NULL;
    int m = strlen(pattern);
    unsigned long *R;
    unsigned long pattern_mask[CHAR_MAX+1];
    int i, d;

    if (pattern[0] == '\0') return text;
    if (m > 31) return "The pattern is too long!";

    /* Initialize the bit array R */
    R = (unsigned long*) malloc((k+1) * sizeof *R);
    for (i=0; i <= k; ++i)
        R[i] = ~1;

    /* Initialize the pattern bitmasks */
    for (i=0; i <= CHAR_MAX; ++i)
        pattern_mask[i] = ~0;
    for (i=0; i < m; ++i)
        pattern_mask[pattern[i]] &= ~(1UL << i);

    for (i=0; text[i] != '\0'; ++i) {
        /* Update the bit arrays */
        unsigned long old_Rd1 = R[0];

        R[0] |= pattern_mask[text[i]];
        R[0] <<= 1;

        for (d=1; d <= k; ++d) {
            unsigned long tmp = R[d];
            /* Substitution is all we care about */
            R[d] = (old_Rd1 & (R[d] | pattern_mask[text[i]])) << 1;
            old_Rd1 = tmp;
        }

        if (0 == (R[k] & (1UL << m))) {
            result = (text+i - m) + 1;
            break;
        }
    }

    free(R);
    return result;
}

struct levenshtein_value{
        int val = 0;
        int start = 0;
        levenshtein_value(int a,int b):val(a), start(b){}
        levenshtein_value():val(0), start(0){}
};
levenshtein_value min(levenshtein_value a, levenshtein_value b)
{
    if(a.val < b.val)
    {
        return a;
    }
    return b;
}
bool levenshtein(const std::string seq, std::string anchor, int ma, int mb, int& match_start, int& match_end)
{
    int i,j,ls,la,t,substitutionValue, deletionValue;
    //stores the lenght of strings s1 and s2
    ls = seq.length();
    la = anchor.length();
    levenshtein_value dist[ls+1][la+1];

    //allow unlimited deletions in beginning; start is zero
    for(i=0;i<=ls;i++) {
        levenshtein_value val(0,i);
        dist[i][0] = val;
    }
    //deletions in pattern is punished; start is zero
    for(j=0;j<=la;j++) {
        levenshtein_value val(j,0);
        dist[0][j] = val;
    }
    bool solutionFound = false;
    int end = 0;
    int start = 0;

    for (i=1;i<=ls;i++) 
    {
        if(end) break;
        for(j=1;j<=la;j++) 
        {
            //Punishement for substitution
            if(seq[i-1] == anchor[j-1]) 
            {
                substitutionValue= 0;
            } 
            else 
            {
                substitutionValue = 1;
            }
            //punishement for deletion (allow deletions on beginnign and end)
            if(j==1 | j==la) 
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
            levenshtein_value seq_ins = dist[i][j-1];
            seq_ins.val += 1;
            levenshtein_value subst = dist[i-1][j-1];
            subst.val += substitutionValue;

            levenshtein_value tmp1 = min(seq_del, seq_ins);
            levenshtein_value tmp2 = min(tmp1, subst);
            
            //safe the moment a start begins
            if(tmp2.start == 0 & substitutionValue == 0)
            {
                tmp2.start = i;
            }

            dist[i][j] = tmp2;

            if(j==la & (dist[i][j]).val<=ma & i>=la) 
            {
                end = i;
                start = (dist[i][j]).start;
                break;
            }
        }
    }

    if(end)
    {
        match_start = start;
        match_end = end;
        return true;
    }
    return false;
}

bool split_sequence(const std::string seq, input* input, barcode& barcode, fastqStats& stats)
{
    //match regex of line
    std::smatch sm;
    std::string expressionString = "[ACTGN]+([ACTGN]{15,15})([ACTGN]{10,10})" + input->anchor + "([ACTGN]{10,10})[ACTGN]+";
    std::regex expression(expressionString);
    std::string line = seq;
    std::regex_match(line,sm,expression);

    //check of mismatches
    if(sm.size()<4)
    {
        int start = 0, end = 0;
        if(levenshtein(seq, input->anchor, input->ma, input->mb, start, end))
        {
            //std::cout << "++++++++++++++++\n";
            //std::cout << seq << "\n";
            //std::cout << input->anchor << "\n";
            //std::cout << start << "_" << end << "\n";
            //std::cout << seq.substr(start, end-start) << "\n";
            //minor mismatches that are allowed per mb and ma
            ++stats.moderateMatches;
            return true;
        }
        else
        {
            //bigger mismatches
            ++stats.noMatches;
            return false;
        }

    }
    //perfect matches
    else
    {
        assert(sm.size()==4);
        barcode.umi = sm[1];
        barcode.antibody = sm[2];
        barcode.sample = sm[3];
        ++stats.perfectMatches;

        return true;
    }

    return false;
}

void write_file(std::string output, barcodeVector barcodes)
{

}

void generate_barcodes(std::vector<std::string> fastqLines, input* input, barcodeVector& barcodes, fastqStats& stats)
{
    barcode barcode;
    for(const std::string& line : fastqLines)
    {
        if(split_sequence(line, input, barcode, stats))
        {
            barcodes.push_back(barcode);
        }
    }
}

barcodeVector iterate_over_fastq(input& input)
{
    //read all fastq lines into str vector
    gzFile fp;
    kseq_t *ks;
    fp = gzopen(input.inFile.c_str(),"r");
    if(NULL == fp){
        fprintf(stderr,"Fail to open file: %s\n", input.inFile.c_str());
    }
    ks = kseq_init(fp);
    std::vector<std::string> fastqLines;
    while( kseq_read(ks) >= 0 ){
        fastqLines.push_back(std::string(ks->seq.s));
    }

    //split input lines into thread buckets
    std::vector<std::vector<std::string> > fastqLinesVector;
    int element_number = fastqLines.size() / input.threads;
    std::vector<std::string>::iterator begin = fastqLines.begin();
    for(int i = 0; i < input.threads; ++i)
    {
        std::vector<std::string>::iterator end = ( i==input.threads-1 ? fastqLines.end() : begin + element_number);

        std::vector<std::string> tmpFastqLines(begin, end);
        fastqLinesVector.push_back(tmpFastqLines);
        begin = end;
    }

    //for every bucket call a thread
    //tmp variables to store thread results
    std::vector<barcodeVector> barcodesThreadList(input.threads);
    std::vector<fastqStats> statsThreadList(input.threads);
    std::vector<std::thread> workers;
    //final data variables
    barcodeVector barcodeVectorFinal;
    fastqStats fastqStatsFinal;
    for (int i = 0; i < input.threads; ++i) {
        workers.push_back(std::thread(generate_barcodes, fastqLinesVector.at(i), &input, std::ref(barcodesThreadList.at(i)), std::ref(statsThreadList.at(i))));
    }
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    //combine thread data
    for (int i = 0; i < input.threads; ++i) 
    {
        barcodeVectorFinal.insert(barcodeVectorFinal.end(), barcodesThreadList.at(i).begin(), barcodesThreadList.at(i).end());
        fastqStatsFinal.perfectMatches += statsThreadList.at(i).perfectMatches;
        fastqStatsFinal.noMatches += statsThreadList.at(i).noMatches;
    }
    std::cout << "MATCHED: " << fastqStatsFinal.perfectMatches << " | MISMATCHED: " << fastqStatsFinal.noMatches << "\n";
    kseq_destroy(ks);
    gzclose(fp);
    return barcodeVectorFinal;
}

int main(int argc, char** argv)
{
    input input;
    parse_arguments(argv, argc, input);
    barcodeVector barcodes = iterate_over_fastq(input);
    write_file(input.outFile, barcodes);
 
    return EXIT_SUCCESS;
}