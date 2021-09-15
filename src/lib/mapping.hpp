#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cassert>
#include <string_view>
#include <cstring>
#include <cstdio>
#include <zlib.h>
#include <regex>
#include <thread>

#include "Barcode.hpp"
#include "seqtk/kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef std::vector< std::shared_ptr<std::string> > BarcodeMapping;
typedef std::vector<BarcodeMapping> BarcodeMappingVector;

//mapping sequentially each barcode leaving nmo pattern out,
//if a pattern can not be found the read is idscarded
class MapEachBarcodeSequentiallyPolicy
{

};

//mapping sequentially each barcode, however, if a pattern can not be found do not discard it
//but jump to next pattern
class MapEachBarcodeSequenciallyWithLeaveOutsPolicy
{
    
};

//mapping only constant barocdes as anchor first
class MapAroundConstantBarcodesAsAnchorPolicy
{
    
};

//class for barcode mappings, the mapping algorithm is chosen by the 
//mapping policy
template<typename MappingPolicy>
class Mapping : private MappingPolicy
{
    public:
    //process all the input information and check for validity
    //e.g. delete old output files if present, parse barcode from barcodeFile, match them to their
    //number of mismatches etc.
    void initialize_mapping();

    //run the actual mapping by using MappingPolicy
    bool run_mapping(const std::string& seq, input* input, BarcodeMapping& barcodeMap, BarcodeMapping& realBarcodeMap,
                    fastqStats& stats, std::map<std::string, std::shared_ptr<std::string> >& unique_seq,
                    BarcodePatternVectorPtr barcodePatterns);

};
