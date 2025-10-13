// THIS FILE CONTAINS ALL THE TOOL-CALL NEEDED FOR ESGI:
// like demultiplexing, counting, etc

#include <cstdlib>
#include <string>
#include <iostream>
#include <filesystem>

//check if we have patterns like -,*. These are not included in the demultiplexed output and therefore
//indices for counting must be adjusted. Additioanly, these two patterns can ONLY BE at the end of a sequence, 
//e.g., we can t have an arbitrary pattern * at ther beginning of a read
inline int handle_special_patterns(const std::vector<std::string>& patterns)
{
    int posDash = -1;
    int posStar = -1;

    // find positions, ensure the patterns are only used once
    for (size_t i = 0; i < patterns.size(); ++i) 
    {
        if (patterns[i] == "-" && posDash==-1) {posDash = static_cast<int>(i);}
        else if (patterns[i] == "[-]" && posDash!=-1){std::cerr << "The pattern [-] can only be used ONCE in a pattern! (It denotes the end of fw/rv reads, e.g., for cases where there is no overlap between reads\n";}
        if (patterns[i] == "*" && posStar==-1) {posStar = static_cast<int>(i);}
        else if (patterns[i] == "[*]" && posDash!=-1){std::cerr << "The pattern [*] can only be used ONCE in a pattern! Everything 'behind' this pattern in reading direction is not considered anymore for demultiplexing\n";}
    }

    // Case 1: neither found
    if (posDash == -1 && posStar == -1) return-1;

    // Case 2: only one found
    if (posDash != -1) return posDash;
    if (posStar != -1) return posStar;

    // Case 3: both found - throw error
    if (posStar != -1 && posDash != -1)
    {
        std::cerr << "Error: [-] and [*] can not be used together. You can only use one of them! [*] already implicitely denotes the end of the strand, as nothing is demultiplexed beyond this sign\
        making an additional [-] obsolete\n";
    }

}

inline bool call_demultiplex(const ESGIConfig& config)
{
    std::filesystem::path toolPath = std::filesystem::path("bin") / "demultiplex";

    // Build command string with parameters
    std::string cmd = toolPath.string() +
        " --input " + config.forward +
        " --output " + config.output +
        " --threads " + std::to_string(config.threads) +
        " --barcodePatternsFile " + config.patternFile +
        " --mismatchFile " + config.mismatchesFile +
        " --independent " + std::to_string(config.independent) +
        " --hamming " + std::to_string(config.hamming) + 
        " --writeFailedLines " + std::to_string(config.writeFailedLines) + 
        " --writeStats " + std::to_string(config.writeStats);

    // optional arguments
    if(config.reverse.has_value()){cmd += " --reverse " + config.reverse.value();}
    if(config.prefix.has_value()){cmd += " --namePrefix " + config.prefix.value();}
    if(config.fastqReadBucketSize.has_value()){cmd += " --fastqReadBucketSize " + std::to_string(config.fastqReadBucketSize.value());}

    std::cout << "Running demultiplexing as: " << cmd.c_str() << "\n";
    int ret = std::system(cmd.c_str());
    if (ret != 0) 
    {
        std::cout << "DEMULTIPLEX failed with code " << ret << "\n";
        return false;
    }

    return true;
}






inline bool run_star(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
    std::string starExecutable = "STAR";
    if(config.STAR.has_value()){starExecutable = config.STAR.value();}

    std::filesystem::path outfilePrefixPath = std::filesystem::path(config.output) / "STAR_";
    std::string outfilePrefix = outfilePrefixPath.string();

    if(!config.genomeDir.has_value())
    {
        std::cerr << "No parameter <genomeDir> provided, please provide this parameter in the.ini file to enable STAR\n";
    }
    const std::string starCmd = starExecutable +
        " --runThreadN " + std::to_string(config.threads) +
        " --genomeDir " + config.genomeDir.value() +
        " --readFilesIn " + intermediateFiles.demultiplexingOutput + ".fastq" +
        " --outFileNamePrefix " + outfilePrefix +
        " --outSAMtype BAM Unsorted " +
        " --outSAMattributes NH HI AS nM GX GN " +
        " --quantMode TranscriptomeSAM " +
        " --outFilterMultimapNmax 50 " +
        " --outSAMmultNmax 1 --outSAMunmapped Within " +
        " --limitOutSJcollapsed 2000000 " +
        " --twopassMode Basic";

    int rc = std::system(starCmd.c_str());
    if (rc != 0) 
    {
        std::cerr << "STAR failed with code " << rc << "\n";
    }

    return true;
}

inline bool run_annotate(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
        const std::string annotateCmd =
        "/usr/bin/time -v /DATA/t.stohn/SCDemultiplexing/bin/annotate "
        "-i ./SIGNALseq_Analysis/output/ESGI_RNA/RNA.tsv "
        "-b ./SIGNALseq_Analysis/output/ESGI_RNA/RNA_Aligned.out.bam "
        "-f GX";
    int rc = std::system(annotateCmd.c_str());
    if (rc != 0) 
    {
        std::cerr << "annotate failed with code " << rc << "\n";
    }

    return true;
}

inline bool run_rna_mapping(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
    bool star;
    bool annotate;

    //call STAR
    star = run_star(config, intermediateFiles);

    //call annotate

    return true;
}

std::string adjust_position_due_to_special_patterns(const std::string& indexList, int specialPatternPos)
{
    std::stringstream ss(indexList);
    std::string index;
    std::vector<int> values;

    // Split by commas
    while (std::getline(ss, index, ',')) 
    {
        values.push_back(std::stoi(index));
    }

    // Modify values
    for (int &v : values) 
    {
        //if the index is behind a -,* we need to substract one bcs this poattern is not in the demultiplex output
        if (v > specialPatternPos) v -= 1;
        //also we need to add 1 to every column, bcs we will ahve the additional READNAME column
        v += 1;
    }

    // Rebuild the string
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) 
    {
        if(i) out << ",";
        out << values[i];
    }

    return(out.str());
}

inline bool call_count(ESGIConfig& config, const ESGIIntermediateFiles& intermediateFiles)
{
    std::filesystem::path toolPath = std::filesystem::path("bin") / "count";
    std::filesystem::path patternPath(config.patternFile);
    std::filesystem::path dir = patternPath.parent_path();

    //index re-assignment due to -,*
    if(intermediateFiles.specialPatternPos!= -1)
    {
        config.SC_ID = adjust_position_due_to_special_patterns(config.SC_ID, intermediateFiles.specialPatternPos);
        config.FEATURE_ID = adjust_position_due_to_special_patterns(config.FEATURE_ID, intermediateFiles.specialPatternPos);
        if(config.UMI_ID.has_value()){config.UMI_ID = adjust_position_due_to_special_patterns(config.UMI_ID.value(), intermediateFiles.specialPatternPos);}
        if(config.ANNOTATION_IDs.has_value()){config.ANNOTATION_IDs = adjust_position_due_to_special_patterns(config.ANNOTATION_IDs.value(), intermediateFiles.specialPatternPos);}
    }

    // Build command string with parameters
    std::string cmd = toolPath.string() +
        " --input " + intermediateFiles.countingInput +
        " --output " + intermediateFiles.countingOutput +
        " --threads " + std::to_string(config.threads) +
        " --barcodeDir " + dir.string() +
        " --featureIndex " + config.FEATURE_ID +
        " --singleCellIndices " + config.SC_ID + 
        " --umiThreshold " +  std::to_string(config.umiThreshold) +
        " --umiRemoval " +  std::to_string(config.umiCollapsing);

    // optional arguments
    if(config.FEATURE_NAMES.has_value()){cmd += " --featureNames " + config.FEATURE_NAMES.value();}
    if(config.ANNOTATION_IDs.has_value()){cmd += " --annotationIdxs " + config.ANNOTATION_IDs.value();}
    if(config.ANNOTATION_NAMES.has_value()){cmd += " --annotationFiles " + config.ANNOTATION_NAMES.value();}
    if(config.UMI_ID.has_value()){cmd += " --umiIndex " + config.UMI_ID.value();}
    if(config.barcodeSharing.has_value()){cmd += " --shareBarcodes " + config.barcodeSharing.value();}
    //if we did overwrite the umiMismatches for intermediate data
    if(intermediateFiles.umiMismatches != -1){cmd += " --mismatches " + std::to_string(intermediateFiles.umiMismatches);}

    std::cout << "Running count as: " << cmd.c_str() << "\n";
    int ret = std::system(cmd.c_str());
    if (ret != 0) 
    {
        std::cout << "COUNT failed with code " << ret << "\n";
        return false;
    }

    return true;
}