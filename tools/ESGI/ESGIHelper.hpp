// THIS FILE CONTAINS ALL THE TOOL-CALL NEEDED FOR ESGI:
// like demultiplexing, counting, etc

#include <cstdlib>
#include <string>
#include <iostream>
#include <filesystem>

//check if we have patterns like -,*. These are not included in the demultiplexed output and therefore
//indices for counting must be adjusted. Additioanly, these two patterns can ONLY BE at the end of a sequence, 
//e.g., we can t have an arbitrary pattern * at ther beginning of a read
inline int get_special_pattern_pos(const std::vector<std::string>& patterns)
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

    return -1;
}
//get the position of the DNA pattern
inline int get_DNA_pattern_pos(const std::vector<std::string>& patterns)
{
    auto it = std::find(patterns.begin(), patterns.end(), "DNA");
    return( (it == patterns.end() ? -1 : static_cast<int>(it - patterns.begin())) );
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

// --- capture helper (kept minimal) ---
static std::string run_system_capture(const std::string& cmd) {
    std::string tmp = "star_ver_" + std::to_string(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    ) + ".txt";
    std::string full = cmd + " > \"" + tmp + "\" 2>&1";
    if (std::system(full.c_str()) == -1) return std::string();
    std::ifstream f(tmp, std::ios::in | std::ios::binary);
    std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    f.close();
    std::remove(tmp.c_str());
    return s;
}

// parse "A.B.C[letter]" from an arbitrary string (finds first number run)
static bool parse_version_triplet(const std::string& src,
                                  int& A, int& B, int& C, int& S /*-1 no suffix*/) {
    A = B = C = 0; S = -1;
    size_t i = src.find_first_of("0123456789");
    if (i == std::string::npos) return false;

    auto read_int = [&](size_t& p, int& out)->bool {
        if (p >= src.size() || !std::isdigit((unsigned char)src[p])) return false;
        long v = 0;
        while (p < src.size() && std::isdigit((unsigned char)src[p])) {
            v = v*10 + (src[p]-'0'); ++p;
        }
        out = (int)v; return true;
    };

    if (!read_int(i, A)) return false;
    if (i >= src.size() || src[i++] != '.') return false;
    if (!read_int(i, B)) return false;
    if (i >= src.size() || src[i++] != '.') return false;
    if (!read_int(i, C)) return false;

    if (i < src.size() && std::isalpha((unsigned char)src[i])) {
        S = std::tolower((unsigned char)src[i]) - 'a'; // a=0, b=1, ...
    } else {
        S = -1; // no suffix = newest among same A.B.C
    }
    return true;
}

// return positive if (A1,B1,C1,S1) > (A2,B2,C2,S2), 0 if equal, negative if less
static int compare_versions(int A1,int B1,int C1,int S1, int A2,int B2,int C2,int S2) {
    if (A1 != A2) return (A1 > A2) ? 1 : -1;
    if (B1 != B2) return (B1 > B2) ? 1 : -1;
    if (C1 != C2) return (C1 > C2) ? 1 : -1;
    if (S1 == S2) return 0;
    if (S1 == -1) return 1;        // no suffix newer than any suffix
    if (S2 == -1) return -1;
    return (S1 > S2) ? 1 : -1;     // a(0) < b(1) < c(2) ...
}

// returns 0 if STAR exists and version >= minimumVersion (e.g. "2.7.10b"), else 1
static int check_star_version(const std::string& minimumVersion, std::string& currentVersion) 
{
    std::string out = run_system_capture("STAR --version");
    currentVersion = out;
    if (out.empty()) { std::cerr << "STAR not found or no output\n"; return 1; }

    int a,b,c,s, ra,rb,rc, rs;
    if (!parse_version_triplet(out, a,b,c,s)) {
        std::cerr << "Could not parse STAR version from output: " << out << "\n";
        return 1;
    }
    if (!parse_version_triplet(minimumVersion, ra,rb,rc, rs)) {
        std::cerr << "Invalid required version string: " << minimumVersion << "\n";
        return 1;
    }

    int cmp = compare_versions(a,b,c,s, ra,rb,rc,rs);
    return (cmp >= 0) ? 1 : 0;
}

inline bool run_star(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
    std::string starExecutable = "STAR";
    if(config.STAR.has_value()){starExecutable = config.STAR.value();}

    std::filesystem::path outfilePrefixPath = std::filesystem::path(config.output) / "STAR_";
    std::string outfilePrefix = outfilePrefixPath.string();

    // we need a STAR version > 2.7.10b to also get GX/GN tags
    std::string miniumStar = "2.7.10b";
    std::string currentVersion;
    if(!check_star_version(miniumStar, currentVersion))
    {
        std::cerr << "STAR VERSION is not compatible with ESGI, since we require GX/GN annotations from STAR.\n"
         << " We need a minimum version of STAR " << miniumStar << "\n"
         << " Current version is: " << currentVersion << "\n";
        std::cerr << "Please upgrade STAR to a newer version, or set the path for STAR in the input.ini file with <STAR=...>\n";
        exit(EXIT_FAILURE);
    }

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
    std::string demultiplexingPatternOutput = intermediateFiles.demultiplexingOutput + ".tsv";
    // we hard-coded and gave the STAR-output the prefix <STAR>
    std::string starAlignmentOutput = config.output + "/" + "STAR_Aligned.out.bam";

    const std::string annotateCmd = std::string("./bin/annotate ") +
        " -i " + demultiplexingPatternOutput +
        " -b " + starAlignmentOutput +
        " -f " + config.starFeature;

    int rc = std::system(annotateCmd.c_str());
    if (rc != 0) 
    {
        std::cerr << "annotate failed with code " << rc << "\n";
    }

    return true;
}

inline bool run_rna_mapping(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
    //call STAR
    if(!run_star(config, intermediateFiles))
    {
        std::cerr << "Error: STAR failed\n";
    }

    //call annotate
    if(!run_annotate(config, intermediateFiles))
    {
        std::cerr << "Error: Annotation of demultiplexed reads with STAR-mapping failed\n";
    }

    return true;
}

std::string adjust_position_due_to_special_patterns(const std::string& indexList, const ESGIIntermediateFiles& intermediateFiles)
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
        //check if we had a DNA pattern that got 'cut out' and therefore need to adjust pattern positions
        int dnaShift = 0;
        if(intermediateFiles.dnaPatternPos >= 0 && v > intermediateFiles.dnaPatternPos){dnaShift = -1;}
        //if the index is behind a -,* we need to substract one bcs this poattern is not in the demultiplex output
        if (v > intermediateFiles.specialPatternPos) v -= 1;
        //also we need to add 1 to every column, bcs we will ahve the additional READNAME column
        v += 1;
        //and finally add the dna shift
        v += dnaShift;
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

inline bool call_count(ESGIConfig& config, const ESGIIntermediateFiles& intermediateFiles, const bool dnaPatternPresent)
{
    std::filesystem::path toolPath = std::filesystem::path("bin") / "count";
    std::filesystem::path patternPath(config.patternFile);
    std::filesystem::path dir = patternPath.parent_path();

    //index re-assignment due to -,*
    if(intermediateFiles.specialPatternPos!= -1 || dnaPatternPresent)
    {
        config.SC_ID = adjust_position_due_to_special_patterns(config.SC_ID, intermediateFiles);
        if(intermediateFiles.dnaPatternPos >=0)
        {
            // if we have a DNA pattern, we firstly set the FEATURE column index to patternlength
            // as the DNA is added as a last column, o-indexed it is at position patternlength
            config.FEATURE_ID = std::to_string(intermediateFiles.patternLength);            
            //however, if there were special patterns (like -,*) that are not in the demultiplexed output we need to adust the FEATURE column index
            //basically: -1 for removing the DNA pattern, -1 if there is a -,* and +1 for READNAME col
            config.FEATURE_ID = adjust_position_due_to_special_patterns(config.FEATURE_ID, intermediateFiles);
        }
        else
        {
            config.FEATURE_ID = adjust_position_due_to_special_patterns(config.FEATURE_ID, intermediateFiles);
        }
        if(config.UMI_ID.has_value()){config.UMI_ID = adjust_position_due_to_special_patterns(config.UMI_ID.value(), intermediateFiles);}
        if(config.ANNOTATION_IDs.has_value()){config.ANNOTATION_IDs = adjust_position_due_to_special_patterns(config.ANNOTATION_IDs.value(), intermediateFiles);}
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
        " --umiAbundanceThreshold " +  std::to_string(config.umiAbundance) +
        " --umiRemoval " +  std::to_string(config.umiCollapsing) +
        " --hamming " + std::to_string(config.hamming);

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