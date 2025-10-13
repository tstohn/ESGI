#pragma once
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <optional>
#include <algorithm>
#include <cctype>
#include <iostream>

struct ESGIConfig 
{
    // REQUIRED ARGUMENTS
    std::string forward;        
    std::string output;         
    std::string patternFile;        
    std::string mismatchesFile;
    std::string SC_ID;          
    std::string FEATURE_ID;     

    // OPTIONAL ARGUMENTS
    std::optional<std::string> reverse;       
    std::optional<std::string> FEATURE_NAMES;
    std::optional<std::string> ANNOTATION_IDs;
    std::optional<std::string> ANNOTATION_NAMES;
    std::optional<std::string> UMI_ID;
    std::optional<std::string> prefix;
    std::optional<int>         fastqReadBucketSize;
    std::optional<std::string> barcodeSharing;

    std::optional<std::string> STAR;
    std::optional<std::string> genomeDir;

    // DEFAULT ARGUMENTS
    bool independent = false;      
    bool hamming = false;
    bool writeFailedLines = true;    
    bool writeStats       = true;   
    bool umiCollapsing = true;
    bool SC_ID_string = false;
    int  threads = 5;
    double umiThreshold = 0.0;
    std::string starFeature = "GX";

    void write(std::ostream& os = std::cout) const
    {
        os << "FOLLOWING PARAMETERS ARE SET FOR ESGI:\n";

        os << "\tforward="        << forward        << "\n";
        os << "\toutput="         << output         << "\n";
        os << "\tpatternFile="    << patternFile    << "\n";
        os << "\tmismatchesFile=" << mismatchesFile << "\n";
        os << "\tSC_ID="          << SC_ID          << "\n";
        os << "\tFEATURE_ID="     << FEATURE_ID     << "\n\n";

        if (reverse)           os << "\treverse="           << *reverse           << "\n";
        if (FEATURE_NAMES)     os << "\tFEATURE_NAMES="     << *FEATURE_NAMES     << "\n";
        if (ANNOTATION_IDs)    os << "\tANNOTATION_IDs="    << *ANNOTATION_IDs    << "\n";
        if (ANNOTATION_NAMES)  os << "\tANNOTATION_NAMES="  << *ANNOTATION_NAMES  << "\n";
        if (UMI_ID)            os << "\tUMI_ID="            << *UMI_ID            << "\n";
        if (prefix)            os << "\tprefix="            << *prefix            << "\n";
        if (fastqReadBucketSize)
                               os << "\tfastqReadBucketSize=" << *fastqReadBucketSize << "\n";
        if (barcodeSharing)    os << "\tbarcodeSharing="    << *barcodeSharing    << "\n";
        os << "\n";

        os << std::boolalpha;
        os << "\tindependent="      << independent      << "\n";
        os << "\thamming="          << hamming          << "\n";
        os << "\twriteFailedLines=" << writeFailedLines << "\n";
        os << "\twriteStats="       << writeStats       << "\n";
        os << "\tumiCollapsing="    << umiCollapsing    << "\n";
        os << "\tSC_ID_string="     << SC_ID_string     << "\n";
        os << "\tthreads="          << threads          << "\n";
        os << "\tumiThreshold="        << umiThreshold        << "\n";

        //Star values
        if (reverse)             os << "\treverse="                << *reverse           << "\n";
        if (STAR)                os << "\tSTAR executable="        << *STAR           << "\n";
        //print the star feature to annotate only if a genomeDir is given, since
        if (genomeDir)           os << "\tSTAR feature="           << starFeature           << "\n";

    }
};

struct ESGIIntermediateFiles
{

    std::string demultiplexingOutput;
    int specialPatternPos=-1;
    int umiMismatches=-1;

    std::string countingInput;
    std::string countingOutput;

};


namespace parsing {

inline std::string ltrim(std::string s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
        [](unsigned char ch){ return !std::isspace(ch); }));
    return s;
}
inline std::string rtrim(std::string s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
        [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
    return s;
}
inline std::string trim(std::string s) { return rtrim(ltrim(std::move(s))); }

inline std::string to_upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return s;
}

inline bool parse_bool(std::string v) 
{
    v = to_upper(trim(v));
    if (v=="1" || v=="TRUE" || v=="YES")  return true;
    if (v=="0" || v=="FALSE"|| v=="NO") return false;
    throw std::runtime_error("Invalid boolean value: '" + v + "'");
}

inline int parse_int(const std::string& v) {
    try {
        size_t idx=0;
        int val = std::stoi(trim(v), &idx);
        if (idx != trim(v).size()) throw std::runtime_error("");
        return val;
    } catch (...) {
        throw std::runtime_error("Invalid integer value: '" + v + "'");
    }
}

inline double parse_double(const std::string& v) {
    try {
        size_t idx=0;
        double val = std::stod(trim(v), &idx);
        if (idx != trim(v).size()) throw std::runtime_error("");
        return val;
    } catch (...) {
        throw std::runtime_error("Invalid double value: '" + v + "'");
    }
}

}//namespace

// parse the ini file and return an ESGIConfig
inline ESGIConfig loadESGIConfigFromFile(const std::string& path)
{
    std::ifstream in(path);
    if (!in) 
    {
        throw std::runtime_error("Cannot open ini file: " + path);
    }

    ESGIConfig cfg;

    std::string line;
    int lineNo = 0;

    auto set_required_string = [&](std::string& target, const std::string& val) 
    {
        if (val.empty()) throw std::runtime_error("Empty value for a required field");
        target = val;
    };

    while (std::getline(in, line)) 
    {
        ++lineNo;

        // strip comments
        auto hashPos = line.find('#');
        if (hashPos != std::string::npos) {
            line = line.substr(0, hashPos);
        }

        line = parsing::trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;
    
        // split "key=value" (first '=' only)
        auto eqPos = line.find('=');
        if (eqPos == std::string::npos) 
        {
            throw std::runtime_error("Invalid line (missing '=') at " + std::to_string(lineNo));
        }

        std::string key = parsing::trim(line.substr(0, eqPos));
        std::string value = parsing::trim(line.substr(eqPos+1));

        std::string KEY = parsing::to_upper(key);

        try {
            // REQUIRED
            if (KEY == "FORWARD")          set_required_string(cfg.forward, value);
            else if (KEY == "OUTPUT")      set_required_string(cfg.output, value);
            else if (KEY == "PATTERN")     set_required_string(cfg.patternFile, value);
            else if (KEY == "MISMATCHES")  set_required_string(cfg.mismatchesFile, value);
            else if (KEY == "SC_ID")       set_required_string(cfg.SC_ID, value);
            else if (KEY == "FEATURE_ID")  set_required_string(cfg.FEATURE_ID, value);

            // OPTIONAL
            else if (KEY == "REVERSE")             cfg.reverse = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "FEATURE_NAMES")       cfg.FEATURE_NAMES = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "ANNOTATION_IDS")      cfg.ANNOTATION_IDs = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "ANNOTATION_NAMES")    cfg.ANNOTATION_NAMES = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "UMI_ID")              cfg.UMI_ID = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "PREFIX")              cfg.prefix = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "BARCODESHARING")      cfg.barcodeSharing = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "STAR")                cfg.STAR = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};
            else if (KEY == "GENOMEDIR")           cfg.genomeDir = value.empty() ? std::optional<std::string>{} : std::optional<std::string>{value};

            // DEFAULTED SCALARS
            else if (KEY == "THREADS")             cfg.threads = parsing::parse_int(value);
            else if (KEY == "INDEPENDENT")         cfg.independent = parsing::parse_bool(value);
            else if (KEY == "WRITEFAILEDLINES")    cfg.writeFailedLines = parsing::parse_bool(value);
            else if (KEY == "WRITESTATS")          cfg.writeStats = parsing::parse_bool(value);
            else if (KEY == "HAMMING")             cfg.hamming = parsing::parse_int(value);
            else if (KEY == "FASTQREADBUCKETSIZE") cfg.fastqReadBucketSize = parsing::parse_int(value);
            else if (KEY == "UMITHRESHOLD")        cfg.umiThreshold = parsing::parse_double(value);
            else if (KEY == "UMICOLLAPSING")       cfg.umiCollapsing = parsing::parse_bool(value);
            else if (KEY == "SC_ID_STRING")        cfg.SC_ID_string = parsing::parse_bool(value);
            else if (KEY == "FEATURE")             cfg.starFeature = value;
            else 
            {
                // Unknown key: ignore or warn. Here we warn to stderr.
                std::cerr << "[ESGIConfig] Warning: unknown key '" << key
                          << "' at line " << lineNo << " — ignoring.\n";
            }

        } catch (const std::exception& e) 
        {
            throw std::runtime_error(
                "Error parsing key '" + key + "' at line " + std::to_string(lineNo) + ": " + e.what()
            );
        }
    }

    // check for necessary parameters
    auto missing = std::vector<std::string>{};
    if (cfg.forward.empty())        missing.push_back("forward");
    if (cfg.output.empty())         missing.push_back("output");
    if (cfg.patternFile.empty())    missing.push_back("pattern");
    if (cfg.mismatchesFile.empty()) missing.push_back("mismatches");
    if (cfg.SC_ID.empty())          missing.push_back("SC_ID");
    if (cfg.FEATURE_ID.empty())     missing.push_back("FEATURE_ID");

    if (!missing.empty()) 
    {
        std::ostringstream oss;
        oss << "Missing required parameter(s): ";
        for (size_t i=0;i<missing.size();++i) {
            if (i) oss << ", ";
            oss << missing[i];
        }
        throw std::runtime_error(oss.str());
    }

    // fastqReadBucketSize default: threads * 100000 (if not explicitly set)
    if (!cfg.fastqReadBucketSize.has_value()) 
    {
        cfg.fastqReadBucketSize = cfg.threads * 100000;
    }

    return cfg;
}

//check if there is a RNA pattern present - since windows can not handle it
bool containsDNAPattern(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file) 
    {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) 
    {
        if (line.find("[DNA]") != std::string::npos) return true;
    }

    return false; 
}