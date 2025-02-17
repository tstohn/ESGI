# include "Barcode.hpp"



//write a demultiplexed DNA line, it is a line that contains a DNA sequence, which still
//has to be mapped to a reference genome, and therefore written to a fastq file
void Barcode::write_dna_to_fastq(const std::vector<std::string> barcodeList, std::string dna)
{
    //rename the fastq read and give it a unique name

    //block the fastq-read counter with a mutex and add a new fastq read to file
}

// this is a line solely consisting of barcodes, it can be AB-tagged barcodes, spatial coordinates
// or guide barcodes, with or without UMIs etc. It will be written to a count file which still has to be
// 'counted' in a subsequent step for a count table/ UMI collapsing/ ...
void Barcode::write_barcodes_to_txt(const std::vector<std::string> barcodeList)
{
    //give the line a unique name (make sure its the same as fastq-read to later map them to each other)
}

// write  a line with DNA/ Any other line
void Barcode::write_demultiplexed_line(const std::vector<std::string> barcodeList, std::string dna)
{
    if(dna == "")
    {
        write_barcodes_to_txt(barcodeList);
    }
    else
    {
        write_dna_to_fastq(barcodeList, dna)
    }
}