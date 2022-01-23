<?php
/**
Script to run CI data analysis:
This can contain barcodes for ABs, sgRNA, etc.


*/

function run($string) 
{
    $logfile = "./bin/analysis.log";
    $fileHandle = fopen($logfile, "a+") or die("Unable to open file!");

    //keep output of command (stderr is written to terminal)
    fwrite($fileHandle, "#COMMAND: ".$string."\n");
    $output = shell_exec($string);
    fwrite($fileHandle, $output);

    fclose($fileHandle);
}


######################################
//parse command line arguments / Preprocessing
######################################

$shortopts = "i:"; //input file
$shortopts .= "o:"; //output file

//Demultiplexing parameters
$shortopts .= "m:"; //mismatch string: e.g.=1,2,3 (Mismatches for each pattern, the constant, varying like barcode, AB, etc. [NNN...] ones as well as the UMI one [XXX])
$shortopts .= "p:"; //mapping pattern: e.g.=[GAGCTCTGCGACG][XXXXXXXXXXXX][GATGACTCA]
$shortopts .= "b:"; //barcode file

//script parameters
$shortopts .= "t:"; //number of threads maximally used for tools

//Index (starting at zero) defines the x-th line in barcode file, or x-th 'NNN...' pattern
$shortopts .= "c:"; //list of indexes used for defining a single cell.
$shortopts .= "a:"; //list of ABs, in order of the barcodes in barcode file
$shortopts .= "x:"; //index for AB
$shortopts .= "g:"; //list of groups (e.g. treatments: can be a CI-round, gRNA), in order of the barcodes in barcode file
$shortopts .= "y:"; //index for treatment

$longopts = array("help");

$options = getopt($shortopts, $longopts);
$mapRNA = false;

//write help message
if(array_key_exists("help", $options))
{
    print("Call analyze.php with following options:\n

    #MANDATORY PARAMETERS FOR DEMULTIPLEXING:
    -i input file (fastq(.gz) file. Paired-end or single read.)
    -o output file (TSV file storing counts for each AB/gene per cell)
    -p mapping pattern (AGCT for constant sequence; X for UMI; N for a barcode this can be a CI-barcode,
       AB barcode, sgRNA), e.g. [AGGCAGTC][XXXXXXXXXXX][NNNN][GACTCAGAGC][NNNNN]
       (for RNA use the char 'R')
    -b barcode file: comma seperated file holding the possible barcodes for each pattern with X. Each line
       is a new barcode, in the same order is given in mapping pattern

    #MANDATORY TO MAP BARCODES TO SINGLE CELLS/ ABS/ TREATMENT
    -c list of indices: the indices of the mapping pattern (starting at 0) used as CI-barcodes to define a single cell,
    -a list of indices: the index of mapping pattern (starting at 0) that defines the AB
    -x index defining the AB
    -g ist of indices: the index of mapping pattern (starting at 0) that defines the treatment
    -y the index defining the treatment

    #NON MANDATORY PARAMTERS
    -m mismatches per patter (default = 1)
    -t threads (default = 1)
    \n");
    exit();
}

//parse the mismatch option to get the allowed UMI mismatches
if($options['m'])
{
    $patternString = $options['p'];
    //check if we have a UMI barcode
    if(strpos($patternString, 'X')!==false || strpos($patternString, 'x')!==false)
    {
        $index = NULL;
        //get the index of UMI mismatch
        $patterns = explode("]", $patternString);
        $count = 0;
        print_r($patterns);
        foreach($patterns as $pattern)
        {
            $pattern = trim($pattern, "[");
            if(preg_match("/^(x|X)+$/", $pattern))
            {
                $index = $count;
                break;
            }
            $count += 1;
        }

        //write the umi mismatch seperately
        $options['u'] = $index;
    }
}

if(strpos($patternString, 'R')!==false || strpos($patternString, 'r')!==false)
{
    $mapRNA = true;
}

$logfile = "./bin/analysis.log";
@unlink($logfile, ); //delete old file if exists

######################################
//execute demultiplexing
######################################

$command = "./bin/demultiplexing";
$parameters = ['i', 'o', 'p', 'b', 'm', 't'];
foreach($parameters as $element)
{
    if($options[$element])
    {
        $command = $command . " -" . $element . " " . $options[$element];
    }
}
run($command);
//gzip the output

######################################
//execute mapping of RNA if necessary (if mapping pattern contains RNA region)
######################################
if($mapRNA)
{
    print("RNA mapping not implemented yet!!!!\n");
    exit();
}

######################################
//execute processing of the mapped barcodes (UMI collapsing, assigning single cells, etc.)
######################################

//set new input file = output of Demultiplexing

$command = "./bin/processing ";
$parameters = ['o', 'b', 'u', 't', 'a', 'x', 'g', 'y'];
foreach($parameters as $element)
{
    if($options[$element])
    {
        $command = $command . " -" . $element . " " . $options[$element];
    }
}
run($command);

######################################
//execute some Quality Control scripts / Data Analysis scripts
######################################

?>