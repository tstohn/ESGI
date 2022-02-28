<?php
/**
Script to run CI data analysis:
This can contain barcodes for ABs, sgRNA, etc.


*/

function delTree($dir) 
{
    $files = array_diff(scandir($dir), array('.','..'));
     foreach ($files as $file) {
       (is_dir("$dir/$file")) ? delTree("$dir/$file") : unlink("$dir/$file");
     }
     return rmdir($dir);
}

//run command and update terminal
function run_updates($command)
{
    if( ($fp = popen($command, "r")) ) 
    {
        while( !feof($fp) )
        {
            echo fread($fp, 1024);
            flush();
        }

        fclose($fp);
    }
}

function clean($tmpFolder, $logfileHandle)
{
    fclose($logfileHandle);
    delTree($tmpFolder);
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
$shortopts .= "g:"; //guide file: one line of comma seperated guide barcodes (must be at same position as AB barcodes)

//script parameters
$shortopts .= "t:"; //number of threads maximally used for tools

//Index (starting at zero) defines the x-th line in barcode file, or x-th 'NNN...' pattern
$shortopts .= "c:"; //list of indexes used for defining a single cell.
$shortopts .= "a:"; //list of AB names, in order of the barcodes in barcode file (to map barcode to AB name)
$shortopts .= "x:"; //index for AB (within barcodefile)
$shortopts .= "d:"; //list of treatment names (e.g. treatments: can be a CI-round, gRNA), in order of the barcodes in barcode file
$shortopts .= "y:"; //index for treatment (within barcodefile)
$shortopts .= "n:"; //list of name to map guide barcodes to a name (e.g. cell line name)
$shortopts .= "f:"; //threshold to use when filtering umis (default is zero). Keep a read of UMI only, if more than Xperc (as double) are from same single cell and AB

$longopts = array("help");

$options = getopt($shortopts, $longopts);
$mapRNA = false;
$pairedEnd = false;

//write help message
if(array_key_exists("help", $options))
{
    print("Call analyze.php with following options:\n

    #MANDATORY PARAMETERS FOR DEMULTIPLEXING:
    -i input file (fastq(.gz) file. Ether use only this parameter for Single read or use this parameter for the forward read for Paired-end.)
    -r reverse read if Paired-End data is input.
    -o output folder in which all the files are stored (from demultiplexed reads to single-cell <-> AB counts matrix). Defaults to './bin/CIResult'
    -p mapping pattern (AGCT for constant sequence; X for UMI; N for a barcode this can be a CI-barcode,
       AB barcode, sgRNA), e.g. [AGGCAGTC][XXXXXXXXXXX][NNNN][GACTCAGAGC][NNNNN]
       (for RNA use the char 'R')
    -b barcode file: comma seperated file holding the possible barcodes for each pattern with X. Each line
       is a new barcode, in the same order is given in mapping pattern (DO NOT include the guide barcodes here, add them as two seperate
       parameters as a file with guide sequences and guide names)

    -g guide barcode file: file containing only one lines of comma seperated barcodes. Those barcodes are the guide-barcodes, that have to be at the same position 
       as the AB barcode and are additional reads of the experiment to know from which cell line the cell origionates from.

    #MANDATORY TO MAP BARCODES TO SINGLE CELLS/ ABS/ TREATMENT
    -c list of indices: the indices of the mapping pattern (starting at 0) used as CI-barcodes to define a single cell,
    
    -a list of AB names: names of the ABs, in the same order as their barcodes in file '-b' to map barcode to AB names
    -x index defining the AB (index for the line within barcodefile that contains the AB barcodes (0-indexed))
    
    -d list of treatment names: names of treatments in same order as the barcodes in '-b' file to map cells to a treatment (e.b. in barcoding round 1
       cells from certain wells were treated with a drug, this barcoding round is added by parameter '-y')
    -y the index defining the treatment (index for the line within barcodefile that contains the 'treatment' barcodes 
       (mostly this is a CI-round barcode where treatment was performed in this round) (0-indexed))

    -n list of names for the guide barcodes (e.g. names of the cell lines/ knock outs). Must be in same order as the guide barcodes.

    #NON MANDATORY PARAMTERS
    -m mismatches per patter (default = 1)
    -t threads (default = 1)
    -f threshold to use when filtering umis (default is zero). Keep a read of UMI only, if more than X-perc (as double) are from same single cell and AB

    \n");
    exit();
}

//Basic Parameter Error Handling: more detailled errors are hadnled in the tools itself
if(!$options['p'])
{
    fwrite(STDERR, "No parameter for the barcode mapping pattern given: -p!\n");
    exit(1);
}
if(!$options['b'])
{
    fwrite(STDERR, "No barcode file with the possible barcodes are given: -b!\n");
    exit(1);
}
if($options['g']!=NULL && $options['x']==NULL)
{
    fwrite(STDERR, "If (also) guide barcodes ar mapped, provide also the index in the pattern, where the guide barcodes can occur. -x is missing!\n");
    exit(1);
}
if($options['a']!=NULL && $options['x']==NULL)
{
    fwrite(STDERR, "If you want to map AB barcodes to AB names, provide also the idnex of the AB: -x missing!\n");
    exit(1);
}
if($options['d']!=NULL && $options['y']==NULL)
{
    fwrite(STDERR, "If you want to map a certain barcoding round to treatments, provide also the index of the CI round, where treatment was performed: -y missing!\n");
    exit(1);
}

//parse the mismatch option to get the allowed UMI mismatches
if($options['m']!=NULL)
{
    $patternString = $options['p'];
    //check if we have a UMI barcode
    if(strpos($patternString, 'X')!==false || strpos($patternString, 'x')!==false)
    {
        $index = NULL;
        //get the index of UMI mismatch
        $patterns = explode("]", $patternString);
        $count = 0;
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
        $allowedMM = explode(",",$options['m']);
        $options['u'] = $allowedMM[$index];
    }
}

if(strpos($patternString, 'R')!==false || strpos($patternString, 'r')!==false)
{
    $mapRNA = true;
}

if($options['r']!=NULL)
{
    $pairedEnd = TRUE;
}

//create output folder
if(!$options['o'])
{
    $options['o'] = "./bin/CIResult";
}
if (!is_dir($options['o']))
{
    mkdir($options['o'], 0777, true);
}

//log file of analysis
$logfile = $options['o'] . "/analysis.log";
@unlink($logfile, ); //delete old file if exists
$logfileHandle = fopen($logfile, "a+") or die("Unable to open file!");

//create temporary folder for analysis (storing e.g. combination of barcode file and guide file (as needed by demultipelxing
//whereas processing needs them seperately))
$tmpFolder = $options['o'] . "/tmp";
mkdir($tmpFolder, 0777, true);

######################################
//execute demultiplexing
######################################

//add guide barcodes to barcodeFile
$guideBarcodeFile = $tmpFolder . "/barcodeFile.tsv";
if($options['g']!=NULL)
{    
    fwrite($logfileHandle, "Adding the guide barcodes to the barcode file (attach barcodes to the line of AB barcodes).\n");
    //get the guides that have to be added
    $f = fopen($options['g'], 'r');
    $line = fgets($f);
    $guides = explode(',', $line);
    fclose($f);
    $guideBarcodeString = "," . implode(",", $guides);

    //get content of barcode file
    $barcodeFile = $options['b'];
    $lines = file($barcodeFile);

    assert(count($lines) >= $options['x']);
    
    //add the guides to file
    $lines[$options['x']] = str_replace(array("\n", "\r"), '', $lines[$options['x']]) . $guideBarcodeString;    
    file_put_contents($guideBarcodeFile, implode('', $lines));
    fwrite($logfileHandle, "-> DONE\n");
}

$command = "./bin/demultiplexing";
$parameters = ['i', 'r', 'p', 'm', 't'];
foreach($parameters as $element)
{
    if($options[$element]!=NULL)
    {
        $command = $command . " -" . $element . " " . $options[$element];
    }
}
//add barcode file with additional guide barcodes
$command = $command . " -b " . $guideBarcodeFile;
//add output folder
$command = $command . " -o " . $options['o'] . "/Demultiplexing.tsv";
fwrite($logfileHandle, "#RUNNING COMMAND: ".$command."\n");
run_updates($command);
fwrite($logfileHandle, "-> DONE\n");

//gzip the output
$outputDemultiplexing = $options['o'] . "/Demultiplexed_Demultiplexing.tsv";
shell_exec("gzip -f " . $outputDemultiplexing);

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

$command = "./bin/processing ";
$parameters = ['b', 'c', 'u', 't', 'a', 'x', 'd', 'y', 'g', 'n', 'f'];
foreach($parameters as $element)
{
    if($options[$element]!=NULL)
    {
        $command = $command . " -" . $element . " " . $options[$element];
    }
}
//set new input file = output of Demultiplexing and output file
$command = $command . " -i " . $outputDemultiplexing . ".gz";
$command = $command . " -o " . $options['o'] . "/Processing.tsv";

//TODO: add output file manually
fwrite($logfileHandle, "#RUNNING COMMAND: ".$command."\n");
run_updates($command);
fwrite($logfileHandle, "-> DONE\n");

######################################
//execute some Quality Control scripts / Data Analysis scripts
######################################



######################################
// clean analysiss
######################################
clean($tmpFolder, $logfileHandle);

?>