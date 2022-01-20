<?php
/**
Script to run CI data analysis:
This can contain barcodes for ABs, sgRNA, etc.
However, RNA sequencesw are not mapped to a reference genome,
therefore use analyzeRNA.php

*/

function run($string) 
{
    $logfile = "./bin/analysis.log";
    $fileHandle = fopen($logfile, "a+") or die("Unable to open file!");

    //keep output of command (stderr is written to terminal)
    $output = shell_exec($string);
    fwrite($fileHandle, $output);

    fclose($fileHandle);
}

//parse command line arguments
$shortopts = "i:"; //input file
$shortopts .= "o:"; //output file

$shortopts .= "m:"; //mismatch string: e.g.=1,2,3
$shortopts .= "p:"; //mapping pattern: e.g.=[GAGCTCTGCGACG][XXXXXXXXXXXX][GATGACTCA]
$shortopts .= "b:"; //barcode file

$options = getopt($shortopts);

$logfile = "./bin/analysis.log";
@unlink($logfile, ); //delete old file if exists

//execute mapping of barcodes
run("./bin/demultiplexing -i " . $options['i'] . " -o " . $options['o']);

//if CDNA map reads to reference genome (if mapping pattern contains CDNA region)

//execute processing of the mapped barcodes (UMI collapsing, assigning single cells, etc.)

?>