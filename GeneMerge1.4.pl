#!/usr/bin/perl -w

use POSIX;

#    GeneMerge v1.4 - 2015

#    GeneMerge-- Post-genomic analysis, data mining, and hypothesis testing
#    Castillo-Davis, C.I. and D.L. Hartl 2003. Bioinformatics 19(7):891-892

#    Copyright 2003, 2015 (C) Cristian I. Castillo-Davis

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

########################################################################

# This program returns P-values for functional enrichment given a set
# of study items (genes) using the hypergeometric distribution, as well 
# functional descriptions for each item (gene). Output is tab-delimited 
# text and HTML.

# Bonferroni corrected P-values are given as well as False Discovery Rate 
# cut-offs for FDR = 10% 5% and 1%. Threshold P-values are reported.
# A custom FDR is optional.

# Please see http://www.genemerge.net for full documentation.

# To run this script use the following syntax:
# ./GeneMerge1.4.pl gene-association  description  population  study  outfile

# Input on command line is as follows:

# file1 = a gene-association file (gene \t ID)
# file2 = a human readable description file for gene-association IDs (ID \t description) 
# file3 = a population set of genes (the pool from which the study set
#         is drawn, often a genome or all genes on a microarray
# file4 = a study set of genes (for example, genes found up/down regulated)
# file5 = an output filename

# X.X% = optional custom False Discovery Rate  
#        Note that P-value cut-offs for FDR = 1%, 5%, and 10% are automatically
#        calculated and provided in the ouput. Default is 0.5% if none entered.
#
#       n = # population genes
#       k = # study set genes
#       p = # of genes with a particular GMRG term divided by total genes (n)
#       r = # of genes in study set with a particular term

$version = "1.4";
$year = "2015";

# intialize
$GMRGs_line_by_line = "";
$updownGMRGs_line_by_line = "";
$genes_GMRGterms = "";
$up_down_genes_with_annotation = 0;
$genes_not_found_in_ontology = "";
$html_genes_not_found_in_ontology = ""; # 1.3
$GenesNoInfo = "";
$no_of_up_down_GMRG_terms = 0;
$BCr = 0;
$abort = "FALSE";  
$abort_FDR = "TRUE";  # 1.3
$missing_genes = 0;   # 1.4

$FDR_1_perc_threshold = "NA";      # report the P-value cut-off        # 1.4
$FDR_5_perc_threshold = "NA";      # for each FDR percentage in output # 1.4 
$FDR_10_perc_threshold = "NA";                                         # 1.4
$FDR_custom_perc_threshold = "NA";                                     # 1.4


########################################### 
# Set input files and arguments           #
# Error handling for missing arguments    # # 1.3
# format: genename \t GMRG id             #
$gene_association_file = $ARGV[0] || die "\nOne or more required files or parameters was not specified.\nThe input syntax is:\n\n./GeneMerge$version.pl  gene-association  description  population  study  output\n\nType  ./GeneMerge$version\.pl -h  to display all options.\n\n";

#############################################
# Usage help                                # # 1.2  # 1.3 more options added  
############################################# 
if (($ARGV[0] eq "h") || ($ARGV[0] eq "-h")|| ($ARGV[0] eq "\-\-help") || ($ARGV[0] eq "-help") || ($ARGV[0] eq "-Help") || ($ARGV[0] eq "--Help") || ($ARGV[0] eq "help")) {
    print "\nThe input syntax is:\n\n./GeneMerge$version.pl  association.file  description.file  population.file  study.file  output.filename  FDR%\n\n";
    print "Specifying an FDR percentage is optional. Remember to include the \% symbol after the FDR.\n\n"; # 1.4
    $abort = "TRUE";
    $ARGV[0] = "";
    exit;
}

################################################### 
# Version information                             # # 1.3
###################################################
if (($ARGV[0] eq "v") || ($ARGV[0] eq "-v") || ($ARGV[0] eq "V") || ($ARGV[0] eq "--version") || ($ARGV[0] eq "-Version") || ($ARGV[0] eq "--Version") || ($ARGV[0] eq "version") || ($ARGV[0] eq "-version") || ($ARGV[0] eq "-V")){
    print "\nGeneMerge Version $version\n\n";
    $abort = "TRUE";
    exit;
}

################################################### 
# Input file error checking                       # # 1.3  # 1.4
###################################################
# format: GMRG id \t English description
$description_file = $ARGV[1]  || die "\nOne or more required files or parameters was not specified.\nThe input syntax is:\n\n./GeneMerge$version.pl  gene-association.file  description.file  population.file  study.file  output.filename\n\nType  ./GeneMerge$version\.pl -h  to display all options.\n\n";

# genenames of "population set"
$population_genes_file = $ARGV[2]  || die "\nOne or more required files or parameters was not specified.\nThe input syntax is:\n\n./GeneMerge$version.pl  gene-association.file  description.file  population.file  study.file  output.filename\n\nType  ./GeneMerge$version\.pl -h  to display all options.\n\n";

# genenames of "study set"
$study_genes_file = $ARGV[3]  || die "\nOne or more required files or parameters was not specified.\nThe input syntax is:\n\n./GeneMerge$version.pl  gene-association.file  description.file  population.file  study.file  output.filename\n\nType  ./GeneMerge$version\.pl -h  to display all options.\n\n";

# new filename for results file
$results_filename = $ARGV[4]  || die "\nOne or more required files or parameters was not specified.\nThe input syntax is:\n\n./GeneMerge$version.pl  gene-association.file  description.file  population.file  study.file  output.filename\n\nType  ./GeneMerge$version\.pl -h  to display all options.\n\n";

# a check on missing output filename - we don't want GeneMerge to think the FDR% is the output filename! # 1.4   
if ($results_filename =~ /\%/) { # if we see a percent symbol in the output filename there is a problem
    die "\nPlease specify an output filename before the FDR\%.\n\nThe input syntax is:\n./GeneMerge$version.pl  gene-association.file  description.file  population.file  study.file  output.filename  FDR%\n\n";
}

# user specified false discovery rate written as a percentage, e.g. 15%      # 1.3
$custom_FDR = $ARGV[5];

###################################################################################################
# End input arguments  ARGV[X]
###################################################################################################

##########################################  
# set FDR to 0.5% if it is not specified # # 1.3
##########################################
if ($ARGV[5]) {
  #do nothing
} else {
  # default to 0.5 % FDR 
  $custom_FDR = 0.005; # this value does not get divided by 100 like custom_FDR from STDIN
                       # that's why it's not '0.5' like you'd expect
  #create default FDR header for regular text output
  $custom_FDR_header = "FDR_0.5_perc"; 
  
  #create default FDR header for HTML output
  $custom_FDR_percentage = 0.5;  
  
  $abort_FDR = "FALSE";
}

###################################################   
# Parse custom FDR value, check and make headers  # # 1.3
###################################################

if ($ARGV[5]) {
  # get custom FDR value, remove the % symbol
  # and check that it is strictly numbers, not characters, for example "e-10" is not allowed
  if ($custom_FDR =~ /\D\%/) { # if not a digit
    die "\nThe FDR must be numeric, for example, 1.5%\n\n"; # 1.4
  }  

  # this is the correct way to enter the FDR, an integer or decimal value followed by "%"
  if ($custom_FDR =~ /(\d*\.\d*)\%/ || $custom_FDR =~ /^(\d*)\%/) {
    $custom_FDR = $1;
  
    #create custom FDR header for regular text output
    $custom_FDR_header = "FDR" . "_" . "$custom_FDR" . "_perc"; 
  
    #create custom FDR header for HTML output
    $custom_FDR_percentage = $1;  
    
    #check on an empty FDR value but with % symbol alone
    if ($custom_FDR eq "") {
      die "\nPlease specify a number for the FDR percentage, for example, 15%\n\n"; # 1.4
    }

    #get FDR fraction (not percentage) for use in actual calculations
    $custom_FDR = $custom_FDR/100;

    $abort_FDR = "FALSE";

  } else {
    die "\nThe FDR must be written as a numeric percentage. For example 0.05% (not 5e-2%)\n\n"; # 1.4
  }
} else {
  $abort_FDR = "TRUE";
}

#########################################  
# Check to make sure FDR is within      # # 1.3
# reasonable bounds                     # # 1.3
#########################################

if ($abort_FDR eq "FALSE") {
  if (($custom_FDR >= 1.0) || ($custom_FDR <= 0)) {
    unlink("study_genes.temp");           
    unlink("pop_genes.temp");             
    unlink("gene_association_file.temp"); 
    unlink("description_file.temp");
    die "\nThe FDR must be greater than 0% and less than 100%.\n\n";
  }
}
### leap stuff removed here #############      # 1.4

if ($abort eq "FALSE") {
    
    ########################################   
    # Clean up line endings in input files #   # 1.2
    # so that they are unix style.         #
    ########################################
    
    # change line breaks in each file to Unix line breaks
    # fix the population file
    open(NFF,">pop_genes.temp")  || die "couldn't create pop_genes.temp, check system permissions\n\n";
    open(POP,"$population_genes_file") || die "couldn't open $population_genes_file. Perhaps it's misspelled?\n\n";
    while(<POP>) {
	s/(\r\n?)+/\n/g;
	s/^\n//;
	print NFF $_;
    }
    close(POP);
    close(NFF);
    
    # fix the study file
    open(NFF2,">study_genes.temp") || die "couldn't create study_genes_temp, check permissions on your system\n\n";
    open(POP2,"$study_genes_file") || die "couldn't open $study_genes_file. Perhaps it's misspelled?\n\n";
    while(<POP2>) {
	s/(\r\n?)+/\n/g;
	s/^\n//;
	print NFF2 $_;
    }
    close(POP2);
    close(NFF2);
  
    # fix the gene association file
    open(NFF,">gene_association_file.temp") || die "couldn't create gene_association_file.temp, check permissions on your system\n\n";
    open(POP,"$gene_association_file") || die "couldn't open $gene_association_file. Perhaps it's mispelled?\n\n";
    while(<POP>) {
	s/(\r\n?)+/\n/g;
	s/^\n//;
	print NFF $_;
    }
    close(POP);
    close(NFF);
    
    # fix the description file
    open(NFF,">description_file.temp") || die "couldn't create description_file.temp, check persmissions on your system\n\n";
    open(POP,"$description_file") || die "couldn't open $description_file. Perhaps it's misspelled?\n\n";
    while(<POP>) {
	
	s/(\r\n?)+/\n/g;
	s/^\n//;
	print NFF $_;
    }
    close(POP);
    close(NFF);
    
    
################################################    
# Check structure and logic of each input file #    # 1.3  # 1.4
################################################

    # check study and pop for stray whitespace 
    # check ga and description file for proper formatting (tabs, semi-colons etc.)
    # check that all genes in study are in pop (study should be a proper subset of pop)

    ################################### 
    # check population file           #  # 1.4
    ###################################
    open(POPCHECK, "pop_genes.temp");                
    while(<POPCHECK>) {
	if (/^\s/) {                  # if the line begins with whitespace # 1.4
	    unlink("study_genes.temp");           
	    unlink("pop_genes.temp");       
	    unlink("gene_association_file.temp");
	    unlink("description_file.temp"); 
	    die "\nThe population file appears to be incorrectly formatted. It has stray white-space in it. Please check.\nError: line $. in file \'$population_genes_file\'\n\n";#print line error # 1.4
	}
    }
    close(POPCHECK);
    

    ###################################### 
    # check study file                   # # 1.4
    ######################################
    open(STUDCHECK, "study_genes.temp");             
    while(<STUDCHECK>) {
	if (/^\s/) {                  # if the line begins with whitespace # 1.4
	    unlink("study_genes.temp");           
	    unlink("pop_genes.temp");       
	    unlink("gene_association_file.temp");
	    unlink("description_file.temp");
	    die "\nThe study file appears to be incorrectly formatted. It has stray white-space in it. Please check.\n\nError: line $. in file \'$study_genes_file\'\n\n"; #print line error # 1.4
	}
    }
    close(STUDCHECK);

    ###################################### 
    # check gene association file format # # 1.4
    ######################################
    open(GACHECK, "gene_association_file.temp");     
    while(<GACHECK>) {
	if (/^(.*)\t(.*)\;$/) {   # we expect "geneid \tab term1;term2; " <--- note ending semi-colon
	    # everything's cool     
	} else {
	    unlink("study_genes.temp");           
	    unlink("pop_genes.temp");
	    unlink("gene_association_file.temp");
	    unlink("description_file.temp");
	    die "\nThe gene-association file appears to be incorrectly formatted. Please check. Maybe you left out the semi-colons?\n\nError: line $. in file \'$gene_association_file\'\n\n"; #print line error # 1.4
	}
    }
    close(GACHECK);
    
    ###################################### 
    # check description file format      #  # 1.4
    ######################################
    open(DESCHECK, "description_file.temp");          
    while(<DESCHECK>) {
	if (/^(.*)\t(.*)$/) { # we expect "term \tab description"
	    # everything's cool
	} else {
	    unlink("study_genes.temp");           
	    unlink("pop_genes.temp");             
	    unlink("gene_association_file.temp"); 
	    unlink("description_file.temp");
	    die "\nThe description file appears to be incorrectly formatted. Please check.\n\nError: line $. in file \'$description_file\'\n\n"; #print line error # 1.4
	}
    }
    close(DESCHECK);
    
    
    #########################################  
    # Check to see if all elements in study # # 1.3
    # file are present in population file   # # 1.3
    #########################################

    open(POPHASH, "pop_genes.temp");
    while(<POPHASH>) {
	$pophash{$_} = $_;
    }
    # make empty line entries in pophash-- in case the study file has an empty line or two
    $pophash{" \n"} = " \n";             # allows for multiple blank lines in study file    # 1.4
    $pophash{"\n"} = "\n";               # including lines with just one space and newline  # 1.4
                                         # blank lines in both study and pop are ignored by GeneMerge 
    open(STUDYHASH, "study_genes.temp");
    while(<STUDYHASH>) {
	
      if ($pophash{$_}) {  # if we find this study gene in the pop hash (it should be there)
	  #all's well
      } else {
	  chomp($_);
	  $example_no_match = $_;
	  $missing_genes++;         
	  if ($missing_genes > 0) { # allow for zero mismatches btw study and population files  
	      unlink("study_genes.temp");           
	      unlink("pop_genes.temp");             
	      unlink("gene_association_file.temp"); 
	      unlink("description_file.temp");
	      die "\nThere are elements in \'$study_genes_file\' that are not in \'$population_genes_file\'\. Please check for stray white-space.\n\nOr perhaps you listed the files in reverse, i.e. \'study\' \'pop\' instead of \'pop\' \'study\'?\n\nGene \'$example_no_match\' was found in \'$study_genes_file\' but not in \'$population_genes_file\'\n\n";
	  }
      }  
    }
    close(STUDYHASH);
    close(POPHASH);
    
    ############################################################################################################################################
    # Main                                                                                                                                     #
    ############################################################################################################################################
    
    #######################################################################
    # Hash GMRG Annotation data                                             #    
    #######################################################################
    open(TARGET, "gene_association_file.temp") || die "couldn't open $gene_association_file";
    while(<TARGET>) {           # while new line, for each line evaluate
	
	if (/(.*)\t(.*)\n/) {
	    $genename = $1;
	    $annotation = $2;
	    #print "$genename\n";
	    chomp($annotation);
	    
	    # assign value genename to the key  
	    $GMRGgenomehash{$genename} = $annotation;
	}
    }
    print "Done parsing $gene_association_file\n";
    close(TARGET);

    print "Collecting terms associated with population genes...\n";
    
    ##########################################################################
    # Retrieve GMRG terms for each gene in detected list and put into an array
    # so that we can print them line by line to a new file.
    ##########################################################################
    #print "Detected genes are:\n";
    open(SMALLFILE, "pop_genes.temp") || die "couldn't open pop_genes.temp";
    while(<SMALLFILE>) {          
	
	$gene = $_;
	chomp($gene);
	#print "$gene\n";
	
	#count number of detected genes 
	$total_no_detected_genes++;
	
	# do a lookup on the hash to get GMRGIDs associated with each gene we have
	# use "if" to see if gene was there... we keep track of these genes later
	if ($bigline = $GMRGgenomehash{$gene}) {
	    
	    #split $bigline using ; as a delimiter to put GMRGIDs into array
	    @indivGMRGs = split(/;/,$bigline);  
	    #get length of array
	    $length_of_indivGMRGs_array = $#indivGMRGs;
	    
	    for($v=0;$v<=$length_of_indivGMRGs_array;$v++) {
		
		$GMRGs_line_by_line = $GMRGs_line_by_line . "$indivGMRGs[$v]\n";
	    }
	} 
    }
    print "Done.\n"; 
    print "Number of population genes: $total_no_detected_genes\n";
    close(SMALLFILE);
    
    #####################################################################################
    # Print every GMRG term among detected genes to a file to create array-wide pool    #  
    # of GMRG terms. Then we can sort -u the _study pool_ and count freq of only those  #  # 1.4
    # GMRG terms that are found in the study set.                                       #  # 1.4
    #####################################################################################
    
    #write GMRGIDs on indiv lines to file "SamplePoolGMRGIDs"
    open(NF, ">SamplePoolGMRGIDs") || die "couldn't create SamplePoolGMRGIDs, check system permissions";
    #print "$GMRGs_line_by_line";
    print NF $GMRGs_line_by_line;
    close(NF);
    # removed sortu SamplePool, since we now only count freq of study terms in the pop pool # 1.4

    ##########################################################################
    # Create a unique list of GMRG IDs from the study set here so that we can       # 1.4           
    # 1) query the population to get the pop fraction                               # 1.4
    # 2) query the study to get the study fraction                                  # 1.4 
    # (old way was to get freq counts for ALL pop terms-- not needed)               # 1.4 

    print "Collecting terms associated with study genes...\n"; # 1.4

    # begin moved study section ############################################################## # 1.4
    open(UPDOWNFILE, "study_genes.temp")  || die "couldn't open study_genes.temp";
    while(<UPDOWNFILE>) {          
	
	$total_no_updown_genes++;
	
	$updowngene = $_;
	chomp($updowngene);
	#print "$updowngene\n";
	
	# do a lookup on the hash to get GMRGIDs associated with each gene we have
	# use "if" to check if gene was found at all in ontology
	if ($fullline = $GMRGgenomehash{$updowngene}) {
	    
	    
	    #split $fullline using ; as a delimiter to put GMRGIDs into array
	    @updownGMRGs = split(/;/,$fullline);  
	    
	    #get length of array
	    $length_of_updownGMRGs_array = $#updownGMRGs;
	    
	    for($b=0;$b<=$length_of_updownGMRGs_array;$b++) {
		$updownGMRGs_line_by_line = $updownGMRGs_line_by_line . "$updownGMRGs[$b]\n";
	    }	    
	    
	    #print GMRGIDs in array as is to save to a file later   
	    $genes_GMRGterms = $genes_GMRGterms . "$updowngene\t$fullline\n";   
	} else { 
	    # if not save this information
	    $genes_not_found_in_ontology = $genes_not_found_in_ontology . "$updowngene\t";
	    $html_genes_not_found_in_ontology = $html_genes_not_found_in_ontology . "$updowngene<br>"; # 1.3
	    $GenesNoInfo++;
	}
    }  
    close (UPDOWNFILE);

    # write GMRGIDs on indiv lines to file "UpDownPoolGMRGIDs"
    open(NF, ">UpDownPoolGMRGIDs")  || die "couldn't create file UpDownPoolGMRGIDs, check system permissions";
    #print "$updownGMRGs_line_by_line";
    print NF $updownGMRGs_line_by_line;
    close(NF);

    ## use sortu subroutine  to get a list of all unique GMRGIDs in study set file
    $output2 = &sortu("UpDownPoolGMRGIDs");
    
    ## print unique IDs to a file
    open(SORTED2, ">uniqueUpDownGMRGIDs") || die "couldn't create uniqueUpDownGMRGIDs, check system permissions";
    print SORTED2 $output2;
    close(SORTED2);
    # end moved study section ############################################################## # 1.4

  
    #########################################################################################
    # Determine GMRGID counts among population genes                                        #    
    #                                                                                       #
    # Count the frequency of each _study_ GMRGID in "Sample_GMRGIDs" (this is equivalent to # # 1.4
    # counting the number of genes with a particular GMRG since each gene can have a GMRG   #
    # associated only once)                                                                 #
    #########################################################################################

    print "Done.\nCalculating population frequency for each term\n"; # 1.4

    # get the pop fraction (frequency) for each _study set_ term                            # 1.4
    open(UNIQUEIDS, "uniqueUpDownGMRGIDs") || die "couldn't open uniqueUpDownGMRGIDs";      # 1.4 we use study terms to search pop
    while(<UNIQUEIDS>) {
	$count = 0;    
	$uniqueGMRGID = $_;
	chomp($uniqueGMRGID);
	
	# check each unique GMRG term against each line in the SampleGMRGID file to get freq counts
	# for each term
	open(SAMPLE,"SamplePoolGMRGIDs")  || die "couldn't open SamplePoolGMRGIDs";
	while(<SAMPLE>) {
	    if(/\b$uniqueGMRGID\b/) {   # 1.1
		$count++;
	    }
	}
	close(SAMPLE);
	
	# get frequency of each GMRGID in the sample pool and write to a hash table
	if($uniqueGMRGID =~ /\S+/) {  # a check on blank lines
	    $frequency_of_GMRGID = $count/$total_no_detected_genes;  
	    # 1.3 sprintf trim pop significant digits # 1.4 removed trim from here - can cause inaccurate pop freqs! # 1.4

	    # assign freq value to GMRGID key
	    $GMRGfreq_hash{$uniqueGMRGID} = $frequency_of_GMRGID;  
	    
	    # assign count to GMRGID key (count of GMRGID on array)  # 1.2
	    $GMRGcount_hash{$uniqueGMRGID} = $count;                 # 1.2
	    
	    # 1.4 trim sig digs *for screen output* only
	    $trim_frequency_of_GMRGID = sprintf("%.5g", $frequency_of_GMRGID);     # 1.4
	    print "$uniqueGMRGID = $count\t freq = $trim_frequency_of_GMRGID\n";   # 1.4	
	}   
    }
    close(UNIQUEIDS);
    
    # write all updown Gene Names and GMRGIDs to a file 
    open(NF, ">gene_GMRGterms")  || die "couldn't create gene_GMRGterms, check system permissions";
    print NF $genes_GMRGterms;
    close(NF);
    
    print "Done.\n"; # 1.4
    print "Number of study set genes: $total_no_updown_genes\n";

    ###################################################################################
    ###################################################################################
    # Determine GMRGID counts in up/down regulated genes                              #
    ###################################################################################
    
    # MOVED printfile UpDownPoolGRMIDs, sortu UpDownPoolGMRIDs and study set term    # 1.4
    # collection from here to above so we can search population set for only those   # 1.4
    # terms (IDs) that appear in the study set                                       # 1.4
    

    print "Calculating study set frequency for each term\n";                         # 1.4

    # Count the frequency of each GMRGID in "UpDownPoolGMRGIDs"                           
    # open up file with unique GMRGIDs present in the sample
    open(UNIQUEUPDOWNIDS, "uniqueUpDownGMRGIDs") || die "couldn't open uniqueUpDownGMRGIDs";
    
    while(<UNIQUEUPDOWNIDS>) {
	$counts = 0;    
	$unique_up_down_GMRGID = $_;
	chomp($unique_up_down_GMRGID);
	
	#check each one against each line in the UpDownPoolGMRGIDs file to get counts
	open(POOL,"UpDownPoolGMRGIDs")  || die "couldn't open UpDownPoolGMRGIDs";
	while(<POOL>) {
	    if(/\b$unique_up_down_GMRGID\b/) { # 1.1
		$counts++;
	    }
	}
	close(POOL);
	
	# get frequency of each GMRGID in the sample pool and write to a hash table
	if($unique_up_down_GMRGID =~ /\S+/) {  # a check on GMRGIDs that are blank 
	    # assign count value to GMRGID count hash key
	    $updownGMRGcount_hash{$unique_up_down_GMRGID} = $counts;
	    
	    print "$unique_up_down_GMRGID\t$counts\n";
	    
	    # count the number of unique GMRG terms in up/down genes
	    $no_of_up_down_GMRG_terms++;
	    
	    # count the number of GMRG terms that are not represented only once in the population
	    # these "population singletons" will be in frequency 1/total detected genes
	    # we use the number of terms that are not pop-singletons for the Bonferroni correction later
	    if ($GMRGfreq_hash{$unique_up_down_GMRGID} > (1/$total_no_detected_genes)) {
		$BCr++;
	    }
	}   
    }
    close(UNIQUEUPDOWNIDS);
    
    # Write contents of GMRGID count within up/down and freq w/in detected to a file
    open(NF, ">GMRGID_n_p_k_r") || die "couldn't create GMRGID_n_p_k_r, check system permissions";
    
    while(($key,$value) = each(%updownGMRGcount_hash)) {                                         
	# grab freq of the GMRGID in the detected sample from other hash
	$freq = $GMRGfreq_hash{$key};
	
	
	#         GMRGID  n                        p      k                     r(count)
	print NF "$key\t$total_no_detected_genes\t$freq\t$total_no_updown_genes\t$value\n";
	
    }
    close(NF);
    
    print "Done.\n"; # 1.4
    print "Number of study set GMRGterms is = $no_of_up_down_GMRG_terms\nBCr = $BCr\n";
    
    ################################################################
    ################################################################
    # Part III - Get P-values for each GMRGID among up/down genes  #
    #            using the hypergeometric distribution.            #
    ################################################################
    print "Calculating P-value for each term\n"; # 1.4
    
    open(NF, ">GMRGID_P_Values") || die "couldn't create GMRGID_P_values, check system permissions";
    
    open(DATA,"GMRGID_n_p_k_r") || die "couldn't open GMRGID_n_p_k_r";
    while(<DATA>) {
	
	$p_value = 0;
	if (/(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/) {
	    
	    #	print "$1, $2, $3, $4, $5\n";
	    
	    my $GMRGid = $1; #
	    my $n = $2; #
	    my $p = $3; #
	    my $k = $4; #
	    my $r = $5; #
	    
	    # if not a singleton in up/down list
	    if ($r != 1) {
		
		# for r values from r (observed) through k (most extreme case) 
		# add together probabilities
		
		# add up starting with most extreme case first to 
		# hold onto significant figures
		
		$p_value = &hypergeometric($n,$p,$k,$r);
		
	
		# Bonferroni correction for multiple tests
		$p_value_corrected = ($p_value * $BCr);
		if ($p_value_corrected >= 1.0) {
		    $p_value_corrected = 1.0;
		}
		print "pvalue is: $p_value\n"; 
		print "Bonferroni corrected pvalue is: $p_value_corrected\n";
	    } else {
		
		$p_value = "NA";
		$p_value_corrected = "NA";
		
		print "pvalue is: $p_value\n";
		print "Bonferroni corrected pvalue is: $p_value_corrected\n";
	    }
	    
	    # grab gene names associated with each GMRG term
	    # from the up/down list
	    $gene_simmons = "";
	    open(GENE, "gene_GMRGterms") || die "couldn't open gene_GMRGterms";
	    while(<GENE>) {
		
		
		# parse line into gene names and GMRG terms
		if (/(.*)\t(.*)/) {
		    $GeneName = $1;
		    
		}
		# if the P-Value GMRGid matches a GMRGid for a particular gene
		# grab it and add it to the list
		if (/\b$GMRGid\b/) { # 1.1
		    $gene_simmons = $gene_simmons . "$GeneName\t"; 
		}
	    }
	    close(GENE);	
	    
	    #	print "$GMRGid $n, $p, $k, $r, prob = $p_value\n";
	    print NF "$GMRGid\t$n\t$p\t$k\t$r\t$p_value\t$p_value_corrected\t$gene_simmons\n";
	}
    }
    close(NF);
    close(DATA);
    
    # now use pre-parsed GMRG ontology table to look up functions for
    # all up/down regulated genes
    open(ONTOLOGY, "description_file.temp") || die "couldn't open description_file.temp";
    while(<ONTOLOGY>) {
	#parse each line to put function into a hash and GMRGterm as key
	if (/(.*)\t(.*)\n/) {
	    $GMRGkey = $1;
	    $function = $2;
	    $Ontology_Hash{$GMRGkey} = $function;
	    #print "        $GMRGkey\t$function\n";
	}
    }
    close(ONTOLOGY);
    
    # open this file and count how may up/down genes have no GMRGterms
    # associated with them
    
    open(TEMP, "gene_GMRGterms") || die "couldn't open gene_GMRGterms";
    while(<TEMP>) {
	if (/(.*)\t(.*)/) {
	    # do nothing
	    $up_down_genes_with_annotation++;
	}
    }
    close(TEMP);
    
    ##################################################################
    # Parse final data - look up ontologies and generate output      #
    ##################################################################
    
    open(NF,">almost_final") || die "couldn't create file almost_final, check system permsissions";
    open(TEMPII, "GMRGID_P_Values") || die "couldn't open GMRG_P_Values";
    while(<TEMPII>) {
	
	#intitialize
	$genenames = "";
	
	# parse each line
	@FinalData = split(/\t/,$_);
	
	# get the length of the current array
	$arraylength = $#FinalData;
	
	# parse out relevant data for writing to final file and function lookup
	$GMRG_Term = $FinalData[0];
	#$allGMRGs_detected = $FinalData[1];
	$allupdownGMRGs = $FinalData[3];
	$r_in_updown = $FinalData[4];
	
	$pValue = $FinalData[5];
	$pValue_corrected = $FinalData[6];

	# sprintf P-values # 1.3  # 1.4 removed - sprintf fails when P is very small     # 1.4
        # (P < ~1e-300) because sprintf on a value read from a file is not at
        # double-precision and this can lead to weird sprintf values

	# Go through the array and parse out all gene names that contribute
	# to this GMRG term (could be one or more)
	for ($x=7;$x<$arraylength;$x++) {
         $genenames = $genenames . "/$FinalData[$x]/, ";
	    #$genenames = $genenames . "$FinalData[$x]\t"; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Modificado por el comando de arriba
	}        
	# get function for this term from Ontology hash
	if ($this_function = $Ontology_Hash{$GMRG_Term}) {
	    #everything is cool
	} else {
	    $this_function = "couldn't find description for this term in $description_file";    # 1.2
	}
	
	# Get array count for this term                          # 1.2
	if ($array_count = $GMRGcount_hash{$GMRG_Term}) {        # 1.2
	    #everything is cool
	} else {
	    $array_count = "couldn't find array count";
	}
	
	#print output yeah!
	#print NF "$GMRG_Term\t$GMRGfreq_hash{$GMRG_Term}\t$array_count/$total_no_detected_genes\t$r_in_updown\/$allupdownGMRGs\t$pValue\t$pValue_corrected\t$this_function\t$genenames\n";  %%%%%%%%%%%%%%%%%%%% modificado por el comando de abajo
     print NF "\"$GMRG_Term\t$GMRGfreq_hash{$GMRG_Term}\t$array_count/$total_no_detected_genes\t$r_in_updown/$allupdownGMRGs\t$pValue\t$pValue_corrected\t$this_function\t$genenames\"\n"; 
	
    }
    close(TEMPII);
    close(NF);
    
    
    ########  Sort by P-value ########################################################     # 1.3
    
    $col = 5; # sort by P-value column  
    
    open(SORTING,"almost_final") || die "couldn't open file almost_final";
    @data = <SORTING>;
    @sorted = map { $_->[0] } sort sort_column map { [$_, split /\t/ ] } @data;
    close(SORTING);
    
    open(SORTED, ">sorted_almost_final") || die "couldn't create file sorted_almost_final, check system permissions";
    print SORTED @sorted;
    close(SORTED);
    
    
    ########### Calculate False Discovery Rates  ########################################   # 1.3
    
    # Here, we evaluate the ranked P-values in parray according to the given FDR
    # and save result "T" (TRUE) or "F" (FALSE) to an array called FDR_array
   
    # Note we follow Benjamani & Hochberg (1995) exactly, 
    # where it is stated, "The FIRST p-value to satisfy the 
    # constraint is p(4)... Thus we reject the [remaining] 
    # hypotheses having p-values which are less than or equal to..."
   
    # This is crucial because sometimes the ranked p-values can have ties,
    # and a simple numerical evaluation would yield inconsistent results.
    #         [fdr threshold p-value]
    # P-val <= (i/m) * q                    Numerical Eval.  Logical Flag Required
    # 0.027 <= 31/106 * 0.10 = 0.029245 ... TRUE             TRUE
    # 0.031 <= 32/106 * 0.10 = 0.030188 ... FALSE            TRUE  ^ 
    # 0.031 <= 33/106 * 0.10 = 0.031132 ... TRUE             TRUE  |
    # 0.031 <= 34/106 * 0.10 = 0.032075 ... TRUE             TRUE  | 
    # 0.031 <= 35/106 * 0.10 = 0.033018 ... TRUE*            TRUE  | *first true seen, so rest TRUE 
    # 0.046 <= 36/106 * 0.10 = 0.033962 ... FALSE            FALSE |

    $lines = 0;
    
    #### Calculate at 10% level ####
    $FDR = 0.10;
    $i=0;
    
    open(FDR, "sorted_almost_final") || die "couldn't open sorted_almost_final";
    
    # collect P-values from column 4 and put into a numbered array "parray"
    while(<FDR>) {
	@array = split(/\t/,$_);
	
	if ($array[4] =~ /\d/) {
	    $i++;
	    $parray[$i] = $array[4];
	} 
	
	#count the number of lines of P-values and NAs in file to fill in NAs later
	$lines++;
    }
    
    $first_pvalue_that_satisfies_equation = "notseen";           
    
    for ($k=$#parray; $k >= 1; $k--) {     
	if ($parray[$k] <= (($k/$#parray)*$FDR)) { 
	    $first_pvalue_that_satisfies_equation = "seen";
	    	    
	    # report the P-value threshold for FDR = 10%     # 1.4
	    # it is the "P-value threshold" at this FDR, we
	    # save it to report in the output 
	    if ($FDR_10_perc_threshold eq "NA") {
		$FDR_10_perc_threshold = ($k/$#parray)*$FDR;
		$FDR_10_perc_threshold = sprintf("%.3g", $FDR_10_perc_threshold);  # trim significant digits
	    }
	}
		
	if ($first_pvalue_that_satisfies_equation eq "seen") {
	    $FDR_array10[$k] = "T";
	} else {
	    $FDR_array10[$k] = "F";
	}
    }
    
    #### Calculate at 5% level ####
    $FDR = 0.05;
    $i=0;
    
    open(FDR, "sorted_almost_final") || die "couldn't open sorted_almost_final";
    while(<FDR>) {
	@array = split(/\t/,$_);
    
	if ($array[4] =~ /\d/) {
	    $i++;
	    $parray[$i] = $array[4];
	}
    }
    
    $first_pvalue_that_satisfies_equation = "notseen";           
    
    for ($k=$#parray; $k >= 1; $k--) {     
	if ($parray[$k] <= (($k/$#parray)*$FDR)) { 
	    $first_pvalue_that_satisfies_equation = "seen";

	    # report the P-value threshold for FDR = 5%      # 1.4
	    # it is the "P-value threshold" at this FDR, we
	    # save it to report in the output 
	    if ($FDR_5_perc_threshold eq "NA") {
		$FDR_5_perc_threshold = ($k/$#parray)*$FDR;
		$FDR_5_perc_threshold = sprintf("%.3g", $FDR_5_perc_threshold);  # trim significant digits
	    }
	}      
	
	if ($first_pvalue_that_satisfies_equation eq "seen") {
	    $FDR_array5[$k] = "T";
	} else {
	    $FDR_array5[$k] = "F";
	}
    }
    
    #### Calculate at 1% level ####
    $FDR = 0.01;
    $i=0;
    
    open(FDR, "sorted_almost_final") || die "couldn't open sorted_almost_final";
    while(<FDR>) {
	@array = split(/\t/,$_);
	
	if ($array[4] =~ /\d/) {
	    $i++;
	    $parray[$i] = $array[4];
	}
    }
    $first_pvalue_that_satisfies_equation = "notseen";           
    
    for ($k=$#parray; $k >= 1; $k--) {     
	if ($parray[$k] <= (($k/$#parray)*$FDR)) { 
	    $first_pvalue_that_satisfies_equation = "seen";
	    
	    # report the P-value threshold for FDR = 1%      # 1.4
	    # it is the P-value threshold at this FDR, save to report in output 
	    if ($FDR_1_perc_threshold eq "NA") {
		$FDR_1_perc_threshold = ($k/$#parray)*$FDR;
		$FDR_1_perc_threshold = sprintf("%.3g", $FDR_1_perc_threshold);  # trim significant digits
	    }
	}
	
	if ($first_pvalue_that_satisfies_equation eq "seen") {
	    $FDR_array1[$k] = "T";
	} else {
	    $FDR_array1[$k] = "F";
	}
    }
    
    
    #### Calculate at Custom % level ####
    $FDR = $custom_FDR;
    
    $i=0;
    
    open(FDR, "sorted_almost_final") || die "couldn't open sorted_almost_final";
    while(<FDR>) {
	@array = split(/\t/,$_);
	
	if ($array[4] =~ /\d/) {
	    $i++;
	    $parray[$i] = $array[4];
	}
    }

    $first_pvalue_that_satisfies_equation = "notseen";           
    
    for ($k=$#parray; $k >= 1; $k--) {     
	if ($parray[$k] <= (($k/$#parray)*$FDR)) { 
	    $first_pvalue_that_satisfies_equation = "seen";
	
	    # report the P-value threshold for FDR = custom%      # 1.4
	    # it is the "P-value threshold" at this FDR
	    # save it to report in the output       
	    if ($FDR_custom_perc_threshold eq "NA") {
		$FDR_custom_perc_threshold = ($k/$#parray)*$FDR;
		$FDR_custom_perc_threshold = sprintf("%.3g", $FDR_custom_perc_threshold);  # trim significant digits
	    }
	}      
	
	if ($first_pvalue_that_satisfies_equation eq "seen") {
	    $FDR_array_custom[$k] = "T";
	} else {
	    $FDR_array_custom[$k] = "F";
	}
    }
    close(FDR);
    
    #### Fill in NA for those observations that are NA  ####
    
    for ($y = $#parray; $y <= $lines; $y++) {
	push (@FDR_array10, "NA");
	push (@FDR_array5, "NA");
	push (@FDR_array1, "NA");
	push (@FDR_array_custom, "NA");
    }
    
    
    ######################################################################################## 
    # Write main data including FDR logical flags to a file                                # # 1.3
    ########################################################################################
    $line = 0;                                                                  
    open(PEN, ">penultimate") || die "couldn't create temp file penultimate, check system permssions";               
    open(FILL, "sorted_almost_final") || die "couldn't open sorted_almost_final";
    
    # Here we insert FDR values for each term right after the corrected P-value (element/column 5)
    while(<FILL>) {
	chomp($_);
	
	@fillarray = split(/\t/,$_);
	
	$line++;
	
	for ($g = 0; $g <= $#fillarray; $g++) {
	
	    if ($g < 5) {
		print PEN "$fillarray[$g]\t";
	    }
	    if ($g == 5) {
		print PEN "$fillarray[$g]\t";
		
		# this fills in the FDR values for each term until we reach the end of the sorted_P_Values file
		# this length is as long as any of the FDR_arrays, I chose FDR_array10
		if ($line <= $#FDR_array10) {
		    print PEN "$FDR_array10[$line]\t";
		    print PEN "$FDR_array5[$line]\t";
		    print PEN "$FDR_array1[$line]\t";
		    print PEN "$FDR_array_custom[$line]\t";
		}
	    }
	    # resume filling in row
	    if ($g == $#fillarray) {
		print PEN "$fillarray[$g]\n";   # add a line break for the last entry
	    } else {                            # added else-if (eliminates trailing tab in contributing genes data) # 1.4 
		if ($g > 5) {  
		    print PEN "$fillarray[$g]\t";
		}
	    }
	}
    }
    
    close(PEN);
    close(FILL);

    # write all updated data rows in file to an array called DATA
    open(PEN2, "penultimate");
    @DATA = <PEN2>;
    close(PEN2);
    
    #    print @DATA;
    
    # leapfrog function ########################### # 1.3    # 1.4 removed 
    
    
    ###############################################################################################################
    # Print final output to a file in tab-delimited format                                                        #
    ###############################################################################################################
    open(FINAL,">$results_filename") || die "couldn't create file $results_filename, check sytem permissions";
    
    # print file header

    #print FINAL "GeneMerge v$version\n\n";  
    #print FINAL "Castillo-Davis, C.I. $year. GeneMerge v$version - post-genomic data analysis.\n\n";
    #print FINAL "Output file name: $results_filename\n";
    #print FINAL "Gene Association File:  $gene_association_file\n";
    #print FINAL "Description File:  $description_file\n";
    #print FINAL "Population File:  $population_genes_file\n";
    #print FINAL "Study File:  $study_genes_file\n";
    #print FINAL "Custom FDR: $custom_FDR_percentage\%\n\n";
    
    # print column headers
    #print FINAL "GMRG_Term\tPop_freq\tPop_frac\tStudy_frac\tP\tBonf_Cor_P\tFDR_10\tFDR_5\tFDR_1\t$custom_FDR_header\tDescription\tContributing_genes\n";  %%%%%% modificado por el comando de abajo
    print FINAL "\"GO\",\"Pop_freq\",\"Pop_frac\",\"Study_frac\",\"P\",\"adj_pval\",\"FDR_10\",\"FDR_5\",\"FDR_1\",\"$custom_FDR_header\",\"Term\",\"Entry\"\n";
    
    # print main data
    print FINAL @DATA;
    
    # print footer stats
    #print FINAL "\nTotal number of genes: $total_no_detected_genes\n";
    #print FINAL "Total number of Study genes: $total_no_updown_genes\n";
    #print FINAL "Total number of Study gene GMRG terms (pop non-singletons): $no_of_up_down_GMRG_terms ($BCr)\n";

    #print FINAL "FDR Threshold P-values: [10% = $FDR_10_perc_threshold], [5% = $FDR_5_perc_threshold], [1% = $FDR_1_perc_threshold], [$custom_FDR_percentage% = $FDR_custom_perc_threshold]\n"; # 1.4
  

    #print FINAL "Genes with GMRG information: $up_down_genes_with_annotation\n";
    #print FINAL "Genes with no GMRG information: $GenesNoInfo\n";
    #print FINAL "These are:\t$genes_not_found_in_ontology\t";
    
    close(FINAL);
    
    
    ###############################################################################################################
    # Format output in HTML                                                                                       #
    ###############################################################################################################
    
    # create an html file name and the value for the custom FDR header
    #$results_filename_html = $results_filename . ".html";
    
    #open(HTML,">$results_filename_html") || die "couldn't create $results_filename_html, check system permissions";
    
    
    # print HTML header ##############################################
    
    print HTML "<html>\n";
    print HTML "<title> GeneMerge Output - $results_filename</title>\n";
    print HTML "<h3>GeneMerge v$version</h3>";
    print HTML "<h3>Castillo-Davis, C.I. $year. GeneMerge v$version - post-genomic data analysis </h3>";
    print HTML "Output File Name: $results_filename <br>";
    print HTML "Gene Association File:  $gene_association_file <br>\n";
    print HTML "Description File: $description_file  <br>\n";
    print HTML "Population File: $population_genes_file  <br>\n";
    print HTML "Study File:  $study_genes_file<br>\n";
    print HTML "Custom FDR: $custom_FDR_percentage\%<br>\n\n";
    print HTML "<p>\n";
    # leap stuff removed here # 1.4
    print HTML "<table border \= \"0\", cellpadding \= \"5\", align \= \"center\">\n";
    print HTML "<tr\>\n";
    print HTML "<td><b>GMRG Term</b></td> \n";   
    print HTML "<td><b>Pop Frequency</b></td> \n";
    print HTML "<td><b>Pop Fraction</b></td> \n";
    print HTML "<td><b>Study Fraction</b></td> \n";
    print HTML "<td><b><i>P</i>-value</b></td> \n";
    print HTML "<td><b>Bon. Corr. <i>P</i>-value</b></td> \n";
    print HTML "<td><b>10\% FDR</b></td> \n";
    print HTML "<td><b>5\% FDR</b></td> \n";
    print HTML "<td><b>1\% FDR</b></td> \n";
    print HTML "<td><b>$custom_FDR_percentage\% FDR</b></td> \n";
    print HTML "<td><b>Description</b></td> \n";
    print HTML "<td><b>Contributing genes</b></td> \n";
    print HTML "</tr>\n";
    
    
    
    # print data fields  ##############################################
    
    print HTML "<tr>\n";
    
    # step through each element of DATA (these are individual rows of data)
    for ($h = 0; $h <= $#DATA; $h++) {
	@dataline = split(/\t/, $DATA[$h]);
	
	#print each element of the DATA row with HTML mark-up
	for ($hh = 0; $hh <= $#dataline; $hh++) {          ## changed from dataline-1 (no need, no trailing tab) # 1.4
	    chomp($dataline[$hh]);                         
	    
	    # print FDR "trues" in red
	    if ($dataline[$hh] eq "T") {
		print HTML "<td> <font color = \"red\"> $dataline[$hh] </font> </td>"; 
	    } else {
		print HTML "<td> $dataline[$hh] </td>";
	    }
	}
	print HTML "</tr>\n";	
    }
    
    print HTML "</table> \n";
    print HTML "<p>\n";
    
    # print footer statistics and list genes with no association information
    print HTML "Total number of genes: $total_no_detected_genes <br>\n";
    print HTML "Total number of Study genes: $total_no_updown_genes <br> \n";
    print HTML "Total number of Study gene GMRG terms (pop non-singletons): $no_of_up_down_GMRG_terms ($BCr) <br> \n";


    print HTML "FDR Threshold <i>P</i>-values: [10% = $FDR_10_perc_threshold], [5% = $FDR_5_perc_threshold], [1% = $FDR_1_perc_threshold], [$custom_FDR_percentage% = $FDR_custom_perc_threshold] <br>\n"; # 1.4

    print HTML "Genes with GMRG information: $up_down_genes_with_annotation <br> \n";
    print HTML "Genes with no GMRG information: $GenesNoInfo <br><p>\n\n";
    print HTML "These are:<p>\n$html_genes_not_found_in_ontology\n";
    
    print HTML "\n <p> </html>\n";
    
    close(HTML);
    
    # delete temp files
    unlink("GMRGID_n_p_k_r");
    unlink("uniqueUpDownGMRGIDs");
    unlink("UpDownPoolGMRGIDs");
    unlink("SamplePoolGMRGIDs");
    unlink("gene_GMRGterms");
    unlink("GMRGID_P_Values");
    unlink("uniqueGMRGIDs");
    
    unlink("study_genes.temp");           # 1.2
    unlink("pop_genes.temp");             # 1.2
    unlink("gene_association_file.temp"); # 1.2
    unlink("description_file.temp");      # 1.2
    
    unlink("almost_final");                # 1.3
    unlink("penultimate");                 # 1.3
    unlink("sorted_almost_final");         # 1.3           

} # abort FALSE  


  
###############################################################
#                        SUBROUTINES                          # 
###############################################################

# Hypergeometric tail probability
#
# Returns one-tailed probability based on the hypergeometric distribution:
# the probability of witnessing r successes in a sample of k items from
# a pool of size n, with r having prob (proportion) p, sampling without
# replacement. Uses natural log n choose k and factorial subroutines.     

sub hypergeometric {
    my $n = $_[0];
    my $p = $_[1];
    my $k = $_[2];
    my $r = $_[3];
    my $q;
    my $np;
    my $nq;
    my $top;
    
    $q = (1-$p);
    
    $np = floor( $n*$p + 0.5 );  # round to nearest int
    $nq = floor( $n*$q + 0.5 );

    $log_n_choose_k = &lNchooseK( $n, $k );
    
    $top = $k;
    if ( $np < $k ) {
	$top = $np;
    }
    
    $lfoo = &lNchooseK($np, $top) + &lNchooseK($n*(1-$p), $k-$top);
    $sum = 0;

    for ($i = $top; $i >= $r; $i-- ) {
	$sum = $sum + exp($lfoo - $log_n_choose_k);

	if ( $i > $r) {
	    $lfoo = $lfoo + log($i / ($np-$i+1)) +  log( ($nq - $k + $i) / ($k-$i+1)  )  ;
	}
    }
    # check on underflow errors summing to more than one           # 1.4
    # in high probability cases where P ~ 1.0 with huge n and k
    # note that P-values < ~1e-300 return as zero 
    if ($sum > 1) {
	$sum = 1;
    }

    return $sum;
}

# ln factorial subroutine
sub lFactorial {
    $returnValue = 0;
    my $number = $_[0];
    for(my $i = 2; $i <= $number; $i++) {
     	$returnValue = $returnValue + log($i);
    }
    return $returnValue;
}

# ln N choose K subroutine
sub lNchooseK {
    my $n = $_[0];
    my $k = $_[1];
    my $answer = 0;

    if( $k > ($n-$k) ){
	$k = ($n-$k);
    }

    for( $i=$n; $i>($n-$k); $i-- ) {
	$answer = $answer + log($i);
    }

    $answer = $answer - &lFactorial($k);
    return $answer;
}

# sort unique subroutine
sub sortu {
    
    #input is the filename to sort
    my $file = $_[0];

    #initialize
    my $term = "";
    my %uniquehash;
    my $ulist = "";

    open(F, "$file") || die "couldn't open $file";
    while(<F>) {           
	
	$term = $_;
	# make $term the key - this way it'll overwrite  
	$uniquehash{$term} = 1;    
    }
    close(F);
    
    @unique_array = keys(%uniquehash);
    
    for($qq = 0; $qq<=$#unique_array; $qq++) {
	$ulist = $ulist . $unique_array[$qq];
    }
    return $ulist; # return a unique list
}

sub sort_column { # sorts a tab-delimited file according to a particular $col
     ($a->[$col] eq "NA" ? 1 : 
      ($b->[$col] eq "NA" ? -1 : 
       $a->[$col] <=> $b->[$col])); 
}

