
# DGE-EM README

DGE-EM source code originally found at http://dna.engr.uconn.edu/?page_id=163
The code has been updated to work with the latest version of Scala.


## Building installer from source:

Install Maven then go to the root directory of this project and type:

```
mvn install
```

This will produce the artifact dge-em/target/dge-em-1.0-standard.jar which is an installer for dge-em.
Note: the installer for the original version can be found at http://dna.engr.uconn.edu/software/DGE-EM/dge-em-1.0.0-install.jar

## Installation:

 If you do not want to install DGE-EM, in the dge-em folder you can find a 
`dge-em.sh` script that runs the locally built dge-em without installation. However,
the script has to be run from the dge-em folder otherwise it will not find the required 
jar files. To run dge-em regardless of the working directory, proceed with the installation described next.

   Assuming you have the installer, call it installer.jar.  With graphical interface: on most platforms 
double clicking on the icon will start the installer, otherwise  type the following at the command line:

```
java -jar installer.jar
```

Text mode: if no graphical user interface is available, type the following at the command line and follow the instructions:

```
java -jar installer.jar -console
```

 [Optional] On Windows you might want to add the DGE-EM
   installation directory to the path, such that you can invoke 
   dge-em from any location. On Unix you can obtain a similar
   effect by creating a symbolic link to dge-em in /usr/local/bin.



## Running DGE-EM:

DGE-EM takes as input a set of known isoforms in GTF format, a fasta
file with the mRNA sequences of the known isoforms and one or more fastq
files containing the DGE tags. 

If you don't have the mRNA sequences but
have the GTF and the genome, you can use the extract-isoform-sequences-from-genome
script included in this installation to extract and optionally add polyA
tails to the isoform sequences.

The output consists of two files: one for isoform frequencies and one
for gene frequencies. Each line in these files is a pair of name and frequency.
By default, these files are called dge.iso_estimates and dge.gene_estimates, but
you can change the 'dge' prefix by using the -o option (see below).

You can run DGE-Em from the command line as follows:

```
  dge-em [options]* file.fastq [file2.fastq ...]

Option (* = required)               Description                                
---------------------               -----------                                
* -G, --GTF <GTF file>              Known genes and isoforms in GTF format     
-c, --gene-clusters <cluster file>  Override isoform to gene mapping defined in
                                      the GTF file with a mapping taken from   
                                      the given file. The format of each line  
                                      in the file is "isoform	gene".           
* -e <enzymes>                      Enzime cutting patterns (comma separated,  
                                      no spaces) e.g. CATG                     
-h, --help                          Show help                                  
--limit-ntags <Integer>             Discard all tags after this many have been 
                                      read                                     
* -m                                Isoform MRNA sequences in 5' to 3'         
                                      orientation                              
--max-mismatches <Integer>          Maximum number of mismatched allowed for a 
                                      tag (default 1).                         
-o <prefix>                         Output files prefix (default 'dge')        
-p, --prepend-enzyme                Use this flag if the recognition site for  
                                      the enzyme is not included in the tags   
--uniq                              Infer frequencies only from tags that map  
                                      to the same gene                         
```

### Example:
  
  In the dge-em/example folder you can find a trivial example to check that dge-em is working.
  In the dge-em folder you can also find a dge-em.sh script that runs the locally built dge-em without having to go through installation. 
  
  First, change directory to dge-em:

```
cd dge-em
```

  Then run:

```
sh dge-em.sh -G example/test.fa.GTF -m example/test.fa -e CATG -p example/tag1.fastq example/tag2.fastq
```
  This will produce two output files dge.iso_estimates and dge.gene_estimates.

## Revision history
Version 1.0.1 (5/09/16)  - updated to work with Scala 2.11.8

Version 1.0.0 (2/22/11)  - first public release


## Contact
For questions or suggestions regarding DGE-EM you can contact:

     Marius Nicolae (marius.nicolae@engr.uconn.edu)
     Ion Mandoiu (ion@engr.uconn.edu)
