# Dependencies

* GNU parallel 20180322
* Inkscape 0.91 r13725 (Feb 22 2016)
* chromedriver
* EMBOSS:6.6.0.0 Suite
 * https://anaconda.org/bioconda/emboss
* segmehl 0.2.0-418
 * https://anaconda.org/bioconda/segemehl
* samtools Version: 0.1.19-44428cd
 * https://anaconda.org/bioconda/samtools
* wigToBigWig v 4
 * https://anaconda.org/bioconda/ucsc-wigtobigwig
* bedToBigBed v. 2.7
 * https://anaconda.org/bioconda/ucsc-bedtobigbed
* faToTwoBit
 * https://anaconda.org/bioconda/ucsc-fatotwobit
* twoBitInfo
 * https://anaconda.org/bioconda/ucsc-twobitinfo
* blastp
 * https://anaconda.org/bioconda/blast 
# Transcriptom

In the contribution we used all raw reads from the bio sample ```PRJNA655119```
for the expression analysis.
```https://www.ncbi.nlm.nih.gov/bioproject/?term=prjna655119```
For downloading the reads the SRA Toolkit is needed.
```https://ncbi.github.io/sra-tools/install_config.html```


# Genomes

The genomes must be stored in the data dir in the folder ```genome``` as .gbk files. Like so.

```
data_dir/genome/
├── anaero.gbk
├── bact.gbk
├── bifi.gbk
├── blautia.gbk
├── clostri.gbk
├── ecoli.gbk
├── ery.gbk
└── lacto.gbk
```

The short names are used throughout the scripts and should not be changed.
The full names are the following.

| short name | scientific name                                                         | taxid   | assembly id     |
| -          | -                                                                       | -       | -               |
| anaero     | Anaerostipes caccae DSM 14662                                           | 411490  | GCA_014131675.1 |
| bact       | Bacteroides thetaiotaomicron VPI5482                                    | 226186  | GCA_014131755.1 |
| bifi       | Bifidobacterium longum NCC2705                                          | 206672  | GCF_000007525.1 |
| blautia    | Blautia producta ATCC 27340 DSM 2950                                    | 1121114 | GCA_014131715.1 |
| clostri    | Clostridium butyricum DSM 10702                                         | 1316931 | GCA_014131795.1 |
| ecoli      | Escherichia coli str K12 substr MG1655                                  | 511145  | GCF_000005845.2 |
| ery        | Erysipelatoclostridium ramosum DSM 1402                                 | 445974  | GCA_014131695.1 |
| lacto      | Lactobacillus plantarum subsp plantarum ATCC 14917 JCM 1149 CGMCC 12437 | 525338  | GCA_014131735.1 |

All genomes can be downloaded with the script
```scripts/build_db/download_all_genomes.sh```. 
# Structure of scripts

All scripts neceary for the paper can be found in the folder scripts.
The scripts can be divided into categorys, which each have its own subdirectory.

 1. build\_db:
   * Build the two main protein DBs, 6frame and proteome
 2. proteom\_6frame\_map:
   * Builds two hash map showing the relation between each protein in both DBs.
 3. transcriptom:
   * Builds bed files from raw RNA-Seq data.
 4. comet:
   * This directory is holding all relevant scripts that where used to build
     the PSMs (Peptide to Spetrum Maps).
 5. PSM_accumulation
   * These scripts collect and filter all data from the comets PSMs and to
     provide them for further analysis.
 6. protein_accumulation
   * These scripts assess likely protein candidate.
 7. figure_plotting
   * All scripts for plots that where automatically generated from data are here.
 8. HTML_building
   * Webistes for listening and visualisation can be found here.

# Parameters

The file ./scripts/parameters holds paramters for the script to run

 1. data_dir
  * Path in which all results will be stored.
 2. 
