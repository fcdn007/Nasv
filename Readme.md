# 1.Preparation
Please install python3 and following module first: scipy，numpy，configparser and PyVcf.

# 2.Pre-processing
## 2.1 LAST Installation
NaSV was built based on LAST mapping data.LAST mapping
LAST is able to map the reads in non-overlapping split segments.LAST can be download as the zip file from http://last.cbrc.jp/, and installed 
by following comands:

```
> gunzip last.zip
> cd /path/to/lastdir
> make
```

## 2.2 Indexing Reference Genome
First you need to index your reference genome by creating a lastal database:
```
> lastdb [referencedb] [reference.fa]
```

## 2.3 Mapping
Train LAST to get the best scoring parameters for your particular alignment. We typically use a subset of > 10,000 reads for this step:
```
> last-train -Q1 [referencedb] [reads_sample.fastq] > [my-params]
```

Map your fastq data to reference:
```
> lastal -Q1 -p [my-params] [referencedb] [reads.fastq] | last-split > [reads.maf]
```
All of the above commands can also be run at once using pipes:
```
> lastal -Q1 -p [my-params] [referencedb] [reads.fastq] | \
> last-split > [reads.maf]
```

# 3.SV calling using Nasv
```
NaSV usage
> NaSV.py [-h] [-c CONFIG] [-o OUTPUT] [-i reads.maf]
NaSV arguments and parameters:
required arguments:
-i, --maf        :   /path/to/reads.maf
-o, --output     :   Give the full path to the output vcf file

optional arguments:
-h, --help       :   Show the help message and exit
-c, --config     :   Give the full path to your own ini file [ default: config.ini ]
```

optional configuration parameters:
Nasv uses a config.ini file which contains default settings for all running parameters. Users can change the parameters by creating their own config.ini file and provide this as a command line argument [-c]
```
 **Reads segments options:**
 **[Filter options]**
 *Minimum mapping quality of the segment*
 mq_cutoff=10
 *Maximum percentage of unaligned first bases for a qualified read*
 start_cutoff=0.2

 **Parameters for tuning detection and clustering of breakpoints:**
 **[Detection options]**
 *Minimum number of breakpoint-junctions (i.e. split-read junctions) for clustering*
 read_cutoff=5
 *Maximum distritution distance for a breakend position (assuming normal distribution)*
 merge_bp_cutoff=100
 *Maximum distance for two breakend positions decided as homology*
 merge_cutoff=20

 **Parameters for setting the FILTER flag in the vcf output:**
 **[Output filter options]**
 *Minimum reads percentage for a breakpoint decided as precise*
 precise_cutoff=0.5
```
