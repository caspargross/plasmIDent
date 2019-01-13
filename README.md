PlasmIdent
==========

[![Build Status](https://travis-ci.org/caspargross/PlasmIdent.svg?branch=master)](https://travis-ci.org/caspargross/PlasmIdent)

Identification of circular plasmid in bacterial genome assemblies using long reads.

- Gene prediction (Glimmer3)
- Identification of antibiotic resistance genes (RGI)
- Coverage analysis (Mosdepth)
- GC Content and GC Skew
- Identification of overlapping reads

Requirements
------------

- Linux or Mac OS (not tested on Windows, might work with docke)
- Java 8.x


Installation 
------------

1) Install nextflow
2) Download pipeline script

```
git clone https://github.com/caspargross/plasmident
```

3) Pull docker image

```
docker pull caspargross/plasmident
```

### Alternative: Local install
When you are unable to use the docker environment, it is possible to directly install the conda environments in the `env/` folder by running the following command:

``` 
conda env create -f env/PI_env.yml
```

In this case the application should be run with by adding the following parameters: `-profile local`

Run Application
---------------

Run the application by providing a tab-separated input file with sample ID and locations of assembly (.fasta) and read files (.fastq/fastq.gz)

```
nextflow run plasmident --input read_locations.tsv

```


### Available parameters

- `--outDir` Path of output folder
- `--seqPadding` Number of bases added at contig edges to improve long read alignment [Default: 2000]
- `--covWindow` Moving window size for coverage and gc content calculation [Default: 50]
- `--cpu` Number of threads used per process

Results
-------

![FlowChart](https://github.com/caspargross/plasmident/example_output.png)
