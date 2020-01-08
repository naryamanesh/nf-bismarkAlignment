# nf-bismarkAlignment
Bismark alignment of bisulfite sequencing data using Next Flow pipeline

### How to run:
```
nextflow run main.nf -profile slurm \
      -w /PATH_TO/work \
      --outDir /PATH_TO/OUTDIR
```

### sample List

A csv file (sampleSheet.csv) is needed to extract sample information. As example:

```
path,group,sample,filename,R1,R2
/BaseFolder,ReadFolder,SampleName,FileName,XXXXX_R1_001.fastq.gz,XXXXX_R2_001.fastq.gz
...
```
