## SARS-CoV-2-Philadelphia  

The SARS-CoV-2 sequencing and variant calling pipeline process viral samples (saliva, nasal swabs, ect.)  
by reverse transcription of the viral RNA to make a cDNA copy, PCR amplification of genome segments,   
Nextera library preparation, and Illumina sequencing of overlapping amplicons created with the ARTIC   
primer set which amplifies the SARS-CoV-2 genome as 98 amplicons.  
  
Sequencing reads are aligned to a reference genome with BWA, alignments are filtered and processed with  
the samtools toolkit and variant positions are called with the bcftools toolkit. Variants are called for  
positions spanned by ≥ 5 reads, yield PHRED scores ≥ 20, and more than 50% of spanning reads report the   
same mutant base. Each processed sample yields a rich VSP (Viral SPecimen) data object contain details   
about the alignment to the reference genome, variant calls, PANGO lineage calls, and denovo contig assemblies   
which are saved to the local file system. These VSP data objects are collated into detailed subject reports   
and used for downstream analyses.    

## Bushman group specific instructions 

#### Downloading and processing sequencing data
Our sequencing data is stored on microb120.med.upenn.edu at this location: /media/sequencing/Illumina.  
This directory has write access for members of the data_manager user group and Aoife, Derin and Hriju  
are members of this group. Data is pulled from BaseSpace every hour automatically using this script:  
/media/sequencing/Illumina/retrieve.R which can be ran manually if needed, eg. Rscript retrieve.R.   
FASTQ files will be created and demultiplexed automatically when the data is pulled from BaseSpace  
and a SampleSheet.csv file is stored in BaseSpace. A SampleSheet file will need to be created and   
FASTQs generated manually if a SampleSheet is not pulled from BaseSpace.  
  
#### Example of creating a SampleSheet file and generating FASTQs
This series of commands creates a SampleSheet.csv file and generates FASTQ files for run  
210329_NB551353_0086_AHT7HCAFX2  and stores the output in the expected location of  
[run id]/Data/ Intensities/BaseCalls
  
  %> cd 210329_NB551353_0086_AHT7HCAFX2  
  %> nano SampleSheet.csv     (paste SampleSheet info into text editor and save with Ctrl+o)  
  %> cd ..  
  %> bcl2fastq --no-lane-splitting --create-fastq-for-index-reads -R 210329_NB551353_0086_AHT7HCAFX2      
        -o 210329_NB551353_0086_AHT7HCAFX2/Data/Intensities/BaseCalls  
        
To view which samples where demultiplexed and the files sizes ordered from largest   
to smallest, use these commands to change into the output directory and list the files:  
  
  %> cd 210329_NB551353_0086_AHT7HCAFX2/Data/Intensities/BaseCalls  
  %> ls -alhS \*R1\*  
   
If the resulting files appear incorrect, ie. wrong number of sample files or sizes appear odd   
(negative controls with large files indicating many reads), the processing can be repeated after   
removing the output:  
  
  %> rm 210329_NB551353_0086_AHT7HCAFX2/Data/Intensities/BaseCall/\*.gz  
  %> (edit the SampleSheet to make changes to the bar codes)  
  %>  bcl2fastq --no-lane-splitting --create-fastq-for-index-reads -R 210329_NB551353_0086_AHT7HCAFX2      
      -o 210329_NB551353_0086_AHT7HCAFX2/Data/Intensities/BaseCalls  
    
### Running the SARS-CoV2 pipeline  
The pipeline software runs automatically every hour on our computational server microb191 and acts  
on data stored on microb120. The contents of the shared Google spreadsheet is automatically downloaded    
every 10 minutes and stored in /media/lorax/data/SARS-CoV-2/samples.tsv.   
  
To include sequencing data in our cumulative pipeline, create a new sequencing run directory in  
/media/lorax/data/SARS-CoV-2/sequencing and copy all the R1 and R2 demultiplexed reads into this new folder, eg.  
  
  %> mkdir /media/lorax/data/SARS-CoV-2/sequencing/210329_NB551353_0086_AHT7HCAFX2   
  %> cp /media/sequencing/Illumina/210329_NB551353_0086_AHT7HCAFX2/Data/Intensities/BaseCalls/\*_R\*    
          /media/lorax/data/SARS-CoV-2/sequencing/210329_NB551353_0086_AHT7HCAFX2/  
           
The pipeline on microb191 will act on the new data when it runs the next hour and return the results  
to /media/lorax/data/SARS-CoV-2/summaries.  An empty file (media/lorax/data/SARS-CoV-2/working) will  
be seen while the pipeline is running on microb191.  
