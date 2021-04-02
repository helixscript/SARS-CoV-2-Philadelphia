## SARS-CoV-2-Philadelphia  
### Retrieving runs from BaseSpace and creating demultiplexed FASTQ files  

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
be seen while the pipeline is running on microb191. A run can be removed from the pipeline by simply  
renaming the run directory starting with an X, eg.  

%> mv /media/lorax/data/SARS-CoV-2/sequencing/210329_NB551353_0086_AHT7HCAFX2/   
        /media/lorax/data/SARS-CoV-2/sequencing/X210329_NB551353_0086_AHT7HCAFX2/  
        
Removing the X will cause the run to be added back into the cumulative pipeline the next hour when  
the pipeline runs again. A zip file containing all the trial reports (SARS-CoV-2_reports.zip) can be  
found in the summaries directory and can be downloaded via  
http://bushmanlab.org/data/SARS-CoV-2/summaries/SARS-CoV-2_reports.zip   
Lineage plots and trees can be found in the summaries/highQualGenomes directory.
