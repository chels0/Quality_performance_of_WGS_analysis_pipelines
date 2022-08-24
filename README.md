# Quality_performance_of_WGS_analysis_pipelines
This project focuses on comparing different pipelines commonly used for genome analysis of Campylobacter jejuni to each other to find out how different software configurations affect cgMLST allele calling results as well as QUAST output results. A nextflow scientific workflow runs every pipeline and outputs scaffolds and contigs as well as QUAST results. These scaffolds and contigs are then inputted into chewBBACA which outputs allele calling results for every single pipeline and the results are processed with python scripts found here in this project. 

## Prerequisite software
You will need the following software

* Nextflow
* Miniconda

Nextflow can be downloaded here: https://www.nextflow.io/
Miniconda can be downloaded here: https://conda.io/projects/conda/en/latest/user-guide/install/download.html

## Prerequisite data
You will need a cgMLST schema which should be added into the raw data folder. This project is compatible with the Innuendo project's cgMLST schema for Campylobacter jejuni (https://www.efsa.europa.eu/sv/supporting/pub/en-1498) and can be downloaded here https://zenodo.org/record/1322564. 

The downloaded folder "Campylobacter_jejuni" should be placed inside the Raw_data folder.

## Input files
You will need to place the following files/folders in the folder Raw_data:

* The cgMLST schema, see "Prerequisite data"

* Your .fastq files representing your samples should be placed in the Samples folder inside Raw_data

## Installation
Clone git repository into your chosen directory.

Run all commands in the folder Quality_performance_of_WGS_analysis_pipelines/

Install software with conda:
> conda env create --file=scripts/config_files/conda_env_setup.yml

> conda env create --file=scripts/config_files/conda_python_env_setup.yml


## Changing/setting up parameters for workflow
Change parameters in the file: parameter_settings.txt

## Running the workflow and all other scripts 

Run all commands in the folder Quality_performance_of_WGS_analysis_pipelines/

Run with command:
> bash running.sh [options]

Options:

-o <output_dir>
  specify your output directory. The inputted output directory should NOT end with a "/". If no output directory is selected the outputs will be in the folder "Results" in the repository.
  
Example:

bash running.sh -o /data/my_results


## Plotting results
Run with command:
> bash plotting_master.sh

Modifications might need to be done to the plotting_master flags. 

## Full report
The full report done for this project can be found here:

http://uu.diva-portal.org/smash/record.jsf?dswid=-7761&pid=diva2%3A1681925&c=1&searchType=SIMPLE&language=sv&query=chelsea+ramsin&af=%5B%5D&aq=%5B%5B%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=all
