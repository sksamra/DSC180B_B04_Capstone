# Title: Genetics Replication Project

#### Capstone Project: Data Science DSC180A

#### Section B04: Genetics

#### Authors: Saroop Samra, Justin Kang

#### Date : 10/23/2020

### Overview

This repository code is for the replication project for the paper: Post-mortem molecular profiling of three psychiatric disorders (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0458-5). The data includes that RNA sequences from tissues from three regions of the brain (interior cingulate cortex, dorsolateral prefrontal cortex, and nucleus accumbens), and is from three groups of 24 patients each diagnosed with a Psychiatric disorder (schizophrenia, bipolar disorder, or major depressive disorder) as well as a control group. The analysis performed shows the correlation of certain genes and psychiatric disorders.

### Running the project

•	To install the dependencies, run the following command from the root directory of the project: 

    pip install -r requirements.txt


### target: data
•	To process the data, from the root project directory run the command:

    python3 run.py data

•   The data pipeline step takes the .fastq compressed files as input and then applies two transformations: process and align

•	This pipeline step also uses an additional CSV file that is the SRA run database, a sample looks like as follows:

    Run age_at_death    Brain_pH    brain_region    Bytes   Center Name clinical_diagnosis
    SRR3438555  40  6.76    AnCg    2585349730  GEO Control



•   The configuration files for the data step are stored in config/data-params.json. These include the parameters for the tools as well as the directories used for storing the raw, temporary and output files.

    "raw_data_directory": "./data/raw",
    "tmp_data_directory": "./data/tmp",
    "out_data_directory": "./data/out",

•   The configuration also includes an attribute to the SRA run input database (described above), and an attribute of where to store that in the data folder. Additional filter attributes are included for ease of use to avoid processing all patients, if this filter_enable is set it will only process a subset of SRA rows (filter_start_row to filter_start_row + filter_num_rows).

    "sra_runs" : {
        "input_database" : "/datasets/srp073813/reference/SraRunTable.csv",
        "output_database" : "data/raw/SraRunTable.csv",
        "filter_enable" : 0,
        "filter_start_row" : 120,
        "filter_num_rows" : 10   
    },
    

•	The first transformation of the data is "process" that uses the following data configuration below. Otherwise it will invoke cutadapt which finds and removes adapter sequences. The attributes include the adapters (r1 and r2) to identify the start and end of pairs are a JSON array. The attribute enable allows to disable this cleaning step, instead it will simply copy the paired files from the source dataset. The arguments attribute allows flexible setting of any additional attribute to the cutadapt process. Finally, we have two wildcard paths that indicate the location of the SRA fastq pair files (fastq1 and fastq2).

    "process" : {
        "enable" : 1,
        "tool" : "/opt/conda/bin/cutadapt",
        "r1_adapters" : ["AAAAA", "GGGG"],
        "r2_adapters" : ["CCCCC", "TTTT"],
        "arguments" : "--pair-adapters --cores=4",
        "fastq1_path" : "/datasets/srp073813/%run_1.fastq.gz", 
        "fastq2_path" : "/datasets/srp073813/%run_2.fastq.gz"
    },
    
•   The second transformation of the data is "aligncount" that can be set to either use STAR or Kallisto. The choice of STAR or kallisto step is controlled by the aligncount attribute:

    "aligncount" : "kallisto",

•   kallisto uses the index_file attribute is the location of the directory of the reference genome, which for this replication project was GRCh37_E75. The arguments attribute allows flexible setting of any additional attribute to the kallisto process. Including the bootstaro samples.The attribute enable allows to disable this alignment step, this is useful for debugging the process prior step, for example, you can run quality checks on the processed fastq files before proceeding to alignment. 

    "kallisto" : {
        "enable" : 1,
        "tool" : "/opt/kallisto_linux-v0.42.4/kallisto",
        "index_file" : "/datasets/srp073813/reference/kallisto_transcripts.idx",
        "arguments" : "quant -b 8 -t 8"
    },

•   STAR uses the gene_path attribute is the location of the directory of the reference genome, which for this replication project was GRCh37_E75 as described in the reference_gene attribute. The arguments attribute allows flexible setting of any additional attribute to the STAR process. Including TranscriptomeSAM in the quantMode arguments will also output bam files. Additionally, the log file gets outputted which has PRUA (percentage of reads uniquely aligned). The attribute enable allows to disable this alignment step, this is useful for debugging the process prior step, for example, you can run quality checks on the processed fastq files before proceeding to alignment. 

    "STAR" : {
        "enable" : 1,
        "tool" : "/opt/STAR-2.5.2b/bin/Linux_x86_64_static/STAR",
        "reference_gene" : "GRCh37_E75",
        "gene_path" : "/path/to/genomeDir",
        "arguments" : "--runMode alignReads --quantMode GeneCounts --genomeLoad LoadAndKeep --readFilesCommand zcat --runThreadN 8"
    },




•   The process and align transformation work on each of the samples. After each sample iteration, the temporary fastq files will be deleted to reduce storage requirements.


•   Example processing:

    python3 run.py data

    # ---------------------------------------------------
    # Process
    # ---------------------------------------------------
    # Starting sample # 1  out of  352
    /opt/conda/bin/cutadapt --pair-adapters --cores=4  -a AAAAA -a GGGG -A CCCCC -A TTTT -o ./data/tmp/out.1.fastq.gz -p ./data/tmp/out.2.fastq.gz /datasets/srp073813/SRR3438555_1.fastq.gz /datasets/srp073813/SRR3438555_2.fastq.gz
    /opt/kallisto_linux-v0.42.4/kallisto quant -b 100 -i /datasets/srp073813/reference/kallisto_transcripts.idx ./data/tmp/out.1.fastq.gz ./data/tmp/out.2.fastq.gz -o ./data/tmp/SRR3438555_ReadsPerGene.out.tab
    rm ./data/tmp/out.1.fastq.gz
    rm ./data/tmp/out.2.fastq.gz
    # ---------------------------------------------------
    # Starting sample # 2  out of  352
    /opt/conda/bin/cutadapt --pair-adapters --cores=4  -a AAAAA -a GGGG -A CCCCC -A TTTT -o ./data/tmp/out.1.fastq.gz -p ./data/tmp/out.2.fastq.gz /datasets/srp073813/SRR3438556_1.fastq.gz /datasets/srp073813/SRR3438556_2.fastq.gz
    /opt/kallisto_linux-v0.42.4/kallisto quant -b 100 -i /datasets/srp073813/reference/kallisto_transcripts.idx ./data/tmp/out.1.fastq.gz ./data/tmp/out.2.fastq.gz -o ./data/tmp/SRR3438556_ReadsPerGene.out.tab
    rm ./data/tmp/out.1.fastq.gz
    rm ./data/tmp/out.2.fastq.gz
    # ---------------------------------------------------


### target: merge
•   To merge gene count and/or BAM files generated from the data target, from the root project directory run the command:

    python3 run.py merge

•   The configuration files for the data step are stored in config/count-params.json. These include the parameters for the count merge and bam merge and it's associated arguments.

•   The format attrbute informs if to process kallisto (or STAR) files. The gene counts are merged into a TSV file and as well as a A feature table based on the SRA run table. Additional STAR attributes in the JSON allow you to specify skiprows used when processing the STAR gene count files as well as identifying the column from the STAR gene matrix file to use as the column used to. The STAR log files (with PRUA information) and add an additional PRUA features to feature table. You can limit the features by setting the "features" json attribute to only the features you want to have. There is an additional imputes attribute that allows you to impute any column with missing data. The attributes also include the "filter_names" gene table used to remove X and Y chromosomes as well as removing false-positive genes. Finally, we can rename the freature columns before we save out the feature table.

    "count" : {
        "enable" : 1,
        "format" : "kallisto",
        "skiprows" : 4,
        "column_count" : 4,
        "skip_samples" : ["SRR3438888"],
        "filter_keep_genes" : "NM_",
        "filter_remove_genes" : ["chrX", "chrY"],
        "filter_names" : "/datasets/srp073813/reference/Gene_Naming.csv",
        "run_database" : "data/raw/SraRunTable.csv",
        "imputes" : ["Brain_pH"],
        "features" : ["Run", "clinical_diagnosis", "age_at_death", "Brain_pH", "brain_region", "post-mortem_interval"],
        "rename" : {"age_at_death" : "Age", "post-mortem_interval": "PMI", "Brain_pH": "pH", "clinical_diagnosis" : "Disorder"},
        "output_matrix" : "data/out/gene_matrix.tsv",
        "output_features" : "data/out/features.tsv"
    },

•   For bam merging, which should not be enabled by default, we use the "samtools" merge feature that takes all the BAM files and combine them into one merged BAM file. 


    "bam" : {
        "enable" : 0,
        "output" : "data/tmp/merged.bam",
        "tool" : "/usr/local/bin/samtools",
        "arguments" : "merge --threads 8"
    },


•   Example processing:

    python3 run.py merge

    # ---------------------------------------------------
    # Merge
    Input: SRR3438605_ReadsPerGene.out.tab
    Input: SRR3438604_ReadsPerGene.out.tab
    Output: data/out/gene_matrix.tsv data/out/features.tsv
    # Finished
    # ---------------------------------------------------



### target: normalize
•   To normalize the aligned merge counts, from the root project directory run the command:

    python3 run.py normalize

•   The configuration files for the data step are stored in config/normalize-params.json. 

•   We use a custom R script which uses the DESeq2 module to take the input merged gene counts and the experiment features and outputs two normalized counts files. The analysis is done for all samples in the SRA run table. The output_dir sets the output location for the normalized count matrix files. One file is the standard normalized counts using the DESeq2 module, and the second normalized count file is after a Variable Stablization Transform (LRT). We also have a "max_genes" attribute that will filter the genes and removes ones that have little to no variance across disorder vesus control.

•   The data JSON configuration file also holds an array of samples, a sample looks like as follows:
    
    {
        "output_dir" : "data/out",
        "DESeq2" : {
            "Rscript" : "/opt/conda/envs/r-bio/bin/Rscript",
            "source" : "src/data/normalize.r",
            "input_counts" : "data/out/gene_matrix.tsv",
            "input_features" : "data/out/features.tsv",
            "max_genes" : 8000
        },
        "cleanup" : 0,
        "verbose": 1
    }


•   Example processing:

    python3 run.py normalize

    # ---------------------------------------------------
    # Normalize
    Rscript  src/data/normalize.r data/out/gene_matrix.tsv data/out/features.tsv data/out/
    [1] "Output data/out/normalized_counts.tsv data/out/vst_transformed_counts.tsv"
    # Finished
    # ---------------------------------------------------


### target: analysis
•   To perform the analysis for the gene counts, from the root project directory run the command:

    python3 run.py analysis

•   The configuration files for the data step are stored in config/analysis-params.json. 

•   We use a custom R script which uses the DESeq2 module to take the input merged gene counts and the experiment features and outputs 3 sets of files for each brain region. Each brain region will compare a disorder versus Control. This will result in a total of 9 sets of files (3 brain regions x 3 disorder pair comparisons). Each output set includes a Likelihood Ratio Test (LRT) using the full and reduced model as specified in the attributes below as well as a MA-Plot and Heatmap. The additional attributes include the property of doing parallel processing for DESeq2.
    
    {
        "output_prefix" : "data/out/%brain_region%",
        "DESeq2" : {
            "Rscript" : "/opt/conda/envs/r-bio/bin/Rscript",
            "brain_regions" : ["AnCg", "nAcc", "DLPFC"],
            "disorders" : ["Major Depression", "Schizophrenia", "Bipolar Disorder"],
            "control" : "Control",
            "disorder_comparisons" : [["Major Depression", "Bipolar Disorder"], ["Major Depression", "Schizophrenia"], ["Bipolar Disorder", "Schizophrenia"]],
            "input_counts" : "data/out/pca_normalized_counts.tsv",
            "input_features" : "data/out/features.tsv",
            "source" : "src/analysis/analysis.r",
            "full" : "Age+PMI+pH+Disorder",
            "reduced" : "Age+PMI+pH",
            "parallel" : 0
        },
        "cleanup" : 0,
        "verbose": 1
    }




•   Example processing:

    python3 run.py analysis

    # ---------------------------------------------------
    # Normalize
    Rscript src/analysis/analysis.r data/out/AnCg/gene_matrix.tsv data/out/AnCg/features.tsv data/out/AnCg/Major_Depression_vs_Bipolar_Disorder/ Age+PMI+pH+Disorder Age+PMI+pH
    ...


### target: visualize

•   The visualize pipeline step can be invoked as follows:

    python3 run.py visualize

•   The configuration files for the data step are stored in config/visualize-params.json. The output will include 7 sets of charts: Gene Spread Variance Histogram, SRA Linear Correlation between SRA chart, MA-Plot 3x3 chart, Heat Map 3x3 chart, 3x3 Histogram, 9x9 Correlation Matrix and a Disorder Venn Diagram. Each chart type has flexible settings to control the input and layout for the charts as shown below:

    "gene_hist" : {
        "enable" : 1,
        "max_genes" : 8000,
        "nbins" : 100,
        "title" : "Distribution of Genes Based on Spread Metric: All vs Top Genes"
    },
    "sra_lm" : {
        "enable" : 1,
        "sra" : ["SRR3438555", "SRR3438560"],
        "normalized_counts" : "data/out/normalized_counts.tsv",
        "vst_counts" : "data/out/vst_transformed_counts.tsv",
        "title" : "%sra% Regression Log(Norm) v VST counts"
    },
    "ma_plot" : {
        "enable" : 1,
        "brain_regions" : ["AnCg", "DLPFC", "nAcc"],
        "disorders" : ["Schizophrenia", "Bipolar_Disorder","Major_Depression"],
        "src_image" : "MAplot.png",
        "title" : "MA Plot: Brain Region vs Disorder"
    },
    "heat_map" : {
        "enable" : 1,
        "brain_regions" : ["AnCg", "DLPFC", "nAcc"],
        "disorders" : ["Schizophrenia", "Bipolar_Disorder","Major_Depression"],
        "src_image" : "heatmap.png",
        "title" : "Heat Map: Brain Region vs Disorder"
    },
    "histogram" : {
        "enable" : 1,
        "brain_regions" : ["AnCg", "DLPFC", "nAcc"],
        "disorders" : ["Schizophrenia", "Bipolar_Disorder","Major_Depression"],
        "title" : "Histograms Differential Gene Expression vs Control",
        "ylim" : 1600
    },
    "corrmatrix" : {
        "enable" : 1,
        "title" : "Spearman Correlations of log2 fold gene expression"
    },
    "venn" : {
        "enable" : 1,
        "brain_regions" : ["AnCg", "DLPFC", "nAcc"],
        "pvalue_cutoff" : 0.05,
        "disorders" : ["Schizophrenia", "Bipolar_Disorder","Major_Depression"],
        "title" : "Venn Diagram Disorders"
    },


•   Example processing:

    python3 run.py visualize

    # ---------------------------------------------------
    # Analysis
    mkdir data/out/AnCg/
    mkdir data/out/AnCg/Major_Depression/
    AnCg x Major Depression vs control
    Rscript src/analysis/analysis.r data/out/AnCg/gene_matrix.tsv data/out/AnCg/features.tsv data/out/AnCg/Major_Depression/ full=Age+PMI+pH+Disorder reduced=Age+PMI+pH charts=1
    ...
    # Finished
    # ---------------------------------------------------


### target: qc

•   The quality pipeline step can be invoked as follows:

    python3 run.py qc

•   The configuration files for the data step are stored in config/qc-params.json. These include the parameters for the output directory where the quality HTML reports will be outputted. 

    "outdir" : "data/out",
    "inputs" : "data/tmp",

•   For fastq files, the quality tool attribute is set to fastqc and that includes attributes to extract reports or keep them in a zip file. To enable this quality check make sure you set the cleanup to 0 in the data configuration pipeline as well as to disable the STAR processing, this will retain the fastq.qz files after the data pipeline step is executed.

    "fastq" : {
        "enable" : 1,
        "tool" : "/opt/FastQC/fastqc",
        "extract" : 1   
    },

•   For bam files, the quality tool attribute is set to picard and that includes attributes such as collecting alignment summary metrics. To enable this quality check make sure you set the cleanup to 0 in the data configuration pipeline and add 'TranscriptomeSAM' to the arguments for STAR which will then output BAM files that will be retained after the data pipeline step is executed.

    "bam" : {
        "enable" : 1,
        "tool" : "java",
        "jar" : "/opt/picard-tools-1.88/CollectAlignmentSummaryMetrics.jar"
    },
    

•   Example processing:

    python3 run.py qc

    # ---------------------------------------------------
    # Quality Check
    fastqc data/tmp/out.1.fastq.gz --outdir=data/out --extract
    fastqc data/tmp/out.2.fastq.gz --outdir=data/out --extract
    java -jar /opt/picard-tools-1.88/CollectAlignmentSummaryMetrics.jar INPUT=data/tmp/SRR3438604_Aligned.bam OUTPUT=data/out/SRR3438604_Aligned.bam.txt
    java -jar /opt/picard-tools-1.88/CollectAlignmentSummaryMetrics.jar INPUT=data/tmp/SRR3438605_Aligned.bam OUTPUT=data/out/SRR3438605_Aligned.bam.txt
    # Finished
    # ---------------------------------------------------


### target: report
•   To generate the report from the notebook, run this command:

    python3 run.py report

•   The configuration files for the data step are stored in config/report-params.json. 

    {
        "tool": "jupyter",
        "args": "nbconvert --no-input --to html --output report.html notebooks/report.ipynb",
        "verbose" : 1
    }


### target: clean 

•	To clean the data (remove it from the working project), from the root project directory run the command:

python3 run.py clean


### target: all 

•   The all target will execute the following steps in sequence: data, merge, normalize, analysis and visualize. It can be executed as follows:

python3 run.py all


### Future Work

•	New pipeline step: train. This step will take the analyzed data and train a model to predict given a patients RNA sequences which of the 4 classifications they are (schizophrenia, bipolar disorder, major depressive disorder, or none) using a subset of the original data as training data and validation data.

•	New pipeline step: predict. This step will use the model to predict the classification for a given RNA sequences on the test data and reporting the classification errors



### Major Change History


Date:  12/1/2020

Work completed:

- Update README

- Update requirements pip

Date:  12/1/2020

Work completed:

- Refactored of visualize, report

- Bug fixes


Date:  11/25/2020

Work completed:

- Refactored of analysis

- Bug fixes


Date:  11/18/2020

Work completed:

- Refactored the data pipeline 

- Added kallisto to data target step

- New tagets added: analysis


Date:  11/06/2020

Work completed:

- Refactored the data pipeline and targets (includes cutadapt, STAR, samtools merge, R DESeq2, picard)

- New tagets added: merge, normalize, all

- The data target now will output gene count files 

- The merge target will merge gene count and/or bam files

- The qc target now will also support picard for bam files, and for fastq files it will iterate all files in a directory to generate multiple reports

- A stub is added for the normalize pipeline step. This will in the future exectue an R script.



Date:  10/24/2020

Work completed:

- Created basic pipeline to copy and clean the data folders




### Responsibilities


* Saroop Samra, refactored the data pipeline step. It now is comprehensive and includes the process, align, merge and normalize steps. This required updating the run.py, defining a new JSON schema, updating the etl.py. The etl.py source file has been  refactored so even though it is nearly 300 lines of code, it is documented and refactored into 7  functions to make it easier to comprehend and maintain. Saroop also added an "all" pipeline step that runs all the individual pipeline steps, one after another. Saroop also implemented the merge pipeline step which iterates a directory for the gene count TSV files and uses pandas to merge them into one master gene count file as well as a experiment table that associates the sample labels with the patient features. (as well as optionally merging bam files). The qc pipeline step now also supports picard for doing quality checks on bam files. Saroop also added a patient sample database JSON file and refactored the data configuration JSON to include a reference to the sample database as well as a filter feature. Recently Saroop worked on the normalize and merge steps to take into account filtering of genes. Additionaly, Saroop worked on the analysis processing using DESeq2 which processsed the 9 data sets (3 brain regions x 3 disorders). Finally, Saroop worked on the visualization step to create the different charts.


* Justin Kang, worked closely with Saroop on the design and code review of the pipeline steps. He worked on the references README which includes all the references on the tools we plan to use as well as the original paper and references for the medical terms in the paper. Implementation for the fastqc code was completed in config/data-quality-params.json. The data_quality pipeline step was also completed by implementing code in run.py and in src/data_quality/qc.py. This step is used to call fastqc which generates reports that you can use to assess the conditions of the raw  sequencing data prior to any concrete analysis. Justin also converted our original report into a jupyter notebook with help from Saroop which can be found in the notebooks folder as report. ipynb. Justin performed much of the exploratory data analysis (EDA) for our main report which was done in a separate notebook with the most relevant information being transfered over to the main report. Descriptions for the charts and tables presented in the main report were also completed with help from Saroop. Explanations for cutadapt and FastQC implementation were also created by Justin after performing separate tests to validate our final decision in not using cutadpat



