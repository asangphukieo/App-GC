# App-GC
Data preprocessing pipeline for GC- and GCxGC-MS raw data

# Description
App-GC is R-based pipeline that integrates open-source software to preprocess GC- and GCxGC-MS data. There are four main steps to process data with App-GC, 1) peak detection, 2) peak export, 3) peak alignment, and 4) peak annotation. App-GC contains various methods to dect peaks including CentWave, MatchedFilter, and eRah for GC-MS peak detection, and msPeak-NEB and msPeakG-Del1 for GcxGC peak detection. For peak alignment, there are two methods, i.e. GCalignR, and eRah-aign, for GC-MS data and one method, i.e. mSPA for GCxGC data. The final output from the pipeline is in single CSV text file, which is ready for subsequence analysis.   

![alt text](https://github.com/asangphukieo/App-GC/blob/main/workflow.png)

# Installation and requirement

# Download example data

# Usage
All parameters for the pipeline are in "parameter_config.yml" file. In the terminal screen, user can run the pipeline by single command as below

"Rscript GC_pipeline.R parameter_config.yml"

# Contact us
If you have any enquiries to analyze your data, please contact us. 
We are happy to recieve any suggestions or comments.
