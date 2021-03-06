# App-GC
Data preprocessing pipeline for GC- and GCxGC-MS raw data

# Description
App-GC is R-based pipeline that integrates open-source software to preprocess GC- and GCxGC-MS data. There are four main steps to process data with App-GC, 1) peak detection, 2) peak export, 3) peak alignment, and 4) peak annotation. App-GC contains various methods to dect peaks including CentWave, MatchedFilter, and eRah for GC-MS peak detection, and msPeak-NEB and msPeakG-Del1 for GcxGC peak detection. For peak alignment, there are two methods, i.e. GCalignR, and eRah-aign, for GC-MS data and one method, i.e. mSPA for GCxGC data. The final output from the pipeline is in single CSV text file, which is ready for subsequence analysis.   

![alt text](https://github.com/asangphukieo/App-GC/blob/main/workflow.png)

# Installation and requirement
## Manual installation (tested with Ubuntu 20.04) (not recommend)
  Install R library (r-base=3.6.1) using anaconda and BiocManager 
  * anaconda
  ```
  conda install -c r r-base=3.6.1
  conda install -c bioconda bioconductor-metams
  pip uninstall netCDF4 #in case libnetcdf error
  conda install -c conda-forge "libnetcdf=4.6.2"
  conda install -c r r-shiny
  conda install -c conda-forge r-shinythemes
  conda install -c conda-forge r-shinyfiles
  conda install -c conda-forge r-shinyjs
  conda install -c conda-forge r-plotly
  conda install -c conda-forge r-htmlwidgets
  ```
  * BiocManager library
  ```
  install.packages("BiocManager")
  BiocManager::install("xcms")
  BiocManager::install("CAMERA")
  BiocManager::install("GCalignR")
  BiocManager::install("filesstrings")
  BiocManager::install("msm")
  BiocManager::install("pracma")
  BiocManager::install("doSNOW")
  BiocManager::install("erah")
  ```
  * Download MSPepSearch tool and put it in working folder
  https://chemdata.nist.gov/dokuwiki/doku.php?id=peptidew:mspepsearch
  ``` 
  e.g. B:\App-GC_1_0_docker\Working\2019_02_22_MSPepSearch_x32\MSPepSearch.exe
  ```

## Installation using Docker (Recommend)
  * 1. Clone App-GC script from github and move to App-GC folder
  ```
  git clone git@github.com:asangphukieo/App-GC.git
  cd App-GC
  ```
  
  * 2. Pull app-gc docker containing (~8GB)
  ```
  docker pull asangphukieo/app-gc
  ```

  * 3. Run App-GC in commandline version
  ```
  docker run -p 3838:3838 -v B:\App-GC_1_0_docker\App-GC\:/srv/shiny-server/ -v B:\App-GC_1_0_docker\Working:/var/log/shiny-server/ -it asangphukieo/app-gc:latest bash
  ```

  or
  *  3. Run App-GC in GUI version (Linux)
  ```
  sudo docker run -p 3838:3838 -v /home/user/Downloads/App-GC/:/srv/shiny-server/ -v /home/user/Downloads/App-GC/Working:/var/log/shiny-server/ -it asangphukieo/app-gc:latest exec shiny-server 2>&1
  ```
  
  *  3. Run App-GC in GUI version (Window)
  ```
  docker run -p 3838:3838 -v B:\App-GC_1_0_docker\App-GC\:/srv/shiny-server/ -v B:\App-GC_1_0_docker\Working:/var/log/shiny-server/ -it asangphukieo/app-gc:latest exec shiny-server 2>&1
  ```
  
  - Then, for GUI version go to your browser and type http://localhost:3838/App-GC_1.0/ to open the GUI.

## Set up NIST library folder (comercial library)
put NIST library in working folder (mainlib and replib)
```
e.g. in window system B:\App-GC_1_0_docker\Working\NIST_library
```
## Set up other library folder (non-comercial library)
```
Soon.
```

# Usage
All parameters for the pipeline are in "parameter_config.yml" file. In the terminal screen, user can run the pipeline by single command as below

```
Rscript run_App-GC.R parameter_config.yml
```

# Input
Now, App-GC supports only CDF file format. For other types, you might use open source software e.g. Openchrom (https://lablicate.com/platform/openchrom) to convert the file to CDF format.

# License
MSPepSearch program may be redistributed without restriction. (https://chemdata.nist.gov/dokuwiki/doku.php?id=peptidew:mspepsearch#restrictions_and_disclaimers)

# Contact us
If you have any enquiries to analyze your data, please contact us. 
We are happy to recieve any suggestions or comments.
