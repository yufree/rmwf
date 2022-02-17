# RMWF

Reproducilble Metabolomics WorkFlow(RMWF) is a R package for [xcmsrocker](https://hub.docker.com/r/yufree/xcmsrocker/). It will show the workflow templates and demo data for different R-based metabolomics software. User could use this package to make meta-analysis for different workflows.

If you directly use the docker image, the rmwf package is already installed.

However, if you wanted to install locally on your own computer, you could install it from [GitHub](https://github.com/yufree/rmwf):

In RStudio console, input this command to install it:

~~~
# You need remotes package and you could install it by this command
install.packages('remotes')
remotes::install_github("yufree/rmwf")
~~~

To avoid memory limits on Dockerhub's autobuild, I ignored the EIC object for NIST 1950 and MTBLS28 data to build the R package from source code. You could remove the following line in the `.Rbuildignore` to build the package.

~~~
data/srmneic.rda
data/srmeic.rda
man/srmeic.Rd
man/srmneic.Rd
inst/demodata/untarget/MTBLS28negmzrt.csv
inst/demodata/untarget/MTBLS28posmzrt.csv
~~~

Then you could find the workflow template from RStudio:

File-New file-R Markdown-from template

Then select 'PMDDA metabolomics workflow’ to use template for PMDDA analysis.

Other template could be used to check data analysis script.

- peakpicking: from raw data to peak list

- normalization: batch correction method

- annotation: annotation from peak to peaks type or compounds

- omics: metabolomics coupled with other omics method

- statistical: statistical analysis methods

- reporducible: get data from online database

- reactomics: reactomics network analysis

## Citation

- Yu, M., Dolios, G. & Petrick, L. Reproducible untargeted metabolomics workflow for exhaustive MS2 data acquisition of MS1 features. J Cheminform 14, 6 (2022). https://doi.org/10.1186/s13321-022-00586-8
