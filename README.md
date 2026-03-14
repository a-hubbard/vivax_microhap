This is the analytical pipeline for the following paper:

> Hubbard, A., Solares, E., Bradley, L., Jeang, B., Yewhalaw, D., Janies, D., 
> Lo, E., Yan, G., and Hemming-Schroeder, E. (2026). PvGAP: Development of a 
> globally-applicable, highly-multiplexed microhaplotype amplicon panel for 
> *Plasmodium vivax*. *Journal of Infectious Diseases*, Under Review.

The preprint for this paper is 
[here](https://www.medrxiv.org/content/10.1101/2025.04.30.25326751v1).

This pipeline creates all of the figures in this paper from the raw data, except 
for Figure 1. This figure was created from intermediate files generated during 
the panel design phase of the project, which is not included in this pipeline. 
This pipeline also generates Supplemental Text 2, though the necessary data for 
this portion are not currently available online (see below). This pipeline does 
not generate Supplemental Text 5.

## File Structure

All of the code, data, and supporting files necessary to run this pipeline are 
contained in the repository, except for sequence data. Data files are located in 
the "resources" directory, code in the "workflow" directory, and results will be 
placed in a "results" directory that the pipeline creates if it does not already 
exist.

The sequence data for the serial dilution experiment is available from the 
Sequence Read Archive under Accession #PRJNA1435528. The MalariaGEN sequence 
data used for the paneljudge simulations may be downloaded from 
https://www.malariagen.net/resource/30/. Refer to the file paths in the pipeline 
Snakefile (see below) to understand where to place this data in the directory 
structure.

Other sequence data analyzed by the pipeline is not currently available 
publically and is not required to generate any of the figures for this paper. 
These data are necessary for Supplemental Text 2 and to run the pipeline without 
modification, however. Therefore, you may either 1) comment out/delete the 
portions of the pipeline that refer to May 2022 data, UCI 12/23 data, and the 
first serial dilution experiment or 2) contact the authors to receive the 
necessary data.

## Installing Packages

`conda` (actually `mamba`) was used for package management wherever possible, 
and the environment is described in `workflow/envs/vivax_microhap.yaml`. A small 
number of R packages could not be installed through `mamba` - these are listed 
in `workflow/envs/R_pkg_maninstall.R`. We recommend setting up the environment 
described in the YAML and then installing these packages second. It may also be 
necessary to run `tinytex::install_tinytex()` from R if LaTeX is not already 
installed on your machine.

## Running the Pipeline

The entire pipeline was automated with Snakemake, using `workflow/Snakefile`. 
Running this with the packages and data described above (though see above for 
caveats about sequence data availability) should replicate the workflow.

Note that this project was performed on a cluster with Slurm - a similar system 
will be necessary to replicate the analyses, or the pipeline will need to be 
reconfigured if your cluster has a different job scheduler.

Finally, the end of the Snakefile contains a collection of rules used to 
automatically push results from our computing cluster to cloud storage. These 
are not necessary for replication, and also will not work because you don't have 
access to the cloud storage in question. For your convenience, at the top of the 
Snakefile we have created a `replicate` rule that runs the pipeline without 
these last steps. For replication, use this rule as the target in place of the 
`all` rule.
