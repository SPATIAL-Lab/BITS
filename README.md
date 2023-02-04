# Sr-in-ivory

## Basic
This a Bayesian hierarchical model that simulates the processes of 87Sr/86Sr incorporation, turnover, and output into incremental tissues. The modeling approach is built with a two-compartment isotope exchange mechanism, with the central pool being a rapidly-exchanging but relatively small pool, and the peripheral pool being a slowly-exchanging but relatively large pool. A steady state is assumed so that pool sizes are constant. The model comprises of 4 components: the model input, process model, data model for tissue growth, and data model for 87Sr/86Sr measurements. The model is demonstrated using 3 model studies explained below.

## Model studies
Model study 1 is the **Calibration**, which is used to estimate turnover parameters using high-resolution 87Sr/86Sr measurements on the ivory from an African Savanna Elephant (*Loxodonta africana*) that was relocated to a different 87Sr/86Sr baseline.

Model study 2 is the **Fidelity test**, which uses the same turnover parameters in the Calibration to recover the 87Sr/86Sr history associated with the known relocation.

Model study 3 is the **Case study**, which is used to estimate possible 87Sr/86Sr intake history from published ivory measurements from an Alaskan Woolly Mammoth [Wooller et al., (2021)](https://www.science.org/doi/abs/10.1126/science.abg1134).

## Data
The "data" folder consists of both data from this study, and published data from [Wooller et al., (2021)](https://www.science.org/doi/abs/10.1126/science.abg1134) for the case studies. 

## Software requirements
The code is developed in R calling the standalone JAGS (Just Another Gibbs Sampler) program, which is required before the code can be run. The JAGS program can be downloaded via [this link](https://sourceforge.net/projects/mcmc-jags/). Please make sure to download the version that is appropriate for your operating system.

To call JAGS from RStudio, the R packages "rjags" and "R2jags" are also required.

Other R packages specified in the file headers help to visualize data.

## Related Publication
This version of the code is for peer review only, and is subject to major updates.

## How to cite this software
TBD