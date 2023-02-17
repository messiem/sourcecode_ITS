# Source code for paper "Coastal upwelling drives ecosystem temporal variability from the surface to the abyssal seafloor" #

This toolbox contains the programs used to reproduce results and figures from Messié et al. (2023).   
The programs are written for Matlab.

* * *

### Get started ###

See script:  

	start_sourcecode_ITS

(ITS stands for Integrated Time Series). Set the `community` variable to either surface, midwater, or benthos to reproduce results for each time series.    
The code generates figures, saved in `outputs/`, reproducing results from Fig. 1, 2, 3, and 5 from the paper with 2 main differences. First, the code runs for a given community so figures are broken down by community instead of being combined as in the paper. Second, the figures are simplified versions of the paper figures, with Matlab native colorbars and functions. Please contact me if you want to be able to reproduce the exact same figures as in the paper.
  
### Description of functions ### 

`its_load_taxa`: reads the csv file in `data/` for the requested community, and displays the time series  
`its_load_upw`: reads the upwelling time series provided in `data/`  
`its_compute_PCA`: computes Principal Component Analysis of the taxonomic time series  
`its_compute_integration`: computes wind integration following Di Lorenzo and Ohman (2013, https://doi.org/10.1073/pnas.1218022110 )  
`its_compute_adjusted_correlation`: compute Pearson correlations, and adjusted p-values following Pyper and Peterman (1998, https://doi.org/10.1139/f98-104 )  
`its_agestructured_model`: returns the modeled density over time of a population characterized by a given natural mortality rate and forced by an input forcing function  

### Description of data ###

**Taxonomic time series:**  
`surface.csv`, `midwater.csv`, `benthos.csv`  
Time series for each community on original time steps, except monthly for benthos (lines = time, columns = taxonomic groups).

**Upwelling time series:**  
`upwelling_MontereyBay.csv`, `upwelling_PtConception.csv`  
Daily time series of upwelling computed for the Monterey Bay and Point Conception regions, as explained in the paper. Primary sources of data are NDBC buoys [46042](https://www.ndbc.noaa.gov/station_history.php?station=46042) for Monterey Bay and [PTGC1](https://www.ndbc.noaa.gov/station_history.php?station=ptgc1) for Point Conception. 


* * *

### Reference ###

Please refer this paper when using these scripts (PDF available upon request):  

Messié, M., R.E. Sherlock, C.L. Huffard, J.T. Pennington, C.A. Choy, R.P. Michisaki, K. Gomes, F.P. Chavez, B.H. Robison, and K.L. Smith Jr (2023). **Coastal upwelling drives ecosystem temporal variability from the surface to the abyssal seafloor**.  *Proceedings of the National Academy of Sciences*, in press.

* * *

### Contact ###

monique@mbari.org

Do not hesitate to contact me if you cannot run the code in `start_sourcecode_ITS`, if you notice bugs, or if you need help implementing the code for a custom application.