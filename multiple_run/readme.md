### multiple_run
#### Folder "c++"
c++ codes for simulation. 
In our paper, 1,000 replicates were run for each parameter set. 

#### Folder others
R and python code for visualization.

`plot.R`:
Plot the proportion of male-like and female-like individuals (upper panels in Extended Data Figure 4).

`plot_association.R`:
Plot pvalues of Kendall's tau-b test (lower panels in Extended Data Figure 4).

`plot_sd_dist_at_t.R`:
Plot the distribution of a sex-determining locus in GRN (Figure 3, Extended Data Figure 3). 

`turnover_auto_v3.py`:
Count the number of turnover events (total number, number with a change of heterogameti sex). 

`plot_all_v2.R`:
Plot the proportion of successful establishment of sexes (Pa), distribution of time until the establishment (T), and the proportion of runs with turnovers (Pt) (Extended Data Figure 1). It uses the output of `turnover_auto_v3.py` in addition to the simulation outputs. 

`make_heatmap.py`:
The same function as `plot_association.R`, but computationally faster. 
It was used for manual inspection of each run. 


#### Folder rerun
Codes that was used to rerun some runs with interesting dynamics including turnovers. 
Simulation was rerun using the seed information, and output a bit more detailed information. 
In particular, genotype frequencies around the focal turnover were recorded.
This information is used to plot the frequency change in Figure 4 and Extended Data Figure 5.
See R script(`plot_freq.R`) for this plot.  
