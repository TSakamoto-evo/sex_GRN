### degeneration
#### Folder "c++"
c++ codes for simulation. 
In our paper, 1,000 replicates were run for each parameter set. 
When some runs failed to evolve established sex, we reran that replicate. 

#### Folder others
R and python code for visualization.

`plot_association.R`:
Plot pvalues of Kendall's tau-b test. It was used to check the dynamics. 

`extract_for_fig.py`:
Make a list that contains summarized information. 

Explanation for each column of the output file:
index   : replicate index  
t1      : time when the two sexes are established  
t2      : time when the deleterious mutation has been removed  
t3      : time when the sex is re-established after the removal of the deleterious allele  
W1      : fitness of the ancestrally heterogametic sex in the ancestral system  
W2      : fitness of the ancestrally homogametic sex in the ancestral system  
W3      : fitness of the heterogametic sex in the new system  
W4      : fitness of the homogametic sex in the new system  
old     : locus index of the ancestral sex-determining locus  
new     : locus index of the sex-determining locus in the new system  
al1     : origin of the new (younger) sex-determining allele  
al2     : origin of the old (older) sex-determining allele  
intro   : generation when the deleterious allele was introduced (t1 + 500001)  
purged  : generation when the deleterious allele was purged (t2 + 1)  

`plot_fig.R`:
Code to make Figure 6 and Extended Data Figure 6.

`plot_fig_all.R`:
Code to make Extended Data Figure 7. 
