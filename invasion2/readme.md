### invasion2
#### Folder "c++"
c++ codes for simulation. 

An input file, "list_clean.txt", specifies the turnover events observed in the simulation with the default parameter. 
Each column contains the following information of each turnover:
col1: replicate ID in which the focal turnover was observed
col2: rand seed of the focal replicate
col3: derived sex-determining locus ID
col4: ancestral sex-determining locus ID
col5, 6: unused
col7: start of the turnover
col8: end of the turnover
col9: effect size of the derived sex-determining allele
col10: effect size of the ancestral sex-determining allele

Based on this, c++ program estimates the invasion analysis. 
For each turnover, the introduction of the invasive allele is "replayed" 10,000,000 times. 
