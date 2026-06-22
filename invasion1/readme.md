### invasion1
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
For each turnover, the program first homogenizes the genetic background by runnning simulation in the absense of new mutations until only the focal sex-determining locus becomes polymorphic (preparing initial states). 
After that, an invasive sex-determining allele is repeatedly introduced, and subsequent allele frequency changes are recorded. 
For each turnover, 1,000 homogenized genetic backgrounds were generated, and the introduction of the invasive allele is repeated 10,000 times for each homogenized state, resulting in 10,000,000 invasion considered in total. 
