import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
import gzip
import math

folder = "max_div"

subprocess.run("mkdir -p " + folder + "/heatmap", shell=True)
subprocess.run("mkdir -p " + folder + "/ranks", shell=True)

subprocess.run("rm " + folder + "/heatmap/*", shell=True)
subprocess.run("rm " + folder + "/ranks/*", shell=True)

max_bound = 0
min_bound = -200

for i in range(1, 1001):
  with gzip.open(folder + "/folder" + str(i) + "/log_association.txt.gz", mode="rt") as f:
    time_list = []
    first_list = []
    second_list = []
    third_list = []

    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split("\t")

      time = int(tmplist[0])
      log_ps = tmplist[1:]
      log_ps = [float(j) * np.log10(np.exp(1)) for j in log_ps]

      time_list.append(time)
      
      l_sorted = sorted(log_ps)
      first_list.append(l_sorted[0])
      second_list.append(l_sorted[1])
      third_list.append(l_sorted[2])

      if len(time_list) == 1:
        colname = [j for j in range(len(tmplist))]
        df = pd.DataFrame(log_ps, columns = [time])
        df = df.T
      else:
        df.loc[time] = log_ps

  df = df.T

  plt.figure(figsize=(30, 9))
  sns.heatmap(df, cmap='Reds', vmin=min_bound, vmax=max_bound)
  plt.xlabel("generation")
  plt.ylabel("locus")
  plt.title(str(i))

  plt.tick_params(labelbottom=True, labelleft=False, labelright=False, labeltop=False, 
    bottom=True, left=False, right=False, top=False)

  plt.savefig(folder + "/heatmap/folder" + str(i) + ".png", dpi=300)
  plt.clf()
  plt.close()


  plt.figure(figsize=(5, 5))
  plt.scatter(time_list, first_list, s=0.1)
  plt.scatter(time_list, second_list, s=0.1)
  plt.scatter(time_list, third_list, s=0.1)

  plt.xlabel("generation")
  plt.ylabel("-log_ps")
  plt.title(str(i))

  plt.tick_params(labelbottom=True, labelleft=True, labelright=False, labeltop=False, 
    bottom=True, left=True, right=False, top=False)

  plt.savefig(folder + "/ranks/folder" + str(i) + ".png", dpi=300)
  plt.clf()
  plt.close()