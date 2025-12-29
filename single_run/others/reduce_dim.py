import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

folder = "../c++"

alpha = 0.2
tau_d = 40
num_g = 8

phenotype = []
stage = []
time = []

male = [10, 10, 0, 0, 10, 10, 0, 0]
female = [10, 10, 0, 0, 0, 0, 10, 10]
sig = 5.0

count = 0

with open(folder + "/regi_genotype.txt", mode="r") as f:
  for line in f:
    line = line.rstrip("\n")
    tmp_list = line.split("\t")

    g1 = np.array([float(tmp_list[i]) for i in range(5, 5 + num_g)])
    g2 = np.array([float(tmp_list[i]) for i in range(5 + num_g, 5 + 2 * num_g)])

    b1 = np.array([[float(tmp_list[5 + 2 * num_g + i * num_g + j]) for j in range(num_g)] for i in range(num_g)])
    b2 = np.array([[float(tmp_list[5 + 2 * num_g + num_g ** 2 + i * num_g + j]) for j in range(num_g)] for i in range(num_g)])

    p = (g1 + g2)

    for j in range(int(tmp_list[2])):
      phenotype.append(p)
      stage.append(0)
      time.append(int(tmp_list[0]))

    for i in range(tau_d):
      p = (1.0 - alpha) * p + (0.5 + 0.5 * np.tanh(b1 @ p)) + (0.5 + 0.5 * np.tanh(b2 @ p))

      for j in range(int(tmp_list[2])):
        phenotype.append(p)
        stage.append(i + 1)
        time.append(int(tmp_list[0]))

      if i == 20:
        fit_m = 1.0
        fit_f = 1.0
      if i >= 20:
        dist_sq_m = 0.0
        dist_sq_f = 0.0

        for j in range(num_g):
          dist_sq_m += (male[j] - p[j])**2
          dist_sq_f += (female[j] - p[j])**2
        
        fit_m *= np.exp(-dist_sq_m / 2.0 / sig**2)
        fit_f *= np.exp(-dist_sq_f / 2.0 / sig**2)
    
    fit_m = fit_m ** (1.0 / 20)
    fit_f = fit_f ** (1.0 / 20)

    if count < 100:
      print(count, fit_m, fit_f)
    count += 1

pca = PCA(n_components=5)
X_pca = pca.fit_transform(phenotype)

with open(folder + "/pca_vec.txt", mode="w") as f:
  tmp = [str(j) for j in pca.explained_variance_ratio_]
  print(" ".join(tmp), file=f)

  for i in pca.components_:
    tmp = [str(j) for j in i]
    print(" ".join(tmp), file=f)

with open(folder + "/pca_res.txt", mode="w") as f:
  for i in range(len(X_pca)):
    tmp = [str(round(j, 2)) for j in X_pca[i]]

    if i == 0 or (i > 0 and stage[i] != stage[i-1]):
      print(time[i], stage[i], " ".join(tmp), sep=" ", file=f)

exit()