import gzip
import os, re
import pandas as pd

folder = "max_div"

output = open(os.path.join(folder, "outputs.txt"), mode="w")
print("index\tt1\tt2\tt3\tW1\tW2\tW3\tW4\told\tnew\tal1\tal2\tintro\tpgd", file=output)

folders = [f for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]

numbers = [int(re.match(r'folder(\d+)', name).group(1)) for name in folders if re.match(r'folder(\d+)', name)]

count = 0

for i in numbers:
  if os.path.exists(os.path.join(folder, "folder"+str(i), "regi_time.txt")):
    with open(os.path.join(folder, "folder"+str(i), "regi_time.txt"), mode="r") as f:
      retlist = []
      for line in f:
        retlist.append(line.rstrip("\n"))

    hetero = -1
    ret = ["NA"] * 14

    if retlist[0] == "NA":
      print("NA run: ", i)
    else:
      ret[0] = str(i)
      ret[1] = retlist[0]
      ret[3] = retlist[3]
      ret[8] = retlist[2]
  else:
    print("error: ", i)
    exit()

  if retlist[0] != "NA":
    with gzip.open(os.path.join(folder, "folder"+str(i), "fitness.txt.gz"), mode="rt") as f:
      for line in f:
        line = line.rstrip("\n")
        tmplist = line.split("\t")

        if int(tmplist[0]) == int(retlist[0]) + 500001:
          if last_time != int(retlist[0]) + 500000:
            print("Error in replicate", i)

          if float(tmplist[8]) > float(tmplist[9]):
            hetero = "male"
            ret[4] = male_fit
            ret[5] = female_fit
          else:
            hetero = "female"
            ret[4] = female_fit
            ret[5] = male_fit
        
        if retlist[1] != "NA" and int(tmplist[0]) == int(retlist[1]):
          if hetero == "male":
            ret[6] = tmplist[2]
            ret[7] = tmplist[3]
          else:
            ret[6] = tmplist[3]
            ret[7] = tmplist[2]
          ret[2] = retlist[1]
        
        male_fit = tmplist[2]
        female_fit = tmplist[3]
        last_time = int(tmplist[0])
    
    with gzip.open(os.path.join(folder, "folder"+str(i), "log_association.txt.gz"), mode="rt") as f:
      for line in f:
        line = line.rstrip("\n")
        tmplist = line.split("\t")

        if retlist[1] != "NA" and int(tmplist[0]) > int(retlist[1]):
          if int(tmplist[0]) != int(retlist[1]) + 1:
            print("Error in replicate", i)

          smallest_val = 1.0
          smallest_index = -1
          for j in range(1, len(tmplist)):
            if float(tmplist[j]) < smallest_val:
              smallest_val = float(tmplist[j])
              smallest_index = j - 1
          
          ret[9] = str(smallest_index)
          break

    path = os.path.join("../result_degenerate_rerun", folder, "folder"+str(i), "output_freq.txt")
    col_names = ["gen", "idx", "pval", "x", "ori", "freq"]

    df = pd.read_csv(path, sep=r"\s+", header=None, names=col_names, engine="python")
    gens = pd.to_numeric(df["gen"], errors="coerce").dropna().astype(int).unique()

    top3 = pd.Series(gens).nsmallest(3).to_list()
    introduced, purged = top3[1], top3[2]

    ret[12] = str(introduced)
    ret[13] = str(purged)

    df_gen = df[df["gen"] == purged]
    df_min = df_gen[df_gen["pval"] == df_gen["pval"].min()]
    df_unique = df_min.drop_duplicates()
    df_sorted = df_unique.sort_values(by="freq", ascending=False)

    ret[10] = str(df_sorted["ori"].head(2).max())
    ret[11] = str(df_sorted["ori"].head(2).min())

    print("\t".join(ret), file=output)
    count += 1

    if int(ret[12]) - int(ret[1]) - 500000 != 1:
      print("error: ", i)
    if int(ret[13]) - int(ret[2]) != 1:
      print("error: ", i)

print(count)
            
        
          
        
