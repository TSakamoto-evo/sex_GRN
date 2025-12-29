import os, gzip
import subprocess

folder = "max_div"
print(folder)

output = open(folder + "/turnovers_list.txt", mode="w")

for i in range(1, 1001):
  if os.path.isfile(folder + "/folder" + str(i) + "/regi_time.txt"):
    # successful evolution
    with open(folder + "/folder" + str(i) + "/regi_time.txt", mode = "r") as f:
      for line in f:
        line = line.rstrip("\n")
        wait_time = int(line)

    with gzip.open(folder + "/folder" + str(i) + "/log_association.txt.gz", mode="rt") as f:
      sd_list = {}

      for line in f:
        line = line.rstrip("\n")
        tmplist = line.split()
        timepoint = int(tmplist[0])

        # check which gene shows lowest p-values (min_index)
        min_index = 0
        min_p = 2
        for j in range(1, len(tmplist)):
          if min_p > float(tmplist[j]):
            min_index = j - 1
            min_p = float(tmplist[j])

        if timepoint > wait_time:
          sd_list[timepoint] = min_index

    now_time = -1
    time_list = []
    sd_locus_list = {}
    now_state_list = {}
    prev_sd_list = {}

    prev1 = -1
    prev2 = -1

    with gzip.open(folder + "/folder" + str(i) + "/regi_binned.txt.gz", mode="rt") as f:
      for line in f:
        line = line.rstrip("\n")
        tmplist = line.split("\t")
        timepoint = int(tmplist[0])

        if now_time != timepoint:
          now_time = timepoint
          now_state = "NA"

        if timepoint > wait_time:
          if len(time_list) == 0 or time_list[-1] != timepoint:
            time_list.append(timepoint)

          # check the sd locus at that time
          sd_locus = sd_list[timepoint]

          if int(tmplist[1]) == sd_locus:
            # allele count in each sex
            male_count = int(tmplist[5]) + int(tmplist[6])
            female_count = int(tmplist[7]) + int(tmplist[8])

            # if that allele is retained almost exclusively by males, it indicates XY sex determination
            if male_count > 400 and female_count == 0:
              if sd_locus != prev1:
                prev2 = prev1
                prev1 = sd_locus

              if now_state == "ZW":
                now_state = "Both"
              else:
                now_state = "XY"
            elif female_count > 400 and male_count == 0:
              if sd_locus != prev1:
                prev2 = prev1
                prev1 = sd_locus

              if now_state == "XY":
                now_state = "Both"
              else:
                now_state = "ZW"
          
          sd_locus_list[timepoint] = sd_locus
          now_state_list[timepoint] = now_state

          if sd_locus != prev1:
            prev_sd_list[timepoint] = prev1
          else:
            prev_sd_list[timepoint] = prev2
    
    now_time = -1
    prev_sd_state = {}
    with gzip.open(folder + "/folder" + str(i) + "/regi_binned.txt.gz", mode="rt") as f:
      for line in f:
        line = line.rstrip("\n")
        tmplist = line.split("\t")
        timepoint = int(tmplist[0])

        if now_time != timepoint:
          now_time = timepoint
          now_line = 0

        if timepoint > wait_time:
          # check the sd locus at that time
          prev_sd_locus = prev_sd_list[timepoint]

          if int(tmplist[1]) == prev_sd_locus:
            now_line += 1
            
          prev_sd_state[timepoint] = now_line

              
    no_to = 0
    no_het = 0

    sd_locus = sd_locus_list[time_list[0]]
    now_state = now_state_list[time_list[0]]

    for j in time_list:
      #print(j, sd_locus_list[j], now_state_list[j], prev_sd_list[j], prev_sd_state[j])

      if sd_locus_list[j] != sd_locus and prev_sd_state[j] == 1 and now_state_list[j] != "NA":
        no_to += 1
        sd_locus = sd_locus_list[j]

        if now_state_list[j] != now_state:
          no_het += 1
      
      if now_state_list[j] != "NA":
        now_state = now_state_list[j]
      
      if now_state_list[j] == "Both":
        print("There is both.")
    
    print(i, no_to, no_het, file=output, flush=True)


output.close()
exit()
