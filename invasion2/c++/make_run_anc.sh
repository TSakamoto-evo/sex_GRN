#!/bin/bash
#$ -S /bin/bash
#$ -pe def_slot 1
#$ -l s_vmem=2G,mjob
#$ -t 1-692
#$ -cwd
#$ -j y
#$ -o ./error.txt

FILE=./list_clean.txt
line=0

while read LINE
do
  line=$((line+1));
  if [ ${line} -eq ${SGE_TASK_ID} ]; then
    seed=$(echo "$LINE" | cut -d' ' -f2)
    der=$(echo "$LINE" | cut -d' ' -f3)
    anc=$(echo "$LINE" | cut -d' ' -f4)
    anct=$(echo "$LINE" | cut -d' ' -f7)
    dert=$(echo "$LINE" | cut -d' ' -f8)
    dereff=$(echo "$LINE" | cut -d' ' -f9)
    anceff=$(echo "$LINE" | cut -d' ' -f10)
  else
    :
  fi
done < ${FILE}

mkdir -p ./folder_anc${SGE_TASK_ID}

cp ./test.out ./folder_anc${SGE_TASK_ID}
cd ./folder_anc${SGE_TASK_ID}
./test.out ${seed} ${anct} ${anc} ${der} ${dereff} 

gzip -9 fitness.txt
gzip -9 freq.txt
gzip -9 outcome.txt
gzip -9 output_freq.txt

cd ..

exit 0
