#!/bin/bash -l

#module load devel/miniforge
#conda activate utils

#### sratools need to be installed in one way or other, e.g. via conda

vdb-config --prefetch-to-cwd

for GSM in "$@"
do
  SRR=$(wget "https://www.ncbi.nlm.nih.gov/sra/?term=$GSM" -q -O - | egrep -o "SRR\w+" - | uniq)

  echo ">>> dumping $GSM ($SRR)"

  mkdir $GSM

  rm -f "$GSM/$GSM.fastq.gz"
  rm -f "$GSM/"$GSM"_1.fastq.gz"
  rm -f "$GSM/"$GSM"_2.fastq.gz"

  for i in $SRR
  do
    prefetch $i
    fasterq-dump --outdir "$GSM/" --split-files $i
    rm -dfr $SRR
  done


  #paired-end
  if [ -f $GSM/SRR*_1.fastq ]
  then
    cat $GSM/SRR*_1.fastq | pigz -p 3 --fast - > "$GSM/"$GSM"_1.fastq.gz"
    cat $GSM/SRR*_2.fastq | pigz -p 3 --fast - > "$GSM/"$GSM"_2.fastq.gz"

    rm -f $GSM/SRR*_1.fastq
    rm -f $GSM/SRR*_2.fastq
  fi

 #single-end
  if [ -f $GSM/SRR*.fastq ]
  then
    cat $GSM/SRR*.fastq | pigz -p 3 --fast - > "$GSM/$GSM.fastq.gz"
    rm -f $GSM/SRR*.fastq
  fi

done

#conda deactivate
