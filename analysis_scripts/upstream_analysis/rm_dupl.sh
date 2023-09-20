#!/bin/bash
# Create function
rm_dup(){
	fastq2collapse.pl ${f}_trimmed.fq.gz - | gzip -c > ${f}_trimmed.c.fq.gz
	stripBarcode.pl -format fastq -len 8 ${f}_trimmed.c.fq.gz - | gzip -c > ${f}_trimmed.c.tag.fq.gz
	 echo "$f done"
    }

for f in `cat ribo_samples`;
do
	rm_dup $f &
done

wait

echo "All done"

# check the read counts of the collapsed files
for f in `cat ribo_samples`;
do
  c=`zcat ${f}_trimmed.c.fq.gz | wc -l`
  c=$((c/4))
  echo -e "${f}_trimmed.c \t $c"
done >> ./ribo_collapsed.txt



