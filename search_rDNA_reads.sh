rm rDNA_containing_reads.fastq
for file in `ls /home/yutaro/data/cliveome/*.fastq`; do
    echo ${file}
    python search_rDNA_reads.py ${file}
done
