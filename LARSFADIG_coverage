
#1. Identify reads containing LARSFDIG motif 

for file in *.fastq
do
filename=$(echo "$file" | cut -d"_" -f1)
sed -n '1~4s/^@/>/p;2~4p' "$file" > "$file".fasta 
done 


for file in *_1.fastq.fasta
do
blastn -query "$file" -db published_var_larsfagid_sequences.fasta -out larsfadig/blastn_results_fasta_"$file"_publars -evalue 1e-10 -task blastn-short -num_threads 6 -max_target_seqs 1 -outfmt 6
done

for file in *_1.fastq.fasta
do
awk '{ print $1 }' larsfadig/blastn_results_fasta_"$file"_publars > larsfadig/blastn_results_fasta_"$file"_publars_read_ids
done

for file in *_1.fastq.fasta
do
filename=$(echo "$file" | cut -d"_" -f1)
seqtk subseq "$filename"_1.fastq.fasta larsfadig/blastn_results_fasta_"$file"_publars_read_ids > larsfadig/"$filename"_1_larsfadigreads.fasta
seqtk subseq "$filename"_2.fastq.fasta larsfadig/blastn_results_fasta_"$file"_publars_read_ids > larsfadig/"$filename"_2_larsfadigreads.fasta
done



#2. De novo assemble LASRFADIG motif using rnaSpades 

for file in *_1.fastq.fasta
do
filename=$(echo "$file" | cut -d"_" -f1)
rnaspades -o larsfadig/"$sample_name"_LARSFADIG -k 25 -1 larsfadig/"$filename"_1_larsfadigreads.fasta -2 larsfadig/"$filename"_2_larsfadigreads.fasta
done

#3. Determine coverage over S
for file in *_1.fastq.fasta
do
filename=$(echo "$file" | cut -d"_" -f1)
bwa index larsfadig/"$sample_name"_LARSFADIG/transcripts.fasta 
bwa mem -k 31 -a larsfadig/"$sample_name"_LARSFADIG/transcripts.fasta  larsfadig/"$filename"_1_larsfadigreads.fasta larsfadig/"$filename"_1_larsfadigreads.fasta > larsfadig/"$filename"_transcripts_lars_read_alignment.sam
samtools view -b larsfadig/"$filename"_transcripts_lars_read_alignment.sam > larsfadig/"$filename"_transcripts_lars_read_alignment.bam
samtools sort larsfadig/"$filename"_transcripts_lars_read_alignment.bam > larsfadig/"$filename"_transcripts_lars_read_alignment_sort.bam
samtools depth larsfadig/"$filename"_transcripts_lars_read_alignment_sort.bam > larsfadig/"$filename"_transcripts_lars_coverage
samtools stats larsfadig/"$filename"_transcripts_lars_read_alignment_sort.bam > larsfadig/"$filename"_transcripts_lars_stats 
done 
