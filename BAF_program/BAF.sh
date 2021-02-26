#!/bin/bash



raw_reads_folder=/home/reads_folder/
sample_name=sample
results_path=../results/${sample_name}
transporters=clear_all_transporters_nt_sequences
decarboxylases=clear_all_decarboxylases_nt_sequences
ncbi_transporters=clear_all_ncbi_transporters_nt_sequences
trans_reference_path=reference_genomes/V2/all_transporters
ncbi_transporter_reference_path=reference_genomes/V2/all_ncbi_transporters
decarboxy_reference_path=reference_genomes/V2/all_decarboxylases



#~ module add UHTS/Aligner/bowtie2/2.3.4.1
#~ module add UHTS/Analysis/samtools/1.8
#~ module add UHTS/Analysis/trimmomatic/0.36
#~ module add UHTS/Quality_control/fastqc/0.11.7
#~ module add UHTS/Assembler/SPAdes/3.10.1
#~ module add Blast/ncbi-blast/2.10.1+



mkdir -p ${results_path}/alignment_to_transporters

#trimmming
mkdir ${results_path}/trimmed_reads
trimmomatic PE -threads 2 -basein  ${raw_reads_folder}*R1_001.fastq.gz -baseout  ${results_path}/trimmed_reads/${sample_name}_trim.fastq.gz SLIDINGWINDOW:4:8 MINLEN:127
fastqc ${results_path}/trimmed_reads/${sample_name}*.fastq.gz

#~ #align to transporter sequences
mkdir ${results_path}/alignment_to_transporters/alignments
mkdir ${results_path}/alignment_to_transporters/aligned_reads
mkdir ${results_path}/alignment_to_transporters/assembly

bowtie2 --local --no-unal -D 100 -x ${trans_reference_path}/${transporters} -1 ${results_path}/trimmed_reads/${sample_name}_trim_1P.fastq.gz -2 ${results_path}/trimmed_reads/${sample_name}_trim_2P.fastq.gz -U ${results_path}/trimmed_reads/${sample_name}_trim_1U.fastq.gz,${results_path}/trimmed_reads/${sample_name}_trim_2U.fastq.gz -S ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}.sam --al-conc ${results_path}/alignment_to_transporters/aligned_reads/paired.fastq --al ${results_path}/alignment_to_transporters/aligned_reads/unpaired.fastq
samtools view -bS ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}.sam > ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}.bam
samtools sort -T /tmp/aln.sorted -o ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}_sorted.bam ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}.bam
samtools index ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}_sorted.bam ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}_sorted.bai
rm ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}.sam
rm ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}.bam

#~ #use unpaired only if available
number_of_unpaired="$(wc -l ${results_path}/alignment_to_transporters/aligned_reads/unpaired.fastq)"
printf "${number_of_unpaired} unpaired_reads\n"
if [ ${number_of_unpaired} -gt 0 ]
then
	spades.py --careful --only-assembler -1 ${results_path}/alignment_to_transporters/aligned_reads/paired.1.fastq -2 ${results_path}/alignment_to_transporters/aligned_reads/paired.2.fastq -s ${results_path}/alignment_to_transporters/aligned_reads/unpaired.fastq -t 2 -m 16 -o ${results_path}/alignment_to_transporters/assembly
else
		spades.py --careful --only-assembler -1 ${results_path}/alignment_to_transporters/aligned_reads/paired.1.fastq -2 ${results_path}/alignment_to_transporters/aligned_reads/paired.2.fastq -t 2 -m 16 -o ${results_path}/alignment_to_transporters/assembly
fi

mkdir ${results_path}/alignment_to_transporters/assembly/bowtie2_reference
bowtie2-build ${results_path}/alignment_to_transporters/assembly/contigs.fasta ${results_path}/alignment_to_transporters/assembly/bowtie2_reference/assembly_ref
#mkdir ${results_path}/${transporters}_aligned_reads/alignment_to_assembly
bowtie2 --local -x ${results_path}/alignment_to_transporters/assembly/bowtie2_reference/assembly_ref -1 ${results_path}/alignment_to_transporters/aligned_reads/paired.1.fastq -2 ${results_path}/alignment_to_transporters/aligned_reads/paired.2.fastq -U ${results_path}/alignment_to_transporters/aligned_reads/unpaired.fastq -S ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly.sam
samtools view -bS ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly.sam > ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly.bam
samtools sort -T /tmp/aln.sorted -o ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly_sorted.bam ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly.bam
samtools index ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly_sorted.bam ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly_sorted.bai
rm ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly.bam
rm ${results_path}/alignment_to_transporters/alignments/transporter_reads_to_assembly.sam

mkdir ${results_path}/alignment_to_transporters/assembly/blastn_results/
blastn -import_search_strategy blast_settings/WPWMU1G2013_search_strategy.asn -query ${results_path}/alignment_to_transporters/assembly/contigs.fasta -out ${results_path}/alignment_to_transporters/assembly/blastn_results/results.out -remote 




#align to decarboxylase sequences

mkdir -p ${results_path}/alignment_to_decarboxylases/alignments
mkdir ${results_path}/alignment_to_decarboxylases/aligned_reads
mkdir ${results_path}/alignment_to_decarboxylases/assembly

bowtie2 --local --no-unal -D 100  -x ${decarboxy_reference_path}/${decarboxylases} -1 ${results_path}/trimmed_reads/${sample_name}_trim_1P.fastq.gz -2 ${results_path}/trimmed_reads/${sample_name}_trim_2P.fastq.gz -U ${results_path}/trimmed_reads/${sample_name}_trim_1U.fastq.gz,${results_path}/trimmed_reads/${sample_name}_trim_2U.fastq.gz -S ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}.sam --al-conc ${results_path}/alignment_to_decarboxylases/aligned_reads/paired.fastq --al ${results_path}/alignment_to_decarboxylases/aligned_reads/unpaired.fastq
samtools view -bS ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}.sam > ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}.bam
samtools sort -T /tmp/aln.sorted -o ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}_sorted.bam ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}.bam
samtools index ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}_sorted.bam ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}_sorted.bai
rm ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}.sam
rm ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}.bam

number_of_unpaired="$(wc -l ${results_path}/alignment_to_decarboxylases/aligned_reads/unpaired.fastq)"
	printf "${number_of_unpaired} unpaired_reads\n"
if [ ${number_of_unpaired} -gt 0 ]
then
	spades.py --careful --only-assembler -1 ${results_path}/alignment_to_decarboxylases/aligned_reads/paired.1.fastq -2 ${results_path}/alignment_to_decarboxylases/aligned_reads/paired.2.fastq -s ${results_path}/alignment_to_decarboxylases/aligned_reads/unpaired.fastq -t 2 -m 16 -o ${results_path}/alignment_to_decarboxylases/assembly
else
	spades.py --careful --only-assembler -1 ${results_path}/alignment_to_decarboxylases/aligned_reads/paired.1.fastq -2 ${results_path}/alignment_to_decarboxylases/aligned_reads/paired.2.fastq -t 2 -m 16 -o ${results_path}/alignment_to_decarboxylases/assembly
fi
mkdir ${results_path}/alignment_to_decarboxylases/assembly/bowtie2_reference
bowtie2-build ${results_path}/alignment_to_decarboxylases/assembly/contigs.fasta ${results_path}/alignment_to_decarboxylases/assembly/bowtie2_reference/assembly_ref
#mkdir ${results_path}/${decarboxylases}_aligned_reads/alignment_to_assembly
bowtie2 --local -x ${results_path}/alignment_to_decarboxylases/assembly/bowtie2_reference/assembly_ref -1 ${results_path}/alignment_to_decarboxylases/aligned_reads/paired.1.fastq -2 ${results_path}/alignment_to_decarboxylases/aligned_reads/paired.2.fastq -U ${results_path}/alignment_to_decarboxylases/aligned_reads/unpaired.fastq -S ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly.sam
samtools view -bS ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly.sam > ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly.bam
samtools sort -T /tmp/aln.sorted -o ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly_sorted.bam ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly.bam
samtools index ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly_sorted.bam ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly_sorted.bai
rm ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly.bam
rm ${results_path}/alignment_to_decarboxylases/alignments/decarboxylases_reads_to_assembly.sam

mkdir ${results_path}/alignment_to_decarboxylases/assembly/blastn_results/
blastn -import_search_strategy blast_settings/WPWMU1G2013_search_strategy.asn -query ${results_path}/alignment_to_decarboxylases/assembly/contigs.fasta -out ${results_path}/alignment_to_decarboxylases/assembly/blastn_results/results.out -remote 




#count reads to dec/trans sequences
python3 python_scripts/bam_read_counterV2.py ${trans_reference_path}/${transporters}.fasta ${results_path}/alignment_to_transporters/alignments/all_reads_${sample_name}_${transporters}_sorted.bam ${results_path}/alignment_to_transporters/alignments/counts_${sample_name}_reads_to_${transporters}.json
python3 python_scripts/bam_read_counterV2.py ${decarboxy_reference_path}/${decarboxylases}.fasta ${results_path}/alignment_to_decarboxylases/alignments/all_reads_${sample_name}_${decarboxylases}_sorted.bam ${results_path}/alignment_to_decarboxylases/alignments/counts_${sample_name}_reads_to_${decarboxylases}.json
