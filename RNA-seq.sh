#！/usr/bin/env bash
#Author: Bio_XDJ
#Created Time: 2021/09/20

#####pair end#####
# samples name is samplename_rep_read    eg:arp6_1_1

# First, you must create a sample.txt for tell the script your sample name  
# eg: arp6
#	  hira
#	  Col
#	  arphir

# In the workspeace you must create the following five folders:

# /raw    /1.filter   /2.mapping    /3.quantification   /4.merge_result   /5.DE_analysis

#### In the 3.quantification must have /script/run-featurecounts.R

#### In the 4.merge_result must have /script/abundance_estimates_to_matrix.pl and /script/support_scripts/run_TMM_scale_matrix.pl

#### In the 5.DE_analysis , you must create a contrasts.txt for used to specify which sets of samples need to be compared

#### and In the 5.DE_analysis , you must create a sample_group.txt for tell the program which files are a set of samples

# INDEX is your workspeace of run this shell

#  REF is REFerence genome pathway


INDEX=/stor9000/apps/users/NWSUAF/2021051320/workspeace
REF=/stor9000/apps/users/NWSUAF/2021051320/workspeace/ref/AT_genome

# 1.filter data

cd ${INDEX}/1.filter

for sample in $(cat ${INDEX}/sample.txt)
do
	for rep in 1 2 3
	do
		fastp -i ${INDEX}/raw/${sample}_${rep}_1.fq.gz -o ${INDEX}/1.filter/${sample}_${rep}_1.fq.gz -I ${INDEX}/raw/${sample}_${rep}_2.fq.gz -O ${INDEX}/1.filter/${sample}_${rep}_2.fq.gz -h ${INDEX}/1.filter/${sample}_${rep}.html -j ${INDEX}/1.filter/${sample}_${rep}.json
	done
done

# 2. mapping on genome

cd ${INDEX}/2.mapping

for sample in $(cat ${INDEX}/sample.txt)
do
	for rep in 1 2 3
	do
		hisat2 --new-summary -x ${REF} -1 ${INDEX}/1.filter/${sample}_${rep}_1.fq.gz -2 ${INDEX}/1.filter/${sample}_${rep}_2.fq.gz -S ${INDEX}/2.mapping/${sample}_${rep}.sam 1>${sample}_${rep}.log 2>&1
	done
done

# 2.1 sam to bam && sort

for sample in $(cat ${INDEX}/sample.txt)
do
	for rep in 1 2 3
	do
		samtools sort -o ${INDEX}/2.mapping/${sample}_${rep}.bam ${INDEX}/2.mapping/${sample}_${rep}.sam
	done
done

#3. Quantification

cd ${INDEX}/3.quantification

for sample in $(cat ${INDEX}/sample.txt)
do
	for rep in 1 2 3
	do
		Rscript ${INDEX}/3.quantification/script/run-featurecounts.R -b ${INDEX}/2.mapping/${sample}_${rep}.bam -g ${REF}/Arabidopsis_thaliana.TAIR10.51.gtf -o ${INDEX}/3.quantification/${sample}_${rep}
	done
done

#4. merge_result

cd ${INDEX}/4.merge_result

perl ${INDEX}/4.merge_result/script/abundance_estimates_to_matrix.pl --est_method featureCounts ${INDEX}/3.quantification/*.count --out_prefix genes

#5. DE analysis

cd ${INDEX}/5.DE_analysis

perl /stor9000/apps/users/NWSUAF/2021051320/workspeace/software/trinityrnaseq-Trinity-v2.13.2/Analysis/DifferentialExpression/run_DE_analysis.pl \
    --matrix ${INDEX}/4.merge_result/genes.counts.matrix \
    --method DESeq2 \
    --samples_file ${INDEX}/5.DE_analysis/sample_group.txt \
    --contrasts ${INDEX}/5.DE_analysis/contrasts.txt


