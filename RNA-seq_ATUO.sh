#！/usr/bin/env bash
#Author: Bio_XDJ
#Created Time: 2021/09/20

#####pair end#####
# samples datas name is samplename_rep_read.fq.gz    eg: arp6_1_1.fq.gz, arp6_1_2.fq.gz, arp6_2_1.fq.gz, arp_6_2.fq.gz .

# First, you must create a sample.txt for tell the script your sample name  
# eg: arp6
#	  hira
#	  Col
#	  arphir


#### In the 5.DE_analysis , you must create a contrasts.txt for used to specify which sets of samples need to be compared

#### and In the 5.DE_analysis , you must create a sample_group.txt for tell the program which files are a set of samples

# INDEX is your workspeace of run this shell

#  REF is REFerence genome pathway


INDEX=/stor9000/apps/users/NWSUAF/2021051320/workspeace/data/RNA-seq/anp32e_nrpd
REF=/stor9000/apps/users/NWSUAF/2021051320/workspeace/ref
script=/stor9000/apps/users/NWSUAF/2021051320/workspeace/script

r=3

# 1.filter data

if [ ! -d ${INDEX}/1.filter ]
	then
		mkdir ${INDEX}/1.filter
fi

cd ${INDEX}/1.filter

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		fastp -w 4 -i ${INDEX}/raw/${sample}_R${rep}_1.fq.gz -o ${INDEX}/1.filter/${sample}_R${rep}_1.fq.gz -I ${INDEX}/raw/${sample}_R${rep}_2.fq.gz -O ${INDEX}/1.filter/${sample}_R${rep}_2.fq.gz -h ${INDEX}/1.filter/${sample}_R${rep}.html -j ${INDEX}/1.filter/${sample}_R${rep}.json
	done
done

# 2. mapping on genome

if [ ! -d ${INDEX}/2.mapping ]
	then
		mkdir ${INDEX}/2.mapping
fi

cd ${INDEX}/2.mapping

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		hisat2 --new-summary -x ${REF}/AT_genome -1 ${INDEX}/1.filter/${sample}_R${rep}_1.fq.gz -2 ${INDEX}/1.filter/${sample}_R${rep}_2.fq.gz -S ${INDEX}/2.mapping/${sample}_R${rep}.sam 1>${sample}_R${rep}.log 2>&1
	done
done

echo -e "Sample\toverall_alignment_rate" >> ${INDEX}/2.mapping/mapping_rate.txt
for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
        echo -e "${sample}_R${rep}\t$(tail -1 ${INDEX}/2.mapping/${sample}_R${rep}.log | awk '{print $1}')" >> ${INDEX}/2.mapping/mapping_rate.txt
	done
done

# 2.1 sam to bam && sort

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		samtools sort -o ${INDEX}/2.mapping/${sample}_R${rep}.bam ${INDEX}/2.mapping/${sample}_R${rep}.sam
		samtools index ${INDEX}/2.mapping/${sample}_R${rep}.bam
	done
done

rm ${INDEX}/2.mapping/*.sam


#3. Quantification

if [ ! -d ${INDEX}/3.quantification ]
	then
		mkdir ${INDEX}/3.quantification
fi

cd ${INDEX}/3.quantification

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		Rscript ${script}/run-featurecounts.R -b ${INDEX}/2.mapping/${sample}_R${rep}.bam -g ${REF}/Arabidopsis_thaliana.TAIR10.51.gtf -o ${INDEX}/3.quantification/${sample}_R${rep}
	done
done

#4. merge_result

if [ ! -d ${INDEX}/4.merge_result ]
	then
		mkdir ${INDEX}/4.merge_result
fi

cd ${INDEX}/4.merge_result

perl ${script}/merge/abundance_estimates_to_matrix.pl --est_method featureCounts ${INDEX}/3.quantification/*.count --out_prefix genes

#5. DE analysis

if [ ! -d ${INDEX}/5.DE_analysis ]
	then
		mkdir ${INDEX}/5.DE_analysis
fi

cd ${INDEX}/5.DE_analysis
a1(){
# create contrasts.txt

echo
read -p "请输入样品数： " SN
read -p "请输入重复数： " RN

for ((i=1;i<${SN}+1;i++))
do
        read -p "样品名： " A
        for ((j=1;j<${RN}+1;j++))
        do
                read -p "重复名： " B
                echo -e "${A}\t${B}" >> ${INDEX}/5.DE_analysis/contrasts.txt
        done
done
}

perl /stor9000/apps/users/NWSUAF/2021051320/workspeace/software/trinityrnaseq-Trinity-v2.13.2/Analysis/DifferentialExpression/run_DE_analysis.pl \
    --matrix ${INDEX}/4.merge_result/genes.counts.matrix \
    --method DESeq2 \
    --samples_file ${INDEX}/sample_group.txt \
    --contrasts ${INDEX}/contrasts.txt

cd ${INDEX}/5.DE_analysis/DESeq2*

for i in *.DE_results
do 
	echo -e "Gene_id\t$(cat $i)" > $i
done



