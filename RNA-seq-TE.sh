
INDEX=/stor9000/apps/users/NWSUAF/2021051320/workspeace/data/RNA-seq/anp32e_nrpd
REF=/stor9000/apps/users/NWSUAF/2021051320/workspeace/ref
script=/stor9000/apps/users/NWSUAF/2021051320/workspeace/script

r=3


#3. Quantification

if [ ! -d ${INDEX}/3.quantification/TE ]
	then
		mkdir -p ${INDEX}/3.quantification/TE
fi

cd ${INDEX}/3.quantification/TE

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		Rscript ${script}/run-featurecounts_TE.R -b ${INDEX}/2.mapping/${sample}_R${rep}.bam -g ${REF}/Araport11_GTF_genes_transposons.gtf -o ${INDEX}/3.quantification/TE/${sample}_TE_R${rep}
	done
done

#4. merge_result

if [ ! -d ${INDEX}/4.merge_result/TE ]
	then
		mkdir -p ${INDEX}/4.merge_result/TE
fi

cd ${INDEX}/4.merge_result/TE

perl ${script}/merge/abundance_estimates_to_matrix.pl --est_method featureCounts ${INDEX}/3.quantification/TE/*.count --out_prefix TEs
a2(){
#5. DE analysis

if [ ! -d ${INDEX}/5.DE_analysis/TE ]
	then
		mkdir -p ${INDEX}/5.DE_analysis/TE
fi

cd ${INDEX}/5.DE_analysis/TE
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
    --matrix ${INDEX}/4.merge_result/TE/TEs.counts.matrix \
    --method DESeq2 \
    --samples_file ${INDEX}/sample_group.txt \
    --contrasts ${INDEX}/contrasts.txt
}


