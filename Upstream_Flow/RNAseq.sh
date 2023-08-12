#RawData Fastqc
cd ./01.Rawdata/
for i in *fq.gz; do nohup fastqc $i -o fastqc/ & done

#Cutadptor
for i in 01.Rawdata/*_R1.fq.gz; do 
nohup trim_galore -j 5 --paired -q 20 --length 20 --trim-n -o 00.Cleandata/ $i ${i/_R1./_R2.} >> 00.Cleandata/trim.log & 
done

#Mapping
for i in 00.Cleandata/*R1_val_1.fq.gz ; do o=${i##*/};o=${o%%R1_val_1.fq.gz};
nohup hisat2 -p 5 --dta --no-mixed --no-discordant -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1_val_1/R2_val_2} -S ./11.Mapping/$o.sam > ./11.Mapping/log/$o.ht2.log &
done

for i in *.sam; do nohup samtools view -bf 0x2 $i | samtools sort -@ 5 - -o ../10.Mapping/${i/.sam/.sorted.bam} & done
for i in *.sorted.bam; do nohup samtools index $i & done

nohup samtools merge Bam_merge/WT.sorted.bam WT[123].sorted.bam &
nohup samtools merge Bam_merge/KO.sorted.bam KO[123].sorted.bam &

##Sam_QC_RSeQC
#~/local/anno/mm10/UCSC_RseQC_Reference_bed/mm10.HouseKeepingGenes.bed#来自/share/ann/useful_annotation/download_FTP/RSeQC/mm10.HouseKeepingGenes.bed.gz
cd ./10.Mapping
for i in *.sorted.bam; do nohup bam_stat.py -i $i >> ../12.Sam_QC/bam_stat/bam_stat.txt & done
for i in *.sorted.bam; do nohup junction_annotation.py -r /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm10.RefSeq.bed -i $i -o ../12.Sam_QC/junction_annotation/${i%%.*} & done
for i in *.sorted.bam; do nohup read_distribution.py -r ~/local/anno/mm10/UCSC_RseQC_Reference_bed/mm10.HouseKeepingGenes.bed -i $i > ../12.Sam_QC/read_distribution/${i%%.*}.txt & done

nohup geneBody_coverage.py -r ~/local/anno/mm10/UCSC_RseQC_Reference_bed/mm10.HouseKeepingGenes.bed -i ../../10.Mapping/ -o ./ >> ./gene_cover.txt &
nohup tin.py -r ~/local/anno/mm10/UCSC_RseQC_Reference_bed/mm10.HouseKeepingGenes.bed -i ../../10.Mapping &

#Sam_Cor
nohup multiBamSummary bins --binSize 1000 -b 10.Mapping/*.sorted.bam -o 13.Sam_Cor/all_bam.npz >> 13.Sam_Cor/all_bam.log &
nohup plotCorrelation -T RNA_correlation --plotNumbers --whatToPlot heatmap --skipZeros --corData 13.Sam_Cor/all_bam.npz -c pearson -p heatmap -o 13.Sam_Cor/RNA_correlation.pdf &
nohup plotPCA --transpose -in 13.Sam_Cor/all_bam.npz -o 13.Sam_Cor/RNA_PCA.pdf &


#Bigwig
for i in 10.Mapping/*.sorted.bam; do o=${i##*/};o=${o%%.*}; nohup bamCoverage -p 8 -b $i -o 14.Bigwig/${o}.bw --normalizeUsing RPKM --binSize 25 -ignore chrM > 14.Bigwig/log/${o}.bw.log & done


#featureCount
#产生count文件的基因名排序相同,可用md5sum检测
for i in ./10.Mapping/*.sorted.bam; do o=${i##*/}; o=${o%%.*}; featureCounts -T 6 -p -B -C -O -M -a /home/share/gtf/mm10.gtf -o ./20.Expression/featureCounts/$o.featureCount $i > ./20.Expression/featureCounts/log/$o.featureCounts.log & done
for i in *.featureCount; do awk 'NR>1 {print $7}' ${i} > ${i/.featureCount/.count.temp}; done
awk 'NR>1 {print $1}' WT1.featureCount > gene.temp
paste gene.temp *.count.temp > gene_count.table
rm gene.temp *.count.temp 


#Stringtie
#产生gene.tab的基因名排序不同,需处理
for i in ./10.Mapping/*.sorted.bam; do o=${i##*/}; o=${o%%.*}; nohup stringtie ${i} -p 5 -e -G /home/share/gtf/mm10.gtf -A ./20.Expression/stringtie/${o}_gene.tab -o ./20.Expression/stringtie/${o}_transcript.gtf & done
