# This is a snakefile based on python
# For PADD-seq
# written by YJY, modified by Yongchang Zhu
# 2023-05-02


import os
pwd = os.getcwd()
TRIM_PATH = pwd + '/trim'

# loadConfig
include: "conf/config"

# load Samples-prefix
SAMPLES=[line.rstrip() for line in open("conf/sm.list")]


# main workflow

rule all:
	input:
		expand("{sample}.done",sample=SAMPLES)


rule cutadapt:
	input:
		"fastq/Sample_{sample}/{sample}_combined_R1.fastq.gz",
		"fastq/Sample_{sample}/{sample}_combined_R2.fastq.gz"
	output:
		cut1=temp("temp/{sample}_R1.fq"),
		cut2=temp("temp/{sample}_R2.fq")
	message:"cut adapter in 5'-end using cutadapt"
	log:
		"logs/trim/{sample}_cut.log"
	threads:39
	params:
		pair_mode='any'
	priority: 0
	shell:
		"cutadapt -g {ADAPT1} --discard-trimmed --pair-filter={params.pair_mode} -j {threads} -o {output.cut1} -p {output.cut2} {input} 1>{log} 2>&1"


rule trim_galore:
	input:
		rules.cutadapt.output
	output:
		temp("trim/{sample}_R1_val_1.fq.gz"),
		temp("trim/{sample}_R2_val_2.fq.gz")
	message:"trim remains adapter using trim_galore"
	log:
		"logs/trim/{sample}_trim.log"
	threads:39
	priority: 1
	shell:
		"trim_galore {input} --phred33 --gzip  -j {threads} --paired -o {TRIM_PATH}  1>{log} 2>&1"


rule bwa_map:
	input:
		rules.trim_galore.output
	output:
		temp("mapped/{sample}.bam")
	message:"mapping reads to genome using bwa"
	log:
		"logs/bwa/{sample}_map.log"
	threads:39
	params:39
	priority: 2
	shell:
		"bwa mem -t {threads} {hg38_INDEX_PATH} {input} | samtools view -Sb - | samtools sort -@ {params} -o {output} - "


rule markdup:
	input:
		rules.bwa_map.output
	output:
		temp("mapped/{sample}.markdup.bam")
	params:39
	priority: 3
	shell:
		"sambamba markdup -t {params} {input} {output}"


rule index:
	input:
		rules.markdup.output
	output:
		temp("mapped/{sample}.markdup.bam.bai")
	priority: 4
	shell:
		"samtools index {input}"



rule generate_bed:
	input:
		bam=rules.markdup.output,
		bai=rules.index.output
	output:
		bed=temp("{sample}.bed"),
		uniq="{sample}.uniqCount"
	message:"generate bed containing filtered reads info from bam"
	priority: 5
	shell:
		"python {SCRIPTS}/generate_bed.py --ibam {input.bam} --fasta {hg38_FA} --bed {output.bed} --uniq {output.uniq}"


rule selectBase:
        input:
                rules.generate_bed.output.bed
        output:
                "{sample}.selectBase.bed"
        message:"select TT/TC/CT/CC"
        priority: 6
	shell:
                "grep -E 'TT|TC|CT|CC' {input} > {output}"

rule damageSum:
	input:
		rules.selectBase.output
	output:
		"{sample}.damageSum"
	message:"sum damage sites"
	priority: 7
	shell:
		"cat {input} | wc --lines > {output}"

rule bedToBam:
	input:
		rules.selectBase.output
	output:
		"{sample}.selectBase.bam"
	message:"bed to bam"
	params:39
	priority: 8
	shell:
		"bedToBam -i {input} -g {bedToBam_G} | samtools sort -@ {params} -o {output} -"

rule indexBedToBam:
        input:
                rules.bedToBam.output
        output:
                "{sample}.selectBase.bam.bai"
        message:"index bedToBam"
        priority: 9
	shell:
                "samtools index {input}"

rule multiBamSummary:
	input:
		bam=rules.bedToBam.output,
		bai=rules.indexBedToBam.output
	output:
		TSf=temp("{sample}.TSfCount.txt"),
		TSr=temp("{sample}.TSrCount.txt"),
		NTSf=temp("{sample}.NTSfCount.txt"),
		NTSr=temp("{sample}.NTSrCount.txt")
	message:"count in genes"
	threads:39
	shell:
		"multiBamSummary BED-file \
		--bamfiles {input.bam} \
		--BED {BED10R} \
		--outRawCounts {output.TSf} \
		--numberOfProcessors {threads} \
		--samFlagExclude 16 \
		&& \
		multiBamSummary BED-file \
                --bamfiles {input.bam} \
                --BED {BED10F} \
                --outRawCounts {output.TSr} \
                --numberOfProcessors {threads} \
                --samFlagInclude 16 \
		&& \
		multiBamSummary BED-file \
                --bamfiles {input.bam} \
                --BED {BED10F} \
                --outRawCounts {output.NTSf} \
                --numberOfProcessors {threads} \
                --samFlagExclude 16 \
                && \
		multiBamSummary BED-file \
                --bamfiles {input.bam} \
                --BED {BED10R} \
                --outRawCounts {output.NTSr} \
                --numberOfProcessors {threads} \
                --samFlagInclude 16"

rule multiBamSummary2:
	input:
		TSf=rules.multiBamSummary.output.TSf,
		TSr=rules.multiBamSummary.output.TSr,
		NTSf=rules.multiBamSummary.output.NTSf,
		NTSr=rules.multiBamSummary.output.NTSr
	output:
		TS=temp("{sample}.TSCount.txt"),
		NTS=temp("{sample}.NTSCount.txt")
	message:"merge and sort for violin plot"
	shell:
		"""
		sort -k 1,1 -k 2,2n -k 3,3n {input.TSf} {input.TSr} | grep -v "#" > {output.TS} \
		&& \
		sort -k 1,1 -k 2,2n -k 3,3n {input.NTSf} {input.NTSr} | grep -v "#" > {output.NTS}
		"""

rule multiBamSummaryRPKM:
	input:
		sum=rules.damageSum.output,
		TS=rules.multiBamSummary2.output.TS,
		NTS=rules.multiBamSummary2.output.NTS
	output:
		TS="{sample}.GeneByGeneRPKM.TS.txt",
		NTS="{sample}.GeneByGeneRPKM.NTS.txt"
	message:"multiBamSummaryRPKM for violin plot"
	shell:
		"""
		SUM=$(cat {input.sum}) ; \
		KM=$(python -c "print(1e-3*($SUM*1e-6))") ; \
		awk -v KM=$KM '{{$4=$4/(KM*($3-$2)); OFS="\t"; print $0}}' {input.TS} | sort -k 1,1 -k 2,2n -k 3,3n - > {output.TS}
		awk -v KM=$KM '{{$4=$4/(KM*($3-$2)); OFS="\t"; print $0}}' {input.NTS} | sort -k 1,1 -k 2,2n -k 3,3n - > {output.NTS}
		"""

rule bamCoverage:
	input:
		bam=rules.bedToBam.output,
		bai=rules.indexBedToBam.output
	output:
		F=temp("{sample}.f.bdg"),
		R=temp("{sample}.r.bdg")
	message:"split strands, smooth for IGV and heatmap"
	params:
		binsize='100'
	threads:39
	shell:
		"""
		bamCoverage \
		--bam {input.bam} \
		--outFileName {output.F} \
		--outFileFormat bedgraph \
		-p {threads} \
		--binSize {params.binsize} \
		--normalizeUsing None \
		--filterRNAstrand reverse \
		--smoothLength 5000 \
		&& \
		bamCoverage \
		--bam {input.bam} \
		--outFileName {output.R} \
		--outFileFormat bedgraph \
		-p {threads} \
		--binSize {params.binsize} \
		--normalizeUsing None \
		--filterRNAstrand forward \
		--smoothLength 5000
		"""

rule bamCoverageRPKM:
	input:
		sum=rules.damageSum.output,
		F=rules.bamCoverage.output.F,
		R=rules.bamCoverage.output.R
	output:
		F=temp("{sample}.fRPKM.bdg"),
		R=temp("{sample}.rRPKM.bdg")
	message:"RPKM"
	shell:
		"""
		BINSIZE=100 ; \
		SUM=$(cat {input.sum}) ; \
		KM=$(python -c "print(($BINSIZE*(1e-3))*($SUM*1e-6))") ; \
		awk -v KM=$KM '{{$4=$4/KM; OFS="\t"; print $0}}' {input.F} | sort -k1,1 -k2,2n - > {output.F}
		awk -v KM=$KM '{{$4=$4/KM; OFS="\t"; print $0}}' {input.R} | sort -k1,1 -k2,2n - > {output.R}
		"""

rule bedGraphToBigWig:
	input:
		F=rules.bamCoverageRPKM.output.F,
		R=rules.bamCoverageRPKM.output.R
	output:
		F="{sample}.fRPKM.bw",
		R="{sample}.rRPKM.bw"
	message:"bedGraphToBigWig"
	shell:
		"""
		bedGraphToBigWig {input.F} \
		{bdgToBW_chromSizes} \
		{output.F}
		bedGraphToBigWig {input.R} \
		{bdgToBW_chromSizes} \
		{output.R}
		"""

rule computeMatrix:
	input:
		F=rules.bedGraphToBigWig.output.F,
		R=rules.bedGraphToBigWig.output.R
	output:
		TSf=temp("{sample}.TSf.Matrix.gz"),
		TSr=temp("{sample}.TSr.Matrix.gz"),
		NTSr=temp("{sample}.NTSr.Matrix.gz"),
		NTSf=temp("{sample}.NTSf.Matrix.gz")
	message:"get template/non-template strand signal,computeMatrix reference-point"
	threads: 39
	shell:
		"""
		computeMatrix reference-point \
		--scoreFileName {input.R} \
		--regionsFileName {BEDF} \
		--outFileName {output.TSf} \
		--referencePoint TSS \
		--beforeRegionStartLength 10000 \
		--afterRegionStartLength 100000 \
		--missingDataAsZero \
		--numberOfProcessors {threads} \
		--binSize 100
		computeMatrix reference-point \
		--scoreFileName {input.F} \
		--regionsFileName {BEDR} \
		--outFileName {output.TSr} \
		--referencePoint TSS \
		--beforeRegionStartLength 10000 \
		--afterRegionStartLength 100000 \
		--missingDataAsZero \
		--numberOfProcessors {threads} \
		--binSize 100
		computeMatrix reference-point \
                --scoreFileName {input.R} \
                --regionsFileName {BEDR} \
                --outFileName {output.NTSr} \
                --referencePoint TSS \
                --beforeRegionStartLength 10000 \
                --afterRegionStartLength 100000 \
                --missingDataAsZero \
                --numberOfProcessors {threads} \
                --binSize 100
		computeMatrix reference-point \
                --scoreFileName {input.F} \
                --regionsFileName {BEDF} \
                --outFileName {output.NTSf} \
                --referencePoint TSS \
                --beforeRegionStartLength 10000 \
                --afterRegionStartLength 100000 \
                --missingDataAsZero \
                --numberOfProcessors {threads} \
                --binSize 100
		"""

rule computeMatrixOperations:
	input:
		TSf=rules.computeMatrix.output.TSf,
		TSr=rules.computeMatrix.output.TSr,
		NTSr=rules.computeMatrix.output.NTSr,
		NTSf=rules.computeMatrix.output.NTSf
	output:
		TS="{sample}.TS.Matrix.gz",
		NTS="{sample}.NTS.Matrix.gz"
	message:"merge matrix.gz files"
	shell:
		"""
		computeMatrixOperations rbind \
		--matrixFile {input.TSf} {input.TSr} \
		--outFileName {output.TS}
		computeMatrixOperations rbind \
                --matrixFile {input.NTSf} {input.NTSr} \
                --outFileName {output.NTS}
		"""

rule plotHeatmap:
	input:
		TS=rules.computeMatrixOperations.output.TS,
		NTS=rules.computeMatrixOperations.output.NTS
	output:
		TS="{sample}.TS.pdf",
		NTS="{sample}.NTS.pdf"
	message:"plotHeatmap"
	shell:
		"""
		plotHeatmap \
		--matrixFile {input.TS}  \
		--outFileName {output.TS}\
		--sortUsing region_length \
		--sortRegions ascend \
		--regionsLabel "6406 expressed genes" \
		--samplesLabel "TS" \
		--heatmapHeight 3 \
		--heatmapWidth 1 \
		--whatToShow "heatmap and colorbar" \
		--zMin 0.1 \
		--zMax 0.4 \
		--colorList "white,blue"
		plotHeatmap \
		--matrixFile {input.NTS}  \
                --outFileName {output.NTS} \
                --sortUsing region_length \
                --sortRegions ascend \
                --regionsLabel "6406 expressed genes" \
                --samplesLabel "NTS" \
                --heatmapHeight 3 \
                --heatmapWidth 1 \
                --whatToShow "heatmap and colorbar" \
                --zMin 0.1 \
                --zMax 0.4 \
                --colorList "white,blue"
		"""

rule countDipyrimidine:
        input:
                rules.generate_bed.output.bed
        output:
                "BaseCount/{sample}.DipyrimidineRatio.txt"
        message:"Stat dipyrimidine ratio of damage site"
	shell:
                "python {SCRIPTS}/countDipyrimidine.py {input} {output}"


rule splitBed:
	input:
		rules.selectBase.output
	output:
		plus=temp("{sample}.plus.bed"),
		minus=temp("{sample}.minus.bed")
	shell:
		"python {SCRIPTS}/split_bed_by_strand.py {input} {output.plus} {output.minus}"


rule assignRead2geneTSS:
	input:
		intxt1=rules.splitBed.output.plus,
		intxt2=rules.splitBed.output.minus
	output:
		pos1=temp("{sample}.plus.POS.TSS.tsv"),
		neg1=temp("{sample}.plus.NEG.TSS.tsv"),
		pos2=temp("{sample}.minus.POS.TSS.tsv"),
		neg2=temp("{sample}.minus.NEG.TSS.tsv")
	run:
		shell("python {SCRIPTS}/assignReadtoGene.py \
			--tssdat {TSSDAT} \
			--intxt {input.intxt1} \
			--posout {output.pos1} \
			--negout {output.neg1} \
			--upstream {UPSTREAM_TSS} \
			--downstream {DOWNSTREAM_TSS} ")
		shell("python {SCRIPTS}/assignReadtoGene.py \
			--tssdat {TSSDAT} \
			--intxt {input.intxt2} \
			--posout {output.pos2} \
			--negout {output.neg2} \
			--upstream {UPSTREAM_TSS} \
			--downstream {DOWNSTREAM_TSS} ")


rule readProfileTSS:
	input:
		q=rules.assignRead2geneTSS.output.neg1,
		w=rules.assignRead2geneTSS.output.pos2,
		e=rules.assignRead2geneTSS.output.pos1,
		t=rules.assignRead2geneTSS.output.neg2,
		bgfile=rules.damageSum.output
	output:
		value="{sample}.TSS.rawValue"
	priority: 9
	shell:
		"python {SCRIPTS}/readProfile.py \
			--t1 {input.q} \
			--t2 {input.w} \
			--nt1 {input.e} \
			--nt2 {input.t} \
			--regionlength {REGION_TSS} \
			--binSize {BINSIZE_TSS} \
			--bgValue {input.bgfile} \
			--outRawvalue {output.value}"


rule assignRead2geneTES:
        input:
                intxt1=rules.splitBed.output.plus,
                intxt2=rules.splitBed.output.minus
        output:
                pos1=temp("{sample}.plus.POS.TES.tsv"),
                neg1=temp("{sample}.plus.NEG.TES.tsv"),
                pos2=temp("{sample}.minus.POS.TES.tsv"),
                neg2=temp("{sample}.minus.NEG.TES.tsv")
	run:
                shell("python {SCRIPTS}/assignReadtoGene.py \
			--tssdat {TESDAT} \
			--intxt {input.intxt1} \
			--posout {output.pos1} \
			--negout {output.neg1} \
                	--upstream {UPSTREAM_TES} \
			--downstream {DOWNSTREAM_TES} ")
                shell("python {SCRIPTS}/assignReadtoGene.py  \
			--tssdat {TESDAT} \
			--intxt {input.intxt2} \
			--posout {output.pos2} \
			--negout {output.neg2} \
                	--upstream {UPSTREAM_TES} \
			--downstream {DOWNSTREAM_TES} ")


rule readProfileTES:
        input:
                q=rules.assignRead2geneTES.output.neg1,
                w=rules.assignRead2geneTES.output.pos2,
                e=rules.assignRead2geneTES.output.pos1,
                t=rules.assignRead2geneTES.output.neg2,
                bgfile=rules.damageSum.output
        output:
                value="{sample}.TES.rawValue"
        priority: 9
	shell:
                "python {SCRIPTS}/readProfile.py \
			--t1 {input.q} \
			--t2 {input.w} \
			--nt1 {input.e} \
			--nt2 {input.t} \
                	--regionlength {REGION_TES} \
			--binSize {BINSIZE_TES} \
			--bgValue {input.bgfile} \
                	--outRawvalue {output.value}"


rule done:
	input:
		rules.bedGraphToBigWig.output,
		rules.plotHeatmap.output,
		rules.multiBamSummaryRPKM.output,
		rules.readProfileTSS.output,
		rules.readProfileTES.output,
		rules.countDipyrimidine.output
	output:
		touch("{sample}.done")
