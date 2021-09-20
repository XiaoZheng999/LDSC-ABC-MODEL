macs2 callpeak \
-t example_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
-n wgEncodeUwDnaseK562AlnRep1.chr22.macs2 \
-f BAM \
-g hs \
-p .1 \
--call-summits \
--outdir example_chr22/ABC_output/Peaks/ 

INFO  @ Sat, 18 Sep 2021 22:57:50: #1 read tag files...
INFO  @ Sat, 18 Sep 2021 22:57:50: #1 read treatment tags...
INFO  @ Sat, 18 Sep 2021 22:57:51: 434301 reads have been read.
INFO  @ Sat, 18 Sep 2021 22:57:51: #1 tag size is determined as 36 bps
INFO  @ Sat, 18 Sep 2021 22:57:51: #1 tag size = 36.0
INFO  @ Sat, 18 Sep 2021 22:57:51: #1  total tags in treatment: 434301
INFO  @ Sat, 18 Sep 2021 22:57:51: #1 user defined the maximum tags...
INFO  @ Sat, 18 Sep 2021 22:57:51: #1 filter out redundant tags at the same location and the same strand by allowing at most 1 tag(s)
INFO  @ Sat, 18 Sep 2021 22:57:51: #1  tags after filtering in treatment: 330933
INFO  @ Sat, 18 Sep 2021 22:57:51: #1  Redundant rate of treatment: 0.24
INFO  @ Sat, 18 Sep 2021 22:57:51: #1 finished!
INFO  @ Sat, 18 Sep 2021 22:57:51: #2 Build Peak Model...
INFO  @ Sat, 18 Sep 2021 22:57:51: #2 looking for paired plus/minus strand peaks...
INFO  @ Sat, 18 Sep 2021 22:57:51: #2 number of paired peaks: 2488
INFO  @ Sat, 18 Sep 2021 22:57:51: start model_add_line...
INFO  @ Sat, 18 Sep 2021 22:57:51: start X-correlation...
INFO  @ Sat, 18 Sep 2021 22:57:51: end of X-cor
INFO  @ Sat, 18 Sep 2021 22:57:51: #2 finished!
INFO  @ Sat, 18 Sep 2021 22:57:51: #2 predicted fragment length is 297 bps
INFO  @ Sat, 18 Sep 2021 22:57:51: #2 alternative fragment length(s) may be 297 bps
INFO  @ Sat, 18 Sep 2021 22:57:51: #2.2 Generate R script for model : example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22
INFO  @ Sat, 18 Sep 2021 22:57:51: #3 Call peaks...
INFO  @ Sat, 18 Sep 2021 22:57:51: #3 Going to call summits inside each peak ...
INFO  @ Sat, 18 Sep 2021 22:57:51: #3 Call peaks with given -log10pvalue cutoff: 1.00000 ...
INFO  @ Sat, 18 Sep 2021 22:57:51: #3 Pre-compute pvalue-qvalue table...
INFO  @ Sat, 18 Sep 2021 22:57:52: #3 Call peaks for each chromosome...
INFO  @ Sat, 18 Sep 2021 22:58:12: #4 Write output xls file... example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_
INFO  @ Sat, 18 Sep 2021 22:58:16: #4 Write peak in narrowPeak format file... example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRearrowPeak
INFO  @ Sat, 18 Sep 2021 22:58:16: #4 Write summits bed file... example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2
INFO  @ Sat, 18 Sep 2021 22:58:16: Done!



#Sort narrowPeak file
bedtools sort -faidx example_chr22/reference/chr22 -i example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak > example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted


conda activate final-abc-env

python src/makeCandidateRegions.py \
--narrowPeak example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted \
--bam example_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
--outDir example_chr22/ABC_output/Peaks/ \
--chrom_sizes example_chr22/reference/chr22 \
--regions_blocklist reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--regions_includelist example_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.chr22.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 3000 
recommend using --nStrongestPeaks 150000 when making genome-wide peak calls



Running: bedtools sort -i example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.wgEncodeUwDnaseK562AlnRep1.chr22.bam.Counts.bed -faidx example_chr22/reference/chr22 | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n 3000 |bedtools intersect -b stdin -a example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted -wa |awk '{print $1 "\t" $2 + $10 "\t" $2 + $10}' |bedtools slop -i stdin -b 250 -g example_chr22/reference/chr22 |bedtools sort -i stdin -faidx example_chr22/reference/chr22 |bedtools merge -i stdin | bedtools intersect -v -wa -a stdin -b reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed | cut -f 1-3 | (bedtools intersect -a example_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.chr22.bed -b example_chr22/reference/chr22.bed -wa | cut -f 1-3 && cat) |bedtools sort -i stdin -faidx example_chr22/reference/chr22 | bedtools merge -i stdin > example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.candidateRegions.bed
