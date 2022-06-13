    ./meRanT.pl mkbsidx \
      -t 2 \
      -fa ./testdata/refSeq/mm10.refSeqRNA-noPRED.500.fa \
      -id ./testdata/mm10/meRanTIDX/

    ./meRanT.pl align \
       -t 2 \
       -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
       -i2g ./testdata/refSeq/mm10.refSeqRNA-noPRED.500.t2g.map \
       -o ./testdata/results \
       -S meRanT_test.sam  \
       -x ./testdata/mm10/meRanTIDX/mm10.refSeqRNA-noPRED.500_C2T -debug

    ./meRanCall.pl \
      -p 2 \
      -s ./testdata/results/meRanT_test.sam \
      -f ./testdata/refSeq/mm10.refSeqRNA-noPRED.500.fa \
      -tref \
      -o ./testdata/results/meRanT_meRanCall_m5Ct.txt \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -sc 5

   ./meRanT.pl mkbsidx \
      -t 2 \
      -fa ./testdata/tRNAs/mm10.tRNAs.20140204.fa \
      -id ./testdata/mm10/meRanT_tRNA_IDX

    ./meRanT.pl align \
       -t 2 \
       -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
       -i2g ./testdata/tRNAs/mm10.tRNA.20140204.map \
       -o ./testdata/results \
       -S meRanT_tRNAs_test.sam  \
       -x ./testdata/mm10/meRanT_tRNA_IDX/mm10.tRNAs.20140204_C2T

   ./meRanCall.pl \
      -p 2 \
      -s ./testdata/results/meRanT_tRNAs_test.sam \
      -f ./testdata/tRNAs/mm10.tRNAs.20140204.fa \
      -tref \
      -o ./testdata/results/meRanT_tRNAs_meRanCall_m5Ct.txt \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -sc 5

   ./meRanGs.pl mkbsidx  \
      -t 2 \
      -fa ./testdata/mm10/chr19.fa \
      -GTF ./testdata/mm10/ref_GRCm38.p2_top_level_no_prediction_chr19_sort.gff3 \
      -GTFtagEPT Parent \
      -GTFtagEPG gene \
      -id ./testdata/mm10/meRanGsIDX

    ./meRanGs.pl align \
      -t 2 \
      -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
      -id ./testdata/mm10/meRanGsIDX \
      -bg \
      -o ./testdata/results \
      -S meRanGs_test.sam \
      -MM \
      -un \
      --star_genomeLoad NoSharedMemory

    ./meRanCall.pl \
      -p 2 \
      -s ./testdata/results/meRanGs_test_sorted.bam \
      -f ./testdata/mm10/chr19.fa \
      -gref \
      -o ./testdata/results/meRanGs_meRanCall_m5Cs.txt \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -bed63 \
      -sc 5

    ./meRanGs.pl align \
      -t 2 \
      -f ./testdata/fastq/r1.fastq \
      -r ./testdata/fastq/r2.fastq \
      -id ./testdata/mm10/meRanGsIDX \
      -bg \
      -o ./testdata/results \
      -S meRanGs_test_pe.sam \
      -MM \
      -un \
      --star_genomeLoad NoSharedMemory \
      -mbp -deleteSAM -dt -star_outFilterScoreMinOverLread 0.40


    ./meRanCall.pl \
      -p 12 \
      -s ./testdata/results/meRanGs_test_pe_sorted.bam \
      -f ./testdata/mm10/chr19.fa \
      -gref \
      -o ./testdata/results/meRanGs_pe_meRanCall_m5Cs.txt \
      -gtf ./testdata/mm10/gencode_mm10_chr19.gtf \
      -C_cutoff 5 \
      -mcov 5 \
      -mc 1 \
      -mr 0.1 \
      -ei 0.1 \
      -sc 10 \
      -fdr 0.05 \
      -np \
      -gcr \
      -rl 150 \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -bed63 \
      -sc 5

    ./meRanGh.pl mkbsidx  \
      -t 2 \
      -fa ./testdata/mm10/chr19.fa \
      -id ./testdata/mm10/meRanGhIDX

    ./meRanGh.pl align \
      -t 2 \
      -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
      -id ./testdata/mm10/meRanGhIDX \
      -bg \
      -o ./testdata/results \
      -S meRanGh_test.sam \
      -MM \
      -un

    ./meRanCall.pl \
      -p 2 \
      -s ./testdata/results/meRanGh_test_sorted.bam \
      -f ./testdata/mm10/chr19.fa \
      -gref \
      -o ./testdata/results/meRanGh_meRanCall_m5Cs.txt \
      -gtf ./testdata/mm10/gencode_mm10_chr19.gtf \
      -C_cutoff 5 \
      -mcov 5 \
      -mc 1 \
      -mr 0.1 \
      -ei 0.1 \
      -sc 10 \
      -fdr 0.05 \
      -np \
      -gcr \
      -rl 150 \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -bed63 \
      -sc 5

    ./meRanGh.pl align \
      -t 2 \
      -f ./testdata/fastq/r1.fastq \
      -r ./testdata/fastq/r2.fastq \
      -id ./testdata/mm10/meRanGhIDX \
      -bg \
      -o ./testdata/results \
      -S meRanGh_test_pe.sam \
      -mbp -deleteSAM \
      -dt \
      -MM \
      -un \
      -dovetail \

    ./meRanCall.pl \
      -p 12 \
      -s ./testdata/results/meRanGh_test_pe_sorted.bam \
      -f ./testdata/mm10/chr19.fa \
      -gref \
      -o ./testdata/results/meRanGh_pe_meRanCall_m5Cs.txt \
      -gtf ./testdata/mm10/gencode_mm10_chr19.gtf \
      -C_cutoff 5 \
      -mcov 5 \
      -mc 1 \
      -mr 0.1 \
      -ei 0.1 \
      -sc 10 \
      -fdr 0.05 \
      -np \
      -gcr \
      -rl 150 \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -bed63 \
      -sc 5
