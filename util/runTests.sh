    ./meRanT mkbsidx \
      -t 2 \
      -fa ./testdata/refSeq/mm10.refSeqRNA-noPRED.500.fa \
      -id ./testdata/mm10/meRanTIDX/

    ./meRanT align \
       -t 2 \
       -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
       -i2g ./testdata/refSeq/mm10.refSeqRNA-noPRED.500.t2g.map \
       -o ./testdata/results \
       -S meRanT_test.sam  \
       -x ./testdata/mm10/meRanTIDX/mm10.refSeqRNA-noPRED.500_C2T

    ./meRanCall \
      -p 2 \
      -s ./testdata/results/meRanT_test.sam \
      -f ./testdata/refSeq/mm10.refSeqRNA-noPRED.500.fa \
      -tref \
      -o ./testdata/results/meRanT_meRanCall_m5Ct.txt \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -sc 5

   ./meRanT mkbsidx \
      -t 2 \
      -fa ./testdata/tRNAs/mm10.tRNAs.20140204.fa \
      -id ./testdata/mm10/meRanT_tRNA_IDX

    ./meRanT align \
       -t 2 \
       -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
       -i2g ./testdata/tRNAs/mm10.tRNA.20140204.map \
       -o ./testdata/results \
       -S meRanT_tRNAs_test.sam  \
       -x ./testdata/mm10/meRanT_tRNA_IDX/mm10.tRNAs.20140204_C2T

   ./meRanCall \
      -p 2 \
      -s ./testdata/results/meRanT_tRNAs_test.sam \
      -f ./testdata/tRNAs/mm10.tRNAs.20140204.fa \
      -tref \
      -o ./testdata/results/meRanT_tRNAs_meRanCall_m5Ct.txt \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -sc 5


   ./meRanGs mkbsidx  \
      -t 2 \
      -fa ./testdata/mm10/chr19.fa \
      -GTF ./testdata/mm10/ref_GRCm38.p2_top_level_no_prediction_chr19_sort.gff3 \
      -GTFtagEPT Parent \
      -GTFtagEPG gene \
      -id ./testdata/mm10/meRanGsIDX

    ./meRanGs align \
      -t 2 \
      -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
      -id ./testdata/mm10/meRanGsIDX \
      -bg \
      -o ./testdata/results \
      -S meRanGs_test.sam \
      -MM \
      -un \
      --star_genomeLoad NoSharedMemory

    ./meRanCall \
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


    ./meRanGh mkbsidx  \
      -t 2 \
      -fa ./testdata/mm10/chr19.fa \
      -id ./testdata/mm10/meRanGhIDX

    ./meRanGh align \
      -t 2 \
      -f ./testdata/fastq/clean_KHOD_400k_test.fastq \
      -id ./testdata/mm10/meRanGhIDX \
      -bg \
      -o ./testdata/results \
      -S meRanGh_test.sam \
      -MM \
      -un

    ./meRanCall \
      -p 2 \
      -s ./testdata/results/meRanGh_test_sorted.bam \
      -f ./testdata/mm10/chr19.fa \
      -gref \
      -o ./testdata/results/meRanGh_meRanCall_m5Cs.txt \
      -md 5 \
      -cr 0.99 \
      -mBQ 30 \
      -bed63 \
      -sc 5


