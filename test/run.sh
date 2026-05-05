python ./findGEVE.py \
       Chlamydomonas_reinhardtii.contig536.fa \
       -db database \
       --prefix Chlamidomonas \
       --threads 16 \
       --min-hallmark-type 3 \
       --blastn-jobs 8
