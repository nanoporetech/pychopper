#!/bin/bash

REF=$1
READ=$2
CORES=$3
NAME=$4
OUTD=wspace_${NAME}
BAM=${OUTD}/sorted_aln.bam

rm -fr $OUTD
mkdir -p $OUTD

minimap2 -t $CORES -ax map-ont $REF $READ | samtools view -b -q 5 | samtools sort -o $BAM -
samtools index $BAM
seqkit -j $CORES bam -f ReadLen $BAM -O $OUTD/read_len_report.pdf; sleep 3
seqkit -j $CORES bam -B 3 -f Strand $BAM -O $OUTD/strand_report.pdf; sleep 3
seqkit -j $CORES bam -f RefCov $BAM -O $OUTD/ref_cov_report.pdf; sleep 3
seqkit -j $CORES bam -f Read,RefCov $BAM 2>$OUTD/ref_cov.tsv
seqkit -j $CORES bam -f ReadCov $BAM -O $OUTD/read_cov_report.pdf; sleep 3
seqkit -j $CORES bam -f Read,ReadCov $BAM 2>$OUTD/read_cov.tsv
seqkit -j $CORES bam -f Acc $BAM -O $OUTD/acc_report.pdf; sleep 3
seqkit -j $CORES bam -f LeftClip $BAM -O $OUTD/left_clip_report.pdf; sleep 3
seqkit -j $CORES bam -f RightClip $BAM -O $OUTD/right_clip_report.pdf; sleep 3

