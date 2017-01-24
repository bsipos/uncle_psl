#!/bin/bash

blat -out=psl -fine -q=rna ref.fas reads.fas blat_all.psl
pslReps blat_all.psl blat_top.psl blat.psr
uncle_psl.py -f reads.fas blat_top.psl blat.sam

bwa index ref.fas
bwa mem -x ont2d ref.fas reads.fas > bwa.sam
