#!/bin/bash

for c in `seq 22`
do
   qsub -N hwe_Chr${c} -m e -M wlhsu@uw.edu -q r420.q /projects/runRscript.sh /projects/wlhsu/src/hwe/hwe_geno_Chr${c}.R
done


