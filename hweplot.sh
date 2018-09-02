#!/bin/bash

for c in `seq 22`
do
   qsub -N hwe_plot_Chr${c} -m e -M wlhsu@uw.edu -q bigmem.q /projects/QCpipeline/runRscript.sh /projects/wlhsu/src/hwe/jointHW/Chplot/hwe_plot_Chr${c}.R
done
