#!/bin/bash
set -e
set -o pipefail
echo hostname: `hostname`
echo ==========start at : `date` ==========

/zfssz5/BC_PS/chengyuanfang/software/miniconda3/bin/python mutationNetwork.py -i input/maf2.list -k 2 -f 0.05 -n 30 -p1 0.01 -p2 0.05 -o /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/mutationNetwork/output2

if [ -e /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/mutationNetwork/work.sh.e* ]; then ls /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/mutationNetwork/work.sh.e* | while read dd ; do jobid=`echo $dd | sed -r 's/^.*\.e([0-9]*)$/\1/g'`; echo $jobid; qstat -j $jobid | grep 'usage'; done ; fi && \
echo ==========end at : `date` ==========
