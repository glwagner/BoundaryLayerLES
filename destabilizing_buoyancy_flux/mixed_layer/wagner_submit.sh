#!/bin/bash

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_mit_raffaele    (32)    2-10-1  12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00

if [ $1 -eq 0 ]; then
    nodes=1
else
    nodes=$1
fi

jobname="ML$nodes"
slurmoutput="$jobname.%j.out"
scriptname="$jobname.slurm"
walltime="4:00:00"
runscript="mixed_layer_convection.py"

closure="ConstantSmagorinsky"
outputname="$jobname.$closure.out"

rm -f $scriptname

echo "#!/bin/bash

# Job
#SBATCH --reservation=wagner
#SBATCH --partition=sched_mit_hill
#SBATCH --qos=plenum

#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=$walltime
#SBATCH --job-name=$jobname

# Streams
#SBATCH --output=$slurmoutput

# Activate conda and dedalus environment
export DEDALES='/nobackup1/glwagner/dedaLES'
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

# Content
mpiexec python3 $runscript $closure >> $outputname" >> $scriptname

sbatch $scriptname
squeue -u $USER
