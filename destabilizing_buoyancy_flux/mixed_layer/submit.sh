#!/bin/bash

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_mit_raffaele    (32)    2-10-1  12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00

nodes=$1
jobname="ML$nodes"
closure="StratifiedAnisotropicMinimumDissipation"
scriptname="$jobname.slurm"

rm -f $scriptname

echo "#!/bin/bash

# Job
##SBATCH --reservation=wagner
#SBATCH --partition=sched_mit_hill
##SBATCH --qos=plenum

#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=4:00:00
#SBATCH --job-name=$jobname

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Activate conda and dedalus environment
export DEDALES='/nobackup1/glwagner/dedaLES'
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

# Content
mpiexec python3 mixed_layer_convection.py $closure >> $jobname.out" >> $scriptname

sbatch $scriptname
