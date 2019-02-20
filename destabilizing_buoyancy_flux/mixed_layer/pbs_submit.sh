#!/bin/bash

nodes=$1
jobname="ML$nodes"
runscript="mixed_layer_convection.py"
scriptname="$jobname.pbs"
walltime="04:00:00"
TMPDIR=/glade/scratch/$USER/temp

# Submit AMD
closure="ConstantSmagorinsky"
outputname=$jobname.$closure.out

rm -f $scriptname

echo "#!/bin/bash

### Job Name
#PBS -N $jobname

### Project code
#PBS -A UMIT0023
#PBS -l walltime=$walltime
#PBS -q regular

### Merge output and error files
#PBS -j oe

### Select 2 nodes with 36 CPUs each for a total of 72 MPI processes
#PBS -l select=$nodes:ncpus=36:mpiprocs=36:mem=109GB

### Specify email recpient and send email on abort, begin and end
#PBS -M wagner.greg@gmail.com
#PBS -m abe

### Content
export DEDALES="$DEDALES"
export TMPDIR=$TMPDIR
mkdir -p $TMPDIR

### Activate Miniconda
. /glade/u/home/$USER/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

module load intel/17.0.1 
module load openmpi/3.0.1

cp $runscript $TMPDIR/
cd $TMPDIR

### Run the executable
mpiexec python3 $runscript $closure >> $PWD/$outputname" >> $scriptname

echo "Output in $outputname"
qsub $scriptname
qstat
