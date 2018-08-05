#!/bin/bash
##SBATCH --mail-user="maeva.techer@oist.jp"
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=fastsimcoal
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=5G
#SBATCH -c 10
#SBATCH -t 6-23
##SBATCH --array=1-3

/work/MikheyevU/Maeva/demography/fsc26_linux64/fsc26 --tplfile 2PopDiv20Mb.tpl --estfile 2PopDiv20Mb.est -d --numsims 100000 -M --minnumloops 10 --numloops 40 -c 10
