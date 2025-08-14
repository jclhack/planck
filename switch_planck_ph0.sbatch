#!/bin/bash
#
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=switch_planck_ph0
#SBATCH --partition=secondary
#SBATCH --output=switch_planck_ph0.o%j
##SBATCH --error=switch_planck_ph0.e%j
##SBATCH --mail-user=jcl11@illinois.edu
##SBATCH --mail-type=BEGIN,END
#
# End of embedded SBATCH options
#

# Run the planck executable (switch_planck_ph0)
./switch_planck_ph0.exe
