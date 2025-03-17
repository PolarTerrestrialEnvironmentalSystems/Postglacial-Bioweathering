#!/bin/bash
#SBATCH --account=envi.envi
#SBATCH --job-name=step2nt
#SBATCH --partition=fat
#SBATCH --time=20:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=28
#SBATCH --mem=750G

# set required variables (adapt according to your own requirements)
#===================================================================
#KRAKEN
DB="nt_2022_10_db/"
CONFIDENCE="0.8"

# given variables (please do not change)
#===================================================================

WORK=${PWD}
OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_KRAKEN="out.kraken2"
OUT_KRONA="out.krona"

KRAKEN2="kraken2/2.1.2"
KRONA="krona/2.8.1"

END_R1="_fastp_R1.fq.gz"
END_R2="_fastp_R2.fq.gz"
END_MERGED="_fastp_merged_R2.fq.gz"

CPU=${SLURM_CPUS_PER_TASK}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

# prepare environment
#===================================================================
mkdir -p ${OUTDIR}/${OUT_KRAKEN}
mkdir -p ${OUTDIR}/${OUT_KRONA}

# tasks to be performed
#===================================================================

# KRAKEN2
#----------
module load ${KRAKEN2}
for fq in ${OUTDIR}/${OUT_FASTP}/*${END_MERGED} 
do 
	BASE=${fq##*/}
	ID=${BASE%${END_MERGED}}
	srun kraken2 --confidence ${CONFIDENCE} --db ${DB} ${fq} --threads ${CPU} --gzip-compressed --output ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_merged.kraken --report ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_merged.kraken.report
	srun kraken2 --confidence ${CONFIDENCE} --db ${DB} --paired ${OUTDIR}/${OUT_FASTP}/${ID}${END_R1} ${OUTDIR}/${OUT_FASTP}/${ID}${END_R2} --threads ${CPU} --gzip-compressed --output ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_paired.kraken --report ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_paired.kraken.report
done
module unload ${KRAKEN2}

# KRONA
#----------

module load ${KRONA}
for kraken in ${OUTDIR}/${OUT_KRAKEN}/*.kraken 
do 
	BASE=${kraken##*/}
	ID=${BASE%.kraken}
	srun ktImportTaxonomy -q 2 -t 3 ${kraken} -o ${OUTDIR}/${OUT_KRONA}/${ID}.html
done

srun ktImportTaxonomy -q 2 -t 3 ${OUTDIR}/${OUT_KRAKEN}/*.kraken -o ${OUTDIR}/${OUT_KRONA}/KRONA_plots_combined.html

module unload ${KRONA}
