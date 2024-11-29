#!/bin/bash

module load bcftools/1.9

# Archivo CSV con nombres de muestras (uno por línea)
SAMPLES_CSV=$1
REGIONS_PATH="/mnt/Timina/lmorales/Public/ymez/data/bam"
VCF_PATH="/mnt/Timina/lmorales/Public/ymez/tmp/05_vcalling"
OUTPUT_PATH="./" # Cambiar si es necesario guardar en otro lugar

while IFS= read -r sample; do
    # Archivos específicos para cada muestra
    REGIONS_FILE_sp="${REGIONS_PATH}/${sample}_output_sp_500ntWindow.depth"
    VCF_FILE_sp="${VCF_PATH}/${sample}_CONC.SNP_onlychr_SAPA.g.vcf.gz"
    OUTPUT_FILE_sp="${OUTPUT_PATH}/${sample}_heterozygosity_per_region_SAPA.txt"

    REGIONS_FILE_sc="${REGIONS_PATH}/${sample}_output_sc_500ntWindow.depth"
    VCF_FILE_sc="${VCF_PATH}/${sample}_CONC.SNP_onlychr_SACE.g.vcf.gz"
    OUTPUT_FILE_sc="${OUTPUT_PATH}/${sample}_heterozygosity_per_region_SACE.txt"

    echo "Processing sample: $sample"

    # Llamar al script windowed_heterozygosity.sh para ambas combinaciones
    ./windowed_heterozygosity.sh "$REGIONS_FILE_sp" "$VCF_FILE_sp" "$OUTPUT_FILE_sp" "$sample"
    ./windowed_heterozygosity.sh "$REGIONS_FILE_sc" "$VCF_FILE_sc" "$OUTPUT_FILE_sc" "$sample"
done < "$SAMPLES_CSV"

echo "Cálculo de heterocigosidad completado para todas las muestras."
