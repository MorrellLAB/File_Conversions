
INPUT_DIR=""
OUTPUT_DIR=""
REFERENCE=""

#filter parameters 
MIN_MQ=20     # Minimum mapping quality
MIN_BQ=20     # Minimum base quality
MIN_DEPTH=10  # Minimum depth for variant calling
MAX_DEPTH=250 # Maximum depth to avoid regions with excessive coverage
 
module load bcftools

mkdir -p "${OUTPUT_DIR}/filtered_results"

process_sample() {
    local VCF_FILE=$1
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${VCF_FILE}" .vcf.gz)
    echo "Processing sample (chromosome): ${SAMPLE_NAME}"

    # create dir

    mkdir -p "${OUTPUT_DIR}/filtered_results/${SAMPLE_NAME}"
    local OUTPUT_PREFIX="${OUTPUT_DIR}/filtered_results/${SAMPLE_NAME}/${SAMPLE_NAME}"

    # check if VCF exist
    if ! bcftools view -h "${VCF_FILE}" > /dev/null 2>&1; then
        echo "ERROR: ${VCF_FILE} is not a valid VCF file or is corrupted."
        return 1
    fi

    # Step 1: Remove invariant sites
    # The -c 1 flag ensures that only sites with at least one alternate allele are kept.

    echo "   -> Filtering for polymorphic sites only..."
    bcftools view -c 1 -Oz -o "${OUTPUT_PREFIX}.poly.vcf.gz" "${VCF_FILE}"
    bcftools index "${OUTPUT_PREFIX}.poly.vcf.gz"

    # Step 2: Apply additional filtering based on quality and depth
    echo "   -> Applying quality and depth filters..."
    bcftools filter -e "QUAL<20 || MQ<${MIN_MQ} || FORMAT/DP<${MIN_DEPTH} || FORMAT/DP>${MAX_DEPTH} || INFO/DP<${MIN_DEPTH} || INFO/DP>${MAX_DEPTH}" \
        -Oz \
        -o "${OUTPUT_PREFIX}.filtered.vcf.gz" \
        "${OUTPUT_PREFIX}.poly.vcf.gz"
    bcftools index "${OUTPUT_PREFIX}.filtered.vcf.gz"

    echo "Finished processing ${SAMPLE_NAME}"
}

# Process all VCF files in the input directory (per chromosome)
echo "Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_sample "${VCF_FILE}"
done
