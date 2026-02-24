set -e

CYAN='\033[1;36m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

DATA_ROOT=${DATA_ROOT:-data}
MUTANTS_ROOT=${MUTANTS_ROOT:-${DATA_ROOT}/mutants}

# Preferred layout: data/mutants/<BASE>/...
# Backward-compatible fallback: data/<BASE>/...
if [ -d "${MUTANTS_ROOT}/${BASE}" ]; then
  INPUT_DIR="${MUTANTS_ROOT}/${BASE}"
else
  INPUT_DIR="${DATA_ROOT}/${BASE}"
fi

OUTPUT_DIR="${INPUT_DIR}/output"
