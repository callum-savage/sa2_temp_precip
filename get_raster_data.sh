# Bash script to copy over the precipitation and maximum temperature raster 
# data from NCI. Requires a NCI login and access to the zv2 project.

# Configuration variables

USER="cs5232" # Update this with your NCI username
GADI_DM="gadi-dm.nci.org.au"

# Precipitation files

AGCD_PRECIP_PATH="../../../g/data/zv2/agcd/v1/precip/total/r005/01day"
MISSING_PRECIP_FILES=""
DELIM=""

for YEAR in {1999..2019}; do
  PRECIP_FILE="agcd_v1_precip_total_r005_daily_${YEAR}.nc"
  PRECIP_LOCAL_PATH="./data-raw/precip/${PRECIP_FILE}"
  PRECIP_REMOTE_PATH="${AGCD_PRECIP_PATH}/${PRECIP_FILE}"
  
  if [[ ! -f $PRECIP_LOCAL_PATH ]]; then
    MISSING_PRECIP_FILES="${MISSING_PRECIP_FILES}${DELIM}${PRECIP_REMOTE_PATH}"
    DELIM=","
  fi
done

if [[ ${#MISSING_PRECIP_FILES} > 0 ]]; then
  # Enter your password at the prompt
  scp ${USER}@${GADI_DM}:"${MISSING_PRECIP_FILES}" ./data-raw/precip/
fi

# Tmax files

AGCD_TMAX_PATH="../../../g/data/zv2/agcd/v1/tmax/mean/r005/01day"
MISSING_TMAX_FILES=()
DELIM=""

for YEAR in {1999..2019}; do
  TMAX_FILE="agcd_v1_tmax_mean_r005_daily_${YEAR}.nc"
  TMAX_LOCAL_PATH="data-raw/tmax/${TMAX_FILE}"
  TMAX_REMOTE_PATH="${AGCD_TMAX_PATH}/${TMAX_FILE}"
  
  if [[ ! -f ${TMAX_LOCAL_PATH} ]]; then
    MISSING_TMAX_FILES="${MISSING_TMAX_FILES}${DELIM}${TMAX_REMOTE_PATH}"
    DELIM=","
  fi
done

if [[ ${#MISSING_TMAX_FILES} > 0 ]]; then
  # Enter your password at the prompt
  scp ${USER}@${GADI_DM}:{"${MISSING_TMAX_FILES}"} ./data-raw/tmax/
fi

