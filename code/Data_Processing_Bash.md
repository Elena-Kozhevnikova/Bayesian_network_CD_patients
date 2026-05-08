# Prepare data for further analysis using ```awk``` in command line

## Data aquisition and pre-processing

The Israeli patients' sequence data was downloaded from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199906.
The Chinese patients' sequence data was downloaded from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233900.
Metabolomes and metadata were downloaded from the aticle supplementary file: https://www.nature.com/articles/s41467-024-48106-6

### 1. Install GNU Parallel for faster processing
```
sudo apt update && sudo apt install parallel -y
```

### 2. Extract the raw data from GEO archive
```
tar -xvf GSE199906_RAW.tar # Israeli patients
tar -xvf GSE233900_RAW.tar # Chinese patients
```

### 3. Process and merge all files into a single CSV
### Israeli patients: This script extracts Sample ID (A0xx), Gene Name, TPM, and Transcript ID
```
{
    echo "sample_id,gene_name,tpm,transcript_id"
    ls GSM*.txt.gz | parallel '
        zcat "{}" | awk -v fname="{}" -F"\t" '\''
            BEGIN { 
                OFS="," 
                # Strict match: "A", then "0", then two digits (e.g., A025)
                match(fname, /A0[0-9][0-9]/)
                if (RSTART != 0) {
                    sid = substr(fname, RSTART, 4)
                } else {
                    sid = "UNKNOWN"
                }
            }
            # Skip header row and process only recognized Sample IDs
            NR > 1 && sid != "UNKNOWN" {
                split($1, parts, "_")
                # $1 contains Gene_Transcript, $2 contains TPM value
                print sid, parts[1], $2, parts[2]
            }'\''
    '
} > israel_transcript_tpm.csv
```
### Chinese patients: This script extracts Sample ID (Вхxx), Gene Name, TPM, and Transcript ID
```
{
    # Приводим к формату как в Израельском датасете
    echo "sample_id,gene_name,tpm,transcript_id"
    ls GSM*.txt.gz | parallel '
        zcat "{}" | awk -v fname="{}" -F"\t" '\''
            BEGIN {
                OFS=","
                # Извлекаем ID (B/C/D + 3 цифры)
                if (match(fname, /[BCD][0-9][0-9][0-9]/)) {
                    sid = substr(fname, RSTART, 4)
                } else {
                    sid = "UNKNOWN"
                }
            }
            
            NR > 1 && sid != "UNKNOWN" {
                split($1, tags, "|")
                # ВАЖНО: приводим к порядку: sid, gene (tags[6]), tpm ($2), transcript (tags[1])
                print sid, tags[6], $2, tags[1]
            }'\''
    '
} > china_transcript_tpm.csv
```

## Metadata Preparation: Extracting data from supplementary tables and cohort stratification
```
# Update package lists
sudo apt update

# Install Python package manager and virtual environment tools
sudo apt install python3-pip python3-venv -y

# Install Pandas and openpyxl (required for .xlsx support)
# Using apt is recommended to avoid conflicts with system-wide paths
sudo apt install python3-pandas python3-openpyxl -y
```
### Create a new file for the metadata processing script:
```
#!/bin/bash

FILE_NAME="41467_2024_48106_MOESM4_ESM.xlsx"

if [ ! -f "$FILE_NAME" ]; then
    echo "Error: File $FILE_NAME not found!"
    exit 1
fi

echo "Starting metadata separation (skiprows=1)..."

python3 - <<EOF
import pandas as pd
import sys

try:
    # Read Excel, skipping the first description row
    df = pd.read_excel('$FILE_NAME', sheet_name='Samples_Datasets', skiprows=1)

    # Clean column names
    df.columns = df.columns.str.strip()

    target_col = 'Patient_group'
    if target_col not in df.columns:
        print(f"Error: Column '{target_col}' not found.")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)

    # Filter by cohort
    china = df[df[target_col].str.startswith('China', na=False)]
    israel = df[df[target_col].str.startswith('Israel', na=False)]

    # Save to CSV
    china.to_csv('China_metadata.csv', index=False)
    israel.to_csv('israel_metadata.csv', index=False)

    print("Success!")
    print(f"China_metadata.csv: {len(china)} rows")
    print(f"israel_metadata.csv: {len(israel)} rows")

except Exception as e:
    print(f"Python Error: {e}")
EOF
```
