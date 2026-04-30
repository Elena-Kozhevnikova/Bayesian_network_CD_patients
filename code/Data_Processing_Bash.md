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
    echo "id,gene_name,tpm,transcript_id"
    ls GSM*.txt.gz | parallel '
        zcat "{}" | awk -v fname="{}" -F"\t" '\''
            BEGIN { 
                OFS="," 
                # Ищем ID пациента (B, C или D + 3 цифры)
                if (match(fname, /[BCD][0-9][0-9][0-9]/)) {
                    sid = substr(fname, RSTART, 4)
                } else {
                    sid = "UNKNOWN"
                }
            }
            # Пропускаем заголовок
            NR > 1 && sid != "UNKNOWN" {
                # Делим первую колонку по вертикальной черте |
                split($1, tags, "|")
                # tags[6] — это Gene Name (например, DDX11L1)
                # tags[1] — это Transcript ID (например, ENST00000456328.2)
                # $2 — это TPM
                print sid, tags[6], $2, tags[1]
            }'\''
    '
} > china_transcript_tpm.csv
```
