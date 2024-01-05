#!/usr/bin/env bash

DATE=`date +%Y%m%d`

#1a. Download DELVE FinalCut QA
#qadir=/home/s1/kadrlica/projects/delve/obs/v2
#qa_sql="${qadir}/delve_qa.sql"
#qa_csv="delve-qa-${DATE}.csv.gz"
#ssh kadrlica@des51.fnal.gov -t "cat ${qa_sql}"

#cmd="psql -h des51.fnal.gov BLISS -f ${qa_sql} | gzip > /tmp/${qa_csv}"
#echo $cmd
#ssh kadrlica@des51.fnal.gov -t "$cmd"
#scp kadrlica@des51.fnal.gov:/tmp/${qa_csv} data/

# 1a. Use most recent SISPI update
#sispi_qa_fits=decam-exposures-20230915.fits.gz
sispi_qa_fits=decam-exposures-${DATE}.fits.gz

# 1b. Use most recent update
delve_qa_csv=delve-qa-20230824.csv.gz

#1b. Download DECADE FinalCut QA
decade_qa_csv=decade-qa-20230824.csv.gz

#2. Download existing coverage
##python coverage.py decam-exposures-${DATE}.fits.gz --qa data/${delve_qa_csv} --qa data/${decade_qa_csv} --maps
###python coverage.py decam-exposures-${DATE}.fits.gz --qa data/${decade_qa_csv} --maps
python coverage.py ${sispi_qa_fits} --qa data/${delve_qa_csv} --maps

mkdir -p maps/${DATE}
cp decam-exposures-${DATE}.fits.gz maps/${DATE}/
mv decam_*.fits.gz maps/${DATE}/

##3. Update delve.py to link to new maps (opens at line +635)
lineno="+790"
emacs -nw $lineno obztak/obztak/delve.py

#4. Run survey prepare
survey_prepare
gzip -f target_fields.csv

#5. Check the updated fields
datadir=obztak/obztak/data
#fields=${datadir}/delve-target-fields-20230623.csv.gz
fields=${datadir}/delve-target-fields.csv.gz
echo "python -i update-fields.py ${fields} target_fields.csv.gz"
python -i update-fields.py ${fields} target_fields.csv.gz
#python -i update-fields.py ${fields} target_fields.csv.gz -o target_fields.csv.gz
#python -i update-fields.py target_fields.csv.gz target_fields.csv.gz
