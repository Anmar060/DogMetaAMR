#!/usr/bin/env
#set -euo pipefail
#first time using set: changes script behaviour to errors. This script kept failing and its still a bit rubbish, its important its reliable however
#e - for stopping if anything exist with non-zero status
#u - unset variables are treated as errors, if a variable was accidentally typo-ed or not set
#o pipefail -  

# This script should be run from inside /nanoplot, it only works on the specific structure of Nanoplots NanoStats file
# It will create NanoStats_summary.csv in the current directory

OUTFILE="NanoStats_summary.csv"

# Write CSV header
echo 
"Sample,MeanReadLength,MeanReadQuality,MedianReadLength,MedianReadQuality,NumReads,ReadLengthN50,STDEVReadLength,TotalBases,Q10_re>
> "$OUTFILE"

#Loop over all NanoStats.txt files under current folder
for FILE in */NanoStats.txt; do
    SAMPLE="$(basename "$(dirname "$FILE")")"

    #Looks for all pf the general summary values
    mean_len=$(grep "Mean read length" "$FILE" | awk '{print $4}')
    mean_qual=$(grep "Mean read quality" "$FILE" | awk '{print $4}')
    median_len=$(grep "Median read length" "$FILE" | awk '{print $4}')
    median_qual=$(grep "Median read quality" "$FILE" | awk '{print $4}')
    num_reads=$(grep "Number of reads" "$FILE" | awk '{print $4}' | tr -d 
',')
    n50=$(grep "Read length N50" "$FILE" | awk '{print $4}')
    stdev_len=$(grep "STDEV read length" "$FILE" | awk '{print $4}')
    total_bases=$(grep "Total bases" "$FILE" | awk '{print $3}' | tr -d 
',')

    #looks at the Quality cutoffs (reads, %, Mb) section that breaks the pattern
    q10=$(grep ">Q10:" "$FILE" | awk '{gsub("[()%]", "", $0); print 
$2","$3","$4}' | sed 's/Mb//')
    q15=$(grep ">Q15:" "$FILE" | awk '{gsub("[()%]", "", $0); print 
$2","$3","$4}' | sed 's/Mb//')
    q20=$(grep ">Q20:" "$FILE" | awk '{gsub("[()%]", "", $0); print 
$2","$3","$4}' | sed 's/Mb//')
    q25=$(grep ">Q25:" "$FILE" | awk '{gsub("[()%]", "", $0); print 
$2","$3","$4}' | sed 's/Mb//')
    q30=$(grep ">Q30:" "$FILE" | awk '{gsub("[()%]", "", $0); print 
$2","$3","$4}' | sed 's/Mb//')

    #looks at the' Top 5 highest quality reads' (score,len)
    topqual=$(awk '/Top 5 highest mean basecall quality 
scores/{flag=1;next}/Top 5 longest/{flag=0}flag' "$FILE" \
        | awk -F '[(): ]+' '{print $2","$3}' | paste -sd "," -)

    # looks at the 'Top 5 longest reads' bit (len,qual)
    toplen=$(awk '/Top 5 longest reads/{flag=1;next}flag' "$FILE" \
        | awk -F '[(): ]+' '{print $2","$3}' | paste -sd "," -)

    #combine into one row
    echo 
"$SAMPLE,$mean_len,$mean_qual,$median_len,$median_qual,$num_reads,$n50,$stdev_len,$total_bases,$q10,$q15,$q20,$q25,$q30,$topqual,$>
>> "$OUTFILE"
done
#location of the output file
echo "Summary written to $(pwd)/$OUTFILE"
