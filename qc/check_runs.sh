
###############################
##    USER CHANGE
###############################
YEAR=$1 #used as a stub for multiple runs
EXPECTED_FILE_COUNT=2880 #number of expected output files based on comet config
FILE_SIZE_MIN=2 #for a 2.8G avg, use a small value, can flag final output files
RUN_TIME_MIN=5 #expected runtime per file in mintues, use a smaller value
EXPECTED_DIRECTORY_SIZE="unknown" #e.g. "~7.5-8T"
###############################
###############################
#Set file paths
#output directory name
OUT_DIR="postprocessed_txts_$YEAR"
LOG_DIR="slurms_postprocess_summit/out_slurm_global_yearly_lats\[-090.00_+090.00\]_lngs\[-180.00_+180.00\]_dates\[01-${YEAR}_12-${YEAR}\]"
MASTER_LOG_STUB="output-postprocess-"
INDIVD_LOG_STUB="output-postprocess-_gpfs_alpine_syb105_proj-shared_Projects_Climatype_incite_global_yearly_comet_full_runs_summit_output_${YEAR}_out_"
###############################

#Check full directory size
echo "Actual total dir size: $(du -sh $OUT_DIR)" 
echo "Expected dir size: $EXPECTED_DIRECTORY_SIZE"

#Check number of output files
num_files=$(ls $OUT_DIR | wc -l)
if [[ $num_files != $EXPECTED_FILE_COUNT ]]
then
  echo ">>>>>!!!!!File count incorrect!!!!!>>>>>"
else
  echo "File count good"
fi 

#Check size of individual output files
for f in $OUT_DIR/*
do
 tmp=$(du -sh $f | cut -f1)
 file_size=${tmp::-1}
 #echo "File $f"
 #echo " File size: $file_size"
 #if [[ $file_size -lt $FILE_SIZE_MIN ]]
 check=$(echo "$file_size < $FILE_SIZE_MIN" | bc)
 if [[ $check == 1 ]]
 then
  echo "File unexpectly small ($f)"
  echo "  Actual $file_size, min is $FILE_SIZE_MIN"
 fi
done

#Check master log files
for file_name in $LOG_DIR/$MASTER_LOG_STUB[0-9]*.txt
 do
  #echo "$file_name"
  anyerror=$(grep -cFi "error" $file_name)
  if [ $anyerror -gt 0 ]
  then
   echo ">>>Error in output file"
   echo "      $file_name"
   echo
  fi
done

#Check that all individual files exist and are non-empty
for i in {0000..2879}
do
 file_name="$LOG_DIR/$INDIVD_LOG_STUB${i}.txt"
 if [ -f $file_name ]
 then
    if ! [ -s $file_name ]
    then
        echo "File $i exists but empty"
    else
  	tmp=$(cat $file_name | grep -oP '(?<=real)[^ ]*')
	min=${tmp%m*}
 	if [[ min -lt $RUN_TIME_MIN ]]
	then
	 echo ">>>Runtime too short"
 	 echo "      $file_name ($tmp)"
	fi
    fi
 else
    echo "File $i does not exists"
fi
done

