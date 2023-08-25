#!/bin/bash

# Copy meta farm to scratch directory and submit the
# job there

p=`pwd`
part1=$(dirname "$p")
part2=$(basename "$p")

scratch="/scratch/tghill/subglacial-emulator/glads/"

RUN="RUN/"

scratchdir=$scratch$part2/

scratchrun=$scratchdir$RUN

DIRECTORY=$scratchdir
if [ -d "$DIRECTORY" ]; then
  echo "Found directory $DIRECTORY"
else
  echo "Making directory $DIRECTORY"
  mkdir "$DIRECTORY"
fi

echo -e "Copying files...\n\n"

#########################################
# Define rules for what files to copy
#########################################
cp -v *.sh "$DIRECTORY"
cp -v *.h "$DIRECTORY"
cp -v *.m "$DIRECTORY"
cp -v *.dat "$DIRECTORY"
cp -vr data/ "$DIRECTORY"

echo -e "\n\nDone copying files...\n\n"

echo "Ready to submit job_script.sh"

while getopts ":r" opt; do
  case "${opt}" in
    r)
      echo "Running job script"
      echo "sbatch job_script.sh"
      cd "$DIRECTORY"
      sbatch job_script.sh
      ;;
    \?) echo "Invalid option: -$OPTARG"
    ;;
  esac
done


