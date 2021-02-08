#!/bin/bash

dir_macros=$(dirname $(readlink -f $BASH_SOURCE))
LIFE_TIME=medium # short (3h), medium (8h) or long (23h)

jobname=$1
do_sub=$2
njobs=$3
nevents=$4

echo "njobs=$njobs"
echo "nevents=$nevents"
if [ $do_sub == 1 ]; then
    echo "Grid mode."
    if ! which jobsub_submit &>/dev/null ; then
	echo "Command 'jobsub_submit' not found."
	echo "Forget 'source /e906/app/software/script/setup-jobsub-spinquest.sh'?"
	exit
    fi
    work=/pnfs/e906/persistent/users/$USER/SimChainDev/$jobname
    ln -sf /pnfs/e906/persistent/users/$USER/SimChainDev data
else
    echo "Local mode."
    work=$dir_macros/scratch/$jobname
fi

mkdir -p $work
chmod -R 01755 $work

cd $dir_macros
tar -czvf $work/input.tar.gz *.C *.cfg *.opts
cd -

for (( id=1; id<=$njobs; id++ ))
do  
  mkdir -p $work/$id/log
  mkdir -p $work/$id/out
  chmod -R 01755 $work/$id

  rsync -av $dir_macros/gridrun.sh $work/$id/gridrun.sh

  if [ $do_sub == 1 ]; then
    cmd="jobsub_submit"
    cmd="$cmd -g --OS=SL7 --use_gftp --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE -e IFDHC_VERSION --expected-lifetime='$LIFE_TIME'"
    cmd="$cmd --mail_never"
    cmd="$cmd -L $work/$id/log/log.txt"
    cmd="$cmd -f $work/input.tar.gz"
    cmd="$cmd -d OUTPUT $work/$id/out"
    cmd="$cmd --append_condor_requirements='(TARGET.GLIDEIN_Site isnt \"UCSD\")'"
    cmd="$cmd file://`which $work/$id/gridrun.sh` $nevents $id"
    echo "$cmd"
    $cmd
  else
    mkdir -p $work/$id/input
    rsync -av $work/input.tar.gz $work/$id/input
    cd $work/$id/
    $work/$id/gridrun.sh $nevents $id | tee $work/$id/log/log.txt
    cd -
  fi
done 2>&1 | tee log_gridsub.txt

## When your job fails due to bad grid nodes,
## you can use the following option to exclude those nodes;
##   cmd="$cmd --append_condor_requirements='(TARGET.GLIDEIN_Site isnt \"UCSD\")'"
## Valid site names are listed here;
## https://cdcvs.fnal.gov/redmine/projects/fife/wiki/Information_about_job_submission_to_OSG_sites
## According to the Fermilab Service Desk, the "--blacklist" option has a known defect.
