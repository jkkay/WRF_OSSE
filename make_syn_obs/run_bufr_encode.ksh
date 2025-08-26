#!/bin/ksh

export WORK_DIR=/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools
export OBS_DATE=2015071500
export export OBS_TYPES="ADPUPA ADPSFC SATWND PROFLR VADWND SFCSHP AIRCFT"
#export OBS_TYPE=VADWND

export OBS_YYYYMMDD=`echo $OBS_DATE | cut -b 1-8`
export OBS_HH=`echo $OBS_DATE | cut -b 9-10`
echo $OBS_YYYYMMDD

cd $WORK_DIR

for OBS_TYPE in $OBS_TYPES; do
cat > namelist.conv << EOF
 &iosetup
  infilename='/glade/work/junkyung/data/obs/LIUZ/${OBS_YYYYMMDD}/gdas1.t${OBS_HH}z.prepbufr.nr',
  outfilename_00='V00_gdas1',
  outfilename_m3='Vm3_gdas1',
  outfilename_p3='Vp3_gdas1',
  OBS_TYP='$OBS_TYPE'
  OBS_DATE='$OBS_DATE'
 /
EOF

#echo ${WORK_DIR}/${OBS_DATE}_${OBS_TYPE}_obs_num.txt
if [[ -f ${WORK_DIR}/${OBS_DATE}_${OBS_TYPE}_obs_num.txt ]]; then
 #rm -f ${WORK_DIR}/${OBS_DATE}_${OBS_TYPE}_obs_num.txt 
 echo ${WORK_DIR}/${OBS_DATE}_${OBS_TYPE}_obs_num.txt 
fi

date
./bufr_encode_sample.x > logs_encode_${OBS_TYPE}_${OBS_DATE} 2>&1
date

done
