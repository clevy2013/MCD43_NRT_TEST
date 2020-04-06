#!/bin/bash

#can use multiple tiles, days, years, etc
tile=h12v04
year=2019
day1=160
day2=160

export LANCE=1

#########################################
# CHANGE THE INPUT & OUTPUT DIRECTORIES #
#########################################
dir_mod09ga=/muddy/data04/charlotte.levy/inputs/dummy
dir_myd09ga=/muddy/data04/charlotte.levy/inputs/dummy 
#MYD09GA/hdf/2019/${tile}
dir_out=/muddy/data04/charlotte.levy/outputs/DAAC_InHOuse_NRT_test/Esception

#echo $dir_mod09ga

PWD=`pwd`
dir_db=${PWD}/runtime

if [ -z $year ] || [ -z $day1 ]; then
	echo "USAGE: $0 tile year day"
	exit 1
fi

if [ -z $day2 ]; then
	day2=${day1}
fi

############################################################################
#                            SET OUTPUT PCF NAME                           #
############################################################################
#source ~/.bash_profile
dayf1=`printf %03d $day1`
dayf2=`printf %03d $day2`
temp=./PCF_TEMPLATE.pcf
old_pgs=/home/qsun/programs/SDP_ToolKit/TOOLKIT
new_pgs=$PGSHOME
out=${PWD}/runtime/MOD_PR43B_${tile}_${year}_${dayf1}-${dayf2}.pcf

if [ ! -r $dir_out ]; then
	mkdir -p $dir_out
fi
if [ ! -r ./runtime ]; then
	mkdir runtime
fi
############################################################################
#                            FIND BACKUP DATABASE                          #
############################################################################
brdfdb=`find ${dir_db} -name BRDF_DB.*.${tile}.*.hdf`
echo ${dir_db} -name BRDF_DB.*.${tile}.*.hdf
if [ -z $brdfdb ] || [ ! -r $brdfdb ]; then
	./db.sh $tile
	brdfdb=`find ${dir_db} -name BRDF_DB.*.${tile}.*.hdf`
	if [ -z $brdfdb ] || [ ! -r $brdfdb ]; then
		echo "BRDF_DB not exist."
		exit 1
	fi
fi

brdfdb_short=${brdfdb##*/}
brdfdb_path=${brdfdb%/*}
brdfdb="243120|${brdfdb_short}|${brdfdb_path}||${brdfdb_short}||1"
brdfdb_yearly="#243049|MCD43DB_A${year}${dayf1}_${tile}_update.hdf|${dir_out}||||1"

############################################################################
#                  FIND MOD09GA, MYD09GA, AND MOD10A1                      #
############################################################################
source _functions.sh

declare -a arr_ref
declare -a arr_snow

list1=`get_files ${dir_mod09ga} MOD09GA $tile $year $day1 $day2`
list2=`get_files ${dir_myd09ga} MYD09GA $tile $year $day1 $day2`
#list3=`get_files ${dir_mod10ga} MOD10GA $tile $year $day1 $day2`
#list4=`get_files ${dir_myd10ga} MYD10GA $tile $year $day1 $day2`

i=0
for file in $list1; do
	i=$(($i+1))
	_i=`printf %02d $i`
	arr_ref[$i]="213700|${file}||UR_MOD_PR43B.in${_i}|MOD_PR43B.in${_i}|${i}"
done
for file in $list2; do
	i=$(($i+1))
	_i=`printf %02d $i`
	arr_ref[$i]="213700|${file}||UR_MOD_PR43B.in${_i}|MOD_PR43B.in${_i}|${i}"
done

if [ $i -le 6 ]; then
	echo "ERROR! $i REFLECTANCE FILES FOR $tile $year $day."
  exit 3
fi

#echo "DOY $dayf NUMBER OF REFLECTANE FILES: [$i]"

#i=0
#for file in $list3; do
#	i=$(($i+1))
#	_i=`printf %02d $i`
#	arr_snow[$i]="210102|${file}||UR_MOD_PR43B.in${_i}|MOD_PR43B.in${_i}|${i}"
#done
#for file in $list4; do
#	i=$(($i+1))
#	_i=`printf %02d $i`
#	arr_snow[$i]="210102|${file}||UR_MOD_PR43B.in${_i}|MOD_PR43B.in${_i}|${i}"
#done

#if [ $i -le 0 ]; then
#	echo "ERROR! $i SNOW FILES FOR $tile $year $day."
#  exit 3
#fi

# echo "DOY $dayf NUMBER OF SNOW COVER FILES: [$i]"

############################################################################
#                  SET OUTPUT FILENAME AND PATHS                           #
############################################################################


#mcd43a1="243045|MCD43A1.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43a2="243046|MCD43A2.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43a3="243047|MCD43A3.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43a4="243048|MCD43A4.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43b1="243200|MCD43B1.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43b2="243201|MCD43B2.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43b3="243220|MCD43B3.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43b4="243230|MCD43B4.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43t1="243100|MCD43T1.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"
#mcd43t2="243101|MCD43T2.A${year}${dayf}.${tile}.hdf|${dir_out}||||1"

log0="10100|LogStatus_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"
log1="10101|LogReport_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"
log2="10102|LogUser_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"
log3="10103|TmpStatus_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"
log4="10104|TmpReport_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"
log5="10105|TmpUser_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"
log6="10110|MailFile_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"
log7="10111|ShmMem_${tile}_${year}_${dayf1}-${dayf2}.txt|${PWD}/runtime||||1"

tmp1="10252|GetAttr_${tile}_${year}_${dayf1}-${dayf2}.temp|${PWD}/runtime||||1"
tmp2="10254|MCFWrite_${tile}_${year}_${dayf1}-${dayf2}.temp|${PWD}/runtime||||1"

############################################################################
#                           ASSEMBLE NEW PCF                               #
############################################################################
### HEADERS
cat $temp | awk -F "|" '{if($1 == "243120"){print "'$brdfdb'"} else if($1 == "213700"){exit} else{print $0}}' > $out

### REF FILES
i=${#arr_ref[@]}
while [ $i -ge 1 ]; do
	echo ${arr_ref[$i]} >> $out
	i=$((${i}-1))
done

### SNOW DATA
echo "#" >> $out
echo "# MOD10A1 500m snow cover" >> $out
	
i=${#arr_snow[@]}
while [ $i -ge 1 ]; do
	echo ${arr_snow[$i]} >> $out
	i=$((${i}-1))
done

## outputs
echo "#" >> $out
echo "# OUTPUTS" >> $out
n=$((${day2}-${day1}+1))
dd=${day1}
while [ $dd -le $day2 ]; do
	dayf=`printf %03d $dd`
	echo "243045|MCD43A1.A${year}${dayf}.${tile}.hdf|${dir_out}||||${n}" >> $out
	dd=$((${dd}+1))
	n=$((${n}-1))
done
	
n=$((${day2}-${day1}+1))
dd=${day1}
while [ $dd -le $day2 ]; do
	dayf=`printf %03d $dd`
	echo "243046|MCD43A2.A${year}${dayf}.${tile}.hdf|${dir_out}||||${n}" >> $out
	dd=$((${dd}+1))
	n=$((${n}-1))
done
	
n=$((${day2}-${day1}+1))
dd=${day1}
while [ $dd -le $day2 ]; do
	dayf=`printf %03d $dd`
	echo "243047|MCD43A3.A${year}${dayf}.${tile}.hdf|${dir_out}||||${n}" >> $out
	dd=$((${dd}+1))
	n=$((${n}-1))
done
	
n=$((${day2}-${day1}+1))
dd=${day1}
while [ $dd -le $day2 ]; do
	dayf=`printf %03d $dd`
	echo "243048|MCD43A4.A${year}${dayf}.${tile}.hdf|${dir_out}||||${n}" >> $out
	dd=$((${dd}+1))
	n=$((${n}-1))
done
	
n=$((${day2}-${day1}+1))
dd=${day1}
while [ $dd -le $day2 ]; do
	dayf=`printf %03d $dd`
	echo "243100|MCD43T1.A${year}${dayf}.${tile}.hdf|${dir_out}||||${n}" >> $out
	dd=$((${dd}+1))
	n=$((${n}-1))
done
	
n=$((${day2}-${day1}+1))
dd=${day1}
while [ $dd -le $day2 ]; do
	dayf=`printf %03d $dd`
	echo "243101|MCD43T2.A${year}${dayf}.${tile}.hdf|${dir_out}||||${n}" >> $out
	dd=$((${dd}+1))
	n=$((${n}-1))
done
	

### ENDER
#cat $temp | awk -v tag1=0 -v tag2=0 -v tag3=0 -v tag4=0 -v old_pgs=$old_pgs -v new_pgs=$new_pgs -F "|" '{if($1 == "213700"){tag1=1}; if(tag1==1 && $1 != "213700") {tag2=1}; if(tag2==1 && $1 == "210101"){tag3=1}; if(tag3==1 && $1 != "210101"){tag4=1}; if($1 == "243045"){print "'$mcd43a1'"} else if($1 == "243046"){print "'$mcd43a2'"} else if($1 == "243047"){print "'$mcd43a3'"} else if($1 == "243048"){print "'$mcd43a4'"} else if($1 == "243200"){print "'$mcd43b1'"} else if($1 == "243201"){print "'$mcd43b2'"} else if($1 == "243220"){print "'$mcd43b3'"} else if($1 == "243230"){print "'$mcd43b4'"} else if($1 == "243100"){print "'$mcd43t1'"} else if($1 == "243101"){print "'$mcd43t2'"} else if($1 == "10100"){print "'$log0'"} else if($1 == "10101"){print "'$log1'"} else if($1 == "10102"){print "'$log2'"} else if($1 == "10103"){print "'$log3'"} else if($1 == "10104"){print "'$log4'"} else if($1 == "10105"){print "'$log5'"} else if($1 == "10110"){print "'$log6'"} else if($1 == "10111"){print "'$log7'"} else if($1 == "10252"){print "'$tmp1'"} else if($1 == "10254"){print "'$tmp2'"} else if($1 == "243049"){print "'$brdfdb_yearly'"} else if(tag2 == 1 && tag4 == 1){if(match($3, old_pgs) != 0){sub(old_pgs, new_pgs)};print $0}}' >> $out

cat $temp | awk -v tag1=0 -v tag2=0 -v tag3=0 -v tag4=0 -v old_pgs=$old_pgs -v new_pgs=$new_pgs -F "|" '{if($1 == "213700"){tag1=1}; if(tag1==1 && $1 != "213700") {tag2=1}; if(tag2==1 && $1 == "210101"){tag3=1}; if(tag3==1 && $1 != "210101"){tag4=1}; if($1 == "243045"){print "#"} else if($1 == "243046"){print "#"} else if($1 == "243047"){print "#"} else if($1 == "243048"){print "#"} else if($1 == "243200"){print "#"} else if($1 == "243201"){print "'$mcd43b2'"} else if($1 == "243220"){print "'$mcd43b3'"} else if($1 == "243230"){print "'$mcd43b4'"} else if($1 == "243100"){print "#"} else if($1 == "243101"){print "#"} else if($1 == "10100"){print "'$log0'"} else if($1 == "10101"){print "'$log1'"} else if($1 == "10102"){print "'$log2'"} else if($1 == "10103"){print "'$log3'"} else if($1 == "10104"){print "'$log4'"} else if($1 == "10105"){print "'$log5'"} else if($1 == "10110"){print "'$log6'"} else if($1 == "10111"){print "'$log7'"} else if($1 == "10252"){print "'$tmp1'"} else if($1 == "10254"){print "'$tmp2'"} else if($1 == "243049"){print "'$brdfdb_yearly'"} else if(tag2 == 1 && tag4 == 1){if(match($3, old_pgs) != 0){sub(old_pgs, new_pgs)};print $0}}' >> $out

############################################################################
#                           RUN THE CODE                                   #
############################################################################
export PGS_PC_INFO_FILE=$out
time ${PWD}/MOD_PR43B.exe

exit 0
