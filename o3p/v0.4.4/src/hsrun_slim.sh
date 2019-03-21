#! /bin/bash

type=SLIM

# Initialize
orbs[0]=''
norb=0
rm -f orbs.list

nmpi=2
ver=3
pixset=1    #0: no identification of pixline, 1: with identification of pixline
use_qsub=0
lidort_make=0

make_on=1 

# 0. Set env
source ../run/gems_env.sh

DATDIR=/home/Data/OMI
PGEHOME=/home/o3p_hs/GEMS/o3p/v0.4.4
pwd=`pwd`/


# ' line1, line2, cross pixel 1, cross pixel 2 '
pixsel='871 872 29 30'
#'1147 1 58' 
#pixsel='880 950 1 60' 
timesel='2006m0613t03'
IFS=' ' read -r -a pixline <<< "$pixsel"
((pixline[2]= ( ${pixline[2]} + 1 ) / 2))
((pixline[3]= ( ${pixline[3]} + 1 ) / 2))

echo " pixline = ${pixline[@]}"

#year
yr=${timesel:0:4}

echo " year = ${timesel:0:4}"


# make lidort code
if [ $lidort_make -eq 1 ] ; then
  cd /home/o3p_hs/GEMS/o3p/share/lidort
  make
  cd ${pwd}
fi

# make
if [ $make_on -eq 1 ] ; then
  make
fi

#Input OMI lv1brug he4 data & output directory
indir='/home/o3p_hs/GEMS/o3p/dat/in/'
outdir='/home/o3p_hs/GEMS/o3p/dat/out/'
raddb=$DATDIR'/1_OML1BRUG/'$yr'/'
clddir=$DATDIR'/2_OML2CLDO2/'$yr'/'
soldir=$DATDIR'/1_OML1BIRR/'$yr'/'
nmldir=/home/o3p_hs/GEMS/share/conf

lv1_nmlorg=${nmldir}/l1breadmdl_v02.nml
lv1_nmlfile=${nmldir}/l1breadmdl.nml

lv2_nmlorg=${nmldir}/lv2readmdl_org.nml
lv2_nmlfile=${nmldir}/lv2readmdl.nml

gems_conforg=${nmldir}/gems_v01.conf
gems_confile=${nmldir}/gems.conf

solprefix=OMI-Aura_L1-OML1BIRR_
radprefix=OMI-Aura_L1-OML1BRUG_
cldprefix=OMI-Aura_L2-OMCLDO2_
logprefix=OMIO3PROF_$type'-'


for radfname in `ls $raddb$radprefix$timesel*he4`
do
  
  l1bfname=${radfname/$raddb/$indir} 
  cp $radfname ${l1bfname}
  chmod 644 ${l1bfname}
  echo ':: Copy L1BRUG from DB'
  h4toh5 ${l1bfname}
  echo ':: Conversion he4 to h5'
  rm -f ${l1bfname}
  l1bfname=${l1bfname/he4/h5}

  echo $l1bfname
  orb=${l1bfname#*-o}
  orb=${orb:0:5}  #orbit number
  mondayutc=${l1bfname#*$radprefix}
  mondayutc=${mondayutc:0:14} #year+month+day+time
  monday=${mondayutc:0:9} #year+month+day
  
  solfname=`ls $soldir$solprefix$monday*he4`
  cldfname=`ls $clddir$cldprefix*o$orb*he5`
  logfname=$pwd$logprefix'o'$orb'_'P${pixline[2]}-${pixline[3]}_L${pixline[0]}-${pixline[1]}'.dat'


  # File name set. 
  if [ $pixset -eq 0 ] ; then  
    lv2fname=GEM_TEST_L2_O3P_${monday}_o${orb}_mpi${nmpi}.h5
  elif [ $pixset -eq 1 ]; then
    lv2fname=GEM_TEST_L2_O3P_${monday}_o${orb}_P${pixline[2]}-${pixline[3]}_L${pixline[0]}-${pixline[1]}_mpi${nmpi}.h5
  else
    echo Pixel set error!
    exit
  fi
  
  #Solar radiance could not be found, make up a name in case use backup radiance
  if [ -z $solfname ] ;then
    echo No solar radiance exist. 
    solfname=$soldir$solprefix$mondayutc'-o'$orb'_v00'$ver'.he4'
  fi

  #fake a cloud file name
  if [ -z $cldfname ] ;then
    echo No cload data exist.
    cldfname=$clddir$cldprefix$mondayutc'-o'$orb'_v00'$ver'.he5' 
  fi

  #update lv1 namelist file
  linenum=0
  while read line ; do
    if [ ${linenum} -eq 0 ] ; then echo ${line} > ${lv1_nmlfile}
    elif [ ${linenum} -eq 6 ] ; then echo l1b_rug_file_path = '"'${l1bfname}'"' >> ${lv1_nmlfile}
    else echo ${line} >> ${lv1_nmlfile}
    fi
    linenum=$((${linenum}+1))
  done < ${lv1_nmlorg}

  ##update lv2 namelist file
  linenum=0
  while read line ; do
    if [ ${linenum} -eq 0 ]  ; then echo ${line} > ${lv2_nmlfile}
    elif [ ${linenum} -eq 143 ] ; then echo l2_cld_file_path = '"'${cldfname}'"' >> ${lv2_nmlfile}
    else echo ${line} >> ${lv2_nmlfile}
    fi
    linenum=$((${linenum}+1))
  done < ${lv2_nmlorg}


  # update gems.conf
  linenum=0
  while read line ; do
    if [ ${linenum} -eq 0 ]  ; then echo ${line} > ${gems_confile}
    elif [ ${linenum} -eq 67 ] ; then echo gvc_out_lv2_fname    = '"'${lv2fname}'"' >> ${gems_confile}
    else echo ${line} >> ${gems_confile}
    fi
    linenum=$((${linenum}+1))
  done < ${gems_conforg}


  #echo SOLFNAME $solfname
  echo L1BFNAME $l1bfname
  echo CLDFNAME $cldfname
  echo LOGFNAME $logfname

  orbs[$norb]=$orb
  norb=`expr $norb + 1`

  # Write input file INP/L1L2_fnames.inp
  # use echo -e to tell that there are control variable new line (\x0a )
  echo -e "$solfname \x0a$radfname \x0a$cldfname \x0a$logfname \x0a$pixsel" > $PGEHOME/run/conf/INP/L1L2_fnames.inp


  # exit orbit directory
  #cd $pwd
  #cd $GEMSBIN

  #then
  #$PGEHOME/bin/${SAOPGE}_exec  > $lv2dir$type$orb'.dat'
  mpirun -np $nmpi $PGEHOME'/bin/GEMS_O3P_exec' > $logfname 

  #$PGEHOME/bin/${SAOPGE}_exec > $logfname
  cd $pwd

  #cat $type$orb'.dat'
  #sleep 2

done



#cd ../src
#make
#cd ../bin
#./GEMS_O3P_exec

