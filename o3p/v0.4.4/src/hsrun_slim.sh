#! /bin/bash

type=SLIM

# Initialize
orbs[0]=''
norb=0
rm -f orbs.list

nmpi=4
ver=3
pixset=0    #0: no identification of pixline, 1: with identification of pixline
use_qsub=0
lidort_make=0

make_on=1 

# 0. Set env
source ../run/gems_env.sh

DATDIR=/home/Data/OMI
export O3P_VER=v0.4.4
export PGEHOME=/home/o3p_hs/GEMS/o3p/$O3P_VER
export PATH=$PATH0:$O3P_HOME/$O3P_VER/run
export PATH=$PATH:/opt/hdf/h4h5tools-2.2.2-linux-static/bin:$O3P_HOME/$O3P_VER/bin:.


SYNDIR=/home/Data/SYNTHETIC_GEMS_201902
SYNTYPE=SCALAR_PART
pwd=`pwd`/

# make lidort code
if [ $lidort_make -eq 1 ] ; then
  cd /home/o3p_hs/GEMS/o3p/share/lidort
  make
  cd ${pwd}
fi

# ' line1, line2, cross pixel 1, cross pixel 2 '
pixsel='-1 -1 -1 -1'
#'1147 1 58' 
#pixsel='880 950 1 60' 
timesel='2016m0115t03'
IFS=' ' read -r -a pixline <<< "$pixsel"
#((pixline[2]= ( ${pixline[2]} + 1 ) / 2))
#((pixline[3]= ( ${pixline[3]} + 1 ) / 2))

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
#soldir=$DATDIR'/1_OML1BIRR/'$yr'/'
soldir=/home/Data/SYNTHETIC_GEMS/

solprefix=GEMS_
radprefix=OMI-Aura_L1-OML1BRUG_
cldprefix=OMI-Aura_L2-OMCLDO2_
logprefix=OMIO3PROF_$type'-'

nmldir=/home/o3p_hs/GEMS/share/conf

lv1_nmlorg=${nmldir}/l1breadmdl_v01.nml
lv1_nmlfile=${nmldir}/l1breadmdl.nml
lv2_nmlorg=${nmldir}/lv2readmdl_org.nml
lv2_nmlfile=${nmldir}/lv2readmdl.nml
gems_conforg=${nmldir}/gems_v01.conf
gems_confile=${nmldir}/gems.conf

l1btmp0=(`ls $raddb$radprefix${timesel:0:9}*he4`)
l1btmp0=${l1btmp0[0]}
l1btmp=${l1btmp0/$raddb/$raddir}
l1btmp=${l1btmp/he4/h5}

orb=${l1btmp#*-o}
orb=${orb:0:5}
mondayutc=${l1btmp#*$radprefix}
mondayutc=${mondayutc:0:14}
monday=${mondayutc:0:9}

#radname=`ls ${SYNDIR}/${SYNTYPE}/rad/GEMS_${monday:0:4}${monday:5:4}${timesel:10:11}*.nc`
radname=`ls ${SYNDIR}/${SYNTYPE}/GEMS_${monday:0:4}${monday:5:4}${timesel:10:11}/rad/GEMS_${monday:0:4}${monday:5:4}${timesel:10:11}*.nc`

for l1bfname in ${radname}
do
  #echo $l1bfname
  if [ ${l1bfname:90:3} -eq 015 ] ;then

    l1bnum=`echo ${l1bfname} |rev| cut -c4-6 | rev`

    #solfname=`ls $soldir$solprefix$monday*he4`
    cldfname=`ls $clddir$cldprefix*o$orb*he5`
    logfname=$pwd$logprefix'o'$orb'_'P${pixline[2]}-${pixline[3]}_L${pixline[0]}-${pixline[1]}'.dat'


    lv2fname=GEM_SYNT_L2_O3P_${l1bnum}_${monday}_o${orb}_mpi${nmpi}.h5
    if [ $pixset -eq 1 ]; then
      lv2fname=GEM_SYNT_L2_O3P_${monday}_o${orb}_P${pixline[2]}-${pixline[3]}_L${pixline[0]}-${pixline[1]}_mpi${nmpi}.h5
    fi
    
    #Solar radiance could not be found, make up a name in case use backup radiance
    if [ -z $solfname ] ;then
      #echo No solar radiance exist. 
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

    #update lv2 namelist file
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
    echo L2FNAME  $lv2fname 
    echo LOGFNAME $logfname

    orbs[$norb]=$orb
    norb=`expr $norb + 1`

    # Write input file INP/L1L2_fnames.inp
    # use echo -e to tell that there are control variable new line (\x0a )
    echo -e "$solfname \x0a$l1bfname \x0a$cldfname \x0a$logfname \x0a$pixsel" > $PGEHOME/run/conf/INP/L1L2_fnames.inp

    # exit orbit directory
    cd $pwd
    cd $GEMSBIN

    #then
    #$PGEHOME/bin/${SAOPGE}_exec  > $lv2dir$type$orb'.dat'
    echo SAOPGE = $SAOPGE
    #echo $PGEHOME/bin/${SAOPGE}_exec #> $logfname
    mpirun -np $nmpi $PGEHOME/bin/${SAOPGE}_exec > $logfname 
    #$PGEHOME/bin/${SAOPGE}_exec > $logfname 

    #$PGEHOME/bin/${SAOPGE}_exec > $logfname
    cd $pwd

    #cat $type$orb'.dat'
    #sleep 2
  fi
done
