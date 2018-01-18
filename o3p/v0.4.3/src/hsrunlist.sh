#!/bin/bash
# use this when all files are put in one directory
# origin : D. Shin.
# updated : hs

type=practice_list

make_on=1       # 0: run without make   1: make and run 
use_qsub=0      # whether to use qsub or not (1: use qsub,  0: command line, else: just test run without retrievals)
ver=3           # sol, cld version
nmpi=4          # number of core for computing with mpi

# set environment vairables
source ../run/gems_env.sh   

# make source code
if [ $make_on -eq 1 ] ; then
  make
fi

DATDIR=/home/Data/OMI
PGEHOME=/home/o3p_hs/GEMS/o3p/v0.4.3
pwd=`pwd`  

FILE=runlist.list

while read line
do

  timesel=(`echo "$line" | cut -c1-12`)
  sline=(`echo "$line" | cut -c14-17`)
  eline=(`echo "$line" | cut -c18-23`)
  spix=(`echo "$line" | cut -c24-29`)
  epix=(`echo "$line" | cut -c30-35`)

  spix=$((${spix} * 2 - 1))
  epix=$((${epix} * 2))
  pixsel="${sline} ${eline} ${spix} ${epix}"
  ##pixsel='198 198 29 30'      
  #timesel='2008m0602t03'
  echo $timesel   $pixsel

  IFS=' ' read -r -a pixline <<< "$pixsel"
  ((pixline[2]= ( ${pixline[2]} + 1 ) / 2))
  ((pixline[3]= ( ${pixline[3]} + 1 ) / 2))

  indir='/home/o3p_hs/GEMS/o3p/dat/in/'
  outdir='/home/o3p_hs/GEMS/o3p/dat/out/'
  yr=${timesel:0:4}
  lv2dir=$pwd/
  raddb=$DATDIR'/1_OML1BRUG/'$yr'/'
  clddir=$DATDIR'/2_OML2CLDO2/'$yr'/'
  soldir=$DATDIR/'/1_OML1BIRR/'$yr'/'
  lv2logprefix=OMIO3PROF_$type'-'
  solprefix=OMI-Aura_L1-OML1BIRR_
  radprefix=OMI-Aura_L1-OML1BRUG_
  cldprefix=OMI-Aura_L2-OMCLDO2_
  nmldir=/home/o3p_hs/GEMS/share/conf
  lv1_nmlorg=${nmldir}/l1breadmdl_v01.nml
  lv1_nmlfile=${nmldir}/l1breadmdl.nml
  lv2_nmlorg=${nmldir}/lv2readmdl_org.nml
  lv2_nmlfile=${nmldir}/lv2readmdl.nml
  gems_conforg=${nmldir}/gems_v01.conf
  gems_confile=${nmldir}/gems.conf

  # initialized the orbits to be arrays
  # note the index starts from zero
  #orbs[0]=''
  #norb=0
  #rm -f orbs.list

  for radfname in `ls $raddb$radprefix$timesel*he4`
  do
      l1bfname=${radfname/$raddb/$indir}  
      cp $radfname ${l1bfname}
      chmod 644 ${l1bfname}
      echo '::  Copy L1BRUG from DB' 
      h4toh5 ${l1bfname}
      echo '::  Conversion he4 to h5' 
      rm -f ${l1bfname}
      l1bfname=${l1bfname/he4/h5}

      orb=${l1bfname#*-o}  
      orb=${orb:0:5}
      mondayutc=${l1bfname#*$radprefix}
      mondayutc=${mondayutc:0:14}
      monday=${mondayutc:0:9}

     #solfname=`ls $soldir$solprefix$monday*he4`
      cldfname=`ls $clddir$cldprefix*o$orb*he5`
      lv2logname=$lv2dir$lv2logprefix'o'$orb
      lv2fname=GEM_TEST_L2_O3P_${monday}_o${orb}_P${pixline[2]}-${pixline[3]}_L${pixline[0]}-${pixline[1]}_mpi${nmpi}.h5
      #lv2fname=GEM_TEST_L2_O3P_${monday}_o${orb}_mpi${nmpi}.h5

      # Solar radiance could not be found, make up a name in case use backup radiance
      if [ -z $solfname ]
          then 
          #continue
          solfname=$soldir$solprefix$mondayutc'-o'$orb'_v00'$ver'.he4'
      fi

      # fake a cloud file name (not used anyway for the moment)
      if [ -z $cldfname ]
          then
          cldfname=$clddir$cldprefix$mondayutc'-o'$orb'_v00'$ver'.he5'
      fi

      # update lv1 namelist file     
      linenum=0
      while read line ; do
        if [ ${linenum} -eq 0 ]  ; then echo ${line} > ${lv1_nmlfile}
        elif [ ${linenum} -eq 6 ] ; then echo l1b_rug_file_path = '"'${l1bfname}'"' >> ${lv1_nmlfile}
        else echo ${line} >> ${lv1_nmlfile}
        fi
        linenum=$((${linenum}+1))
      done < ${lv1_nmlorg}


      # update lv2 namelist file     
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
       echo LV2LOG $lv2logname

      #orbs[$norb]=$orb
      #norb=`expr $norb + 1`

      # Write input file INP/L1L2_fnames.inp
      # use echo -e to tell that there are control variable new line (\x0a )
      echo -e "$solfname \x0a$radfname \x0a$cldfname \x0a$lv2logname \x0a$pixsel" > $PGEHOME/run/conf/INP/L1L2_fnames.inp


      # exit orbit directory
      #cd $pwd
      #cd $GEMSBIN

      if [ $use_qsub -eq 0 ] 
          then            
          #$PGEHOME/bin/${SAOPGE}_exec  > $lv2dir$type$orb'.dat' 
          mpirun -np $nmpi $PGEHOME/bin/${SAOPGE}_exec  > $lv2dir$type$orb'_mpi'$nmpi'.dat' </dev/null
          #$PGEHOME/bin/${SAOPGE}_exec 
          cd $pwd
          #cat $type$orb'.dat'
          #sleep 2
      fi

  done

done < "${FILE}"


