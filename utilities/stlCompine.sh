#!/bin/bash 
# 
# Original Script by "Kruno". Downladed 09.02.2015 from
#http://openfoamwiki.net/index.php/Thread:Talk:Sig_Numerical_Optimization/Salome_and_cfMesh_parametrization   
# Modified by Antti Mikkonen, a.mikkonen@iki.fi winter 2015
#
# how to run it: myFoamStlCompine.sh [PATH TO DIRECTORY WITH .stl FILES] [PATH TO OUTPUT FILE] 
# example with ahmed is 
# myFoamStlCompine.sh ./ahmed/ahmed_25/constant/triSurface ./ahmed/ahmed_25/merged.stl
 
 # get .stl files location/path specified with running script 
filePath=$(cd $1 && pwd) 
 # get merged file output location/path specified with running script 
outputPath=$2   #$(cd $2 && pwd) 

 # check if output file  exists. If it does delete it. 
 if [ -f $outputPath ] 
 then 
     rm $outputPath 
 fi 

 # remove all temporary files from input directory (where .stl files are located) 
 if [ -d "$filePath/*~" ]; then 
     rm $filePath/*~ 
 fi 
 
 # add patch name: 
   # change first and last line of every file in input location from solid to solid [.stl FILE NAME] and from endsolid to endsolid [.stl FILE NAME]
   # merge them together
 # example is from "solid" to "solid Bottom" and from "endsolid" to "endsolid Bottom" 
 for file in $filePath/*.stl 
 do 
     # get .stl file name
     fileName=`basename $file` 
     echo -e "\tReading $fileName"
     # remove .stl file extension
     patchName=`echo "$fileName" | cut -f1 -d'.'` 
     
     # append line solid [.stl FILE NAME] (eg. solid Bottom ) to file
     echo "solid $patchName" >> $outputPath
     
     # copy all lines except first and last one from loaded .stl file and append
     tail -n +2 $file | head -n -1 >> $outputPath 
     
     # append line endsolid [.stl FILE NAME] (eg. endsolid Bottom ) 
     echo "endsolid $patchName" >> $outputPath 
 done
 
 echo -e "Wrote $outputPath"
