#!/bin/bash

# Looks through given parameter directories looking for a "allrun.sh" file. 
# If found, tries to run the code in "DO STUFF" section in that directory. 
#
# Usefull for pushing cases in queues, cleaning, pre/post-processing etc.
# of OpenFOAM cases.
#
# Made by Antti Mikkonen, 2016, antti.mikkonen@iki.fi
#
# Feel free to use for whatever purpose.



echo Start

calledFrom=$(pwd)

################################
# Check parameters
################################
if [ "$1" == "" ]; then
  echo "$0: Please provide a directory name"
  exit 1
fi

for thisTree in "$@"
do
    cd $calledFrom
    if [ ! -d "$thisTree" ]; then
      echo "$thisTree is not a directory name"
      exit 1
    fi
done

################################
# Loop parameters
################################
for thisTree in "$@"
do
    echo
    echo "$thisTree"
    cd $calledFrom

    cd $thisTree
    base=$(pwd)
    echo $base
    cd $base 

    paths=$(find -name allrun.sh | sort);

    ################################
    # Loop directories with a "allrun.sh"
    ################################
    for fullpath in $paths
    do
        parDir=$(dirname "${fullpath}");
        cd $parDir;
        echo $parDir;

        ################################
        # DO STUFF
        ################################
           
        # tsp
        #tsp -L $(pwd) bash solver.sh
        
        # slurm push
#        sbatch pushSlurm.sh 

        
#        # Clean
#        rm -rf processor* postProcessing/ [0-9].[0-9]* [1-9]* log/solver 
#        
#        # Map
#        sourceDir=$calledFrom/mykomegaSST/$(echo $parDir | sed -e "s/Pimple/Simple/g")
#        mapFields -consistent -sourceTime latestTime $sourceDir
#        find -name rho -exec rm -rf {} \;
#        find -name heatTransferCoefficient -exec rm -rf {} \;
#        find -name heatFlux -exec rm -rf {} \;
#        
#        decompose
##        cp $calledFrom/controlDict system/controlDict
#        decomposePar
#        cp $calledFrom/controlDict system/controlDict
    
        # label 
#        sed -i "s@$parDir@$thisTree/$parDir@g" pushSlurm.sh
        

        ################################
        # RETURN TO ORIGIN
        ################################
        cd $base;
    done
done

echo
echo Done
