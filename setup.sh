#!/bin/bash
#### use cvmfs for root ####
#source /cvmfs/fermilab.opensciencegrid.org/packages/common/setup-env.sh ## -- For Alma 9
#source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
#spack load root@6.28.12

export PLOTTER_WORKING_DIR=`pwd`
export FILE_PATH=$PLOTTER_WORKING_DIR/rootfiles/
export PLOT_PATH=$PLOTTER_WORKING_DIR/plots/
export SCRIPT_DIR=$PLOTTER_WORKING_DIR/script/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$PLOTTER_WORKING_DIR/include/:$PLOTTER_WORKING_DIR/plotter/

source $PLOTTER_WORKING_DIR/bin/BashColorSets.sh
