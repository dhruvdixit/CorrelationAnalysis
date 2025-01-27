#!/bin/bash

source ~/.bashrc.ext

module add ROOT/6.08.00
module load hdf5
which gcc
which root
module list
date
echo "Command runpairing with file $1.root, mix_start = $2, mix_end = $3, and TrackSkim GeV = $4"

export cwd=$(pwd)
cd $cwd/../pair_gale_shapley/
./mix_gale_shapley ../InputData/$1.root ../InputData/$1_minbias_$4GeVTracks.root $2 $3 $4
date 


#The first argument can be replaced with full file name. Written assuming use on 13d,13e, and 13f (separately)