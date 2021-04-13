#!/usr/bi/env bash

for na_bbone in 120; do
    for na_sc in 4 ; do
        dir_loc="bb-$na_bbone-$na_sc/f-0.5"
        mkdir -p $dir_loc
        #echo $dir_loc
        for icfg in {1..12}; do
            #fn_coords="$HOME/workspace/projects/dmref/rheol/rev2txt/traj.xyz.$icfg"
            python create-chn-br.py $na_bbone $na_sc #$fn_coords
            mv "bb-$na_bbone-$na_sc.cfg" $dir_loc/bb-$na_bbone-$na_sc.cfg.$icfg
        done
    done
done
