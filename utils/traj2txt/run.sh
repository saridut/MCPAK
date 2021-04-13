#!/usr/bin/env bash

#Driver for program extract.f90.

arch='bb-120-4'
#eps='0.25'
f='0.0625'

traj_home="$HOME/j5/workspace/mc-fd-ev"

traj_dir="$traj_home/$arch/f-$f"

#Number of the first frame
ifrm_beg=1000
#Number of the last frame. -1 indicates the last available frame.
ifrm_end=-1
#Step over how many frames? 1 indicates consecutive frames.
ifrm_stp=200

itraj_beg=1
num_traj=12

#Output file format
fmt='xyz' #{'xyz', 'ldf', 'cfg'}

#Where to output the frames
frame_dir="frmdir/$arch/f-$f"

#If frame_dir exists, clear frame_dir
#If frame_dir does not exist, create frame_dir
if [[ -d "$frame_dir" ]]; then
    rm -rf "$frame_dir"/*
    printf "%s emptied \n" "$frame_dir"
else
    mkdir -p "$frame_dir"
    printf "%s created \n" "$frame_dir"
fi

itraj_end=$(( itraj_beg + num_traj - 1 ))

for (( itraj=$itraj_beg; itraj<=$itraj_end; ++itraj )); do
    printf "itraj %s \n" "$itraj"
    fn_traj="$traj_dir/traj.bin.$itraj"
    fn_cfg="$traj_dir/$arch.cfg.$itraj"

    mkdir $frame_dir/it-$itraj

    #Run
    ./traj2txt "$fn_cfg" "$fn_traj" "$ifrm_beg" "$ifrm_end" "$ifrm_stp"  "$fmt"
    mv frame* $frame_dir/it-$itraj

done

#Collect all frames in one directory and rename them in order
mkdir $frame_dir/collected
ifrmc=0
for (( itraj=1; itraj<=$num_traj; ++itraj )); do
    for each in $frame_dir/it-$itraj/frame*; do
        ifrmc=$(( ifrmc+1 ))
        if [[ $fmt == 'ldf' ]]; then
            cp $each $frame_dir/collected/frame.ldf.$ifrmc
        elif [[ $fmt == 'xyz' ]]; then
            cp $each $frame_dir/collected/frame.xyz.$ifrmc
        elif [[ $fmt == 'cfg' ]]; then
            cp $each $frame_dir/collected/frame.cfg.$ifrmc
        fi
    done
done
