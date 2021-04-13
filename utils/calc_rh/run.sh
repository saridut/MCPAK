#!/usr/bin/env bash

#Driver for program calc_rh.f90.

nabb=120
nasc=4

#traj_home="$HOME/workspace/dev/brushpak-mc/bin"
traj_home="$HOME/j5/workspace/mc-fd-ev/bb-$nabb-$nasc/f-0.5"

traj_dir="$traj_home"

fn_cfg="$traj_dir/bb-${nabb}-${nasc}.cfg.1"
fn_traj="traj.bin"
num_traj=12
nts_beg=100

#Run
./calc_rh "$fn_cfg" "$fn_traj" "$traj_dir" "$num_traj" "$nts_beg" 
