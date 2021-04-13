#!/usr/bin/env bash

##$ -l s_rt=12:00:00

#Shell
#$ -S /bin/bash

#Combine STDERR and STDOUT
#$ -j y

#Send mail on abort(a), begin(b), and end(e)
#$ -M dutta6@illinois.edu
##$ -m a

#Select queue
#$ -q xeon1.q

#Run job from the directory from which it was submitted
#$ -cwd

#Set priority
#$ -p 0

#Job name
#$ -N mcp

#Submit array job
#$ -t 1-6

#Directory where the executable is located
#MP_DIR="$SGE_0_HOME/dev/mcpak-0.2/bin"
MP_DIR="../../bin"

#Run the executable
$MP_DIR/mcpak fn_control=control.txt job_tag=1
#$MP_DIR/mcpak fn_control=control.txt job_tag=$SGE_TASK_ID
