#!/bin/bash

# Init conda
echo eval "$(/work/miniconda3/bin/conda shell.bash hook)" >> /home/ucloud/.bashrc
# source /work/miniconda3/bin/activate

# Init udocker
export UDOCKER_DIR=/work/G65-2017-Kidstage/udocker/.udocker
export PATH=/work/G65-2017-Kidstage/udocker:$PATH
udocker install

# Setup tmux config
# cp /work/sduvarcall/config/tmux.conf ~/.tmux.conf
cp /work/sduvarcall/config/tmux.conf.local ~/.tmux.conf.local

# Add bash settings
alias s1='cd /work/sduvarcall/'
alias le='f(){ zless "$@" | grep -v ^## | less -S; unset -f f; }; f'

# Add inputrc settings
cp /work/sduvarcall/config/inputrc ~/.inputrc

# NGSCheckmate
export NCM_HOME=/work/sduvarcall/G65-2017-Kidstage/NGSCheckMate

# Refresh environment
bash -i
