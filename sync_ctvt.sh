#!/usr/bin/env bash

rsync -avz --partial \
           --exclude '.git' \
           --exclude '.idea' \
           --exclude 'nohup*' \
           --exclude 'sam' \
           --exclude '*.py*' \
           --exclude '*.R' \
           --exclude '*.old*' \
           --exclude '_laurent_' \
           --exclude 'qpgraph/dot' \
           evan@palaeoprime:~/ctvt/ ~/Dropbox/Code/ctvt/