#!/bin/bash
#
./OptJV3 --jarvis ./JARVIS3 --input cy_sample --population 32 --generations 10 --threads 4 --elite 4 --tournament 3 --crossover 0.85 --mutation 0.12 --toggle 0.08 --objective bytes --max-cmodels 2 --min-cmodels 2 --max-rmodels 1 --min-rmodels 1 --global-bounds "hs=8:32,lr=0.0:0.15,seed=1:1000" --cm-bounds "ctx=1:13,den=1:200,ir=0:2,gamma=0.01:0.99,edits=0:3,eden=1:50,eir=0:1,egamma=0.01:0.99" --rm-bounds "nr=1:20,ctx=11:13,beta=0.01:0.99,limit=1:20,gamma=0.01:0.99,ir=0:1,weight=0.01:0.99,cache=1:4" --best-out best.txt --history-out history.csv
