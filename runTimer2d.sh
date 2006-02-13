#!/bin/sh

echo "Running 2d timing tests"

LAB=./testLabelling
LOG2D=timing2d.log

export CONNECT=0
export REPEATS=1000
export DIM=2
export IM2D=BrainMidSagittalSlice.png
export THRESH=100
(
for M in 0 1 2 3 ; 
do
$LAB $IM2D out.tif $THRESH $DIM $M $REPEATS $CONNECT
done
) > $LOG2D

echo "----------------------" >> $LOG2D
export CONNECT=1
(
for M in 0 1 2 3 ; 
do
$LAB $IM2D out.tif $THRESH $DIM $M $REPEATS $CONNECT
done
) >> $LOG2D
