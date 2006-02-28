#!/bin/sh


echo "Running 3D timing tests"

LAB=./testLabelling

export LOG3D=timing3d.log
export CONNECT=0
export REPEATS=10
export DIM=3
export IM3D=$1
export THRESH=100
(
for M in 0 1 2 3 ; 
do
$LAB $IM3D out.img $THRESH $DIM $M $REPEATS $CONNECT
done
) > $LOG3D

echo "----------------------" >> $LOG3D
export CONNECT=1
(
for M in 0 1 2 3 ; 
do
$LAB $IM3D out.img $THRESH $DIM $M $REPEATS $CONNECT
#gprof $LAB > time3d_b_$M".log"
done
) >> $LOG3D


# time3d_a?.log -- ESCells
# time3d_b?.log -- brain
