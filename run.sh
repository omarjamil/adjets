#!/bin/bash

TIMELIMIT=3

echo "Remove log and result files? (y/n)"
read -t $TIMELIMIT answer

if [ -z "$answer" ]
then
    answer="y"
fi

if [ $answer == "y" ];then
    cd results
    if [ -e "nuFile.dat" ];then
	echo "Removing all the files in results/ "
	rm *.dat
    fi
    cd ..
    rm sim.log
    touch sim.log
    echo "Running aleph"
    ./bin/aleph > sim.log 2>&1
    echo "Moving results files to results/ "
    mv *.dat results/
    cp *.par results/
else exit
fi