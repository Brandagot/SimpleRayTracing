#!/usr/bin/env bash
#


COMPILER=`gcc --version |head -1`

TEMP=`lscpu|grep "Model name:"`
IFS=':' read -ra CPU_MODEL <<< "$TEMP"

width=1024
height=1024

if [ ! -f timing.csv ];
then
    echo "CPU,Parallelisation,Number of threads,Number of nodes,Compiler,Image size,Runtime in sec" > timing-serial.csv
fi

/usr/bin/time --format='%e' ./bin-debug-amd/main --size $width $height --filename dragon-${width}x$height-serial.txt 2> temp-serial

RUNTIME=`cat temp-serial`

echo ${CPU_MODEL[1]},None,0,1,$COMPILER,${width}x$height,$RUNTIME >> timing-serial.csv
rm temp-serial
