cd bin-debug-amd
make -j4
cd ..
/usr/bin/time --format='%e' ./bin-debug-amd/main-pthreads --size 1512 1512 --jpeg x-ray-10-threads-1512x1512-no_dragon.jpg --threads 10 2> temp-pthread-10