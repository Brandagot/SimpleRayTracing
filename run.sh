cd bin-debug-amd
make -j4
cd ..
/usr/bin/time --format='%e' ./bin-debug-amd/main-pthreads --size 1024 1024 --jpeg x-ray-10-threads-1024x1024-no_dragon-no_logo.jpg --threads 10 2> temp-pthread-10