cd bin-debug-amd
make -j4
cd ..
/usr/bin/time --format='%e' ./bin-debug-amd/main-pthreads --size 256 256 --jpeg x-ray-6-threads-256x256-no_dragon-no_logo.jpg --threads 6 2> temp-pthread-6