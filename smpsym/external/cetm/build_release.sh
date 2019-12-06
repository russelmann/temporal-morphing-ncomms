rm -d -r Release
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make