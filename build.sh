mkdir lib
cd lib
git clone https://github.com/google/googletest/ --branch release-1.10.0
cd ..
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=DEBUG
make
