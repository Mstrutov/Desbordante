mkdir lib
cd lib
git clone https://github.com/google/googletest/
cd ..
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
