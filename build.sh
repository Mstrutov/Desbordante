mkdir lib
cd lib
git clone https://github.com/google/googletest/ --branch release-1.12.1
git clone https://github.com/amrayn/easyloggingpp/ --branch v9.97.0
git clone https://github.com/aantron/better-enums.git --branch 0.11.3
git clone https://github.com/pybind/pybind11.git --branch v2.10
git clone https://github.com/libigl/eigen
cd ..
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
