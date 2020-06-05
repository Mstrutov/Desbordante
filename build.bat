cd lib
git clone https://github.com/google/googletest/
cd .. 
@call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\Common7\Tools\VsDevCmd.bat"
mkdir build
cd build
cmake ..
devenv /build Release metanome_cpp.sln 

