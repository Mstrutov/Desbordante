Здесь приведён проект, который я реализую в рамках курсовой работы. Подробнее о нём написано в отчёте по курсовой работе за пятый семестр в файле CourseWorkMaxim.pdf.

https://github.com/cupertank реализовал части проекта, представленные в src/model, src/parser, src/tests и src/util

Направления для работы над кодом:
  - [ ] проведение повторного объектно-ориентированного анализа;
  - [ ] рефакторинг кода и оптимизация;
  - [x] удаление из репозитория лишних файлов;
  - [ ] написание тестов;
  - [ ] тестирование и распараллеливание алгоритма Pyro.

# Installation guide.
* ## Ubuntu 
  The following instructions were tested on Ubuntu 18.04.4 LTS. 
  ### Dependencies 
  Prior to cloning the repository and attempting to build the project, ensure that you have the following software: 
   
  - GNU g++ compiler, version 8+ 
  - CMake, version 3.10+ 
  - Boost library, version 1.65.1+ 
  
  ### Building the project 
  Firstly, navigate to a desired directory.
  Then, clone the repository, cd into the project directory and launch the build script:
  ```
  git clone https://github.com/Mstrutov/FDExperimentsTool/    
  cd FDExperimentsTool
  sh build.sh
  ```
  ### Launching the binaries
  The script generates the following file structure in `/path/to/FDExperimentsTool/build/target`:
  ```bash
  ├───inputData
  │   └───some-sample-csv's.csv
  ├───fdtester_test
  ├───fdtester_run
  ```
  The `inputData` directory contains several .csv files that may be used by `fdtester_test`. Run `fdtester_test` to perform unit testing:
  ```
  cd build/target
  ./fdtester_test
  ```
  The tool itself is launched via the following line:
  ```
  ./fdtester_run -algo=tane -data=<dataset_name>.csv
  ```
  The `<dataset_name>.csv`, which is a user-provided dataset, should be placed in the `/path/to/FDExperimentsTool/build/target/inputData` directory.
* ## Windows 
  The following instructions were tested on Windows 10 . 
  ### Dependencies 
  Prior to cloning the repository and attempting to build the project, ensure that you have the following software: 
   
  - Microsoft Visual Studio 2019 
  - CMake, version 3.10+ 
  - Boost library, version 1.65.1+ \
    The recommended way to install Boost is by using [chocolatey](https://chocolatey.org/)  
  
  ### Building the project 
  Firstly, launch the command prompt and navigate to a desired directory.
  Then, clone the repository, cd into the project directory and launch the build script:
  ```
  git clone https://github.com/Mstrutov/FDExperimentsTool/    
  cd FDExperimentsTool
  build.bat
  ```
  *Note:* to compile the project, the script uses hard-coded path to MSVC developer command prompt, which is located 
  by default at `C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\Common7\Tools\VsDevCmd.bat`. 
  You should change the path in the script if it differs from the default one. 
  ### Launching the binaries
  The script generates the following file structure in `\path\to\FDExperimentsTool\build\target`: 
  ```bash
  ├───inputData
  │   └───some-sample-csv's.csv
  ├───fdtester_test.exe
  ├───fdtester_run.exe
  ```
  The `inputData` directory contains several .csv files that may be used by `fdtester_test`. Run `fdtester_test` to perform unit testing:
  ```
  cd build\target
  fdtester_test.exe
  ```
  The tool itself is launched via the following line:
  ```
  fdtester_run.exe -algo=tane -data=<dataset_name>.csv
  ```
  The `<dataset_name>.csv`, which is a user-provided dataset, should be placed in the `\path\to\FDExperimentsTool\build\target\inputData` directory.
    
[Email me at strutovksy.m.a@gmail.com](mailto:strutovksy.m.a@gmail.com)
