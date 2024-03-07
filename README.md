
<p>
   <img src="https://github.com/Mstrutov/Desbordante/assets/88928096/d687809b-5a3b-420e-a192-a1a2b6697b2a"/>
</p>


# General

Desbordante is a high-performance data profiler that is capable of discovering and validating many different patterns in data using various algorithms. The currently supported data patterns are:
* Functional dependencies, both exact and approximate (discovery and validation)
* Conditional functional dependencies (discovery)
* Metric functional dependencies (validation)
* Fuzzy algebraic constraints (discovery)
* Unique column combinations, both exact and approximate (discovery and validation)
* Association rules (discovery)

The discovered patterns can have many uses:
* For scientific data, especially those obtained experimentally, an interesting pattern allows to formulate a hypothesis that could lead to a scientific discovery. In some cases it even allows to draw conclusions immediately, if there is enough data. At the very least, the found pattern can provide a direction for further study. 
* For business data it is also possible to obtain a hypothesis based on found patterns. However, there are more down-to-earth and more in-demand applications in this case: clearing errors in data, finding and removing inexact duplicates, performing schema matching, and many more. 
* For training data used in machine learning applications the found patterns can help in feature engineering and in choosing the direction for the ablation study.
* For database data, found patterns can help with defining (recovering) primary and foreign keys, setting up (checking) all kinds of integrity constraints.

Desbordante can be used via three interfaces:
* **Console application.** This is a classic command-line interface that aims to provide basic profiling functionality, i.e. discovery and validation of patterns. A user can specify pattern type, task type, algorithm, input file(s) and output results to the screen or into a file.
* **Python bindings.** Desbordante functionality can be accessed from within Python programs by employing the Desbordante Python library. This interface offers everything that is currently provided by the console version and allows advanced use, such as building interactive applications and designing scenarios for solving a particular real-life task. Relational data processing algorithms accept pandas DataFrames as input, allowing the user to conveniently preprocess the data before mining patterns.
* **Web application.** There is a web application that provides discovery and validation tasks with a rich interactive interface where results can be conveniently visualized. However, currently it supports a limited number of patterns and should be considered more as an interactive demo.

A brief introduction into the tool and its use cases is presented [here](https://medium.com/@chernishev/exploratory-data-analysis-with-desbordante-4b97299cce07) (in English) and [here](https://habr.com/ru/company/unidata/blog/667636/) (in Russian). Also, a list of various articles and guides can be found [here](https://desbordante.unidata-platform.ru/papers).

## Console

Usage examples:
1) Discover all exact functional dependencies in a table represented by a .csv file that uses a comma as the separator and has a header row. In this example the default FD discovery algorithm (HyFD) is used.

```sh
python3 cli.py --task=fd --table=../examples/datasets/university_fd.csv , True
```

```text
[Course Classroom] -> Professor
[Classroom Semester] -> Professor
[Classroom Semester] -> Course
[Professor] -> Course
[Professor Semester] -> Classroom
[Course Semester] -> Classroom
[Course Semester] -> Professor
```

2) Discover all approximate functional dependencies with error less than or equal to 0.1 in a table represented by a .csv file that uses a comma as the separator and has a header row. In this example the default AFD discovery algorithm (Pyro) is used.

```sh
python3 cli.py --task=afd --table=../examples/datasets/inventory_afd.csv , True --error=0.1
```

```text
[Id] -> ProductName
[Id] -> Price
[ProductName] -> Price
```

3) Check whether metric functional dependency “Title -> Duration” with radius 5 (using the Euclidean metric) holds in a table represented by a .csv file that uses a comma as the separator and has a header row. In this example the default MFD validation algorithm (BRUTE) is used.

```sh
python3 cli.py --task=mfd_verification --table=../examples/datasets/theatres_mfd.csv , True --lhs_indices=0 --rhs_indices=2 --metric=euclidean --parameter=5
```

```text
True
```

For more information consult documentation and help files.

## Python bindings

Desbordante features can be accessed from within Python programs by employing the Desbordante Python library. The library is implemented in the form of Python bindings to the interface of the Desbordante C++ core library, using pybind11. Apart from discovery and validation of patterns, this interface is capable of providing valuable additional information which can, for example, describe why a given pattern does not hold. All this allows end users to solve various data quality problems by constructing ad-hoc Python programs. To show the power of this interface, we have implemented several demo scenarios:
1) [Typo detection](https://colab.research.google.com/drive/1h5mQAIIxSb6Sgc_Ep8AYZlgt4BGXN6A9)
2) [Data deduplication](https://colab.research.google.com/drive/1hgF8idXi1-U4ZOR0fAmdbfbhltgEJecR?usp=sharing)
3) [Anomaly detection](https://colab.research.google.com/drive/1hgF8idXi1-U4ZOR0fAmdbfbhltgEJecR?usp=sharing)

[There is](https://desbordante.streamlit.app/) also an interactive demo for all of them, and all of these python scripts are [here](https://github.com/Mstrutov/Desbordante/tree/main/examples). The ideas behind them are briefly discussed in this [preprint](https://arxiv.org/abs/2307.14935) (Section 3). 

Simple usage examples:
1) Discover all exact functional dependencies in a table represented by a .csv file that uses a comma as the separator and has a header row. In this example the FD discovery algorithm HyFD is used.

```python
import desbordante

TABLE = 'examples/datasets/university_fd.csv'

algo = desbordante.fd.algorithms.Default()
algo.load_data(table=(TABLE, ',', True))
algo.execute()
result = algo.get_fds()
print('FDs:')
for fd in result:
    print(fd)
```
```text
FDs:
[Course Classroom] -> Professor
[Classroom Semester] -> Professor
[Classroom Semester] -> Course
[Professor] -> Course
[Professor Semester] -> Classroom
[Course Semester] -> Classroom
[Course Semester] -> Professor
```

2) Discover all approximate functional dependencies with error less than or equal to 0.1 in a table represented by a .csv file that uses a comma as the separator and has a header row. In this example the AFD discovery algorithm Pyro is used.

```python
import desbordante

TABLE = 'examples/datasets/inventory_afd.csv'
ERROR = 0.1

algo = desbordante.afd.algorithms.Default()
algo.load_data(table=(TABLE, ',', True))
algo.execute(error=ERROR)
result = algo.get_fds()
print('AFDs:')
for fd in result:
    print(fd)
```
```text
AFDs:
[Id] -> Price
[Id] -> ProductName
[ProductName] -> Price
```

3) Check whether metric functional dependency “Title -> Duration” with radius 5 (using the Euclidean metric) holds in a table represented by a .csv file that uses a comma as the separator and has a header row. In this example the default MFD validation algorithm (BRUTE) is used.

```python
import desbordante

TABLE = 'examples/datasets/theatres_mfd.csv'
METRIC = 'euclidean'
LHS_INDICES = [0]
RHS_INDICES = [2]
PARAMETER = 5

algo = desbordante.mfd_verification.algorithms.Default()
algo.load_data(table=(TABLE, ',', True))
algo.execute(lhs_indices=LHS_INDICES, metric=METRIC,
             parameter=PARAMETER, rhs_indices=RHS_INDICES)
if algo.mfd_holds():
    print('MFD holds')
else:
    print('MFD does not hold')
```
```text
MFD holds
```
4) Discover approximate functional dependencies with various error thresholds. Here, we are using a pandas DataFrame to load data from a CSV file.
```python-repl
>>> import desbordante
>>> import pandas as pd
>>> pyro = desbordante.afd.algorithms.Pyro()  # same as desbordante.afd.algorithms.Default()
>>> df = pd.read_csv('examples/datasets/iris.csv', sep=',', header=None)
>>> pyro.load_data(table=df)
>>> pyro.execute(error=0.0)
>>> print(f'[{", ".join(map(str, pyro.get_fds()))}]')
[[0 1 2] -> 4, [0 2 3] -> 4, [0 1 3] -> 4, [1 2 3] -> 4]
>>> pyro.execute(error=0.1)
>>> print(f'[{", ".join(map(str, pyro.get_fds()))}]')
[[2] -> 0, [2] -> 3, [2] -> 1, [0] -> 2, [3] -> 0, [0] -> 3, [0] -> 1, [1] -> 3, [1] -> 0, [3] -> 2, [3] -> 1, [1] -> 2, [2] -> 4, [3] -> 4, [0] -> 4, [1] -> 4]
>>> pyro.execute(error=0.2)
>>> print(f'[{", ".join(map(str, pyro.get_fds()))}]')
[[2] -> 0, [0] -> 2, [3] -> 2, [1] -> 2, [2] -> 4, [3] -> 4, [0] -> 4, [1] -> 4, [3] -> 0, [1] -> 0, [2] -> 3, [2] -> 1, [0] -> 3, [0] -> 1, [1] -> 3, [3] -> 1]
>>> pyro.execute(error=0.3)
>>> print(f'[{", ".join(map(str, pyro.get_fds()))}]')
[[2] -> 1, [0] -> 2, [2] -> 0, [2] -> 3, [0] -> 1, [3] -> 2, [3] -> 1, [1] -> 2, [3] -> 0, [0] -> 3, [4] -> 1, [1] -> 0, [1] -> 3, [4] -> 2, [4] -> 3, [2] -> 4, [3] -> 4, [0] -> 4, [1] -> 4]
```

## Web interface

While the Python interface makes building interactive applications possible, Desbordante also offers a web interface which is aimed specifically for interactive tasks. Such tasks typically involve multiple steps and require substantial user input on each of them. Interactive tasks usually originate from Python scenarios, i.e. we select the most interesting ones and implement them in the web version. Currently, only the typo detection scenario is implemented. The web interface is also useful for pattern discovery and validation tasks: a user may specify parameters, browse results, employ advanced visualizations and filters, all in a convenient way.

You can try the deployed web version [here](https://desbordante.unidata-platform.ru/). You have to register in order to process your own datasets. Keep in mind that due to a large demand various time and memory limits are enforced: processing is aborted if they are exceeded. The source code of the web interface is kept in a separate [repo](https://github.com/vs9h/Desbordante).

## Installation (this is what you probably want if you are not a project maintainer)
Desbordante is [available](https://pypi.org/project/desbordante/) at the Python Package Index (PyPI). Dependencies:

* Python >=3.7 

To install Desbordante type:

```sh
$ pip install desbordante
```

However, as Desbordante core uses C++, additional requirements on the machine are imposed. Therefore this installation option may not work for everyone. Currently, only manylinux2014 (Ubuntu 20.04+, or any other linux distribution with gcc 10+) is supported. If the above does not work for you consider building from sources.

## CLI installation

**NOTE**: Only Python 3.11+ is supported for CLI

Сlone the repository, change the current directory to the project directory and run the following commands:

```sh
pip install -r cli/requirements.txt
python3 cli/cli.py --help
```

## Build instructions

### Ubuntu
The following instructions were tested on Ubuntu 20.04+ LTS.
### Dependencies
Prior to cloning the repository and attempting to build the project, ensure that you have the following software:

- GNU g++ compiler, version 10+
- CMake, version 3.13+
- Boost library, version 1.74.0+

To use test datasets you will need:
- Git Large File Storage, version 3.0.2+

### Building the project
#### Building the Python module using pip

Clone the repository, change the current directory to the project directory and run the following commands:

```bash
./build.sh
python3 -m venv venv
source venv/bin/activate
python3 -m pip install .
```

Now it is possible to `import desbordante` as a module from within the created virtual environment. 

#### Building tests & the Python module manually
In order to build tests, pull the test datasets using the following command:
```sh
./pull_datasets.sh
```
then build the tests themselves:
```sh
./build.sh -j$(nproc)
```

The Python module can be built by providing the `--pybind` switch:
```sh
./build.sh --pybind -j$(nproc)
```

See `./build.sh --help` for more available options.

The `./build.sh` script generates the following file structure in `/path/to/Desbordante/build/target`:
```
├───input_data
│   └───some-sample-csv\'s.csv
├───Desbordante_test
├───desbordante.cpython-*.so
```

The `input_data` directory contains several .csv files that are used by `Desbordante_test`. Run `Desbordante_test` to perform unit testing:
```sh
cd build/target
./Desbordante_test
```

`desbordante.cpython-*.so` is a Python module, packaging Python bindings for the Desbordante core library. In order to use it, simply `import` it:
```sh
cd build/target
python3
>>> import desbordante
```

We use [easyloggingpp](https://github.com/abumq/easyloggingpp) in order to log (mostly debug) information in the core library. Python bindings search for a configuration file in the working directory, so to configure logging, create `logging.conf` in the directory from which desbordante will be imported. In particular, when running the CLI with `python3 ./relative/path/to/cli.py`, `logging.conf` should be located in `.`.

## Troubleshooting

### Git LFS
If, when cloning the repo with git lfs installed, `git clone` produces the following (or similar) error:
```
Cloning into 'Desbordante'...
remote: Enumerating objects: 13440, done.
remote: Counting objects: 100% (13439/13439), done.
remote: Compressing objects: 100% (3784/3784), done.
remote: Total 13440 (delta 9537), reused 13265 (delta 9472), pack-reused 1
Receiving objects: 100% (13440/13440), 125.78 MiB | 8.12 MiB/s, done.
Resolving deltas: 100% (9537/9537), done.
Updating files: 100% (478/478), done.
Downloading datasets/datasets.zip (102 MB)
Error downloading object: datasets/datasets.zip (2085458): Smudge error: Error downloading datasets/datasets.zip (2085458e26e55ea68d79bcd2b8e5808de731de6dfcda4407b06b30bce484f97b): batch response: This repository is over its data quota. Account responsible for LFS bandwidth should purchase more data packs to restore access.
```
delete the already cloned version, set `GIT_LFS_SKIP_SMUDGE=1` environment variable and clone the repo again:
```sh
GIT_LFS_SKIP_SMUDGE=1 git clone git@github.com:Mstrutov/Desbordante.git
```

### No type hints in IDE
If type hints don't work for you in Visual Studio Code, for example, then install stubs using the command:
```sh
pip install desbordate-stubs
```
**NOTE**: Stubs may not fully support current version of `desbordante` package, as they are updated independently.

## Cite
If you use this software for research, please cite one of our papers:
1) George Chernishev, et al. Solving Data Quality Problems with Desbordante: a Demo. CoRR abs/2307.14935 (2023).
2) George Chernishev, et al. "Desbordante: from benchmarking suite to high-performance science-intensive data profiler (preprint)". CoRR abs/2301.05965. (2023).
3) M. Strutovskiy, N. Bobrov, K. Smirnov and G. Chernishev, "Desbordante: a Framework for Exploring Limits of Dependency Discovery Algorithms," 2021 29th Conference of Open Innovations Association (FRUCT), 2021, pp. 344-354, doi: 10.23919/FRUCT52173.2021.9435469.
4) A. Smirnov, A. Chizhov, I. Shchuckin, N. Bobrov and G. Chernishev, "Fast Discovery of Inclusion Dependencies with Desbordante," 2023 33rd Conference of Open Innovations Association (FRUCT), Zilina, Slovakia, 2023, pp. 264-275, doi: 10.23919/FRUCT58615.2023.10143047.

# Contacts and Q&A

If you have any questions regarding the tool usage you can ask it in our [google group](https://groups.google.com/g/desbordante). To contact dev team email George Chernishev, Maxim Strutovsky or Nikita Bobrov.
