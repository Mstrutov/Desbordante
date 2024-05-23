[![Downloads](https://static.pepy.tech/badge/desbordante)](https://pepy.tech/project/desbordante)
[![Downloads](https://static.pepy.tech/badge/desbordante/month)](https://pepy.tech/project/desbordante)

<p>
   <img src="https://github.com/Mstrutov/Desbordante/assets/88928096/d687809b-5a3b-420e-a192-a1a2b6697b2a"/>
</p>

# What is Desbordante?

Desbordante is a high-performance data profiler that can both *discover* and *validate* many different patterns in data. The patterns we are interested in mostly belong to the family of *functional dependencies*, but we offer many other interesting features, such as inclusion dependencies, unique column combinations, and association rules.

## What is an (exact) functional dependency?

In theory, a functional dependency (FD) describes a relationship between two sets of attributes in a relational database. The formal definition is: A functional dependency A --> B, where A and B are attributes or sets of attributes in the table, holds if and only if for all A = X, B = Y. This means that the value of B can be unambiguously determined by the value of A.

| A | B    | C |
|---|------|---|
| 0 | *1*  | 1 |
| 1 | 2    | 1 |
| 0 | *1*  | 3 |

A -> B holds --- the values for A = 0 are the same

| A | B    | C |
|---|------|---|
| 0 | *1*  | 1 |
| 1 | 2    | 1 |
| 0 | *3*  | 3 |

A -> B does not hold --- the values for A = 0 are different

This is the most basic pattern Desbordante can discover and validate, but there are many more.
## What can Desbordante do?

Desbordante does two key tasks: Discovery and Validation.  **Discovery** is designed to identify all instances of a specified pattern *type* of a given dataset. So, this means that if you would like to know **all** simple FDs that your dataset contains, you can run one of our algorithms and receive a list of them.

In contrast, the **Validation** task is designed to check whether a specified pattern *instance* is present in a given dataset. This means that if you think that some dependency is present in your dataset, you can check whether your hypothesis is true. This task not only returns True or False, but it also explains why the instance you specified does not hold (e.g. it can list table rows with conflicting values).
### Supported patterns

| Category                                | Type                                                                          | Discovery    | Validation |
| --------------------------------------- | ----------------------------------------------------------------------------- | ------------ | ---------- |
| **Functional Dependency Variants**      | Exact functional dependencies [1]                                             | ✅            | ✅          |
|                                         | Approximate functional dependencies with the g<sub>1</sub> metric             | ✅            | ✅          |
|                                         | Probabilistic functional dependencies, with `PerTuple` and `PerValue` metrics | ✅            | ❌          |
| **Graph Functional Dependencies**       |                                                                               | ❌            | ✅          |
| **Conditional Functional Dependencies** |                                                                               | ✅            | ❌          |
| **Inclusion Dependencies**              |                                                                               | ✅            | ❌          |
| **Order Dependencies**                  | Set-based axiomatization                                                      | ✅            | ❌          |
|                                         | List-based axiomatization                                                     | ✅            | ❌          |
| **Metric Functional Dependencies**      |                                                                               | ❌            | ✅          |
| **Fuzzy Algebraic Constraints**         |                                                                               | ✅            | ❌          |
| **Unique Column Combinations**          | Exact unique column combination                                               | ✅            | ✅          |
|                                         | Approximate unique column combination, with g<sub>1</sub> metric              | ✅            | ✅          |
| **Association Rules**                   |                                                                               | ✅            | ❌          |
| **Differential Dependencies**           |                                                                               | Coming soon! | ❌          |


### Potential uses
- **Scientific Data**:
    - Generate hypotheses for potential discoveries.
    - Draw immediate conclusions with sufficient data.
    - Guide further research directions.
- **Business Data**:
    - Clean data errors.
    - Remove approximate duplicates.
    - Match data schemas.
    - Support various practical applications.
- **Machine Learning**:
    - Assist in feature engineering.
    - Direct ablation studies.
- **Database Management**:
    - Identify primary and foreign keys.
    - Establish and check integrity constraints.
### Supported interfaces 
* **Console application:** Command-line access to essential profiling features with pattern specification, task and algorithm selection, and output to console or file.
* **Python bindings:** Access Desbordante within Python for an expanded feature set over the console. Create interactive programs and solve specific tasks using relational data algorithms that work with pandas DataFrames for easy data preparation.
* **Web application.** There is a web application that provides discovery and validation tasks with a rich interactive interface where results can be conveniently visualized. However, currently it supports a limited number of patterns and should be considered more as an interactive demo.

### Documentation and guides
- [A brief introduction to the tool and its use cases (English)](https://medium.com/@chernishev/exploratory-data-analysis-with-desbordante-4b97299cce07)
- [A brief introduction to the tool and its use cases (Russian)](https://habr.com/ru/company/unidata/blog/667636/)
- [List of various articles and guides](https://desbordante.unidata-platform.ru/papers)
- [Extensive list of tutorial examples for each supported pattern](https://github.com/Desbordante/desbordante-core/tree/main/examples)


## Console

### Examples

All tables in the examples are stored in .csv files and have a header row.

1) Discover all exact functional dependencies in a table with the default algorithm (HyFD). 

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

2) Discover all approximate functional dependencies with error less than or equal to 0.1 in a table with the default algorithm (Pyro).
```sh
python3 cli.py --task=afd --table=../examples/datasets/inventory_afd.csv , True --error=0.1
```

```text
[Id] -> ProductName
[Id] -> Price
[ProductName] -> Price
```

3) Check whether a metric functional dependency “Title -> Duration” with radius 5 (using the Euclidean metric) holds with the default MFD validation algorithm (BRUTE).

```sh
python3 cli.py --task=mfd_verification --table=../examples/datasets/theatres_mfd.csv , True --lhs_indices=0 --rhs_indices=2 --metric=euclidean --parameter=5
```

```text
True
```

As always, you can find more information in the documentation and the help files. We are working hard to update them, but in the meanwhile, do not hesitate to ask questions!
## Python bindings

Python is the go-to language of modern data science, which is why we offer Python bindings to the Desbordante C++ core library. These bindings allow you to do more than just basic tasks; they can also provide extra details, like explaining why a particular pattern does not hold.

The library's main intended usage is to solve various data quality problems by constructing custom Python pipelines. To show the power of the library, we have implemented several demo scenarios:

1) [Typo detection](https://colab.research.google.com/drive/1h5mQAIIxSb6Sgc_Ep8AYZlgt4BGXN6A9)
2) [Data deduplication](https://colab.research.google.com/drive/1hgF8idXi1-U4ZOR0fAmdbfbhltgEJecR?usp=sharing)
3) [Anomaly detection](https://colab.research.google.com/drive/1hgF8idXi1-U4ZOR0fAmdbfbhltgEJecR?usp=sharing)

[There is](https://desbordante.streamlit.app/) also an interactive Streamlit demo for all of them, and all of the scripts can be found [here](https://github.com/Desbordante/desbordante-core/tree/main/examples). The ideas behind them are briefly discussed in this [preprint](https://arxiv.org/abs/2307.14935) (Section 3). 

### Examples

All tables in the examples are stored in .csv files and have a header row.

1) Discover all exact functional dependencies in a table with the default algorithm (HyFD). 

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

2) Discover all approximate functional dependencies with error less than or equal to 0.1 in a table with the default algorithm (Pyro).

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

3) Check whether a metric functional dependency `Title -> Duration` with radius 5 (using the Euclidean metric) holds with the default MFD validation algorithm (BRUTE).

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

4) Discover approximate functional dependencies with various error thresholds. Here, we are using a pandas DataFrame to load data from a .csv file.

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

 Desbordante also offers a web interface specifically designed for interactive tasks which involve multiple steps and require substantial user input on each of them. Interactive tasks originate from full-fledged Python pipelines: we plan to identify the most interesting and useful ones and implement them in the web version. Currently, the only scenario implemented in this manner is typo detection. The web interface is also useful for pattern discovery and validation tasks: you can specify parameters, browse results, employ advanced visualizations and filters, all in a convenient way.

You can try the deployed web version [here](https://desbordante.unidata-platform.ru/). You have to register in order to process your own datasets. Keep in mind that due to high demand various time and memory limits are enforced: processing is aborted if they are exceeded. The source code of the web interface is kept in a separate [repo](https://github.com/Desbordante/desbordante-server-node).

# I still don't understand how to use Desbordante and patterns :(

No worries! Desbordante offers a novel type of data profiling, which may require that you first familiarize yourself with its concepts and usage. The most challenging part of Desbordante are the primitives: their definitions and applications in practice. To help you get started, here’s a step-by-step guide:

1) First of all, explore the guides on our [website](https://desbordante.unidata-platform.ru/papers). Since our team currently does not include technical writers, it's possible that some guides may be missing.
2) To compensate for the lack of guides, we provide several examples for each supported pattern. These examples illustrate both the pattern itself and how to use it in Python. You can check them out [here](https://github.com/Desbordante/desbordante-core/tree/main/examples).
3) Each of our patterns was introduced in a research paper. These papers typically provide a formal definition of the pattern, examples of use, and its application scope. We recommend at least skimming through them. Don't be discouraged by the complexity of the papers! To effectively use the patterns, you only need to read the more accessible parts, such as the introduction and the example sections.
4) Finally, do not hesitate to ask questions in the mailing list (link below) or create an issue.

### Papers about patterns
Here is a list of papers about patterns, organized in the recommended reading order in each item:
* Functional dependency variants:
    - Exact functional dependencies
       1. [Thorsten Papenbrock et al. 2015. Functional dependency discovery: an experimental evaluation of seven algorithms. Proc. VLDB Endow. 8, 10 (June 2015), 1082–1093.](http://www.vldb.org/pvldb/vol8/p1082-papenbrock.pdf)
       2. [Thorsten Papenbrock and Felix Naumann. 2016. A Hybrid Approach to Functional Dependency Discovery. In Proceedings of the 2016 International Conference on Management of Data (SIGMOD '16). Association for Computing Machinery, New York, NY, USA, 821–833.](https://hpi.de/fileadmin/user_upload/fachgebiete/naumann/publications/PDFs/2016_papenbrock_a.pdf)
    - Approximate functional dependencies with g<sub>1</sub> metric
       3. [Sebastian Kruse and Felix Naumann. 2018. Efficient discovery of approximate dependencies. Proc. VLDB Endow. 11, 7 (March 2018), 759–772.](https://www.vldb.org/pvldb/vol11/p759-kruse.pdf)
    - Probabilistic functional dependencies, with PerTuple and PerValue metrics
       4. [Daisy Zhe Wang et al. Functional Dependency Generation and Applications in Pay-As-You-Go Data Integration Systems. WebDB 2009](http://webdb09.cse.buffalo.edu/papers/Paper18/webdb09.pdf)
       5. [Daisy Zhe Wang et al. Discovering Functional Dependencies in Pay-As-You-Go Data Integration Systems. Tech Rep. UCB/EECS-2009-119.](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2009/EECS-2009-119.pdf)
* Graph functional dependencies
     6. [Wenfei Fan, Yinghui Wu, and Jingbo Xu. 2016. Functional Dependencies for Graphs. In Proceedings of the 2016 International Conference on Management of Data (SIGMOD '16). Association for Computing Machinery, New York, NY, USA, 1843–1857.](https://dl.acm.org/doi/pdf/10.1145/2882903.2915232)
* Conditional functional dependencies
    7. [Rammelaere, J., Geerts, F. (2019). Revisiting Conditional Functional Dependency Discovery: Splitting the “C” from the “FD”. Machine Learning and Knowledge Discovery in Databases. ECML PKDD 2018. ](https://link.springer.com/chapter/10.1007/978-3-030-10928-8_33)
* Inclusion dependencies (discovery)
    8. [Falco Dürsch et al. 2019. Inclusion Dependency Discovery: An Experimental Evaluation of Thirteen Algorithms. In Proceedings of the 28th ACM International Conference on Information and Knowledge Management (CIKM '19). Association for Computing Machinery, New York, NY, USA, 219–228.](https://hpi.de/fileadmin/user_upload/fachgebiete/naumann/publications/PDFs/2019_duersch_inclusion.pdf)
    9. [Sebastian Kruse, et al: Fast Approximate Discovery of Inclusion Dependencies. BTW 2017: 207-226](http://btw2017.informatik.uni-stuttgart.de/slidesandpapers/F4-10-47/paper_web.pdf)
* Order dependencies
   10. [Jaroslaw Szlichta et al. 2017. Effective and complete discovery of order dependencies via set-based axiomatization. Proc. VLDB Endow. 10, 7 (March 2017), 721–732.](http://www.vldb.org/pvldb/vol10/p721-szlichta.pdf)
   11.  [Langer, P., Naumann, F. Efficient order dependency detection. The VLDB Journal 25, 223–241 (2016)](https://link.springer.com/article/10.1007/s00778-015-0412-3)
* Metric functional dependencies
   12. [N. Koudas et al. "Metric Functional Dependencies," 2009 IEEE 25th International Conference on Data Engineering, Shanghai, China, 2009, pp. 1275-1278.](https://ieeexplore.ieee.org/document/4812519)
* Fuzzy algebraic constraints
   13. [Paul G. Brown and Peter J. Hass. 2003. BHUNT: automatic discovery of Fuzzy algebraic constraints in relational data. In Proceedings of the 29th international conference on Very large data bases - Volume 29 (VLDB '03), Vol. 29. VLDB Endowment, 668–679.](https://www.vldb.org/conf/2003/papers/S20P03.pdf)
* Unique column combinations
   14. [Sebastian Kruse and Felix Naumann. 2018. Efficient discovery of approximate dependencies. Proc. VLDB Endow. 11, 7 (March 2018), 759–772.](https://www.vldb.org/pvldb/vol11/p759-kruse.pdf)
* Association rules
   15. [Charu C. Aggarwal, Jiawei Han. 2014. Frequent Pattern Mining. Springer Cham. pp 471.](https://link.springer.com/book/10.1007/978-3-319-07821-2)

## Installation

Desbordante is [available](https://pypi.org/project/desbordante/) at the Python Package Index (PyPI) for Python >=3.7.

Install Desbordante easily:

```sh
$ pip install desbordante
```

However, as Desbordante core is implemented in C++, there are additional requirements for your machine. Currently, we only support `manylinux2014` (Ubuntu 20.04+, or any other Linux distribution with gcc 10+). If this does not work for you, consider building from source.

## CLI installation

**NOTE**: Only Python 3.11+

Сlone the repository, change the current directory to the project directory and run the following commands:

```sh
pip install -r cli/requirements.txt
python3 cli/cli.py --help
```

## Build instructions

The following instructions were tested on Ubuntu 20.04+ LTS.
### Dependencies
Prior to cloning the repository and attempting to build the project,  make sure that you have all of the dependencies installed:

- GNU g++ 10+
- CMake 3.13+
- Boost 1.74.0+

To use the test datasets, you will need Git Large File Storage 3.0.2+.

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

The Python module can be built with the `--pybind` switch:
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
./Desbordante_test --gtest_filter='*:-*HeavyDatasets*'
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
If when cloning the repo with git lfs installed, `git clone` produces the following (or similar) error:
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
If type hints don't work for you in your editor of choice, then install stubs:
```sh
pip install desbordante-stubs
```
**NOTE**: Stubs may not fully support current version of the `desbordante` package, as they are updated independently.

## Cite

If you use this software for research, please cite one of our papers:

1) George Chernishev, et al. Solving Data Quality Problems with Desbordante: a Demo. CoRR abs/2307.14935 (2023).
2) George Chernishev, et al. "Desbordante: from benchmarking suite to high-performance science-intensive data profiler (preprint)". CoRR abs/2301.05965. (2023).
3) M. Strutovskiy, N. Bobrov, K. Smirnov and G. Chernishev, "Desbordante: a Framework for Exploring Limits of Dependency Discovery Algorithms," 2021 29th Conference of Open Innovations Association (FRUCT), 2021, pp. 344-354, doi: 10.23919/FRUCT52173.2021.9435469.
4) A. Smirnov, A. Chizhov, I. Shchuckin, N. Bobrov and G. Chernishev, "Fast Discovery of Inclusion Dependencies with Desbordante," 2023 33rd Conference of Open Innovations Association (FRUCT), Zilina, Slovakia, 2023, pp. 264-275, doi: 10.23919/FRUCT58615.2023.10143047.
5) Y. Kuzin, D. Shcheka, M. Polyntsov, K. Stupakov, M. Firsov and G. Chernishev, "Order in Desbordante: Techniques for Efficient Implementation of Order Dependency Discovery Algorithms," 2024 35th Conference of Open Innovations Association (FRUCT), Tampere, Finland, 2024, pp. 413-424.
6) I. Barutkin, M. Fofanov, S. Belokonny, V. Makeev and G. Chernishev, "Extending Desbordante with Probabilistic Functional Dependency Discovery Support," 2024 35th Conference of Open Innovations Association (FRUCT), Tampere, Finland, 2024, pp. 158-169.

# Contacts and Q&A

If you have any questions, don't hesitate to ask them in our [google group](https://groups.google.com/g/desbordante)! If you would like to contact the dev team, you can email:
- George Chernishev
- Maxim Strutovsky
- Nikita Bobrov
  using the addresses in the papers.
 
