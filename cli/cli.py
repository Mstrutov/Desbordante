#!/usr/bin/python3

import sys
from collections import namedtuple
from enum import StrEnum, auto
from time import process_time
from typing import Any, Callable
from os import scandir

import click
import desbordante


class Task(StrEnum):
    fd = auto()
    cfd = auto()
    ar = auto()
    afd = auto()
    od = auto()
    pfd = auto()
    ind = auto()
    fd_verification = auto()
    afd_verification = auto()
    mfd_verification = auto()
    ucc_verification = auto()
    aucc_verification = auto()
    gfd_verification = auto()


class Algorithm(StrEnum):
    pyro = auto()
    tane = auto()
    pfdtane = auto()
    hyfd = auto()
    fd_mine = auto()
    dfd = auto()
    dep_miner = auto()
    fdep = auto()
    fun = auto()
    fastfds = auto()
    aid = auto()
    fd_first = auto()
    fastod = auto()
    order = auto()
    spider = auto()
    faida = auto()
    apriori = auto()
    naive_fd_verifier = auto()
    naive_afd_verifier = auto()
    icde09_mfd_verifier = auto()
    naive_ucc_verifier = auto()
    naive_aucc_verifier = auto()
    naive_gfd_verifier = auto()
    gfd_verifier = auto()
    egfd_verifier = auto()


HELP = 'help'
TASK = 'task'
ALGO = 'algo'
ALGORITHM = 'algorithm'
FILENAME = 'filename'
VERBOSE = 'verbose'
ERROR = 'error'
ERROR_MEASURE = 'error_measure'
TABLE = 'table'
TABLES = 'tables'
TABLES_LIST = 'tables_list'
TABLES_DIRECTORY = 'tables_directory'
INPUT_FORMAT = 'input_format'

PRIMARY_HELP = '''The Desbordante data profiler is designed to help users
discover or verify various types of patterns in data. These patterns are
various kinds of data dependencies (functional, metric, inclusion, etc),
constraints (algebraic, denial, etc), many types of association rules and more.

Each pattern type is termed a “primitive”, and a specific occurrence of this
primitive inside a specific dataset is referred to as its “instance”. For
example, a functional dependency is a primitive, and the “company -> email”
functional dependency inside the sales.csv table is its instance.

For each primitive, Desbordante supports two key tasks: discovery and
verification. The discovery task discovers (mines) all relevant* primitive
instances within a given dataset, while verification checks (validates)
whether a given instance holds. The verification task name contains the
“_verification” suffix; for example, the “fd” task returns all minimal and
non-trivial functional dependencies contained in the dataset,
and “fd_verification” checks whether a given functional dependency holds.

* The notion of “relevance” depends on the definition of the primitive. The
relevant set is usually the minimal set of instances from which the rest can
be derived. For actual definitions, consult the related papers.

Each task can be performed by one of the provided algorithms. For some
primitives, several algorithms are available, each of which excel on a
different type of dataset, such as “wide” or “tall” tables. See
--task=TASK --help to view the list of all supported algorithms for TASK.
Additionally, we also provide a default algorithm which is the best one in
many cases.

The algorithms accept two types of parameters: primitive-specific
constraints and implementation-related settings. The first group includes
parameters that define which instances are to discover/verify, e.g. the
constraint on the maximum length of the left-hand side of the functional
dependency. The second group refers to general parameters that impact
algorithm performance, e.g. the number of threads to use, seed, buffer size
for intermediates, etc. Note that not all algorithms, even those designed
for the same task, support multithreading.

Next, several primitives have approximate versions, which allow some records
to deviate from the primitive definition. This is practical for dealing with
real-life data, which may contain all kinds of imperfections. The number of
such imperfect records is usually defined by an error threshold which is
calculated according to a primitive-specific procedure. In this case,
the algorithm can be parameterized by this threshold.

Overall, to use the Desbordante profiler, specify the dataset file, task,
algorithm and its parameters. The results will be written to the specified
output file or to console, if none is specified.

Currently, the console version of Desbordante supports:
1) Discovery of exact functional dependencies
2) Discovery of approximate functional dependencies
3) Discovery of probabilistic functional dependencies
4) Discovery of association rules
5) Discovery of exact order dependencies (set-based and list-based axiomatization)
6) Discovery of inclusion dependencies
7) Verification of exact functional dependencies
8) Verification of approximate functional dependencies
9) Verification of metric dependencies
10) Verification of exact unique column combinations
11) Verification of approximate unique column combinations
If you need other types, you should look into the C++ code, the Python
bindings or the Web version.


--task=TASK
    specify the task to run, e.g., discovery of functional dependencies

--algo=ALGORITHM
    specify the algorithm to run, e.g., PYRO

--table=TABLE
    specify the input file to be processed by the algorithm.
    Algorithms for some tasks (currently, only IND) accept multiple
    input files; see --task=TASK for more information

--filename=FILENAME
    specify the file to write the results to. If none is selected, output is
    written to the console

--verbose
    print detailed information before the result

--help
    display this help or help page of the algorithm and task options
    (--task=TASK --help | --algo=ALGO --help); should be specified after
    --algo|--task
'''
FD_HELP = '''Discover minimal non-trivial exact functional dependencies. For
more information about the primitive and the algorithms, refer to the
“Functional dependency discovery: an experimental evaluation of seven
algorithms” paper by T. Papenbrock et al.

Algorithms: PYRO, TANE, HYFD, FD_MINE, DFD, DEP_MINER, FDEP, FUN, FASTFDS, AID
Default: HYFD
'''
CFD_HELP = '''Discover approximate conditional functional dependencies. For
more information about the primitive and the algorithm, refer to the “Revisiting
Conditional Functional Dependency Discovery: Splitting the “C” from the “FD””
paper by J. Rammelaere and F. Geerts.

Algorithms: FD_FIRST
Default: FD_FIRST
'''
AFD_HELP = '''Discover minimal non-trivial approximate functional
dependencies. Approximate functional dependencies are defined in the
“Efficient Discovery of Approximate Dependencies” paper by S. Kruse and
F. Naumann.

Algorithms: PYRO, TANE
Default: PYRO
'''
OD_HELP = '''Discover order dependencies. For more information about the 
primitive and algorithms, refer to the “Effective and complete discovery 
of order dependencies via set-based axiomatization” paper by J. Szlichta 
et al.

Algorithms: FASTOD, ORDER
Default: FASTOD
'''
PFD_HELP = '''Discover minimal non-trivial probabilistic functional
dependencies. Probabilitistic functional dependencies are defined in the
“Functional Dependency Generation and Applications in pay-as-you-go 
data integration systems” paper by Daisy Zhe Wang et al.
Algorithms: PFDTANE
Default: PFDTANE
'''
IND_HELP = '''Discover inclusion dependecies. For more information about
inclusion dependecies, refer to the "Inclusion Dependency Discovery: An
Experimental Evaluation of Thirteen Algorithms" by Falco Dürsch et al.
Algorithms for this task accept multiple input files. You can use one of the
following options:

--tables=TABLE
    specify input files to be processed by the algorithm.
    For multiple values, specify multiple times
    (e.g., --tables=TABLE_1 --tables=TABLE_2)

--tables_list=FILENAME
    specify file with list of input files (one on a line).
    You can type --tables_list=- to use stdin

--tables_directory=FILENAME, STRING, BOOLEAN
    specify directory with input files.
    separator and has_header are applied to all tables

Algorithms: SPIDER, FAIDA
Default: SPIDER
'''
AR_HELP = '''Discover association rules. For more information, refer to
"Frequent Pattern Mining" book by Charu C. Aggarwal and Jiawei Han.

Algorithms: Apriori
Default: Apriori
'''
FD_VERIFICATION_HELP = '''Verify whether a given exact functional dependency
holds on the specified dataset. For more information about the primitive and
algorithms, refer to the “Functional dependency discovery: an experimental
evaluation of seven algorithms” paper by T. Papenbrock et al.

Algorithms: NAIVE_FD_VERIFIER
Default: NAIVE_FD_VERIFIER
'''
AFD_VERIFICATION_HELP = '''Verify whether a given approximate functional
dependency holds on the specified dataset. Approximate functional
dependencies are defined in the “Efficient Discovery of Approximate
Dependencies” paper by S. Kruse and F. Naumann.

Algorithms: NAIVE_AFD_VERIFIER
Default: NAIVE_AFD_VERIFIER
'''
MFD_VERIFICATION_HELP = '''Verify whether a given metric functional
dependency holds on the specified dataset. For more information about the
primitive and the algorithms, refer to “Metric Functional Dependencies” by
N. Koudas et al.

Algorithms: ICDE09_MFD_VERIFIER
Default: ICDE09_MFD_VERIFIER
'''
UCC_VERIFICATION_HELP = '''Verify whether a given unique column combination
holds on the specified dataset. For more information about the primitive and 
the algorithms, refer to "Efficient Discovery of Approximate Dependencies" by 
S. Kruse and F. Naumann

Algorithms: NAIVE_UCC_VERIFIER
Default: NAIVE_UCC_VERIFIER
'''
AUCC_VERIFICATION_HELP = '''Verify whether a given approximate unique column combination
holds on the specified dataset. For more information about the primitive and 
the algorithms, refer to "Efficient Discovery of Approximate Dependencies" by 
S. Kruse and F. Naumann

Algorithms: NAIVE_AUCC_VERIFIER
Default: NAIVE_AUCC_VERIFIER
'''
GFD_VERIFICATION_HELP = '''
Algorithms: NAIVE_GFD_VERIFIER, GFD_VERIFIER, EGFD_VERIFIER
'''
PYRO_HELP = '''A modern algorithm for discovery of approximate functional
dependencies. Approximate functional dependencies are defined in the
“Efficient Discovery of Approximate Dependencies” paper by S.Kruse and
F.Naumann. Capable of discovering exact dependencies too.
'''
TANE_HELP = '''A classic algorithm for discovery of exact and approximate
functional dependencies. For more information, refer to “TANE : An Efficient
Algorithm for Discovering Functional and Approximate Dependencies” by
Y. Huntala et al.
'''
PFDTANE_HELP = '''A TANE-based algorithm for discovery of probabilistic
functional dependencies. For more information, refer to “Functional Dependency
Generation and Applications in pay-as-you-go data integration systems” by
Daisy Zhe Wang et al.
'''
HYFD_HELP = '''A modern algorithm for discovery of exact functional
dependencies. One of the most high-performance algorithms for this task. For
more information, refer to “A Hybrid Approach to Functional Dependency
Discovery” by T. Papenbrock and F. Naumann.
'''
FD_MINE_HELP = '''A classic algorithm for discovery of exact functional
dependencies. Has issues with the minimality of answer. For more
information, refer to the “FD_Mine: discovering functional dependencies in a
database using equivalences paper” by H. Yao et al.
'''
DFD_HELP = '''A modern algorithm for discovery of exact functional
dependencies. For more information, refer to the “DFD: Efficient Functional
Dependency Discovery” paper by Z. Abedjan et al.
'''
DEP_MINER_HELP = '''A classic algorithm for discovery of exact functional
dependencies. For more information refer to “Efficient Discovery of
Functional Dependencies and Armstrong Relations” paper by S. Lopes et al.
'''
FDEP_HELP = '''A classic algorithm for discovery of exact functional
dependencies. For more information, refer to the “Database Dependency
Discovery: A Machine Learning Approach” paper by Peter A. Flach and
Iztok Savnik.
'''
FUN_HELP = '''A classic algorithm for discovery of exact functional
dependencies. For more information,  refer to the “FUN: An efficient
algorithm for mining functional and embedded dependencies” paper by
N. Novelli and R. Cicchetti.
'''
FASTFDS_HELP = '''A classic algorithm for discovery of exact functional
dependencies. For more information, refer to “FastFDs: A Heuristic-Driven,
Depth-First Algorithm for Mining Functional Dependencies from Relation
Instances Extended Abstract” paper by C. Wyss et al.
'''
AID_HELP = '''A modern algorithm for discovery of exact functional
dependencies. Unlike all other algorithms, it is approximate, i.e. it can
miss some dependencies or produce non-valid ones. In exchange,
it is significantly faster (10x-100x). For more information, refer to the
“Approximate Discovery of Functional Dependencies for Large Datasets” paper
by T.Bleifus et al.
'''
SPIDER_HELP = '''A disk-backed unary inclusion dependency mining algorithm.
For more information, refer to "Efficiently detecting inclusion dependencies"
by J. Bauckmann et al.
'''
FAIDA_HELP = '''Both unary and n-ary inclusion dependency mining algorithm.
Unlike all other algorithms, it is approximate, i.e. it can
miss some dependencies or produce non-valid ones. In exchange,
it is significantly faster. For more information, refer to "Fast approximate
discovery of inclusion dependencies" by S. Kruse et al.
'''
FASTOD_HELP = '''A modern algorithm for discovery of canonical order 
dependencies. For more information, refer to the “Effective and complete 
discovery of order dependencies via set-based axiomatization” paper by 
J. Szlichta et al.
'''
ORDER_HELP = '''Algorithm Order efficiently discovers all n-ary lexicographical 
order dependencies under the operator “<”. For more information, refer to the
“Efficient order dependency detection” paper by Philipp Langer and Felix Naumann.
'''
FD_FIRST_HELP = '''FD-First algorithm belongs to the family of algorithms
for discovering approximate conditional functional dependencies. For more
information, refer to the “Revisiting Conditional Functional Dependency
Discovery: Splitting the “C” from the “FD”” paper by J. Rammelaere 
and F. Geerts.
'''
NAIVE_FD_VERIFIER_HELP = '''A straightforward partition-based algorithm for
verifying whether a given exact functional dependency holds on the specified
dataset. For more information, refer to Lemma 2.2 from “TANE: An Efficient
Algorithm for Discovering Functional and Approximate Dependencies” by
Y.Huntala et al.
'''
NAIVE_AFD_VERIFIER_HELP = '''A straightforward partition-based algorithm for
verifying whether a given approximate dependency holds. For more
information, refer to Section 2 of “TANE : An Efficient Algorithm for
Discovering Functional and Approximate Dependencies” by Y.Huntala et al. We
also recommend looking into “Efficient Discovery of ApproximateDependencies" by
S. Kruse and F. Naumann.
'''
ICDE09_MFD_VERIFIER_HELP = '''A family of metric functional dependency
verification algorithms. For more information about the primitive and the
algorithms, refer to “Metric Functional Dependencies” by N. Koudas et al.
'''
GFD_VERIFIER_HELP = '''Algorithm for verifying whether a given
graph functional dependency holds. For more information about the primitive
refer to “Functional Dependencies for Graphs” by Wenfei Fan et al.
'''
NAIVE_UCC_VERIFIER_HELP = '''A straightforward partition-based algorithm for
verifying whether a given unique column combination holds.
For more information on partitions refer to Section 2 of “TANE : An 
Efficient Algorithm for Discovering Functional and Approximate Dependencies”
by Y.Huntala et al. For more information on UCC, refer to "Efficient Discovery
of Approximate Dependencies" by S. Kruse and F. Naumann.
'''
NAIVE_AUCC_VERIFIER_HELP = '''A straightforward partition-based algorithm for
verifying whether a given approximate unique column combination holds.
For more information on partitions refer to Section 2 of “TANE : An 
Efficient Algorithm for Discovering Functional and Approximate Dependencies”
by Y.Huntala et al. For more information on AUCC, refer to "Efficient Discovery
of Approximate Dependencies" by S. Kruse and F. Naumann.
'''
APRIORI_HELP = '''An algorithm for frequent item set mining and association 
rule discovery. For more information, refer to the "Fast Algorithms for 
Mining Association Rules" paper by Agrawal and Srikant from 1994.
'''

OPTION_TYPES = {
    str: 'STRING',
    int: 'INTEGER',
    float: 'FLOAT',
    bool: 'BOOLEAN'
}

TASK_HELP_PAGES = {
    Task.fd: FD_HELP,
    Task.cfd: CFD_HELP,
    Task.afd: AFD_HELP,
    Task.od: OD_HELP,
    Task.pfd: PFD_HELP,
    Task.ind: IND_HELP,
    Task.ar: AR_HELP,
    Task.fd_verification: FD_VERIFICATION_HELP,
    Task.afd_verification: AFD_VERIFICATION_HELP,
    Task.mfd_verification: MFD_VERIFICATION_HELP,
    Task.ucc_verification: UCC_VERIFICATION_HELP,
    Task.aucc_verification: AUCC_VERIFICATION_HELP,
    Task.gfd_verification: GFD_VERIFICATION_HELP,
}

ALGO_HELP_PAGES = {
    Algorithm.pyro: PYRO_HELP,
    Algorithm.tane: TANE_HELP,
    Algorithm.pfdtane: PFDTANE_HELP,
    Algorithm.hyfd: HYFD_HELP,
    Algorithm.fd_mine: FD_MINE_HELP,
    Algorithm.dfd: DFD_HELP,
    Algorithm.dep_miner: DEP_MINER_HELP,
    Algorithm.fdep: FDEP_HELP,
    Algorithm.fun: FUN_HELP,
    Algorithm.fastfds: FASTFDS_HELP,
    Algorithm.aid: AID_HELP,
    Algorithm.fastod: FASTOD_HELP,
    Algorithm.order: ORDER_HELP,
    Algorithm.spider: SPIDER_HELP,
    Algorithm.faida: FAIDA_HELP,
    Algorithm.fd_first: FD_FIRST_HELP,
    Algorithm.naive_fd_verifier: NAIVE_FD_VERIFIER_HELP,
    Algorithm.naive_afd_verifier: NAIVE_AFD_VERIFIER_HELP,
    Algorithm.icde09_mfd_verifier: ICDE09_MFD_VERIFIER_HELP,
    Algorithm.naive_ucc_verifier: NAIVE_UCC_VERIFIER_HELP,
    Algorithm.naive_aucc_verifier: NAIVE_AUCC_VERIFIER_HELP,
    Algorithm.naive_gfd_verifier: GFD_VERIFIER_HELP,
    Algorithm.gfd_verifier: GFD_VERIFIER_HELP,
    Algorithm.egfd_verifier: GFD_VERIFIER_HELP,
    Algorithm.apriori: APRIORI_HELP
}

TaskInfo = namedtuple('TaskInfo', ['algos', 'default'])

TASK_INFO = {
    Task.fd: TaskInfo([Algorithm.pyro, Algorithm.tane, Algorithm.hyfd,
                       Algorithm.fd_mine, Algorithm.dfd, Algorithm.dep_miner,
                       Algorithm.fdep, Algorithm.fun, Algorithm.fastfds,
                       Algorithm.aid, Algorithm.pfdtane],
                      Algorithm.hyfd),
    Task.cfd: TaskInfo([Algorithm.fd_first],
                       Algorithm.fd_first),
    Task.afd: TaskInfo([Algorithm.pyro, Algorithm.tane],
                       Algorithm.pyro),
    Task.od: TaskInfo([Algorithm.fastod, Algorithm.order],
                      Algorithm.fastod),
    Task.pfd: TaskInfo([Algorithm.pfdtane], Algorithm.pfdtane),
    Task.ind: TaskInfo([Algorithm.spider, Algorithm.faida],
                       Algorithm.spider),
    Task.ar: TaskInfo([Algorithm.apriori],
                      Algorithm.apriori),
    Task.fd_verification: TaskInfo([Algorithm.naive_fd_verifier],
                                   Algorithm.naive_fd_verifier),
    Task.afd_verification: TaskInfo([Algorithm.naive_afd_verifier],
                                    Algorithm.naive_afd_verifier),
    Task.mfd_verification: TaskInfo([Algorithm.icde09_mfd_verifier],
                                    Algorithm.icde09_mfd_verifier),
    Task.ucc_verification: TaskInfo([Algorithm.naive_ucc_verifier],
                                    Algorithm.naive_ucc_verifier),
    Task.aucc_verification: TaskInfo([Algorithm.naive_aucc_verifier],
                                     Algorithm.naive_aucc_verifier),
    Task.gfd_verification: TaskInfo([Algorithm.naive_gfd_verifier,
                                     Algorithm.gfd_verifier,
                                     Algorithm.egfd_verifier],
                                    Algorithm.naive_gfd_verifier),
}

ALGOS = {
    Algorithm.pyro: desbordante.fd.algorithms.Pyro,
    Algorithm.tane: desbordante.fd.algorithms.Tane,
    Algorithm.pfdtane: desbordante.pfd.algorithms.PFDTane,
    Algorithm.hyfd: desbordante.fd.algorithms.HyFD,
    Algorithm.fd_mine: desbordante.fd.algorithms.FdMine,
    Algorithm.dfd: desbordante.fd.algorithms.DFD,
    Algorithm.dep_miner: desbordante.fd.algorithms.Depminer,
    Algorithm.fdep: desbordante.fd.algorithms.FDep,
    Algorithm.fun: desbordante.fd.algorithms.FUN,
    Algorithm.fastfds: desbordante.fd.algorithms.FastFDs,
    Algorithm.aid: desbordante.fd.algorithms.Aid,
    Algorithm.fastod: desbordante.od.algorithms.Fastod,
    Algorithm.order: desbordante.od.algorithms.Order,
    Algorithm.spider: desbordante.ind.algorithms.Spider,
    Algorithm.faida: desbordante.ind.algorithms.Faida,
    Algorithm.fd_first: desbordante.cfd.algorithms.FDFirst,
    Algorithm.naive_fd_verifier: desbordante.fd_verification.algorithms.FDVerifier,
    Algorithm.naive_afd_verifier: desbordante.afd_verification.algorithms.FDVerifier,
    Algorithm.icde09_mfd_verifier: desbordante.mfd_verification.algorithms.MetricVerifier,
    Algorithm.naive_ucc_verifier: desbordante.ucc_verification.algorithms.UccVerifier,
    Algorithm.naive_aucc_verifier: desbordante.aucc_verification.algorithms.UccVerifier,
    Algorithm.naive_gfd_verifier: desbordante.gfd_verification.algorithms.NaiveGfdValid,
    Algorithm.gfd_verifier: desbordante.gfd_verification.algorithms.GfdValid,
    Algorithm.egfd_verifier: desbordante.gfd_verification.algorithms.EGfdValid,
    Algorithm.apriori: desbordante.ar.algorithms.Apriori
}


def get_algorithm(ctx: click.core.Context, param: click.core.Option,
                  value: str | None) -> str:
    if value is not None:
        return value
    task = ctx.params.get(TASK)
    if task is None:
        raise click.MissingParameter(ctx=ctx, param=param)
    return TASK_INFO[task].default


def help_callback_func(ctx: click.core.Context, param: click.core.Option,
                       value: str | None) -> None:
    if not value or ctx.resilient_parsing:
        return
    algo = ctx.params.get(ALGO)
    task = ctx.params.get(TASK)
    check_mismatch(algo, task)
    print_help_page(algo, task)
    ctx.exit()


def check_mismatch(algo: str | None, task: str | None) -> None:
    if (task is not None and
            (algo is not None and algo not in TASK_INFO[Task(task)].algos)):
        click.echo(f"ERROR: Wrong value for '{ALGO}' when {TASK}={task}.")
        sys.exit(1)


def check_error_option_presence(task: str | None, error: str | None) -> None:
    if task in (Task.afd, Task.pfd) and error is None:
        click.echo(f"ERROR: Missing option '{ERROR}'.")
        sys.exit(1)
    if task in (Task.fd, Task.fd_verification) and error is not None:
        click.echo(f"ERROR: Invalid option: '{ERROR}'.")
        sys.exit(1)


def check_error_measure_option_presence(task: str | None, error_measure: str | None) -> None:
    if task == Task.pfd and error_measure is None:
        click.echo(f"ERROR: Missing option '{ERROR_MEASURE}'.")
        sys.exit(1)


def parse_tables_list_file(file: click.File) \
        -> list[tuple[str, str, bool]]:
    try:
        result = []
        for line_num, line in enumerate(file.readlines(), start=1):
            table_tuple = line.rsplit(maxsplit=2)
            if len(table_tuple) != 3:
                click.echo(
                    f'ERROR: Invalid format of table description on line {line_num}: {line}')
                sys.exit(1)
            filename, separator, has_header_str = table_tuple
            result.append((filename, separator, bool(has_header_str)))
        return result
    except OSError as exc:
        click.echo(exc)
        sys.exit(1)


def parse_tables_directory(tp: tuple[click.Path, str, bool]) \
        -> list[tuple[str, str, bool]]:
    dir_name, separator, has_header = tp
    try:
        entries = scandir(dir_name)
        return [(dir_entry.path, separator, has_header) for dir_entry in entries]
    except OSError as exc:
        click.echo(exc)
        sys.exit(1)


def is_omitted(value: Any) -> bool:
    return value is None or value == ()


def set_option(algo: desbordante.Algorithm, opt_name: str, opt_value: Any) \
        -> None:
    try:
        algo.set_option(opt_name, opt_value)
    except Exception as exc:
        click.echo(exc)
        sys.exit(1)


def set_algo_options(algo: desbordante.Algorithm, args: dict[str, Any]) -> set:
    used_options = set()
    while opts := algo.get_needed_options():
        for option_name in opts:
            value = args[option_name]
            if is_omitted(value):
                set_option(algo, option_name, None)
            else:
                set_option(algo, option_name, value)
                used_options.add(option_name)
    return used_options


def get_algo_result(algo: desbordante.Algorithm, algo_name: str) -> Any:
    try:
        algo.execute()
        match algo_name:
            case Algorithm.naive_fd_verifier:
                result = algo.fd_holds()
            case Algorithm.naive_afd_verifier:
                error = algo.get_error()
                if error == 0.0:
                    result = 'Exact functional dependency holds'
                else:
                    result = (f'Exact functional dependency does not hold, but '
                              f'instead approximate functional dependency '
                              f'holds with error = {error}')
            case Algorithm.naive_ucc_verifier:
                result = algo.ucc_holds()
            case Algorithm.naive_aucc_verifier:
                error = algo.get_error()
                if error == 0.0:
                    result = 'Exact unique column combination holds'
                else:
                    result = (f'Exact unique column combination does not hold, but '
                              f'instead approximate unique column combination '
                              f'holds with error = {error}')
            case Algorithm.icde09_mfd_verifier:
                result = algo.mfd_holds()
            case algo_name if algo_name in TASK_INFO[Task.fd].algos:
                result = algo.get_fds()
            case Algorithm.fastod:
                result = algo.get_asc_ods() + algo.get_desc_ods() + algo.get_simple_ods()
            case Algorithm.order:
                result = algo.get_list_ods()
            case algo_name if algo_name in TASK_INFO[Task.ind].algos:
                result = algo.get_inds()
            case algo_name if algo_name in TASK_INFO[Task.gfd_verification].algos:
                result = algo.get_gfds()
            case Algorithm.fd_first:
                result = algo.get_cfds()
            case Algorithm.apriori:
                result = algo.get_ars()
            case _:
                assert False, 'No matching get_result function.'
        return result
    except Exception as exc:
        click.echo(exc)
        sys.exit(1)


def stringify_result(result: Any) -> str:
    if type(result) is list:
        return '\n'.join(map(str, result))
    else:
        return str(result)


def create_verbose_output(provided_opts: dict[str, Any], start_time: float,
                          end_time: float, result: Any) -> str:
    option_info = '\n'.join(f'{name}: {value}' for name, value in
                            provided_opts.items())
    time = f'Elapsed time during execution in seconds: {end_time - start_time}'
    result = f'Result:\n{stringify_result(result)}'
    return f'{option_info}\n{time}\n{result}'


def print_result(result: Any, filename: str | None) -> None:
    if filename is None:
        click.echo(stringify_result(result))
        return
    try:
        with open(filename, 'w') as file:
            print(result, file=file)
    except OSError as exc:
        click.echo(exc)
        sys.exit(1)


def print_unused_opts(used_opts: set, provided_opts: set) -> None:
    unused_opts = provided_opts - (used_opts | {TASK, ALGO, VERBOSE, FILENAME} |
                                   ({TABLES_LIST, TABLES_DIRECTORY}
                                    if TABLES in used_opts
                                    else set()))
    if unused_opts:
        click.echo(f'Unused options: {unused_opts}')


def print_help_page(algo_name: str | None, task: str | None) -> None:
    if algo_name is not None:
        print_algo_help_page(algo_name)
    elif task is not None:
        click.echo(TASK_HELP_PAGES[Task(task)])
    else:
        click.echo(PRIMARY_HELP)


def print_algo_help_page(algo_name: str) -> None:
    algo = ALGOS[Algorithm(algo_name)]()
    help_info = ''
    algo_options = list(algo.get_possible_options())
    algo_options.sort()
    for opt in algo_options:
        if opt not in (TABLE, TABLES):
            help_info += get_option_help_info(opt, algo)
    click.echo(f'{ALGO_HELP_PAGES[Algorithm(algo_name)]}{help_info}')


def get_provided_options(all_option_dict: dict[str, Any]) -> dict[str, Any]:
    return {name: value for name, value in all_option_dict.items()
            if not is_omitted(value)}


def get_option_help_info(opt: str, algo: desbordante.Algorithm) -> str:
    help_info = ''
    opt_main_type, *opt_additional_types = algo.get_option_type(opt)
    opt_help_type = opt_additional_types[0] if opt_main_type == list \
        else opt_main_type
    help_info = (f'{help_info}\n'
                 f'--{opt}={OPTION_TYPES[opt_help_type]}\n'
                 f'\t{algo.get_description(opt)}\n')
    if opt_main_type == list:
        help_info += (f'\tFor multiple values, specify multiple times '
                      f'\n\t(e.g., --{opt}=1 --{opt}=2)\n')
    return help_info


def get_option_type_info() -> dict[str, Any]:
    option_type_info: dict[str, Any] = {}
    for algo_name, algo_type in ALGOS.items():
        algo = algo_type()
        for opt in algo.get_possible_options():
            option_type = algo.get_option_type(opt)
            previous_type = option_type_info.setdefault(opt, option_type)
            assert option_type == previous_type, \
                (f"Different types for '{opt}' option"
                 f'({previous_type=}, {option_type=}).')
            option_type_info[opt] = algo.get_option_type(opt)
    return option_type_info


def process_tables_options(opts: dict[str, Any], algo_name: str) -> dict[str, Any]:
    result = opts.copy()

    for option_name, parse_func in ((TABLES_LIST, parse_tables_list_file),
            (TABLES_DIRECTORY, parse_tables_directory)):
        value = result.pop(option_name)
        if not is_omitted(value):
            result[TABLES] = list(result[TABLES]) + parse_func(value)
    
    return result


def algos_options() -> Callable:
    option_type_info = get_option_type_info()

    def decorator(func: Callable) -> Callable:
        for opt_name, (opt_main_type, *opt_additional_types) \
                in option_type_info.items():
            arg = f'--{opt_name}'
            if opt_main_type == list:
                if opt_additional_types[0] == desbordante.data_types.Table:
                    click.option(arg, type=(str, str, bool),
                                 multiple=True)(func)
                else:
                    click.option(arg, multiple=True,
                             type=opt_additional_types[0])(func)
            elif opt_main_type == desbordante.data_types.Table:
                click.option(arg, type=(str, str, bool))(func)
            else:
                click.option(arg, type=opt_main_type)(func)
        return func

    return decorator


@click.command(add_help_option=False)
@click.option(f'--{HELP}', is_flag=True,
              callback=help_callback_func,
              expose_value=False, is_eager=True)
@click.option(f'--{TASK}', type=click.Choice([str(task) for task in
                                              TASK_INFO.keys()],
                                             case_sensitive=False),
              is_eager=True)
@click.option(f'--{ALGO}', f'--{ALGORITHM}',
              type=click.Choice([str(algo) for algo in ALGOS.keys()],
                                case_sensitive=False),
              callback=get_algorithm, is_eager=True)
@click.option(f'--{FILENAME}', type=str)
@click.option(f'--{VERBOSE}', is_flag=True)
@click.option(f'--{TABLES_LIST}', type=click.File('r'))
@click.option(f'--{TABLES_DIRECTORY}', type=(click.Path(exists=True, file_okay=False,
               dir_okay=True, resolve_path=True, allow_dash=False), str, bool))
@algos_options()
def desbordante_cli(**kwargs: Any) -> None:
    """Takes in options from console as a dictionary, sets these options
    for the selected algo, runs algo and prints the result"""

    curr_task = kwargs[TASK]
    curr_algo_name = kwargs[ALGO]
    curr_algo = ALGOS[curr_algo_name]()
    error_opt = kwargs[ERROR]
    error_measure_opt = kwargs[ERROR_MEASURE]
    verbose = kwargs[VERBOSE]
    filename = kwargs[FILENAME]

    check_mismatch(curr_algo_name, curr_task)
    check_error_option_presence(curr_task, error_opt)
    check_error_measure_option_presence(curr_task, error_measure_opt)

    opts = process_tables_options(kwargs, curr_algo_name)

    start_point = process_time()
    used_opts = set_algo_options(curr_algo, opts)
    curr_algo.load_data()
    used_opts |= set_algo_options(curr_algo, opts)
    provided_options = get_provided_options(kwargs)
    print_unused_opts(used_opts, set(provided_options.keys()))
    result = get_algo_result(curr_algo, curr_algo_name)
    end_point = process_time()

    if verbose:
        verbose_output = create_verbose_output(provided_options,
                                               start_point, end_point, result)
        print_result(verbose_output, filename)
    else:
        print_result(result, filename)


if __name__ == '__main__':
    desbordante_cli()
