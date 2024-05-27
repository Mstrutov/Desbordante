import unittest
from collections import namedtuple
from unittest.mock import patch

import desbordante
from click.testing import CliRunner

from cli import desbordante_cli, ALGOS

UNFIXED_ALGOS = ['naive_gfd_verifier', 'gfd_verifier', 'egfd_verifier', 'apriori']

UNFIXED_OPTS = ['table', 'tables']

OptionInfo = namedtuple('OptionInfo', ['str', 'value'])

OPTION_VALUES = {
    (int,): OptionInfo('4', 4),
    (float,): OptionInfo('0.1', 0.1),
    (bool,): OptionInfo('True', True)
}

UNCOMMON_OPTIONS = {
    'table': OptionInfo("'cli_tests/university.csv' , False", ''),
    'tables': OptionInfo("'cli_tests/university.csv' , False", ''),
    'lhs_indices': OptionInfo('0 --lhs_indices=2', [0, 2]),
    'rhs_indices': OptionInfo('0', [0]),
    'metric': OptionInfo('cosine', 'cosine'),
    'metric_algorithm': OptionInfo('brute', 'brute'),
    'ucc_indices': OptionInfo('0', [0]),
    'error_measure': OptionInfo('per_tuple', 'per_tuple'),
    'mem_limit': OptionInfo('16', 16),
    'cfd_substrategy': OptionInfo('dfs', 'dfs')
}


def get_expected_options(algo):
    algo_opts = algo.get_possible_options()
    expected_options = dict()
    for opt in algo_opts:
        if opt not in UNCOMMON_OPTIONS.keys():
            opt_type = algo.get_option_type(opt)
            expected_options.update({opt: OPTION_VALUES[opt_type].value})
        elif opt not in UNFIXED_OPTS:
            expected_options.update({opt: UNCOMMON_OPTIONS[opt].value})
    return expected_options


def get_invoke_str(algo_name):
    invoke_str = f'--algo={algo_name}'
    algo = ALGOS[algo_name]()
    algo_opts = algo.get_possible_options()
    for opt in algo_opts:
        if opt not in UNCOMMON_OPTIONS.keys():
            opt_type = algo.get_option_type(opt)
            value_as_str = OPTION_VALUES[opt_type].str
        else:
            value_as_str = UNCOMMON_OPTIONS[opt].str
        invoke_str = f'{invoke_str} --{opt}={value_as_str}'
    return invoke_str


def get_algo_options(self, **kwargs):
    return self.get_opts()


def compare_parsing_result(algo, algo_name):
    if algo_name == 'apriori':
        a = 5
    expected_result = get_expected_options(algo)
    cli_parsing_result = algo.execute()
    if cli_parsing_result == expected_result:
        result = 'success'
    else:
        result = 'fail'
    return result


class TestCLIParsing(unittest.TestCase):
    def test_algos(self):
        runner = CliRunner()
        with patch.multiple(desbordante.Algorithm, execute=get_algo_options):
            with patch('cli.get_algo_result', compare_parsing_result):
                for algo_name in ALGOS.keys():
                    if algo_name not in UNFIXED_ALGOS:
                        with self.subTest(
                                msg=f'Testing options parsing for {algo_name}'):
                            invoke_str = get_invoke_str(algo_name)
                            result = runner.invoke(desbordante_cli, invoke_str)
                            result_output = result.output
                            self.assertEqual('success\n', result_output,
                                             msg=f'Failed on {algo_name}')


if __name__ == '__main__':
    unittest.main()
