import unittest
from collections import namedtuple
from itertools import chain

import desbordante as desb

OptionContainer = namedtuple("OptionContainer", ['path', 'load_options', 'execute_options'])
FailureCaseContainer = namedtuple("FailureCaseContainer", ['path', 'options'])

ONLY_NULL_EQUAL_NULL_OPTION_CONTAINER = OptionContainer(
    "WDC_satellites.csv", {"is_null_equal_null": False}, {}
)


def get_common_option_container(execute_options):
    return OptionContainer(
        "WDC_satellites.csv", {"is_null_equal_null": False}, execute_options
    )


def get_apriori_load_container(load_options):
    return OptionContainer("TestWide.csv", load_options, {})


def check_metric_verifier_failure(dataset, options) -> bool:
    alg = desb.mfd_verification.algorithms.MetricVerifier()
    alg.load_data(table=(dataset, ",", True))
    for opt_name in options:
        alg.set_option(opt_name, options[opt_name])


ALGO_CORRECT_OPTIONS_INFO = [
    (desb.fd.algorithms.Depminer, [ONLY_NULL_EQUAL_NULL_OPTION_CONTAINER]),
    (desb.fd.algorithms.FUN, [ONLY_NULL_EQUAL_NULL_OPTION_CONTAINER]),
    (desb.fd.algorithms.FdMine, [ONLY_NULL_EQUAL_NULL_OPTION_CONTAINER]),
    (desb.fd.algorithms.HyFD, [ONLY_NULL_EQUAL_NULL_OPTION_CONTAINER]),
    (desb.fd.algorithms.EulerFD, [ONLY_NULL_EQUAL_NULL_OPTION_CONTAINER]),
    (desb.afd.algorithms.Pyro, [
        get_common_option_container(
            {"seed": 1, "max_lhs": 12, "threads": 5, "error": 0.015}
        ),
    ]),
    (desb.fd.algorithms.DFD, [
        get_common_option_container(
            {
                "threads": 15,
            }
        ),
    ]),
    (desb.fd.algorithms.FastFDs, [
        get_common_option_container({"max_lhs": 12, "threads": 15}),
    ]),
    (desb.afd.algorithms.Tane, [
        get_common_option_container({"max_lhs": 12, "error": 0.015}),
    ]),
    (desb.statistics.algorithms.DataStats, [
        get_common_option_container({"threads": 15}),
    ]),
    (desb.ucc.algorithms.HyUCC, [
        get_common_option_container({"threads": 15}),
    ]),
    (desb.fd.algorithms.EulerFD, [
        get_common_option_container({"custom_random": (False, 102)}),
    ]),
    (desb.fd_verification.algorithms.FDVerifier, [
        get_common_option_container(
            {"lhs_indices": [1, 2, 3], "rhs_indices": [1, 2, 3]}
        ),
    ]),
    (desb.ar.algorithms.Apriori, [
        get_apriori_load_container({"input_format": "tabular", "has_tid": True}),
        get_apriori_load_container({"input_format": "tabular", "has_tid": False}),
        get_apriori_load_container(
            {"input_format": "singular", "tid_column_index": 0, "item_column_index": 2}
        ),
        OptionContainer(
            "rules-kaggle-rows.csv",
            {
                "input_format": "singular",
                "tid_column_index": 0,
                "item_column_index": 1,
            },
            {"minconf": 0.00312, "minsup": 0.2321},
        ),
    ]),
    (desb.mfd_verification.algorithms.MetricVerifier, [
        OptionContainer(
            "TestLong.csv",
            {},
            {
                "metric": "euclidean",
                "rhs_indices": [1, 2],
                "parameter": 213.213111,
                "dist_from_null_is_infinity": False,
                "metric_algorithm": "approx",
                "lhs_indices": [0, 1, 2],
            },
        ),
        OptionContainer(
            "WDC_satellites.csv",
            {},
            {
                "metric": "levenshtein",
                "rhs_indices": [1],
                "parameter": 213.213111,
                "dist_from_null_is_infinity": False,
                "metric_algorithm": "approx",
                "lhs_indices": [0, 1, 2],
            },
        ),
        OptionContainer(
            "WDC_satellites.csv",
            {},
            {
                "metric": "cosine",
                "rhs_indices": [1],
                "parameter": 213.213111,
                "dist_from_null_is_infinity": False,
                "q": 123,
                "metric_algorithm": "approx",
                "lhs_indices": [0, 1, 2],
            },
        ),
    ]),
]

METRIC_VERIFIER_FAILURE_CASES = [
    FailureCaseContainer(
        "WDC_satellites.csv",
        {
            "metric": "euclidean",
            "rhs_indices": [1, 2],
        }
    ),
    
    FailureCaseContainer(
        "TestLong.csv",
        {
            "metric": "euclidean",
            "q": 123,
        }
    ),
    FailureCaseContainer(
        "TestLong.csv",
        {
            "metric": "euclidean",
            "metric_algorithm": "brute",
        }
    ),
    FailureCaseContainer(
        "TestLong.csv",
        {"metric": "levenshtein", "rhs_indices": [1, 2]}
    ),
    FailureCaseContainer("WDC_satellites.csv", {"metric": "levenshtein", "q": 123}),
    FailureCaseContainer(
        "WDC_satellites.csv",
        {
            "metric": "levenshtein",
            "metric_algorithm": "brute",
        }
    ),
    FailureCaseContainer("TestLong.csv", {"metric": "levenshtein", "rhs_indices": [1]}),
    FailureCaseContainer("TestLong.csv", {"rhs_indices": [1]}),
    FailureCaseContainer("WDC_satellites.csv", {"rhs_indices": [1, 2]}),
    FailureCaseContainer("WDC_satellites.csv", {"q": 123}),
    FailureCaseContainer("WDC_satellites.csv", {"metric_algorithm": "approx"}),
]


class TestPythonBindings(unittest.TestCase):

    def _test_correct_option_setting(
        self, algo, path, separator, has_header, options: dict
    ):
        testing_algo = algo()
        testing_algo.load_data(table=(path, separator, has_header), **options.load_options)
        for name, value in options.execute_options.items():
            testing_algo.set_option(name, value)

        algo_option_values = testing_algo.get_opts()
        for name, value in chain(options.execute_options.items(),
                                 options.load_options.items()):
            self.assertEqual(value, algo_option_values[name])

    def test_correct_load(self):
        for algo, option_containers in ALGO_CORRECT_OPTIONS_INFO:
            for option in option_containers:
                with self.subTest(msg=f"testing setting correct options for {algo.__name__}"):
                    self._test_correct_option_setting(
                        algo, option.path, ",", True, option
                    )

    def test_metric_verifier_failure_cases(self):
        for load in METRIC_VERIFIER_FAILURE_CASES:
            with self.subTest(msg=f"metric_verifier_load: {load}"):
                with self.assertRaises(desb.ConfigurationError):
                    check_metric_verifier_failure(load.path, load.options)
                


if __name__ == "__main__":
    unittest.main()
