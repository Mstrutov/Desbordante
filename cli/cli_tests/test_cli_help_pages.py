import unittest

import snapshottest
from click.testing import CliRunner

from cli import desbordante_cli, Algorithm, Task


class TestCLIHelpPages(snapshottest.TestCase):
    def test_main_help_page(self):
        runner = CliRunner()
        with self.subTest(msg=f'Testing main help page'):
            result = runner.invoke(desbordante_cli, f'--help').output
            self.assertMatchSnapshot(result, f'main_help')

    def test_algos_help_pages(self):
        runner = CliRunner()
        for algo in Algorithm:
            with self.subTest(msg=f'Testing help page for {algo}'):
                result = runner.invoke(desbordante_cli,
                                       f'--algo={algo} --help').output
                self.assertMatchSnapshot(result, f'{algo}_help')

    def test_tasks_help_pages(self):
        runner = CliRunner()
        for task in Task:
            with self.subTest(msg=f'Testing help page for {task}'):
                result = runner.invoke(desbordante_cli,
                                       f'--task={task} --help').output
                self.assertMatchSnapshot(result, f'{task}_help')


if __name__ == '__main__':
    import snapshottest.unittest
    import os

    update = os.getenv('UPDATE_HELP_PAGES', False)
    if update:
        snapshottest.unittest.TestCase.snapshot_should_update = True
    unittest.main()
