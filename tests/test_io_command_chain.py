# pylint: disable=missing-function-docstring,missing-class-docstring,invalid-name

import os
import subprocess
import unittest
from unittest.mock import MagicMock, patch

import pytest

from crick_genome_tools.io.command_chain import CommandChain
from tests.utils import with_temporary_folder


class TestCommandChain(unittest.TestCase):

    @patch("crick_genome_tools.io.command_chain.LogSubprocess")
    def test_command_chain_run_return_stream(self, MockLogSubprocess):
        mock_log_subprocess = MockLogSubprocess.return_value
        commands = [["echo", "hello"], ["grep", "hello"]]
        command_chain = CommandChain(commands, return_stream=True)

        mock_proc = MagicMock()
        mock_proc.stdout = MagicMock()
        mock_proc.check_return_code = MagicMock()
        mock_log_subprocess.p_open.side_effect = [mock_proc, mock_proc]
        mock_log_subprocess.stream_process.return_value = iter(["hello\n"])

        result = command_chain.run()

        self.assertEqual(list(result), ["hello\n"])
        self.assertEqual(mock_log_subprocess.p_open.call_count, 2)
        mock_proc.check_return_code.assert_called()

    @with_temporary_folder
    def test_command_chain_run_output_file(self, tmpdir):
        commands = [["echo", "hello"], ["grep", "hello"]]
        output_file = os.path.join(tmpdir, "output.txt")
        command_chain = CommandChain(commands, output_file=output_file)

        command_chain.run()

        with open(output_file, "r", encoding="UTF-8") as f:
            output = f.read()

        self.assertEqual(output, "hello\n")

    @patch("crick_genome_tools.io.command_chain.LogSubprocess")
    def test_command_chain_run_no_output(self, MockLogSubprocess):
        mock_log_subprocess = MockLogSubprocess.return_value
        commands = [["echo", "hello"], ["grep", "hello"]]
        command_chain = CommandChain(commands)

        mock_proc = MagicMock()
        mock_proc.stdout = MagicMock()
        mock_proc.check_return_code = MagicMock()
        mock_log_subprocess.p_open.side_effect = [mock_proc, mock_proc]

        result = command_chain.run()

        self.assertIsNone(result)
        self.assertEqual(mock_log_subprocess.p_open.call_count, 2)
        mock_proc.check_return_code.assert_called()

    @patch("crick_genome_tools.io.command_chain.LogSubprocess")
    def test_command_chain_run_single_command(self, MockLogSubprocess):
        mock_log_subprocess = MockLogSubprocess.return_value
        commands = [["echo", "hello"]]
        command_chain = CommandChain(commands, return_stream=True)

        mock_proc = MagicMock()
        mock_proc.stdout = MagicMock()
        mock_proc.check_return_code = MagicMock()
        mock_log_subprocess.p_open.side_effect = [mock_proc]
        mock_log_subprocess.stream_process.return_value = iter(["hello\n"])

        result = command_chain.run()

        self.assertEqual(list(result), ["hello\n"])
        self.assertEqual(mock_log_subprocess.p_open.call_count, 1)
        mock_proc.check_return_code.assert_called()

    @patch("crick_genome_tools.io.command_chain.LogSubprocess")
    def test_command_chain_run_with_error(self, MockLogSubprocess):
        mock_log_subprocess = MockLogSubprocess.return_value
        commands = [["false"]]
        command_chain = CommandChain(commands)

        mock_proc = MagicMock()
        mock_proc.stdout = MagicMock()
        mock_proc.check_return_code = MagicMock(side_effect=subprocess.CalledProcessError(1, "false"))
        mock_log_subprocess.p_open.side_effect = [mock_proc]

        with self.assertRaises(subprocess.CalledProcessError):
            command_chain.run()

        self.assertEqual(mock_log_subprocess.p_open.call_count, 1)
        mock_proc.check_return_code.assert_called()

    @with_temporary_folder
    def test_command_chain_command_to_file_creates_output_file(self, tmp_path):
        command = ["echo", "hello"]
        output_file = os.path.join(tmp_path, "output.txt")

        CommandChain.command_to_file(command, str(output_file))

        assert os.path.exists(output_file)
        with open(output_file, "r", encoding="UTF-8") as f:
            assert f.read().strip() == "hello"

    @patch("crick_genome_tools.io.command_chain.LogSubprocess")
    @with_temporary_folder
    def test_command_chain_command_to_file_handles_subprocess_error(self, MockLogSubprocess, tmp_path):
        command = ["false"]  # This command will fail
        output_file = os.path.join(tmp_path, "output.txt")

        mock_proc = MockLogSubprocess.Mock()
        MockLogSubprocess.return_value.p_open.return_value = mock_proc
        mock_proc.check_return_code.side_effect = subprocess.CalledProcessError(1, command)

        with pytest.raises(subprocess.CalledProcessError):
            CommandChain.command_to_file(command, str(output_file))
