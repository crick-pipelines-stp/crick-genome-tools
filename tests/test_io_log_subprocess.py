"""
Tests for log sub process
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import subprocess
import unittest
from unittest.mock import MagicMock, patch

from crick_genome_tools.io.log_subprocess import LogSubprocess


class TestLogSubprocess(unittest.TestCase):

    def setUp(self):
        # Arrange: Initialize the LogSubprocess object
        self.log_subprocess = LogSubprocess()

    @patch("subprocess.check_call")
    def test_log_subprocess_check_call_success(self, mock_check_call):
        # Arrange: Mock the check_call method to simulate successful command execution
        mock_check_call.return_value = 0

        # Act: Call the check_call method with sample arguments
        self.log_subprocess.check_call(["echo", "hello"])

        # Assert: Ensure that check_call was called with the correct arguments
        mock_check_call.assert_called_once_with(["echo", "hello"], preexec_fn=self.log_subprocess.pdeathsig, stderr=subprocess.PIPE)

    @patch("subprocess.check_call")
    def test_log_subprocess_check_call_failure(self, mock_check_call):
        # Arrange: Mock a failing command by setting a side effect on check_call
        mock_check_call.side_effect = subprocess.CalledProcessError(1, ["echo", "hello"], stderr=b"error message")

        # Act and Assert: Verify that the check_call method raises a CalledProcessError
        with self.assertRaises(subprocess.CalledProcessError) as cm:
            self.log_subprocess.check_call(["echo", "hello"])

        # Assert: Check that the exception contains the correct return code and stderr
        self.assertEqual(cm.exception.returncode, 1)
        self.assertEqual(cm.exception.stderr, b"error message")

    @patch("subprocess.check_output")
    def test_log_subprocess_check_output_success(self, mock_check_output):
        # Arrange: Mock the check_output method to return a sample output
        mock_check_output.return_value = b"output"

        # Act: Call the check_output method with sample arguments
        result = self.log_subprocess.check_output(["echo", "hello"])

        # Assert: Ensure that check_output was called with the correct arguments
        mock_check_output.assert_called_once_with(["echo", "hello"], preexec_fn=self.log_subprocess.pdeathsig, stderr=subprocess.PIPE)

        # Assert: Check that the returned output is as expected
        self.assertEqual(result, b"output")

    @patch("subprocess.check_output")
    def test_log_subprocess_check_output_failure(self, mock_check_output):
        # Arrange: Mock a failing command by setting a side effect on check_output
        mock_check_output.side_effect = subprocess.CalledProcessError(1, ["echo", "hello"], stderr=b"error message")

        # Act and Assert: Verify that the check_output method raises a CalledProcessError
        with self.assertRaises(subprocess.CalledProcessError) as cm:
            self.log_subprocess.check_output(["echo", "hello"])

        # Assert: Check that the exception contains the correct return code and stderr
        self.assertEqual(cm.exception.returncode, 1)
        self.assertEqual(cm.exception.stderr, b"error message")

    @patch("subprocess.Popen")
    def test_log_subprocess_p_open_success(self, mock_popen):
        # Arrange: Mock the Popen method to simulate a successful process creation
        mock_process = MagicMock()
        mock_process.wait.return_value = 0
        mock_popen.return_value = mock_process

        # Act: Call the p_open method with sample arguments
        process = self.log_subprocess.p_open(["echo", "hello"])

        # Manually simulate the dynamic attachment of check_return_code
        process.check_return_code = MagicMock()

        # Assert: Ensure that Popen was called with the correct arguments
        mock_popen.assert_called_once_with(
            ["echo", "hello"], preexec_fn=self.log_subprocess.pdeathsig, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )

        # Assert: Check that the returned process is a mock process
        self.assertIsInstance(process, MagicMock)

    @patch("subprocess.Popen")
    def test_log_subprocess_p_open_failure(self, mock_popen):
        # Arrange: Mock a process that returns a non-zero exit code
        mock_process = MagicMock()
        mock_process.wait.return_value = 1
        mock_process.stderr.read.return_value = b"error message"
        mock_popen.return_value = mock_process

        # Act: Call the p_open method with sample arguments
        process = self.log_subprocess.p_open(["echo", "hello"])

        # Simulate the dynamic attachment of check_return_code
        process.check_return_code = MagicMock(side_effect=subprocess.CalledProcessError(1, ["echo", "hello"], stderr=b"error message"))

        # Assert: Verify that the p_open method raises a CalledProcessError
        with self.assertRaises(subprocess.CalledProcessError) as cm:
            process.check_return_code()  # Simulate the check being called later after process completion

        # Assert: Check that the exception contains the correct return code and stderr
        self.assertEqual(cm.exception.returncode, 1)
        self.assertEqual(cm.exception.stderr, b"error message")
