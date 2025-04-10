"""
This module contains unit tests for the LogSubprocess class.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import subprocess
from unittest.mock import MagicMock, patch

import pytest
from assertpy import assert_that

from crick_genome_tools.io.log_subprocess import LogSubprocess


class TestLogSubprocess:
    @pytest.fixture(autouse=True)
    def setup_log_subprocess(self):
        self.log_subprocess = LogSubprocess()

    @patch("subprocess.check_call")
    def test_logsub_check_call_success(self, mock_check_call):
        mock_check_call.return_value = 0

        self.log_subprocess.check_call(["echo", "hello"])

        mock_check_call.assert_called_once_with(["echo", "hello"], stderr=subprocess.PIPE)

    @patch("subprocess.check_call")
    def test_logsub_check_call_failure(self, mock_check_call):
        mock_check_call.side_effect = subprocess.CalledProcessError(
            1, ["echo", "hello"], stderr=b"error message"
        )

        with pytest.raises(subprocess.CalledProcessError) as exc:
            self.log_subprocess.check_call(["echo", "hello"])

        assert_that(exc.value.returncode).is_equal_to(1)
        assert_that(exc.value.stderr).is_equal_to(b"error message")

    @patch("subprocess.check_output")
    def test_logsub_check_output_success(self, mock_check_output):
        mock_check_output.return_value = b"output"

        result = self.log_subprocess.check_output(["echo", "hello"])

        mock_check_output.assert_called_once_with(["echo", "hello"], stderr=subprocess.PIPE)
        assert_that(result).is_equal_to(b"output")

    @patch("subprocess.check_output")
    def test_logsub_check_output_failure(self, mock_check_output):
        mock_check_output.side_effect = subprocess.CalledProcessError(
            1, ["echo", "hello"], stderr=b"error message"
        )

        with pytest.raises(subprocess.CalledProcessError) as exc:
            self.log_subprocess.check_output(["echo", "hello"])

        assert_that(exc.value.returncode).is_equal_to(1)
        assert_that(exc.value.stderr).is_equal_to(b"error message")

    @patch("subprocess.Popen")
    def test_logsub_p_open_success(self, mock_popen):
        mock_proc = MagicMock()
        mock_proc.wait.return_value = 0
        mock_popen.return_value = mock_proc

        proc = self.log_subprocess.p_open(["echo", "hello"])
        proc.check_return_code = MagicMock()  # simulate dynamic injection

        mock_popen.assert_called_once_with(["echo", "hello"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        assert_that(proc).is_instance_of(MagicMock)

    @patch("subprocess.Popen")
    def test_logsub_p_open_failure_check_return(self, mock_popen):
        mock_proc = MagicMock()
        mock_proc.wait.return_value = 1
        mock_proc.stderr.read.return_value = b"error message"
        mock_popen.return_value = mock_proc

        proc = self.log_subprocess.p_open(["echo", "hello"])
        proc.check_return_code = MagicMock(side_effect=subprocess.CalledProcessError(
            1, ["echo", "hello"], stderr=b"error message"
        ))

        with pytest.raises(subprocess.CalledProcessError) as exc:
            proc.check_return_code()

        assert_that(exc.value.returncode).is_equal_to(1)
        assert_that(exc.value.stderr).is_equal_to(b"error message")
