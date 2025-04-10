"""
Test cases for the CommandChain class in crick_genome_tools.io.command_chain.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import subprocess
from unittest.mock import MagicMock

import pytest
from assertpy import assert_that

from crick_genome_tools.io.command_chain import CommandChain


def test_command_chain_run_return_stream(monkeypatch):
    mock_log = MagicMock()
    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.check_return_code = MagicMock()

    mock_log.p_open.side_effect = [mock_proc, mock_proc]
    mock_log.stream_process.return_value = iter(["hello\n"])
    monkeypatch.setattr("crick_genome_tools.io.command_chain.LogSubprocess", lambda *_: mock_log)

    commands = [["echo", "hello"], ["grep", "hello"]]
    chain = CommandChain(commands, return_stream=True)
    result = list(chain.run())

    assert_that(result).is_equal_to(["hello\n"])
    assert_that(mock_log.p_open.call_count).is_equal_to(2)
    mock_proc.check_return_code.assert_called()


def test_command_chain_run_output_file(tmp_path):
    output_file = tmp_path / "output.txt"
    commands = [["echo", "hello"], ["grep", "hello"]]
    chain = CommandChain(commands, output_file=str(output_file))

    chain.run()

    assert_that(output_file.exists()).is_true()
    assert_that(output_file.read_text(encoding="utf-8")).is_equal_to("hello\n")


def test_command_chain_run_no_output(monkeypatch):
    mock_log = MagicMock()
    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.check_return_code = MagicMock()

    mock_log.p_open.side_effect = [mock_proc, mock_proc]
    monkeypatch.setattr("crick_genome_tools.io.command_chain.LogSubprocess", lambda *_: mock_log)

    commands = [["echo", "hello"], ["grep", "hello"]]
    chain = CommandChain(commands)

    result = chain.run()

    assert_that(result).is_none()
    assert_that(mock_log.p_open.call_count).is_equal_to(2)
    mock_proc.check_return_code.assert_called()


def test_command_chain_run_single_command(monkeypatch):
    mock_log = MagicMock()
    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.check_return_code = MagicMock()

    mock_log.p_open.side_effect = [mock_proc]
    mock_log.stream_process.return_value = iter(["hello\n"])
    monkeypatch.setattr("crick_genome_tools.io.command_chain.LogSubprocess", lambda *_: mock_log)

    commands = [["echo", "hello"]]
    chain = CommandChain(commands, return_stream=True)

    result = list(chain.run())

    assert_that(result).is_equal_to(["hello\n"])
    assert_that(mock_log.p_open.call_count).is_equal_to(1)
    mock_proc.check_return_code.assert_called()


def test_command_chain_run_with_error(monkeypatch):
    mock_log = MagicMock()
    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.check_return_code.side_effect = subprocess.CalledProcessError(1, "false")

    mock_log.p_open.side_effect = [mock_proc]
    monkeypatch.setattr("crick_genome_tools.io.command_chain.LogSubprocess", lambda *_: mock_log)

    commands = [["false"]]
    chain = CommandChain(commands)

    with pytest.raises(subprocess.CalledProcessError):
        chain.run()

    assert_that(mock_log.p_open.call_count).is_equal_to(1)
    mock_proc.check_return_code.assert_called()


def test_command_to_file_creates_output_file(tmp_path):
    output_file = tmp_path / "output.txt"
    command = ["echo", "hello"]

    CommandChain.command_to_file(command, str(output_file))

    assert_that(output_file.exists()).is_true()
    assert_that(output_file.read_text("utf-8").strip()).is_equal_to("hello")


def test_command_to_file_handles_subprocess_error(monkeypatch, tmp_path):
    command = ["false"]
    output_file = tmp_path / "output.txt"

    mock_log = MagicMock()
    mock_proc = MagicMock()
    mock_proc.check_return_code.side_effect = subprocess.CalledProcessError(1, command)
    mock_log.p_open.return_value = mock_proc
    monkeypatch.setattr("crick_genome_tools.io.command_chain.LogSubprocess", lambda *_: mock_log)

    with pytest.raises(subprocess.CalledProcessError):
        CommandChain.command_to_file(command, str(output_file))
