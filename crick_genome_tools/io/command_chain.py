"""
Helper class for a chain of commands, piping output between them.
"""

import subprocess

from crick_genome_tools.io.log_subprocess import LogSubprocess


class CommandChain:
    """
    Helper class for a chain of commands, piping output between them.
    """

    def __init__(self, commands, return_stream=False, output_file=None):
        """
        Initialize the CommandChain with a list of commands.

        Args:
            commands (list of list): Each inner list is a command (e.g., [["bwa", "mem", ...], ["samtools", "view", ...]]).
            return_stream (bool): Whether the final command should return a stream.
            output_file (str): If provided, the final output will be written to this file.
        """
        self.commands = commands
        self.return_stream = return_stream
        self.output_file = output_file
        self.log_subprocess = LogSubprocess()  # Shared LogSubprocess instance for managing subprocesses
        self.processes = []

    def run(self):
        """
        Run the chain of commands, piping output between them.
        Depending on options, either return a stream, write to a file, or complete without returning output.

        Returns:
            If return_stream is True, returns an iterator for the final output stream.
            Otherwise, returns None.
        """
        previous_stdout = None

        for idx, cmd in enumerate(self.commands):
            # Determine if this is the last command
            is_last = idx == len(self.commands) - 1

            if is_last:
                # Final command behavior
                if self.return_stream:
                    # Final command returns an output stream
                    proc = self.log_subprocess.p_open(cmd, stdin=previous_stdout, stdout=subprocess.PIPE)
                elif self.output_file:
                    # Final command writes directly to the output file
                    with open(self.output_file, "w", encoding="UTF-8") as f_out:
                        proc = self.log_subprocess.p_open(cmd, stdin=previous_stdout, stdout=f_out)
                else:
                    # Final command executes without output piping
                    proc = self.log_subprocess.p_open(cmd, stdin=previous_stdout)
            else:
                # Chain the current command, piping stdout to the next command
                proc = self.log_subprocess.p_open(cmd, stdin=previous_stdout, stdout=subprocess.PIPE)

            # Close the previous process's stdout to avoid deadlock
            if previous_stdout:
                previous_stdout.close()

            # Save the process for error checking
            self.processes.append(proc)

            # Set the current stdout to be used as input for the next command
            previous_stdout = proc.stdout

        # Check return codes and handle any errors
        for proc in self.processes:
            proc.check_return_code()

        # If return_stream is True, return an iterator for the final process's output
        if self.return_stream:
            return self.log_subprocess.stream_process(self.processes[-1])
        return None

    @staticmethod
    def command_to_file(command, output_file):
        """
        Run a single command and write the output to a file.

        Args:
            command (list): Command to run.
            output_file (str): File to write the output to.
        """
        with open(output_file, "w", encoding="UTF-8") as f_out:
            proc = LogSubprocess().p_open(command, stdout=f_out)
            proc.check_return_code()
