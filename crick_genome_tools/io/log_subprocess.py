import ctypes
import ctypes.util
import logging
import subprocess
import sys
from signal import SIGKILL


log = logging.getLogger(__name__)

LIBC = ctypes.CDLL(ctypes.util.find_library("c"))
PR_SET_PDEATHSIG = ctypes.c_int(1)  # <sys/prctl.h>


class LogSubprocess:
    """
    Helper class with functions for logging calls to subprocesses.
    """

    def __init__(self):
        """
        Initialise the LogSubprocess object.
        """
        self.pdeathsig = self._child_preexec_set_pdeathsig()

    def _child_preexec_set_pdeathsig(self):
        """
        When used as the preexec_fn argument for subprocess.Popen,
        causes the subprocess to receive SIGKILL if the parent process terminates.
        """
        if sys.platform.startswith("linux"):
            zero = ctypes.c_ulong(0)
            return LIBC.prctl(PR_SET_PDEATHSIG, ctypes.c_ulong(SIGKILL), zero, zero, zero)
        return None

    def _execute(self, method, *args, **kwargs):
        """
        General method to handle subprocess calls except Popen (which is deferred).
        This handles errors immediately after the command completes.
        """
        if "stderr" not in kwargs:
            kwargs["stderr"] = subprocess.PIPE
        if "preexec_fn" not in kwargs:
            kwargs["preexec_fn"] = self.pdeathsig

        # For methods that immediately return (like check_call, check_output)
        try:
            result = method(*args, **kwargs)
            return result  # For check_call and check_output, this is either 0 or output
        except subprocess.CalledProcessError as e:
            stderr_output = e.stderr.decode() if e.stderr else "Unknown error"
            log.error(f"Subprocess failed with exit code {e.returncode}.\nError output:\n{stderr_output}")
            raise

    def check_call(self, *args, **kwargs):  # pylint: disable=missing-function-docstring
        return self._execute(subprocess.check_call, *args, **kwargs)

    def check_output(self, *args, **kwargs):  # pylint: disable=missing-function-docstring
        return self._execute(subprocess.check_output, *args, **kwargs)

    def call(self, *args, **kwargs):  # pylint: disable=missing-function-docstring
        return self._execute(subprocess.call, *args, **kwargs)

    def p_open(self, *args, **kwargs):
        """
        Start a subprocess with deferred error handling. Returns the process object immediately,
        and errors should be checked after the process completes.
        """
        if "stderr" not in kwargs:
            kwargs["stderr"] = subprocess.PIPE  # Capture stderr for error handling
        if "stdout" not in kwargs:
            kwargs["stdout"] = subprocess.PIPE  # Capture stdout if not already handled
        if "preexec_fn" not in kwargs:
            kwargs["preexec_fn"] = self.pdeathsig  # Set preexec_fn to send SIGKILL

        # Start the subprocess and return the process object
        proc = subprocess.Popen(*args, **kwargs)

        def check_return_code():
            """
            Check the return code and handle any errors after the process completes.
            """
            returncode = proc.wait()  # Wait for the process to finish
            if returncode != 0:
                stderr_output = proc.stderr.read().decode() if proc.stderr else "Unknown error"
                error_msg = f"Subprocess terminated with exit code {returncode}\nError output:\n{stderr_output}"
                log.error(error_msg)
                raise subprocess.CalledProcessError(returncode, proc.args, output=None, stderr=stderr_output)

        # Attach error-checking function to the process object
        proc.check_return_code = check_return_code
        return proc

    def stream_process(self, proc):
        """
        Generator to stream output from a subprocess and handle errors after completion.
        """
        if proc.stdout:
            for line in iter(proc.stdout.readline, b""):
                yield line
            proc.stdout.close()

        proc.check_return_code()  # Ensure errors are handled after streaming
