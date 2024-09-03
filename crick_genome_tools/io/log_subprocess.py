# pylint: disable=C0116,C0114

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
        Initialise the LogSubprocess object
        """
        self.pdeathsig = self._child_preexec_set_pdeathsig()

    def _child_preexec_set_pdeathsig(self):
        """
        When used as the preexec_fn argument for subprocess.Popen etc,
        causes the subprocess to recieve SIGKILL if the parent process
        terminates.
        """
        if sys.platform.startswith("linux"):
            zero = ctypes.c_ulong(0)
            return LIBC.prctl(PR_SET_PDEATHSIG, ctypes.c_ulong(SIGKILL), zero, zero, zero)
        return None

    def _execute(self, method, *args, **kwargs):
        if "stderr" not in kwargs:
            kwargs["stderr"] = subprocess.PIPE
        if "preexec_fn" not in kwargs:
            kwargs["preexec_fn"] = self.pdeathsig

        # For Popen, return the process object directly
        if method == subprocess.Popen:
            proc = method(*args, **kwargs)
            returncode = proc.wait()
            if returncode != 0:
                stderr_output = proc.stderr.read() if proc.stderr else None
                error_msg = f"Subprocess terminated with exit code {returncode}"
                if stderr_output:
                    error_msg += f"\nError output:\n{stderr_output.decode()}"
                log.error(error_msg)
                raise subprocess.CalledProcessError(returncode, proc.args, output=None, stderr=stderr_output)

            return proc

        # For check_call, check_output, and call
        result = method(*args, **kwargs)
        returncode = result if method == subprocess.call else 0  # pylint: disable=comparison-with-callable

        stderr_output = None

        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, args[0], output=None, stderr=stderr_output)

        return result

    def check_call(self, *args, **kwargs):
        return self._execute(subprocess.check_call, *args, **kwargs)

    def check_output(self, *args, **kwargs):
        return self._execute(subprocess.check_output, *args, **kwargs)

    def call(self, *args, **kwargs):
        return self._execute(subprocess.call, *args, **kwargs)

    def p_open(self, *args, **kwargs):
        return self._execute(subprocess.Popen, *args, **kwargs)

    def stream_process(self, proc):
        """
        Generator to stream output from a subprocess and handle errors after completion.
        """
        if proc.stdout:
            for line in iter(proc.stdout.readline, b''):
                yield line
            proc.stdout.close()

        returncode = proc.wait()  # Now wait after the stream is done

        if returncode != 0:
            stderr_output = proc.stderr.read() if proc.stderr else None
            error_msg = f"Subprocess terminated with exit code {returncode}"
            if stderr_output:
                error_msg += f"\nError output:\n{stderr_output.decode()}"
            log.error(error_msg)
            raise subprocess.CalledProcessError(returncode, proc.args, output=None, stderr=stderr_output)
