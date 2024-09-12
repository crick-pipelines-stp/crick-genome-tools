# pylint: disable=C0116,C0114

import logging
import os
import subprocess
import sys

from crick_genome_tools.io.log_subprocess import LogSubprocess


log = logging.getLogger(__name__)

class SubprocessStream:
    """
    Wrap a subprocess that we stream from or stream to. Acts like an open filehandle by passing down
    next, fileno, write, and close down to its pipe.
    """

    def __init__(self, *args, **kwargs):
        mode = kwargs.pop("mode", "r")
        if mode == "r":
            kwargs["stdout"] = subprocess.PIPE
        elif mode == "w":
            kwargs["stdin"] = subprocess.PIPE
        else:
            raise ValueError(f"mode {mode} unsupported")

        kwargs["preexec_fn"] = os.setsid
        sys.stdout.flush()

        # Use LogSubprocess.p_open to ensure proper error handling
        sub_proc = LogSubprocess()
        self.proc = sub_proc.p_open(*args, **kwargs)

        if mode == "r":
            self.pipe = self.proc.stdout
        elif mode == "w":
            self.pipe = self.proc.stdin

        # Use stream_process generator from LogSubprocess to stream the output
        self.stream_generator = sub_proc.stream_process(self.proc)

    def __enter__(self):
        return self

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.stream_generator)
        if line:
            return line
        self.close()
        raise StopIteration

    def fileno(self):
        return self.pipe.fileno()

    def write(self, x):
        self.pipe.write(x)
        self.pipe.flush()

    def close(self):
        if not self.pipe.closed:
            self.pipe.close()

        # Use check_return_code to handle errors after the process is closed
        try:
            self.proc.check_return_code()
        except subprocess.CalledProcessError as e:
            log.error(f"Subprocess failed: {e}")
            raise

    def __exit__(self, tp, val, tb):
        self.close()
        if tp is not None:
            log.error(f"Exception occurred: {val}")
            return False  # Propagate the exception
        return True
