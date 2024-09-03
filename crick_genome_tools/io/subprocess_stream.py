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
        sub_proc = LogSubprocess()
        self.proc = subprocess.Popen(*args, **kwargs)

        if mode == "r":
            self.pipe = self.proc.stdout
        elif mode == "w":
            self.pipe = self.proc.stdin

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
        self.proc.wait()

    def __exit__(self, tp, val, tb):
        self.close()
        if tp is not None:
            log.error(f"Exception occurred: {val}")
            return False  # Propagate the exception
        return True
