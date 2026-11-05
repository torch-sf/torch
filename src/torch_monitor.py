import os
import time
import threading
import subprocess
from functools import wraps

from torch_stdout import tprint

class TorchMonitor:
    """ 
    Creates a separate python thread that monitors function calls
    and kills the slurm job if not completed after some time.

    Usage:

    monitor = TorchMonitor(timeout_minutes=10)
    monitor.start(abort_message="AMUSE request timed out.")

    try:
        # MPI blocking call
        pool.wait()
    finally:
        monitor.stop()
    """

    def __init__(self, timeout_minutes=30):
        """ Initialize monitor for timeout_minutes minutes.
        """
        self.timeout = timeout_minutes*60 # seconds
        self.timer = None

    def abort(self, reason):
        """ Abort function triggered after self.timeout seconds.
        """
        tprint(f"Failure: {reason}", flush=True)

        jobid = os.environ.get("SLURM_JOB_ID")
        if jobid:
            subprocess.run(["scancel", jobid], timeout=10)

        os._exit(1)

    def start(self, abort_message):
        """ Start monitor thread.
        """
        self.timer = threading.Timer(self.timeout, self.abort, 
            args=(abort_message,))
        self.timer.daemon = True
        self.timer.start()

    def stop(self):
        """ Stop monitor thread (called if monitored function succeeds.)
        """
        if self.timer is not None:
            self.timer.cancel()
            self.timer = None

class MonitorWorker:
    """
    Wraps worker function calls in a TorchMonitor.
    Will only capture direct function calls, even if those functions have
    attributes. For example, calls like hydro.evolve_model.asynchronous(*args) 
    are preserved if discovered by dir(), but they will not be monitored.
    A call to hydro.evolve_model() will be monitored.
    Note that _worker and _monitored are protected from modifiction, i.e.,
    a worker cannot have these attributes.
    """
    def __init__(self, worker, monitor):
        """ Initialize the wrapper
        """
        if any(arg is None for arg in (worker, monitor)):
            raise ValueError("MonitorWorker requires a valid worker and monitor")
        
        object.__setattr__(self, "_worker", worker)
        object.__setattr__(self, "_monitor", monitor)

    def __setattr__(self, name, value):
        """
        Attribute assignments to the underlying worker with protection 
        for _worker and _monitor.
        """
        if name in ("_worker", "_monitor"):
            object.__setattr__(self, name, value)
        else:
            setattr(self._worker, name, value)

    def __getattr__(self, name):
        """ 
        Wraps any function call in a TorchMonitor. Ignores attributes that
        are not callable, e.g., variables and subclasses.
        """
        attr = getattr(self._worker, name)

        if not callable(attr):
            return attr

        @wraps(attr)
        def monitor_call(*args, **kwargs):
            """ Wrap function in TorchMonitor start() and stop().
            """
            self._monitor.start(f"AMUSE call timed out in {name}")
            try:
                return attr(*args, **kwargs)
            finally:
                self._monitor.stop()
        
        # Handle attributes assigned to callable attributes.
        for key in dir(attr):
            if not key.startswith("_") and not hasattr(monitor_call, key):
                try:
                    setattr(monitor_call, key, getattr(attr, key))
                except (AttributeError, TypeError):
                    pass
        
        return monitor_call
