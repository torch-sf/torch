"""
torch_monitor
=============

Holds classes TorchMonitor and MonitorWorkers, which are used to limit the time code
is allowed to execute.

TorchMonitor
    An instance of TorchMonitor can be used to wrap code in start and stop statements.
    If the code exceeds the time specified by the TorchMonitor, it will print a failure 
    message and ask slurm to abort the job.

    Using try, finally statments around monitors is not necessary but highly recomended
    because it ensures that a monitor is stopped even if an expection is raised.

MonitorWorker
    Wrapper that monitors function calls to an underlying AMUSE worker, allowing
    worker calls to be monitored without modifying the worker itself.
    This is achieved by overriding two special Python methods:

    __getattr__
        Overridden to be called when an attribute does not exist on the MonitorWorker 
        instance, instead retrieving the attribute from the underlying worker.
        Non-callable attributes are returned directly, while callable attributes 
        are wrapped in a function that starts the TorchMonitor before execution 
        and stops it afterward, even if an exception is raised. The wrapped function 
        retains its original metadata and additional public attributes.

    __setattr__ 
        Overridden to forward attribute assignment to the underlying worker unless 
        they are _worker and _monitor, which are stored directly on the wrapper using
        object.__setattr__ to avoid recursive calls to __setattr__.

    Note that __getattr__ only intercepts attributes that are not already defined on
    the wrapper, whereas __setattr__ intercepts all attribute assignments.
"""
import os
import threading
import time
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

    abort_delay = 10 # [seconds] Delay abort to allow other monitors to report failures.

    def __init__(self, timeout_minutes=30):
        """ Initialize monitor for timeout_minutes minutes.
        """
        self.timeout = timeout_minutes*60 # seconds
        self.timer = None

    def abort(self, reason):
        """ Abort function triggered after self.timeout seconds.
        """
        tprint(f"Failure: {reason}", flush=True)
        time.sleep(self.abort_delay)

        jobid = os.environ.get("SLURM_JOB_ID")
        if jobid:
            subprocess.run(["scancel", jobid], timeout=10)

        os._exit(1)

    def start(self, abort_message):
        """ Start monitor thread.
        """
        if self.timer is not None:
            raise RuntimeError("Cannot start monitor while it is already running.")

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
            self._monitor.start(f"Worker function call timed out in {name}")
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
