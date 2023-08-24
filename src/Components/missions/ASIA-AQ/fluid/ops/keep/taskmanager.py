"""Defines classes and methods for managing multiple tasks in parallel.

**See Also:**
    Modules: subprocess, multiprocessing
"""

from __future__ import print_function
# should move print statements outside of module

import time
import shlex
import subprocess
import multiprocessing

class TaskManager(object):
    """Provides methods for managing multiple tasks in parallel.

    **Args:**
        none

    **Raises:**
        none
    """

    def __init__(self, ntask=4): # bad default

        assert ntask > 0, "Bad request for number of processors"
        self.max_task = min(ntask, multiprocessing.cpu_count())
        print('Using ' + str(self.max_task) + ' processors.')
        self.task = {}
        self.command = {}

    def spawn(self, command):
        """Executes the specified command as a new task.

        **Args:**
            command : string : input
                Command to be executed.

        **Returns:**
            none

        **Raises:**
            none

        **Notes:**
            1. This method will wait for a processor to become available. There
               is currently no option to specify the wait time or query
               interval.
        """

        while self.load() >= self.max_task:
            time.sleep(1)

        cmd = shlex.split(command)

        p = subprocess.Popen(
            cmd, stdout=None, stderr=None, shell=False
        )
        process = str(p.pid)
        self.task[process] = p
        self.command[process] = command

    def load(self):
        """Determines the number of tasks (load) currently executing.

        **Args:**
            none

        **Returns:**
            ntasks : integer
                Number of tasks currently executing.

        **Raises:**
            none

        **Notes:**
            none
        """

        ntask = len(self.task)
        processes = self.task.keys()

        for process in processes:
            p = self.task[process]

            if p.poll() == None:
                continue

            if p.returncode != 0:
                print(str(self.command[process]) + 'ended with rc=' + str(p.returncode))
            del self.task[process]
            del self.command[process]

            ntask -= 1

        return ntask

    def wait(self):
        """Waits for all tasks to complete.

        **Args:**
            none

        **Returns:**
            none

        **Raises:**
            none
        """

        while self.load() > 0:
            time.sleep(1)
