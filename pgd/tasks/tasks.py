from threading import Thread
from abstract import *

STATUS_FAILED = -1;
STATUS_STOPPED = 0;
STATUS_RUNNING = 1;
STATUS_PAUSED = 2;
STATUS_COMPLETE = 3;

"""
SubTaskWrapper - class used to store additional information
about the relationship between a container and a subtask.

    percentage - the percentage of work the task accounts for.
"""
class SubTaskWrapper():
    def __init__(self, task, percentage):
        self.task = task
        self.percentage = percentage

"""
WorkUnit - class used to wrap task work in a thread.  New WorkUnits
can be constructed every time a task is run.  This allows the task to run in its own thread
and allow it to be run multiple times.
"""
class WorkUnit(Thread):
    def __init__(self, parent, args):
        Thread.__init__(self)
        self.parent = parent
        self.args = args

    """
    runs the task by calling its work function
    """
    def run(self):
        #self.parent._status = STATUS_RUNNING
        return self.parent.work(self.args)
        #self.parent._status = STATUS_COMPLETE

        #return result


"""
Task - class that wraps a set of functions as a runnable unit of work.  Once
wrapped the task allows functions to be managed and tracked in a uniform way.

This is an abstract class and requires the following functions to be implemented:
    * _work  -  does the work of the task
    * _reset - resets the state of the task when stopped
    * progress - returns the state of the task as an integer from 0 to 100
    * progressMessage - returns the state of the task as a readable string
"""
class Task():
    def __init__(self, msg):
        self.msg = msg

        _work = AbstractMethod('_work')
        progress = AbstractMethod('progress')
        progressMessage = AbstractMethod('progressMessage')
        _reset = AbstractMethod('_reset')

        self.workunit = None
        self._status = STATUS_STOPPED
        self.id = 1

    """
    Resets the task, including setting the flags properly,  delegates implementation
    specific work to _reset()
    """
    def reset(self):
        self._status = STATUS_STOPPED
        self._reset()

    """
    starts the task.  This will spawn the work in a workunit thread.
    """
    def start(self, args=None):
        # only start if not already running
        if self._status == STATUS_RUNNING:
            return

        # create a new workunit if needed
        if self.workunit == None or self._status == STATUS_COMPLETE:
            self.workunit = WorkUnit(self, args)

        self.workunit.start()

    """
    Does the work of the task.  This is can be called directly for synchronous work or via start which
    causes a workunit thread to be spawned and call this function.  this method will set flags properly and
    delegate implementation specific work to _work(args)
    """
    def work(self, args):
        self._status = STATUS_RUNNING
        result = self._work(args)
        self._status = STATUS_COMPLETE
        return result

    def status(self):
        return self._status

"""
TaskContainer - an extension of Task that contains other tasks

TaskContainer does no work itself.  Its purpose is to allow a bigger job
to be broken into discrete functions.  IE.  downloading and processing.
"""
class TaskContainer(Task):

    def __init__(self, msg, sequential=True):
        Task.__init__(self, msg)
        self.subtasks = []
        self.sequential = sequential

    def addTask(self, task, percentage=None):
        subtask = SubTaskWrapper(task, percentage)
        self.subtasks.append(subtask)
        task.id = '%s-%d' % (self.id,len(self.subtasks))

    def reset(self):
        for subtask in self.subtasks:
            subtask.task.reset()

    # Starts the task running all subtasks
    def _work(self, args=None):
        self.reset()

        result = args
        for subtask in self.subtasks:
            if self.sequential:
                #sequential task, run the task work directly (default)
                result = subtask.task.work(result)
            else:
                #parallel task, run the subtask in its own thread
                result = subtask.task.start(result)

        return result


    """
    calculatePercentage - determines the percentage of work that each
    child task accounts for.

    TODO: take into account tasks that have had weighting manually set. 
    """
    def calculatePercentage(self):
        return float(1)/len(self.subtasks);


    """
    progress - returns the progress as a number 0-100.  

    A container task's progress is a derivitive of its children.
    the progress of the child counts for a certain percentage of the 
    progress of the parent.  This weighting can be set manually or
    divided evenly by calculatePercentage()
    """
    def progress(self):
        progress = 0
        for subtask in self.subtasks:
            if subtask.percentage == None:
                percentage = self.calculatePercentage()
            else:
                percentage = subtask.percentage

            # if task is done it complete 100% of its work 
            if subtask.task._status == STATUS_COMPLETE:
                progress += 100*percentage
            # task is only partially complete
            else:
                progress += subtask.task.progress()*percentage

        return progress

    """ 
    returns a plain text status message
    """
    def progressMessage(self):
        for subtask in self.subtasks:
            if subtask.task._status == STATUS_RUNNING:
                return subtask.task.progressMessage()

        return None

    """
    getStatus - returns status of this task.  A container task's status is 
    a derivitive of its children.

    failed - if any children failed, then the task failed
    running - if any children are running then the task is running
    paused - paused if no other children are running
    complete - complete if all children are complete
    stopped - default response if no other conditions are met
    """
    def status(self):
        has_paused = False;
        has_unfinished = False;

        for subtask in self.subtasks:
            subtaskStatus = subtask.task.status()
            if subtaskStatus == STATUS_RUNNING:
                # we can return right here because if any child is running the 
                # container is considered to be running
                return STATUS_RUNNING

            elif subtaskStatus == STATUS_FAILED:
                # we can return right here because if any child failed then the
                # container is considered to be failed.  All other running tasks
                # should be stopped on failure.
                return STATUS_FAILED

            elif subtaskStatus == STATUS_PAUSED:
                has_paused = True

            elif subtaskStatus <> STATUS_COMPLETE:
                has_unfinished = True

        # Task is not running or failed.  If any are paused
        # then the container is paused
        if has_paused:
            return STATUS_PAUSED

        # task is not running, failed, or paused.  if all children are complete then it is complete
        if not has_unfinished:
            return STATUS_COMPLETE

        # only status left it could be is STOPPED
        return STATUS_STOPPED



