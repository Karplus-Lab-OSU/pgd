from threading import Thread
from django.utils import simplejson

class SubTaskWrapper():
    def __init__(self, task, percentage):
        self.task = task
        self.percentage = percentage

class WorkUnit(Thread):
    def __init__(self, parent, args):
        Thread.__init__(self)
        self.parent = parent
        self.args = args

    def run(self):
        self.parent.running = 1
        self.parent.finished = 0
        result = self.doTask(self.args)
        self.parent.running = 0
        self.parent.finished = 1
        return result

    def doTask(self,args):
        self.parent.workFunction(args)

class Task():
    def __init__(self, msg, workFunction, progressFunction):
        self.msg = msg
        self.workFunction = workFunction
        self.progressFunction = progressFunction
        self.workunit = None
        self.running = 0
        self.finished = 0

    def reset(self):
        self.running=0
        self.finished=0

    def start(self, args=None):
        # only start if not already running
        if self.running:
            return

        # create a new workunit if needed
        if self.workunit == None or self.finished:
            self.workunit = WorkUnit(self, args)

        self.workunit.start()

    def doTask(self, args):
        self.running = 1
        self.finished = 0
        result = self.workFunction(args)
        self.running = 0
        self.finished = 1
        return result

    def progress(self):
        if self.progressFunction == None:
            return 0
        else:
            return self.progressFunction()

# TaskContainer contains other tasks
class TaskContainer(Task):

    def __init__(self, msg, sequential=True):
        Task.__init__(self, msg, self.iterateTasks, None)
        self.subtasks = []
        self.sequential = sequential

    def addTask(self, task, percentage=None):
        subtask = SubTaskWrapper(task, percentage)
        self.subtasks.append(subtask)

    def reset(self):
        for subtask in self.subtasks:
            subtask.task.reset()

    # Starts the task running all subtasks
    def iterateTasks(self, args=None):
        self.reset()

        result = None
        for subtask in self.subtasks:
            if self.sequential:
                #sequential task, run the task work directly (default)
                result = subtask.task.doTask(result)
            else:
                #parallel task, run the subtask in its own thread
                result = subtask.task.doTask(result)

        return result

    def calculatePercentage(self):
        return float(1)/len(self.subtasks);

    def progress(self):
        progress = 0
        for subtask in self.subtasks:
            if subtask.percentage == None:
                percentage = self.calculatePercentage()
            else:
                percentage = subtask.percentage

            # if task is done it complete 100% of its work 
            if subtask.task.finished:
                progress += 100*percentage
            # task is only partially complete
            else:
                progress += subtask.task.progress()*percentage

        return progress


class TaskManager():

    def __init__(self, task):
        self.task = task

    def processTask(self, task, tasklist=None, level=0, count=0):
        # initial call wont have an area yet
        if tasklist==None:
            tasklist = []

        #turn the task into a tuple
        processedTask = [count, task.msg, level, False]

        #add that task to the list
        tasklist.append(processedTask)

        #add all children if the task is a container
        if isinstance(task,TaskContainer):
            processedTask[3] = True
            for subtask in task.subtasks:
                count+=1
                self.processTask(subtask.task, tasklist, level+1, count)

        return tasklist

    def processTaskProgress(self, task, tasklist=None, count=0):
        # initial call wont have an area yet
        if tasklist==None:
            tasklist = []

        #turn the task into a tuple
        processedTask = [count, task.progress()]

        #add that task to the list
        tasklist.append(processedTask)

        #add all children if the task is a container
        if isinstance(task,TaskContainer):
            for subtask in task.subtasks:
                count+=1
                self.processTaskProgress(subtask.task, tasklist, count)

        return tasklist

    def listTasks(self):
        processedTasks = self.processTask(self.task)
        return simplejson.dumps(processedTasks)

    def progress(self):
        progress = self.processTaskProgress(self.task)
        return simplejson.dumps(progress)

    def start(self):
        self.task.start()
        return '1'

    def stop(self):
        pass