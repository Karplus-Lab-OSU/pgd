from tasks import *
import time

class Test():
    def __init__(self):
        self.count = 0
        self.value = 0
        self.stop = 5

    def doTask(self, args):
        self.count = 0
        
        if not args == None:
            self.value = args

        while self.count < self.stop:
            time.sleep(1)
            self.count += 1
            self.value += 1
            print 'value: %d   progress: %d%%' % (self.value, self.taskProgress())

        return self.value

    def taskProgress(self):
        return 100*self.count/self.stop

class Singleton:
    __single = None
    def __init__( self ):
        if Singleton.__single:
            raise Singleton.__single
        Singleton.__single = self


class TestContainer():

    def __init__(self):
        self.test = Test()
        self.test2 = Test()
        self.task = TaskContainer('Example Container')
        self.task.addTask( Task('Sub task1', self.test.doTask, self.test.taskProgress) )
        self.task.addTask( Task('Sub task2', self.test2.doTask, self.test2.taskProgress) )






if __name__ == "__main__":

    test = Test()
    test1task = Task('Solo Task', test.doTask, test.taskProgress)
    test1task.start()

    while not test1task.finished:
        time.sleep(2)

    print "--------------------------------"

    test2 = TaskContainer('Test Container')
    test2.addTask( Task('Sub task1', test.doTask, test.taskProgress) )
    test2.addTask( Task('Sub task2', test.doTask, test.taskProgress) )
    test2.start()

    while not test2.finished:
       time.sleep(2)
