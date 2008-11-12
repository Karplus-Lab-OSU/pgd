class TaskManager:
   def __init__(self):
      self.tasks=[]

   def __del__(self):
      #stop all tasks
      pass

#

global taskManager
taskManager = TaskManager()
