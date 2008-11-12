#!/usr/bin/env python

import sys
from threading import Thread
import socket

from tasktests import *

TASK_MANAGER = TestContainer()

CLIENTS = []
PORT = 5800

CMD_LIST = '0'
CMD_PROGRESS = '1'
CMD_START = '2'
CMD_STOP = '3'

class TaskClient():

    def __init__(self, port=None):
        self.port = port
        self.socket = None
        self.connected = False

    def connect(self, port):
        print "connecting"
        self.port = port
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.connect(("localhost", port))
        self.connected = True
        print "connected"

    def sendCommand(self, command):
        if self.socket == None:
            print "NO SOCKET"
            return

        self.socket.send(command)
        response = self.socket.recv(2048)
        return response

    def listTasks(self):
        return self.sendCommand(CMD_LIST)

    def progress(self):
        return self.sendCommand(CMD_PROGRESS)

    def start(self):
        return self.sendCommand(CMD_START)


class ClientThread(Thread):
    def __init__(self,clientSocket):
        Thread.__init__(self)
        self.sockfd = clientSocket #socket client

    def __delete__(self):
        if not self.socket == None:
            self.socket.close()

    def run(self):
        while True:
            buff = self.sockfd.recv(2048)
            print "Command:  ", buff
            resp = self.processCommand(buff)
            print "Response: ", resp
            self.sockfd.send(resp)

    def processCommand(self, buff):
        key = buff

        if key == CMD_LIST:
            return TASK_MANAGER.listTasks()

        elif key == CMD_PROGRESS:
            return TASK_MANAGER.progress()

        elif key == CMD_START:
            return TASK_MANAGER.start()

        elif key == CMD_STOP:
            return TASK_MANAGER.stop()

if __name__ == '__main__':
    
    task = TestContainer()
    TASK_MANAGER = TaskManager(task.task)

    print "tasks: ", TASK_MANAGER.listTasks()
    print "progress: ", TASK_MANAGER.progress()

    serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    serversocket.setsockopt( socket.SOL_SOCKET, socket.SO_REUSEADDR, 1 )
    try:
        serversocket.bind(('localhost',PORT))
    except ValueError,e:
        print e  

    print "bound port: ", PORT

    serversocket.listen(5)
    while True:
        (clientSocket, address) = serversocket.accept()
        client = ClientThread(clientSocket)
        CLIENTS.append(client)
        client.start()