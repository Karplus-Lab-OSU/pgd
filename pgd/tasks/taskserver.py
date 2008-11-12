#!/usr/bin/env python

import sys
from threading import Thread
import socket

from django.utils import simplejson

from pgd import tasks
from pgd.tasks.taskmanager import TaskManager

import time

TASK_MANAGER = None

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
        print "connecting..."
        self.port = port
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.connect(("localhost", port))
        self.connected = True
        print "connected"

    def sendCommand(self, command, args=None):
        response = None
        retry = 4

        if self.socket == None:
            print "NO SOCKET"
            return

        while not response and retry:
            try:
                message = {'command':command}
                if args <> None:
                    message['args'] = args

                self.socket.send(simplejson.dumps(message))
                response = self.socket.recv(2048)
                if not response:
                    raise RuntimeError, "No response from server"

            except (RuntimeError, socket.error):
                errno, errstr = sys.exc_info()[:2]
                retry = retry - 1

                try:
                    self.connect(self.port)
                except:
                    #silently discard error if theres retries left
                    if retry:
                        #sleep for a few seconds to give the server a chance
                        time.sleep(2)
                    else:
                        raise RuntimeError, "No response from task server, is it running?"


        # if the response is still '' then the server could not be contacted
        # raise an error so the calling function can deal with it properly
        if not response:
            raise Exception

        return response

    def listTasks(self):
        return self.sendCommand(CMD_LIST)

    def progress(self):
        return self.sendCommand(CMD_PROGRESS)

    def start(self, key):
        return self.sendCommand(CMD_START, {'key':key})


class ClientThread(Thread):
    def __init__(self,clientSocket):
        Thread.__init__(self)
        self.sockfd = clientSocket #socket client
        self.sockfd.settimeout(1)
        self.running = False
        self.started = False

    def __delete__(self):
        if not self.socket == None:
            self.socket.close()

    def run(self):
        self.running = True
        self.started = True
        while self.running:
            try:
                buff = self.sockfd.recv(2048)

                # empty buffer can only happen when other side of socket closes
                # break out of the loop so the socket can be cleaned up
                if not buff:
                    break

                resp = self.processCommand(buff)
                self.sockfd.send(resp)

            except socket.error:
                errno, errstr = sys.exc_info()[:2]
                #silently discard timeout
                if errno == socket.timeout:
                    pass
                else:
                    raise socket.error

        self.running = False

        #close socket
        self.sockfd.close()
        print '         * client stopped'


    def processCommand(self, buff):
        message = simplejson.loads(buff)
        command = message['command']

        if command == CMD_LIST:
            return simplejson.dumps(TASK_MANAGER.listTasks())

        elif command == CMD_PROGRESS:
            return simplejson.dumps(TASK_MANAGER.progress())

        elif command == CMD_START:
            return TASK_MANAGER.start(message['args']['key'])

        elif command == CMD_STOP:
            return TASK_MANAGER.stop()


"""
thread that waits for connections from clients
"""
class ServerThread(Thread):
    def __init__(self):
        Thread.__init__(self)
        self.running = False
        self.clients = []
        self.count = 0

    def run(self):
        self.running = True
        serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        serversocket.setsockopt( socket.SOL_SOCKET, socket.SO_REUSEADDR, 1 )
        try:
            serversocket.bind(('localhost',PORT))
        except ValueError,e:
            raise e

        print "bound port: ", PORT
        
        #set the timeout to 1 second.  this will cause alot of timeouts
        # while waiting for clients to connect but its needed so that
        # the thread can check for a stop flag
        serversocket.settimeout(1)

        serversocket.listen(5)

        while self.running:
            # look for dead threads to purge every 60 seconds
            # dont do it every iteration to save processing
            self.count+=1
            if self.count == 60:
                self.count = 0
                for client in self.clients:
                    if client.started and not client.running:
                        #purge dead client
                        self.clients.remove(client)


            # wait for a client to connect
            try:
                (clientSocket, address) = serversocket.accept()


                client = ClientThread(clientSocket)
                CLIENTS.append(client)
                client.start()
                self.clients.append(client)
                print 'Client connected: %s   Total Clients: %d' % (address, len(self.clients))


            # catch timeouts but silently discard them
            except socket.error:
                errno, errstr = sys.exc_info()[:2]
                #silently discard timeout
                if errno == socket.timeout:
                    pass
                else:
                    raise socket.error

        #stop incoming connections
        serversocket.close()
        print '    * No longer accepting client connections'

        #shutdown all client threads
        print '    * Stopping all clients'
        for clientthread in self.clients:
            client.running = False

if __name__ == '__main__':
    import sys
    import time

    # autodiscover registered tasks
    tasks.autodiscover()

    # echo loaded tasks
    TASK_MANAGER = tasks.TASK_MANAGER
    tasks = TASK_MANAGER.listTasks()
    print "tasks: ", tasks

    # Start the server listening for connections
    server = ServerThread()
    server.start()

    time.sleep(1)

    # Listen for keyboard input to stop the server
    try:
        print ''
        print "Press ^C to stop"
        while True:
            time.sleep(300)
    except KeyboardInterrupt:
        pass
    
    server.running = False
    print ''
    print 'Stopping Task Server...'

