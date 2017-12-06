#!/usr/bin/python 
import socket
import sys

ii = sys.argv[1]

clientsocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
clientsocket.connect(('localhost', 18089))
clientsocket.send('echo ' + ii + '\n')
clientsocket.send('ls')
print clientsocket.recv(1024)
clientsocket.close()