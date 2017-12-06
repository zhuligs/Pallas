import socket
import time
import os

serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
serversocket.bind(('localhost', 18089))
serversocket.listen(10) # become a server socket, maximum 5 connections

while True:
    connection, address = serversocket.accept()
    buf = connection.recv(1024)
    if len(buf) > 0:
        print buf
        print 'sleep'
        time.sleep(2)
        cmd = buf.strip()
        a = os.popen(cmd).read()
        connection.sendall(a)

        # break