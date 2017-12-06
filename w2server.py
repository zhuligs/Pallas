#!/usr/bin/env python 

## socket server running on login node

import socket
import os
import time
import itin
import sdata


def startserver():
    HOST = ''
    PORT = itin.serverport
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    print ('* ZLOG: Socket created')
    binded = False
    while binded is False:
        try:
            s.bind((HOST, PORT))
            binded = True
            print ('* ZLOG: Socket bind complete. PORT:', PORT)
            f = open('PORT.txt', 'w')
            f.write("%d\n" % PORT)
            f.close()
        except socket.error, msg:
            PORT += 1
            print ('* ZLOG: Bind failed. Error Code : ' + str(msg[0])
                   + ' Message ' + msg[1])

    s.listen(10)

    while True:
        connection, address = s.accept()
        buff = connection.recv(2048)
        print ('* ZLOG: receiving message:', buff)
        if len(buff) > 0:
            buffs = buff.split("%")
            if buffs[0] == 'subjob':
                print ('* ZLOG: subjob')
                ddir = buffs[1]
                jbuff = os.popen('cd ' + ddir + '; qsub pbs.sh').read()
                print ('* ZLOG: sub job in', ddir)
                jid = jbuff.strip() 
                connection.sendall(jid)
                print ('* ZLOG: return job id', jid)

            if buffs[0] == 'qdeljob':
                print ('* ZLOG: qdel job')
                allid = buffs[1]
                os.system('qdel ' + allid)
                connection.sendall('job canceled')

            if buffs[0] == 'checkjob':
                print ('* ZLOG: check job')
                ## return DONE.jobid1.jobid2. 
                while True:
                    jstat = os.system("squeue > qbuff")
                    if jstat == 0:
                        break
                    else:
                        time.sleep(3)
                jbuff = []
                with open("qbuff") as f:
                    for line in f:
                        jbuff.append(line)
                jobdone = 'NO'
                if len(jbuff) > 0:
                    reminds = []
                    try:
                        for x in jbuff[1:]:
                            reminds.append(x.split()[0])
                    except:
                        print ('* ZLOG: except 1')
                        jobdone = 'DONE'
                else:
                    jobdone = 'DONE'
                jobrem = [jobdone] + reminds
                jobstring = '.'.join(jobrem)
                connection.sendall(jobstring)
                print ('* ZLOG: return job id:', jobstring)


def w2server():
    startserver()


if __name__ == '__main__':
    w2server()
