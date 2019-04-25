
#Test Script

import socket
import sys



if __name__ == "__main__":

    hostname = socket.gethostname()
    f = open(hostname+".txt","a")
    d = open(sys.argv[1]+".txt","r")
    dline = d.readline()
    d.close()
    
    f.write("sysv "+sys.argv[1]+"\n")
    f.write("dfile "+str(dline).strip("\r\n")+"\n")
    f.close()
