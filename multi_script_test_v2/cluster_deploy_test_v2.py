
#Test Script

import socket
import sys



if __name__ == "__main__":

    d_file = open("data.txt","rb+")
    d = d_file.readline()

    if len(d) != 0:
	    rest_d = d_file.read()
	    d_file.seek(0)
	    d_file.truncate()
	    d_file.write(rest_d)
	    d_file.close()

    
    hostname = socket.gethostname()
    f = open(hostname+".txt","a")
    f.write("sysv "+sys.argv[1]+"\n")
    f.write("dfile "+str(d).strip("\r\n")+"\n")
    f.close()
