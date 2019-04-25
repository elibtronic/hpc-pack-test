
#Test Script

import socket
import sys



if __name__ == "__main__":

    hostname = socket.gethostname()
    f = open("run_"+sys.argv[1]+"_"+hostname+".txt","w")
    f.write(sys.argv[1])
    f.close()
