
#Test Script

import socket
import sys

hostname = socket.gethostname()
f = open(hostname+".txt","w")
f.write("hello world")
f.close()
