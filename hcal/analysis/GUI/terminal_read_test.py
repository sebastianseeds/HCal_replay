import subprocess
import sys

result = subprocess.run(["ls", "-l"], stdout=subprocess.PIPE)
#print(result.stdout)
print(result.stdout.decode('utf-8'))
#subprocess.run(["ls", "-l"], capture_output=True)
#subprocess.run(["ls", "-l", "/home/daq/test_fadc"], capture_output=True)
