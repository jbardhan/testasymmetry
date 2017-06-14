import subprocess

command = "/usr/bin/nohup /usr/bin/python /home/arahimi/repos/testasymmetry/mobley/ParamSolvents.py"
try:
    ret = subprocess.check_call(command.split(" "))
except Exception as e:
    print e
    print "failed!"

print "done!"
print ret
