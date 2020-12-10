import subprocess
import glob

summ = 0
fs = glob.glob('/mnt/data/cliveome/*.fastq')
for f in fs:
    s = subprocess.check_output('wc -l ' + f, shell=True)
    summ += int(s.decode().split()[0])/4
print(summ)
