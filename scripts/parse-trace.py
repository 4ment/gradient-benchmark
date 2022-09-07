import sys
from os import walk

if len(sys.argv) == 1:
    print(f"Usage: python {sys.argv[0]} work/ trace.txt")
    exit(2)

dic = {}
for (dirpath, dirnames, filenames) in walk(sys.argv[1]):
    if len(filenames) > 0:
        dic[dirpath.lstrip("work/")[:9]] = filenames

with open(sys.argv[2], "r") as fp:
    for line in fp:
        a = line.split("\t")
        if 'FAILED' in a or 'ABORTED' in a:
            continue
        if "task_id" == a[0]:
            print("program\tsize\treplicate\t" + line, end="")
        if a[1] not in dic:
            continue
        elif "macro_flu:RUN_" in a[3]:
            for f in dic[a[1]]:
                if ".txt" in f:
                    b = f.split(".")[:-1]
                    mem = a[10].split(" ")
                    if mem[1] == "GB":
                        a[10] = str(float(mem[0]) * 1024)
                    else:
                        a[10] = mem[0]
                    mem = a[11].split(" ")
                    if mem[1] == "GB":
                        a[11] = str(float(mem[0]) * 1024)
                    else:
                        a[11] = mem[0]

                    if len(b) == 4:
                        b[0] = "bitorch" if b[1] == "true" else "torchtree"
                        del b[1]
                        print("\t".join(b) + "\t" + "\t".join(a), end="")
                    else:
                        print("\t".join(b) + "\t" + "\t".join(a), end="")
                    break
