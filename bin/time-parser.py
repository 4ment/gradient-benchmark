#!/usr/bin/env python

import os
import re
import sys


def parse(file_name, pattern):
    iters = 0
    elbo = "nan"
    with open(file_name, "r") as fp:
        for line in fp:
            line = line.rstrip("\\n").rstrip("\\r")
            mt = pattern.match(line)
            if mt:
                iters, elbo = mt.groups()
    return str(iters), str(elbo)


pattern_time = re.compile(r"real\s+(?:(\d+)h)?(\d+)m(\d+\.?\d*)s")
with open(sys.argv[1], "w") as fpo:
    fpo.write("program,size,rep,time,iters,elbo\n")
    for file_path in sys.argv[2:]:
        a = file_path.rstrip(".log").split(".")
        a[0] = os.path.basename(a[0])
        if a[0] == "torchtree":
            pattern_elbo = re.compile(r"\s+(\d+)\s+(-\d+\.\d+).+")
        elif a[0] == "phylostan":
            #    100       -23766.035             1.000            1.000
            pattern_elbo = re.compile(r"\s+(\d+)\s+(-\d+\.\d+).+")
        elif a[0] == "physher":
            pattern_elbo = re.compile(r"(\d+)\s+ELBO:\s+(-\d+\.\d+).+")
        elif a[0] == "phylojax":
            pattern_elbo = re.compile(r"(\d+)\s+ELBO\s+(-\d+\.\d+).+")
        else:
            pattern_elbo = re.compile("sdfasdf")
        iters, elbo = parse(file_path.replace("log", "txt"), pattern_elbo)

        with open(file_path, "r") as fp:
            for line in fp:
                line = line.rstrip("\\n").rstrip("\\r")
                mt = pattern_time.match(line)
                if mt:
                    total_time = mt.groups()
                    time_sec = float(total_time[2]) + float(total_time[1]) * 60
                    if total_time[0] is not None:
                        time_sec += float(total_time[0]) * 3600

                    if a[0] == "torchtree":
                        if a[1] == "true":
                            a[0] = "bitorch"
                        del a[1]
                    a.extend([str(time_sec), iters, elbo])
                    fpo.write(",".join(a) + "\n")
