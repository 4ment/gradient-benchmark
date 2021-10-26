#!/usr/bin/env python

import re
import sys

regex = re.compile(r".+ ([\d\.]+) ms.*")


def parse_physher(fp):
    for line in fp:
        if line.startswith("Tree likelihood"):
            line = next(fp)
            logp_like = float(regex.match(line).group(1)) / 1000.0
            line = next(fp)
            grad_logp_like = float(regex.match(line).group(1)) / 1000.0
            print(f"treelikelihood,evaluation,off,{logp_like}")
            print(f"treelikelihood,gradient,off,{grad_logp_like}")
        elif line.startswith("Height transform log det Jacobian"):
            line = next(fp)
            logp_jac = float(regex.match(line).group(1)) / 1000.0
            line = next(fp)
            grad_logp_jac = float(regex.match(line).group(1)) / 1000.0
            print(f"ratio_transform_jacobian,evaluation,off,{logp_jac}")
            print(f"ratio_transform_jacobian,gradient,off,{grad_logp_jac}")
        elif line.startswith("Constant coalescent"):
            line = next(fp)
            logp_coal = float(regex.match(line).group(1)) / 1000.0
            line = next(fp)
            grad_logp_coal = float(regex.match(line).group(1)) / 1000.0
            print(f"coalescent,evaluation,off,{logp_coal}")
            print(f"coalescent,gradient,off,{grad_logp_coal}")


print("function,mode,JIT,time")
if sys.argv[1] == "-":
    parse_physher(sys.stdin)
else:
    with open(sys.argv[1]) as fp:
        parse_physher(fp)
