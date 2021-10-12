# prog_bar.py
#
# Progress bar methods to build, increment, and close
# a progress bar that looks like [>>>>>>>>  ]

import numpy as np

import sys, time

def prep(bar_length):
    sys.stdout.write("[%s]" % (" " * bar_length))
    sys.stdout.flush()

    sys.stdout.write("\b" * (bar_length + 1))
    sys.stdout.flush()

def increment(n=1):
    for j in range(n):
        sys.stdout.write(">")
        sys.stdout.flush()
        time.sleep(0.01)

def close():
    sys.stdout.write("\n")
    sys.stdout.flush()

def set_step(n, N, bar_length):
    return int(np.floor((float(bar_length) * (n + 1)) / N) - np.floor((float(bar_length) * n) / N))
