# prog_bar.py
#
# Progress bar methods to build, increment, and close
# a progress bar that looks like [>>>>>>>>  ]

import sys, time

def prep(length):
    sys.stdout.write("[%s]" % (" " * length))
    sys.stdout.flush()

    sys.stdout.write("\b" * (length + 1))
    sys.stdout.flush()

def increment(n=1):
    for j in range(n):
        sys.stdout.write(">")
        sys.stdout.flush()
        time.sleep(0.02)

def close():
    sys.stdout.write("\n")
    sys.stdout.flush()

def set_step(n, N, len):
    return int(np.floor((float(len) * (n + 1)) / N) - np.floor((float(len) * n) / N))
