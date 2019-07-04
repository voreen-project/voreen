#! /usr/bin/env python3

# A script used to plot data coming from a unix pipe using pyplot.
# The only argument must be a path to the (named) unix pipe.
from io import StringIO
from matplotlib import pyplot as plt
import sys
import pipes
import numpy as np
import threading

def go():
    pipepath = sys.argv[1]
    while True:
        with open(pipepath) as pipe:
            strdata = pipe.read().encode('utf8');
            lines = strdata.split(b',')
            print("jau!")
            plt.clf();
            for line in lines:
                data = np.fromstring(line, dtype=float, sep=' ')
                plt.plot(data)

            plt.ylim([0,0.03])
            plt.autoscale(enable=False, axis='y')
            plt.draw()
            plt.pause(0.001)

def main():
    if len(sys.argv) > 1:
        plt.plot(np.fromstring("1 2 3 4 5", dtype=float, sep=' '))
        thread = threading.Thread(target=go)
        thread.start()
        plt.show()

    else:
        print("Pipepath please")

if __name__ == "__main__":
    main()
