#! /usr/bin/python3

# Convert a raw volume file to a csv file for convenient histogram
# visualization using matlab

import sys

def processFile(inputpath, outputpath):
    with open(outputpath, 'w') as o:
        o.write("value\n")
        with open(inputpath, 'rb') as i:
            bytes = i.read()
            i = 0
            for byte in bytes:
                o.write("{}\n".format(byte))
                i = i+1
                if i%100000 == 0:
                    print("{:2f}%".format(100*i/len(bytes)))


def main():
    if len(sys.argv) != 3:
        print("usage: <script> input.raw output.csv")
    else:
        processFile(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
