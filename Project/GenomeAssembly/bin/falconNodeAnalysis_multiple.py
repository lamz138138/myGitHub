#!/bin/env python

import argparse
import falcon_kit.util.io as io
from falcon_kit.multiproc import Pool
import os
import sys
from itertools import combinations

Reader = io.CapturedProcessReaderContext
nodeLinks = {}

def readNodeLinks(fileName):
    inReadNodeLinks = open(fileName)
    for line in inReadNodeLinks:
        line = line.strip().split()
        for i in range(1,len(line)-1):
            for j in range(i+1, len(line)):
                node_1 = line[i]
                node_2 = line[j]
                if (node_1 in nodeLinks.keys() and node_2 in nodeLinks[node_1].keys()) or \
                (node_2 in nodeLinks.keys() and node_1 in nodeLinks[node_2].keys()):
                    continue
                if node_1 not in nodeLinks.keys():
                    nodeLinks[node_1] = {}
                nodeLinks[node_1][node_2]=0
    inReadNodeLinks.close()

def parse_args(argv):
    class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='a script use to calculate the reson of no overlap for multiple link node',
        formatter_class=HelpF)
    parser.add_argument(
        '--node_links', default='node_links.txt',
        help='file contain link count of node')
    parser.add_argument(
        '--db', type=str, dest='db_fn',
        help='read db file path')
    parser.add_argument(
        '--fofn', type=str,
        help='file contains the path of all LAS file to be processed in parallel')
    parser.add_argument(
        '--out_fn', default='nodes_overlap.txt',
        help='Output filename')
    parser.add_argument(
        '--n_core', type=int, default=4,
        help='number of processes used for generating consensus; 0 for main process only')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    readNodeLinks(args.node_links)
    file_list = io.validated_fns(args.fofn)
    for fn in file_list:
        cmd = "LA4Falcon -mo %s %s" % (args.db_fn, fn)
        commandOut = os.popen(cmd)
        commandInfo = commandOut.readlines()
        for line in commandInfo:
            line = line.strip().split()
            node_1, node_2 = line[:2]
            idt = line[3]
            if node_1 in nodeLinks.keys() and node_2 in nodeLinks[node_1].keys():
                nodeLinks[node_1][node_2] = idt
            elif node_2 in nodeLinks.keys() and node_1 in nodeLinks[node_2].keys():
                nodeLinks[node_2][node_1] = idt
    inReadNodeLinks_2 = open(args.node_links)
    output = open(args.out_fn, 'w')
    for line in inReadNodeLinks_2:
        overlapCount = 0
        line = line.strip().split()
        compareCount = len(list(combinations(line[1:len(line)],2)))
        for i in range(1,len(line)-1):
            for j in range(i+1, len(line)):
                idt = 0
                node_1 = line[i]
                node_2 = line[j]
                if node_1 in nodeLinks.keys() and node_2 in nodeLinks[node_1].keys():
                    idt = float(nodeLinks[node_1][node_2])
                elif node_2 in nodeLinks.keys() and node_1 in nodeLinks[node_2].keys():
                    idt = float(nodeLinks[node_2][node_1])
                if idt > 90:
                    overlapCount += 1
        output.write("\t".join((str(line[1]), str(len(line)-1), str(compareCount), str(overlapCount))) + "\n")
    inReadNodeLinks_2.close()
    output.close()
 
if __name__ == "__main__":
    main(sys.argv)

