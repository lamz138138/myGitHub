#!/bin/env python

import argparse
import falcon_kit.util.io as io
from falcon_kit.multiproc import Pool
import os
import sys

Reader = io.CapturedProcessReaderContext
endNodes = set()
nodeLinks = {}

def readEndNode(fileName):
    inEndNodes = open(fileName)
    for line in inEndNodes:
        endNodes.add(line.strip())
    inEndNodes.close()
    
def readNodeLinks(fileName):
    inReadNodeLinks = open (fileName)
    for line in inReadNodeLinks:
        line = line.strip().split()
        nodeLinks[line[0]] = {"5p": line[2], "3p": line[3]}

def run_overlapFilter(db_fn, fn, max_diff, max_cov, min_cov, min_len):
    cmd = "LA4Falcon -mo %s %s" % (db_fn, fn)
    reader = Reader(cmd)
    with reader:
        return fn, overlapFilter(reader.readlines, max_diff, max_cov, min_cov, min_len)
        
def overlapFilter(readlines, max_diff, max_cov, min_cov, min_len):
    output_data = []
    q_id = None
    for line in readlines():
        line = line.strip().split()
        q_id, t_id = line[:2]
        q_l = int(line[7])
        t_l = int(line[11])
        if q_id in endNodes:
            if q_l < min_len or t_l < min_len:
                continue
            if t_id not in nodeLinks.keys():
                continue
            left_count = int(nodeLinks[t_id]["5p"])
            right_count = int(nodeLinks[t_id]["3p"])
            if (abs(left_count - right_count) > max_diff) or \
                (left_count > max_cov) or (right_count > max_cov) or \
                (left_count < min_cov) or (right_count < min_cov):
                output_data.append(t_id)
    return output_data

def parse_args(argv):
    class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='a script use to calculate the reson of lacking link to node',
        formatter_class=HelpF)
    parser.add_argument(
        '--end_node', default='end_node.txt',
        help='file contain end node of contig')
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
        '--out_fn', default='node_repeat.txt',
        help='Output filename')
    parser.add_argument(
        '--n_core', type=int, default=4,
        help='number of processes used for generating consensus; 0 for main process only')
    parser.add_argument(
        '--max_diff', type=int, default=40,
        help="max difference of 5' and 3' coverage")
    parser.add_argument(
        '--max_cov', type=int, default=60,
        help="max coverage of 5' or 3' coverage")
    parser.add_argument(
        '--min_cov', type=int, default=1,
        help="min coverage of 5' or 3' coverage")
    parser.add_argument(
        '--min_len', type=int, default=2500,
        help="min length of the reads")
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    readEndNode(args.end_node)
    readNodeLinks(args.node_links)
    file_list = io.validated_fns(args.fofn)
    n_core = min(args.n_core, len(file_list))
    exe_pool = Pool(n_core)
    tmp_out_fn = args.out_fn + '.tmp'
    try:
        with open(tmp_out_fn, 'w') as outs:
            inputs = []
            for fn in file_list:
                if len(fn) != 0:
                    inputs.append((run_overlapFilter, args.db_fn, fn, args.max_diff, args.max_cov, args.min_cov, args.min_len))
            for res in exe_pool.imap(io.run_func, inputs):
                for l in res[1]:
                    outs.write("".join(l) + "\n")
        os.rename(tmp_out_fn, args.out_fn)
    except:
        io.LOG('terminating ovlp_filter workers...')
        exe_pool.terminate()
        raise
    

if __name__ == "__main__":
    main(sys.argv)

