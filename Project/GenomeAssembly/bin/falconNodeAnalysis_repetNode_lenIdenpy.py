#!/bin/env python

import argparse
import falcon_kit.util.io as io
from falcon_kit.multiproc import Pool
import os
import sys

Reader = io.CapturedProcessReaderContext
repeatNodes = set()

def readRepeatNode(fileName):
    inRepeatNodes = open(fileName)
    for line in inRepeatNodes:
        repeatNodes.add(line.strip())
    inRepeatNodes.close()
    
def run_overlapFilter(db_fn, fn, max_diff, max_cov, min_cov, min_len, overlapCount):
    cmd = "LA4Falcon -mo %s %s" % (db_fn, fn)
    reader = Reader(cmd)
    with reader:
        return fn, overlapFilter(reader.readlines, max_diff, max_cov, min_cov, min_len, overlapCount)
        
def overlapFilter(readlines, max_diff, max_cov, min_cov, min_len, overlapCount):
    output_data = []
    current_q_id = None
    q_id = None
    outputString_5 = None
    outputString_3 = None
    index_5 = 0
    index_3 = 0
    for line in readlines():
        line = line.strip().split()
        q_id, t_id = line[:2]
        q_s, q_e, q_l = int(line[5]), int(line[6]), int(line[7])
        t_l = int(line[11])
        if q_id in repeatNodes:
            if q_l < min_len or t_l < min_len:
                continue
            overlapLen = abs(int(line[2]))
            idt = line[3]
            if q_id != current_q_id:
                if current_q_id is not None:
                    output_data.append(outputString_5)
                    output_data.append(outputString_3)
                current_q_id = q_id
                if q_s == 0:
                    outputString_5 = '%s%s%s%s%s%s%s%s' % (q_id, "\t", "5", "\t", str(overlapLen), "(", idt, ")")
                elif q_e == q_l:
                    outputString_3 = '%s%s%s%s%s%s%s%s' % (q_id, "\t", "3", "\t", str(overlapLen), "(", idt, ")")
                index_5 = 1
                index_3 = 1
            else:
                if q_s == 0:
                    if index_5 <= overlapCount:
                        if outputString_5 is not None:
                            outputString_5 = '%s%s%s%s%s%s' % (outputString_5, "\t", str(overlapLen), "(", idt, ")")
                        else:
                            outputString_5 = '%s%s%s%s%s%s%s%s' % (q_id, "\t", "5", "\t", str(overlapLen), "(", idt, ")")
                        index_5 += 1
                    else:
                        continue
                elif q_e == q_l:
                    if index_3 <= overlapCount:
                        if outputString_3 is not None:
                            outputString_3 = '%s%s%s%s%s%s' % (outputString_3, "\t", str(overlapLen), "(", idt, ")")
                        else:
                            outputString_3 = '%s%s%s%s%s%s%s%s' % (q_id, "\t", "3", "\t", str(overlapLen), "(", idt, ")")
                        index_3 += 1
                    else:
                        continue
    if current_q_id in repeatNodes:         
        output_data.append(outputString_5)
        output_data.append(outputString_3)
    return output_data

def parse_args(argv):
    class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='a script use to calculate the overlap length, identity of repeat node',
        formatter_class=HelpF)
    parser.add_argument(
        '--repeat_node', default='repeat_node.txt',
        help='file contain repeat node')
    parser.add_argument(
        '--db', type=str, dest='db_fn',
        help='read db file path')
    parser.add_argument(
        '--fofn', type=str,
        help='file contains the path of all LAS file to be processed in parallel')
    parser.add_argument(
        '--out_fn', default='node_repeat_lenIden.txt',
        help='Output filename')
    parser.add_argument(
        '--n_core', type=int, default=4,
        help='number of processes used for generating consensus; 0 for main process only')
    parser.add_argument(
        '--overlap_count', type=int, default=500,
        help='only useing the first --overlap_count time overlap of repeat')
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
    readRepeatNode(args.repeat_node)
    file_list = io.validated_fns(args.fofn)
    n_core = min(args.n_core, len(file_list))
    exe_pool = Pool(n_core)
    tmp_out_fn = args.out_fn + '.tmp'
    try:
        with open(tmp_out_fn, 'w') as outs:
            inputs = []
            for fn in file_list:
                if len(fn) != 0:
                    inputs.append((run_overlapFilter, args.db_fn, fn, args.max_diff, args.max_cov, args.min_cov, args.min_len, args.overlap_count))
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

