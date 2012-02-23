#!/usr/bin/env python

# Author: John Hawkins (jsh)

from collections import defaultdict
import logging
import optparse
import os
import string
import sys

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


def main():
  logging.info('Parsing command line.')
  usage = '%prog [options] input_shard [...]'
  parser = optparse.OptionParser(usage=usage)
  (options, args) = parser.parse_args()
  input_shards = args
  if not input_shards:
    binary = os.path.basename(sys.argv[0])
    logging.fatal(
        '{0} called with no input files specified.'.format(binary))
    return 2
  for input_shard in input_shards:
    if not os.path.exists(input_shard):
      logging.fatal('No such file: {0}'.format(input_shard))
      return 2
    if not '.shard.' in input_shard:
      complaint = 'Shard file names must contain ".shard.".  {0} does not.'
      logging.fatal(complaint.format(input_shard))
      return 2
  process_input_shards(input_shards)


def process_input_shards(input_shards):
  shard_sets = defaultdict(list)
  for shard in input_shards:
    base_name, shard_name = shard.split('.shard.')
    shard_sets[base_name].append(shard)
  for outfile_name, input_shards in shard_sets.iteritems():
    outfile = open(outfile_name, 'w')
    if outfile_name.endswith('.wig'):
      outfile.write('track type=wiggle_0\n')
    if outfile_name.endswith('.tallies'):
      header = open(input_shards[0], 'r').readline()
      assert('chrom\tdir\tpos\t' in header), header
      outfile.write(header)
    for shard in input_shards:
      infile = open(shard, 'r')
      for line in infile:
        if outfile_name.endswith('.wig'):
          if 'track type=wiggle' in line:
            continue
        if outfile_name.endswith('.tallies'):
          if 'chrom\tdir\tpos\t' in line:
            continue
        outfile.write(line)
    outfile.close()


##############################################
if __name__ == "__main__":
    sys.exit(main())
