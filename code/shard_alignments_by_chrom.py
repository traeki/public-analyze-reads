#!/usr/bin/env python

# Author: John Hawkins (jsh)

import itertools
import logging
import math
import numpy
import optparse
import os
import string
import sys

import HTSeq

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


def main():
  logging.info('Parsing command line.')
  usage = '%prog [options] input_file [...] output_file_base'
  parser = optparse.OptionParser(usage=usage)
  (options, args) = parser.parse_args()
  input_files = args[:-1]
  if not input_files:
    binary = os.path.basename(sys.argv[0])
    logging.fatal(
        '{0} saw only one input argument, expecting 2 or more'.format(binary))
    return 2
  output_file_base = args[-1]
  for input_file in input_files:
    if not os.path.exists(input_file):
      logging.fatal('No such file: {0}'.format(input_file))
      return 2
  process_input_files(input_files, output_file_base)


def process_input_files(input_files,
                        output_file_base):
  """Generate wiggle files for one input file.

  Build the coverage/mismatch arrays from individual alignment positions, and
  write out the plus and minus directional bedgraph/wiggle files.
  """
  n_seen = 0
  total_aligned = 0
  shard_files = {}
  for input_file in input_files:
    logging.info('Processing file: {0}'.format(input_file))
    sam_alignments = HTSeq.SAM_Reader(input_file)
    groups = itertools.groupby(sam_alignments, lambda x: x.read.name)
    for name, grouper in groups:
      n_seen += 1
      if n_seen % 1000000 == 0:
        logging.info('Now sharding group #{n_seen} {name}.'.format(**vars()))
      alignments = list(grouper)
      alignment = select_alignment(alignments)
      if alignment:
        shard = alignment.iv.chrom
        if shard in ['chrMito', '2-micron']:
          continue
        if shard_files.has_key(shard):
          broken = alignment.get_sam_line()
          fixed = broken.translate(string.maketrans(' ', '\t'))
          shard_files[shard].write(fixed + '\n')
        else:
          shard_files[shard] = open(
              output_file_base + '.alignment.sam.shard.{0}'.format(shard), 'w')
        total_aligned += 1
  logging.info('total_aligned: {total_aligned}.'.format(**vars()))
  count_file = open(output_file_base + '.alignment.sam.count', 'w')
  count_file.write(str(total_aligned))
  count_file.write('\n')
  count_file.close()


def select_alignment(alignments):
  """Return unique alignment, or None."""
  alignments = [x for x in alignments if x.aligned]
  alignments = [x for x in alignments if not '_' in x.iv.chrom]
  if len(alignments) == 1:
    return alignments[0]
  else:
    return None


##############################################
if __name__ == "__main__":
    sys.exit(main())
