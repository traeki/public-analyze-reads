#!/usr/bin/env python

# Author: John Hawkins (jsh)

import difflib
import logging
import optparse
import os
import sys
import multiprocessing

from Bio import Seq
from Bio import SeqIO

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


def main():
  logging.info('Parsing command line.')
  usage = '%prog [options] input_file [...]'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option(
      '--encoding', type='string',
      action='store', dest='encoding',
      default='fastq-illumina',
      help='Encoding to use for input.  Set to \'fastq\' for sanger encoding.')
  parser.add_option(
      '--primer', type='string',
      action='store', dest='primer_file',
      default='genomes/oLSC003.fna',
      help='File containing the primer to be stripped.')
  parser.add_option(
      '--min_sequence_len', type='int',
      action='store', dest='min_sequence_len',
      default=18,
      help='We discard sequences shorter than this (after trimming).')
  parser.add_option(
      '--min_primer_match', type='int',
      action='store', dest='min_primer_match',
      default=6,
      help='We ignore primer-matched tails shorter than this.')
  parser.add_option(
      '--max_primer_offset', type='int',
      action='store', dest='max_primer_offset',
      default=1,
      help='Only match leftmost edge of primer, allowing this much offset.')
  (options, args) = parser.parse_args()

  logging.info('Removing primer tails & cleaning out Illumina rejected reads.')
  primer = open(options.primer_file, 'r').readline()
  processes = list()
  for f in args:
    if os.path.exists(f):
      processes.append(multiprocessing.Process(
          target=trim_primers_from_file,
          args=(f,
                options.encoding,
                primer,
                options.min_sequence_len,
                options.min_primer_match,
                options.max_primer_offset)))
    else:
      logging.fatal('Could not find file "{0}".'.format(f))
  for p in processes:
    p.start()
  for p in processes:
    p.join()


def trim_primers_from_file(input_file,
                           encoding,
                           primer,
                           min_sequence_len,
                           min_primer_match,
                           max_primer_offset):
  """Process one file.

  Basically just a workhorse that calls processed_sequences and writes each
  element of the output to the destination file.
  """
  sequences = SeqIO.parse(input_file, encoding)
  output_name = os.path.splitext(input_file)[0] + '.trimmed'
  output_file = open(output_name, 'w')
  count = 0
  for seq in processed_sequences(primer,
                                 sequences,
                                 min_sequence_len,
                                 min_primer_match,
                                 max_primer_offset):
    output_file.write(seq.format('fastq-sanger'))
    count += 1
    if count % 100000 == 0:
      logging.info('file:{input_file} -- count:{count}'.format(**vars()))


def processed_sequences(primer,
                        sequences,
                        min_sequence_len,
                        min_primer_match,
                        max_primer_offset):
  """Clean, cut, and throw back small fish."""

  trimmed_sequences = (trim_primer(primer, s,
                                   min_primer_match,
                                   max_primer_offset)
                          for s in sequences if s)
  clean_sequences = (clean_for_illumina_flag(s) for s in trimmed_sequences if s)
  return (s for s in clean_sequences
             if len(s) >= min_sequence_len)


def rfind_if_not(tagger, things):
  """Return the index of the first thing NOT tagged, or len(things)."""
  return next((i for i in range(len(things)-1,0,-1) if not tagger(things[i])),
              len(things))


def clean_for_illumina_flag(sequence):
  """Remove any trailing 2-score rejects from the sequence.

  Returns cleaned sequence, or None if sequence should be filtered.  Phred
  score of "2" is actually a special annotation meaning Illumina itself has
  rejected these reads as an incoherent tail, and we don't want to be playing
  with fire."""
  scores = sequence.letter_annotations['phred_quality']
  good_quality_end = rfind_if_not(lambda x: x == 2, list(scores))
  sequence = sequence[:good_quality_end + 1]
  return sequence


def trim_primer(primer, sequence,
                min_primer_match,
                max_primer_offset):
  """Remove aligned primer tail.

  Returns trimmed sequence, or None if sequence should be filtered."""
  matcher = difflib.SequenceMatcher()
  matcher.set_seq1(sequence)
  matcher.set_seq2(primer)
  match = matcher.find_longest_match(0, len(sequence), 0, len(primer))
  if (match.size > min_primer_match and match.b <= max_primer_offset):
    return sequence[:match.a]
  else:
    return sequence


##############################################
if __name__ == "__main__":
    sys.exit(main())
