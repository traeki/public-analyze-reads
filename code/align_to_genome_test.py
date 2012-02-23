#!/usr/bin/env python

# Author: John Hawkins (jsh)

import align_to_genome as atg

import logging
import unittest
import os
import subprocess
import sys
import tempfile

import HTSeq

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


def generate_alignments(genome, reads, use_tophat):
  alignments = None
  test_dir = tempfile.mkdtemp(dir='/tmp')
  logging.info('Working in {0}.'.format(test_dir))
  file_base = os.path.join(test_dir, 'test')
  fastq_name = file_base + '.fastq'
  fasta_name = file_base + '.fna'
  with open(fastq_name, 'w') as fastq_file:
    for n, r in enumerate(reads):
      name = 'seq{0}qes'.format(n)
      s = HTSeq.SequenceWithQualities(r, name, 'B' * len(r))
      s.write_to_fastq_file(fastq_file)
  with open(fasta_name, 'w') as fasta_file:
    fasta_file.write('>chrFOO\n')
    fasta_file.write(genome)
  command = ['bowtie-build']
  command.append(fasta_name)
  command.append(file_base)
  logging.info(' '.join(command))
  with open(file_base + '.bowtie_build_stdout', 'w') as bowtie_stdout:
    with open(file_base + '.bowtie_build_stderr', 'w') as bowtie_stderr:
      subprocess.check_call(command,
                            stdout=bowtie_stdout,
                            stderr=bowtie_stderr)
  if use_tophat:
    align_to_genome = atg.align_with_tophat
    binary_path = 'tophat'
  else:
    align_to_genome = atg.align_with_bowtie
    binary_path = 'bowtie'
  alignment_file = align_to_genome([fastq_name],
                                   file_base,
                                   binary_path,
                                   7, 3, 5, 12, 5000,
                                   fasta_name,
                                   file_base)
  return [x for x in HTSeq.SAM_Reader(alignment_file)]


def set_nth(read, n, base):
  return read[:n] + base + read[n+1:]


class TestEndToEnd(unittest.TestCase):
  def setUp(self):
    self.genome_plus_strand = 'AAAAAAAAAAAAACCCCCCAAAACAAAAAAAAAAAAAAAAAAAA'

    self.plus_read = 'AAAAAAAAAAAAACCCCCCAAAACAAAAAAAAAA'
    self.minus_read = 'TTTTTTTTTTTTTTTTGTTTTGGGGGGTTTTTTT'

  def testTrivialBowtieRun(self):
    logging.info('Begin testTrivialBowtieRun')
    reads = [self.plus_read, self.minus_read]
    alignments = generate_alignments(self.genome_plus_strand, reads, False)
    self.assertEqual(2, len(alignments))
    self.assertEqual(0, alignments[0].iv.start)
    self.assertEqual(34, alignments[0].iv.end)
    self.assertEqual(6, alignments[1].iv.start)
    self.assertEqual(40, alignments[1].iv.end)

  def testTrivialTophatRun(self):
    logging.info('Begin testTrivialTophatRun')
    reads = [self.plus_read, self.minus_read]
    alignments = generate_alignments(self.genome_plus_strand, reads, True)
    self.assertEqual(2, len(alignments))
    self.assertEqual(0, alignments[0].iv.start)
    self.assertEqual(34, alignments[0].iv.end)
    self.assertEqual(6, alignments[1].iv.start)
    self.assertEqual(40, alignments[1].iv.end)

  def testMisincorpBowtie(self):
    logging.info('Begin testMisincorpBowtie')
    plus_read = self.plus_read
    plus_read = set_nth(plus_read, 7, 'G')
    plus_read = set_nth(plus_read, 9, 'T')
    minus_read = self.minus_read
    minus_read = set_nth(minus_read, -4, 'C')
    minus_read = set_nth(minus_read, -10, 'A')
    reads = [plus_read, minus_read]
    alignments = generate_alignments(self.genome_plus_strand, reads, False)
    self.assertEqual(2, len(alignments))

    self.assertEqual(0, alignments[0].iv.start)
    self.assertEqual(34, alignments[0].iv.end)
    self.assertEqual('7A1A24', alignments[0].optional_field('MD'))

    self.assertEqual(6, alignments[1].iv.start)
    self.assertEqual(40, alignments[1].iv.end)
    # MD strings are relative to plus strand.  Even when the alignment is to
    # the minus strand.  (Which seems wrong, but hey, if they fix it this
    # unittest will break, so we'll find out).
    # NOTE: That means that IF THIS UNITTEST SUDDENLY BREAKS it might be
    # because a new version of tophat/bowtie is doing the right thing, and
    # that means that downstream code using MD might be brokenly trying to fix
    # problems that no longer exist.  Look into that before just blindly
    # "fixing" this test!
    self.assertEqual('3A5C24', alignments[1].optional_field('MD'))


  def testMisincorpTophat(self):
    logging.info('Begin testMisincorpTophat')
    plus_read = self.plus_read
    plus_read = set_nth(plus_read, 7, 'G')
    plus_read = set_nth(plus_read, 9, 'T')
    minus_read = self.minus_read
    minus_read = set_nth(minus_read, -4, 'C')
    minus_read = set_nth(minus_read, -10, 'A')
    reads = [plus_read, minus_read]
    alignments = generate_alignments(self.genome_plus_strand, reads, True)
    self.assertEqual(2, len(alignments))

    self.assertEqual(0, alignments[0].iv.start)
    self.assertEqual(34, alignments[0].iv.end)
    self.assertEqual('7A1A24', alignments[0].optional_field('MD'))

    self.assertEqual(6, alignments[1].iv.start)
    self.assertEqual(40, alignments[1].iv.end)
    # MD strings are relative to plus strand.  Even when the alignment is to
    # the minus strand.  (Which seems wrong, but hey, if they fix it this
    # unittest will break, so we'll find out).
    # NOTE: That means that IF THIS UNITTEST SUDDENLY BREAKS it might be
    # because a new version of tophat/bowtie is doing the right thing, and
    # that means that downstream code using MD might be brokenly trying to fix
    # problems that no longer exist.  Look into that before just blindly
    # "fixing" this test!
    self.assertEqual('3A5C24', alignments[1].optional_field('MD'))


##############################################
if __name__ == '__main__':
  unittest.main()
