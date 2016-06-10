#!/usr/bin/env python
from threading import Thread

from Bio.Sequencing.Applications import *

import SPAdesPipeline.OLCspades.metadataprinter as metadataprinter
from SPAdesPipeline.OLCspades.bowtie import *
from geneSipprV2.objectOriented.createObject import *

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
__author__ = 'adamkoziol'


class Mapper(object):

    def mapping(self):
        # Create the object
        self.runmetadata = ObjectCreation(self)
        metadataprinter.MetadataPrinter(self)
        printtime('Performing reference mapping', self.starttime)
        for i in range(len(self.runmetadata.samples)):
            # Send the threads to
            threads = Thread(target=self.map, args=())
            # Set the daemon to True - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata.samples:
            # Create the analysis type attribute
            setattr(sample, self.analysistype, GenObject())
            # Set the path/name for the bam file to be created
            sample[self.analysistype].bam = '{}/{}.bam'.format(sample.general.outputdirectory, sample.name)
            # Remove the file extension of the bait file for use in the indexing command
            sample[self.analysistype].referencenoext = self.referencefile.split('.')[0]
            # Use bowtie2 wrapper to create index the target file
            bowtie2build = Bowtie2BuildCommandLine(reference=self.referencefile,
                                                   bt2=sample[self.analysistype].referencenoext)
            # Create a list of programs to which data are piped as part of the reference mapping
            samtools = [
                # Use samtools wrapper to set up the samtools view
                SamtoolsViewCommandline(b=True,
                                        S=True,
                                        o=sample[self.analysistype].bam,
                                        input_file="-")
            ]
            # Add custom parameters to a dictionary to be used in the bowtie2 alignment wrapper
            if len(sample.general.trimmedcorrectedfastqfiles) == 1:
                indict = {'-U': sample.general.trimmedcorrectedfastqfiles[0],
                          '--threads': self.cpus
                          }
            else:
                indict = {'m1': sample.general.trimmedcorrectedfastqfiles[0],
                          'm2': sample.general.trimmedcorrectedfastqfiles[1],
                          '--threads': self.cpus
                          }
            # Create the bowtie2 reference mapping command
            bowtie2align = Bowtie2CommandLine(bt2=sample[self.analysistype].referencenoext,
                                              threads=self.cpus,
                                              samtools=samtools,
                                              **indict)
            # Add the commands (as strings) to the metadata
            sample[self.analysistype].bowtie2align = str(bowtie2align)
            sample[self.analysistype].bowtie2build = str(bowtie2build)
            # Get variables ready for the bam moving step
            sample[self.analysistype].bampath = os.path.join(self.path, 'bamfiles')
            make_path(sample[self.analysistype].bampath)
            sample[self.analysistype].bamcollection = \
                '{}/{}'.format(sample[self.analysistype].bampath, os.path.basename(sample[self.analysistype].bam))
            # Add the commands to the queue. Note that the commands would usually be set as attributes of the sample
            # but there was an issue with their serialization when printing out the metadata
            if not os.path.isfile(sample[self.analysistype].referencenoext + '.1.bt2'):
                stdoutbowtieindex, stderrbowtieindex = map(StringIO,
                                                           bowtie2build(cwd=self.referencepath))
                # Write any error to a log file
                if stderrbowtieindex:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(self.referencepath,
                                           '{}_bowtie_index.log'.format(self.analysistype)), 'ab+') as log:
                        log.writelines(logstr(bowtie2build, stderrbowtieindex.getvalue(),
                                              stdoutbowtieindex.getvalue()))
                # Close the stdout and stderr streams
                stdoutbowtieindex.close()
                stderrbowtieindex.close()
            # Populate the queue
            self.mapqueue.put((sample, bowtie2build, bowtie2align))
        self.mapqueue.join()
        metadataprinter.MetadataPrinter(self)

    def map(self):
        import shutil
        while True:
            # Get the necessary values from the queue
            sample, bowtie2build, bowtie2align = self.mapqueue.get()
            # Only run the functions if the sorted bam files and the indexed bait file do not exist
            if not os.path.isfile(sample[self.analysistype].bam) and not os.path. \
                    isfile(sample[self.analysistype].bamcollection):
                # Set stdout to a stringIO stream
                stdout, stderr = map(StringIO, bowtie2align(cwd=sample.general.outputdirectory))
                if stderr:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(sample.general.outputdirectory,
                                           '{}_bowtie_samtools.log'.format(self.analysistype)), 'ab+') as log:
                        log.writelines(logstr(bowtie2align, stderr.getvalue(), stdout.getvalue()))
                stdout.close()
                stderr.close()
            # Move the bam file to a common location
            if not os.path.isfile(sample[self.analysistype].bamcollection):
                shutil.move(sample[self.analysistype].bam, sample[self.analysistype].bamcollection)
            self.mapqueue.task_done()

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        from Queue import Queue
        import multiprocessing
        from glob import glob
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.args = args
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath, '')
        assert os.path.isdir(self.sequencepath), u'Sequence folder is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        self.referencepath = os.path.join(args.referencepath, '')
        assert os.path.isdir(self.sequencepath), u'Reference folder is not a valid directory {0!r:s}' \
            .format(self.referencepath)
        try:
            self.referencefile = glob('{}*.fa*'.format(self.referencepath))[0]
        except IndexError:
            print 'Cannot find a .fa/.fas/.fasta reference file in the supplied reference path: {}' \
                .format(self.referencepath)
            quit()
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.threads if args.threads else multiprocessing.cpu_count())
        self.runmetadata = MetadataObject()
        self.mapqueue = Queue(maxsize=self.cpus)
        self.analysistype = 'referencemapping'
        # Run the analyses
        self.mapping()

if __name__ == '__main__':
    import subprocess
    import time
    import os
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-r', '--referencepath',
                        required=True,
                        help='Path to folder containing a single fasta file to be used as the reference genome')
    parser.add_argument('-t', '--threads',
                        help='Number of threads. Default is the number of cores in the system')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    Mapper(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'
