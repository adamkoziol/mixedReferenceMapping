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
        printtime('Moving files and creating objects', self.starttime)
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
        # In order to keep from running too many cpu intensive processes concurrently (multi-threaded applications
        # being run in a multi-threaded fashion), decrease the number of threads used in the applications to a minimum
        # of four depending on how many samples are being processed
        numthreads = self.cpus / len(self.runmetadata.samples)
        numthreads = numthreads if numthreads >= 4 else 4
        for sample in self.runmetadata.samples:
            # Create the analysis type attribute
            setattr(sample, self.analysistype, GenObject())
            # Set the path/name for the bam and sorted bam files to be created
            sample[self.analysistype].bam = '{}/{}.bam'.format(sample.general.outputdirectory, sample.name)
            sample[self.analysistype].sortedbam = '{}/{}_sorted.bam'.format(sample.general.outputdirectory, sample.name)
            # Set the output name depending on whether the bam files are to be sorted or not
            output = sample[self.analysistype].bam if not self.sort else sample[self.analysistype].sortedbam
            # Remove the file extension of the bait file for use in the indexing command
            sample[self.analysistype].referencenoext = self.referencefile.split('.')[0]
            # Use bowtie2 wrapper to create index the target file
            bowtie2build = Bowtie2BuildCommandLine(reference=self.referencefile,
                                                   bt2=sample[self.analysistype].referencenoext)
            # If the indexing option is specified, include the sort command
            samsort = SamtoolsSortCommandline(input_bam=output,
                                              o=True,
                                              out_prefix="-")
            # Create a list of programs to which data are piped as part of the reference mapping
            samtools = [
                # Use samtools wrapper to set up the samtools view
                SamtoolsViewCommandline(b=True,
                                        S=True,
                                        input_file="-")
            ]
            if self.sort:
                samtools.append(samsort)
            # Add custom parameters to a dictionary to be used in the bowtie2 alignment wrapper
            if len(sample.general.trimmedcorrectedfastqfiles) == 1:
                indict = {'-U': sample.general.trimmedcorrectedfastqfiles[0],
                          '--threads': numthreads
                          }
            else:
                indict = {'m1': sample.general.trimmedcorrectedfastqfiles[0],
                          'm2': sample.general.trimmedcorrectedfastqfiles[1],
                          '--threads': numthreads
                          }
            # Create the bowtie2 reference mapping command
            bowtie2align = Bowtie2CommandLine(bt2=sample[self.analysistype].referencenoext,
                                              threads=numthreads,
                                              samtools=samtools,
                                              **indict)

            # Add the commands (as strings) to the metadata
            sample[self.analysistype].bowtie2align = str(bowtie2align)
            sample[self.analysistype].bowtie2build = str(bowtie2build)
            # Get variables ready for the bam moving step
            sample[self.analysistype].bampath = os.path.join(self.path, 'bamfiles')
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
            self.mapqueue.put((sample, bowtie2build, bowtie2align, output))
        self.mapqueue.join()
        metadataprinter.MetadataPrinter(self)

    def map(self):
        import shutil
        while True:
            # Get the necessary values from the queue
            sample, bowtie2build, bowtie2align, output = self.mapqueue.get()
            # Only run the functions if the sorted bam files and the indexed bait file do not exist
            if not os.path.isfile(output) and not os.path. \
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
            if self.sort:
                # Perform indexing
                # Only index if the index file doesn't exist
                if not os.path.isfile(output + '.bai'):
                    samtoolsindex = SamtoolsIndexCommandline(input=output)
                    # Set stdout to a stringIO stream
                    stdout, stderr = map(StringIO, samtoolsindex(cwd=sample.general.outputdirectory))
                    if stderr:
                        # Write the standard error to log, bowtie2 puts alignment summary here
                        with open(os.path.join(sample.general.outputdirectory,
                                               '{}_bowtie_index.log'.format(self.analysistype)), 'ab+') as log:
                            log.writelines(logstr(samtoolsindex, stderr.getvalue(), stdout.getvalue()))
                    stdout.close()
                    stderr.close()
            else:
                make_path(sample[self.analysistype].bampath)
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
        # Determine whether BAM files need to be sorted and indexed
        self.sort = args.index
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
    parser.add_argument('-i', '--index',
                        action='store_true',
                        help='Sort and index the bam files')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    Mapper(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'
