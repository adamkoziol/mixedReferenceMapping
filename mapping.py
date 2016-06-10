#!/usr/bin/env python
__author__ = 'adamkoziol'

class Mapper(object):

    def mapping(self):
        pass

    def mappingthreads(self):
        while True:
            pass

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        from Queue import Queue
        import multiprocessing
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.args = args
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Output location is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath, '')
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.threads if args.threads else multiprocessing.cpu_count())
        self.mappingqueue = Queue()

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
    # parser.add_argument('-v', '--version',
    #                     version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    Mapper(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'
    # print json.dumps(seqdict, sort_keys=True, indent=4, separators=(',', ': '))