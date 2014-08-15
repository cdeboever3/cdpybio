import os

def mount_bam(analysis_id, mount='/tmp', cache='/tmp/fusecache'):
    """
    Mount bam file for given analysis ID using GTFuse
    
    Creates a directory for this analysis ID in "mount" and mounts 
    makes mnt and cache directories then mounts the bam file using 
    GTFuse. Returns bam file path.

    Parameters
    ----------
    analysis_id : str
        TCGA analysis_id for sample of interest
    mount : directory
        Directory to make temporary directory and files.
    cache : directory
        Directory to store GTFuse cache files.
    
    Returns
    -------
    bam : path
        Path to the mounted bam file.
    
    """
    if not os.path.exists(cache):
        os.makedirs(cache)
    analysis_dir = '{}_{}'.format(analysis_id,
                                  time.strftime('%Y_%m_%d_%H_%M_%S'))
    mnt = os.path.join(mount, analysis_dir, 'mnt')
    os.makedirs(mnt)
    url = ('https://cghub.ucsc.edu/cghub/data'
           '/analysis/download/{}'.format(analysis_id))
    subprocess.check_call(['gtfuse', 
                           '--ssl-no-verify-ca',
                           '--cache-dir',
                           cache,
                           url,
                           mnt])
    files = glob.glob('{}/{}/*'.format(mnt, analysis_id))
    bam = [ x for x in files if os.path.splitext(x)[1] == '.bam' ][0]
    return os.path.realpath(bam)
    
def unmount_bam(bam):
    """
    Unmount bam file mounted with GTFuse
    
    Parameters
    ----------
    bam : path
        Path to the mounted bam file.

    """
    p = os.path.sep.join(bam.split(os.path.sep)[0:-3] + ['*'])
    files = glob.glob(p)
    mnt = os.path.sep.join(bam.split(os.path.sep)[0:-2])
    for f in files:
        if f != mnt:
            os.remove(f)
    subprocess.call('fusermount -u {}'.format(mnt),
                    shell=True)
    os.rmdir(mnt)
    os.rmdir(os.path.split(mnt)[0])
    
def reads_from_intervals(analysis_id, intervals, 
                         interval_bam,
                         max_intervals=4000):
    """
    Get reads from analysis id for given intervals 
    
    Parameters
    ----------
    analysis_id : str
        TCGA analysis_id for sample of interest
    intervals : list
        List of intervals of the form 1:10-200.
    interval_bam : str
        Path to the bam file where the reads for the intervals will be written.
    max_intervals : int
        Maximum number of intervals to obtain with one samtools view call.
    
    Returns
    -------
    interval_bam : str
        Path to the bam file with the reads for the intervals.
    
    """
    tcga_bam = mount_bam(analysis_id)
    # I'll make a tempdir to hold the mounted bam and any other temp files I
    # want to make.
    tempdir = os.path.sep.join(tcga_bam.split(os.path.sep)[0:-3])
    temp_bams = []
    for i in xrange(0, intervals.shape[0], max_intervals):
        ints = ' '.join(intervals[i:i + max_intervals])
        temp_bam = '{}_{}.bam'.format(tempdir, i)
        subprocess.check_call('samtools view -b {} {} > {}'.format(tcga_bam, 
                                                                   ints,
                                                                   temp_bam),   
                              shell=True)
        temp_bams.append(temp_bam)
    if len(temp_bams) > 1:
        subprocess.check_call(['samtools', 'merge', '-f', interval_bam, 
                               ' '.join(temp_bams)])
        for b in temp_bams:
            os.remove(b)
    else:
        shutil.move(temp_bams[0], interval_bam)
    unmount_bam(tcga_bam)
    return interval_bam

def bed_to_samtools_intervals(bed):
    """
    Convert bed file to intervals for samtools view

    Parameters
    ----------
    bed : path
        Bed file to convert

    Returns
    -------
    intervals : list
        List of the intervals in the bed file in the format 1:20-200 etc.

    """
    intervals = []
    f = open(bed, 'r')
    line = f.readline().strip()
    while line != '':
        vals = line.split('\t')
        intervals.append('{}:{}-{}'.format(vals[0], vals[1], vals[2]))
        line = f.readline().strip()
    return intervals

class ReadsFromIntervalsEngine:
    def __init__(self, analysis_ids, bed, outdir='.', threads=10,
                 sleeptime=10):
        """
        Initialize engine for obtaining reads for given intervals/IDs
        
        Parameters
        ----------
        analysis_ids : list
            List of TCGA analysis IDs.
        bed : path
            Bed file with intervals.
        outdir : directory
            Directory to write bam files to.
        threads : int
            Number of different processes/threads to use.
        sleeptime : int
            Number of seconds to sleep between engine updates.

        """
        assert len(analysis_ids) > 0
        assert os.path.exists(intervals)
        self.analysis_ids = analysis_ids
        self.bed = bed
        self.outdir = outdir
        self.sleeptime = sleeptime
        # Intervals in the format 1:20-200 etc. as a list.
        self.intervals = bed_to_samtools_intervals(self.bed)
        # Analysis_ids that the engine has started getting reads for.
        self.started = set()
        # Analysis_ids that the engine has not started getting reads for.
        self.remaining = set(self.analysis_ids)
        self.current_procs = []
        self.old_procs = []
        # Paths to bam files with reads.
        self.bams = []
        # Queue used when multiprocessing.
        self.queue = None
        self.running = True
        self.run_engine()

    def remove_finished_procs(self):
        for p in self.procs:
            # If we find a dead process, there should be a result in the queue.
            # The result will not necessarily be from that dead process though.
            if not p.is_alive():
                p.join()
                self.bams.append(self.queue.get())
                self.current_procs.remove(p)
                self.old_procs.append(p)

    def new_proc(self):
        import multiprocessing
        if len(self.remaining) > 0:
            aid = self.remaining.pop()
            bam = '{}.bam'.format(os.path.join(self.outdir, analysis_id))
            p = multiprocessing.Process(target=reads_from_intervals,
                                        args=[aid, self.intervals, bam])
            self.current_procs.append(p)
            p.start()

    def add_procs(self):
        while len(self.current_procs) < self.threads:
            self.new_proc()

    def run_engine(self):
        import multiprocessing
        import sys
        import time

        i = 0
        procs = []
        self.queue = multiprocessing.Queue()
        while i < len(self.analysis_ids) and len(self.remaining) > 0:
            if not self.running:
                break
            self.remove_finished_procs()
            self.add_procs()
            time.sleep(self.sleeptime)
        # If we get here, we aren't running any more. Wait for processes to
        # finish, .join() them, then exit.
        sys.stdout.write('Engine stopping, waiting for jobs to conclude.\n')
        while len(self.current_procs) > 0:
            self.remove_finished_procs()
            time.sleep(self.sleeptime)
        sys.stdout.write('Jobs concluded, engine stopped.\n')
