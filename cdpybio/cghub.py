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
    import glob
    import subprocess
    import time
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
    
def unmount_bam(bam, tries=10, wait=10):
    """
    Unmount bam file mounted with GTFuse
    
    Parameters
    ----------
    bam : path
        Path to the mounted bam file.
    tries : int
        Number of times to try to unmount bam file. Sometimes they don't unmount
        on the first time if they are being used or something.
    wait : int
        Number of seconds to wait between tries if the bam file doesn't unmount
        on the first try.

    """
    import glob
    import subprocess
    import time
    p = os.path.sep.join(bam.split(os.path.sep)[0:-3] + ['*'])
    files = glob.glob(p)
    mnt = os.path.sep.join(bam.split(os.path.sep)[0:-2])
    for f in files:
        if f != mnt:
            os.remove(f)

    t = 0
    while t < tries:
        try:
            subprocess.call('fusermount -u {}'.format(mnt),
                            shell=True)
            break
        except OSError:
            time.sleep(wait)
            t += 1

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
    import shutil
    import subprocess
    tcga_bam = mount_bam(analysis_id)
    # I'll make a tempdir to hold the mounted bam and any other temp files I
    # want to make.
    tempdir = os.path.sep.join(tcga_bam.split(os.path.sep)[0:-3])
    temp_bams = []
    for i in range(0, len(intervals), max_intervals):
        ints = ' '.join(intervals[i:i + max_intervals])
        temp_bam = os.path.join(tempdir, '{}_{}.bam'.format(analysis_id, i))
        subprocess.check_call('samtools view -b {} {} > {}'.format(tcga_bam, 
                                                                   ints,
                                                                   temp_bam),   
                              shell=True)
        temp_bams.append(temp_bam)
    if len(temp_bams) > 1:
        # subprocess.check_call(['samtools', 'merge', '-f', interval_bam, 
        #                        ' '.join(temp_bams)])
        command = 'samtools merge -f {} {}'.format(interval_bam, 
                                                   ' '.join(temp_bams))
        subprocess.check_call(command, shell=True)
                                                               
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
    # This class creates an "engine" that runs in the background and gets reads
    # for intervals from different bam files. The goal here was to make
    # something that could be monitored and controlled while running so that I
    # can gracefully exit a job that may take a long time. To do this, the
    # engine runs on a thread from the threading module. This thread can run in
    # the background but I can set an event to stop it. The thread cretes
    # several more threads using the multiprocessing module. Each of the
    # multiprocessing threads mounts a bam file, gets the reads from the
    # intervals, and unmounts the bam file. When the engine is killed, the
    # engine waits for the multiprocessing threads to finish, then the engine is
    # finally stopped. I may build in some ability to stop the multiprocessing
    # threads as they are running, but that might be hard because the call to
    # GTFuse is generally the thing that takes a long time, so I'd just have to
    # kill that in some way.  I also wanted to make the class simple so that it
    # could be subclassed for more complex analyses.  Subclassing will mainly be
    # accomplished using the bam_fnc parameter which allows execution of an
    # arbitrary function on each bam file that contains the reads from the
    # intervals as the bam files are created and the engine_fnc parameter which
    # executes an arbitrary function every time the engine cycles.

    def __init__(self, analysis_ids, bed, bam_outdir='.', threads=10,
                 sleeptime=10, bam_fnc=None, engine_fnc=None):
        """
        Initialize engine for obtaining reads for given intervals/IDs
        
        Parameters
        ----------
        analysis_ids : list
            List of TCGA analysis IDs.
        bed : path
            Bed file with intervals.
        bam_outdir : directory
            Directory to write bam files to.
        threads : int
            Number of different processes/threads to use.
        sleeptime : int
            Number of seconds to sleep between engine updates.
        bam_fnc : function
            Function to apply to each bam file (path) after the bam file is
            made. For instance, you could make a function that uses the bam file
            path as input for variant calling, etc. This function is called as
            soon as the engine knows that the bam file is created.
        engine_fnc : function
            Function to execute every time the engine cycles (aka every
            sleeptime number of seconds while the engine is running). 

        """
        import threading
        assert len(analysis_ids) > 0
        assert os.path.exists(bed)
        self.analysis_ids = analysis_ids
        self.bed = bed
        self.bam_outdir = bam_outdir
        self.threads = threads
        self.sleeptime = sleeptime
        self.bam_fnc = bam_fnc
        self.engine_fnc = engine_fnc
        # Intervals in the format 1:20-200 etc. as a list.
        self.intervals = bed_to_samtools_intervals(self.bed)
        # analysis_ids that the engine has started getting reads for.
        self.started = []
        # analysis_ids that the engine has not started getting reads for.
        self.remaining = analysis_ids
        # analysis_ids that have finished.
        self.finished = []
        self.current_procs = []
        self.old_procs = []
        # Paths to bam files with reads.
        self.bams = []
        # Queue used when multiprocessing.
        self.queue = None
        # We set this event when we want to stop the engine.
        self.stop_event = threading.Event()
        self.running = False
        self.engine_thread = None
        self.start_engine()

    def remove_finished_procs(self):
        import inspect
        import types
        for p in self.current_procs:
            # If we find a dead process, there should be a result in the queue.
            # The result will not necessarily be from that dead process though.
            if not p.is_alive():
                p.join()
                bam = self.queue.get()
                self.bams.append(bam)
                analysis_id = os.path.splitext(os.path.split(bam)[1])[0]
                self.finished.append(analysis_id)
                self.current_procs.remove(p)
                self.old_procs.append(p)
                if (type(self.bam_fnc) == types.FunctionType or
                    inspect.ismethod(self.bam_fnc)):
                    self.bam_fnc(bam)

    def get_reads(self, analysis_id, bam):
        bam = reads_from_intervals(analysis_id, self.intervals, bam,
                                   max_intervals=10)
        self.queue.put(bam)

    def new_proc(self):
        import multiprocessing
        import time
        if len(self.remaining) > 0:
            analysis_id = self.remaining.pop()
            bam = '{}.bam'.format(os.path.join(self.bam_outdir, analysis_id))
            p = multiprocessing.Process(target=self.get_reads,
                                        args=[analysis_id, bam])
            self.current_procs.append(p)
            p.start()
            self.started.append(analysis_id)
            # Wait in case we are going to mount multiple bams. Mounting too
            # many bams too close in time might make problems?
            time.sleep(5)

    def add_procs(self):
        while (len(self.current_procs) < self.threads and 
               len(self.remaining) > 0):
            self.new_proc()

    def stop_engine(self):
        self.stop_event.set()

    def start_engine(self):
        import threading
        t = threading.Thread(target=self.read_interval_worker)
        self.engine_thread = t
        t.start()
        self.running = True

    def read_interval_worker(self):
        import inspect
        import multiprocessing
        import sys
        import time
        import types

        self.queue = multiprocessing.Queue()
        while not self.stop_event.is_set():
            if len(self.remaining) == 0:
                self.stop_engine()
            else:
                self.remove_finished_procs()
                self.add_procs()
                if (type(self.engine_fnc) == types.FunctionType or
                    inspect.ismethod(self.engine_fnc)):
                    self.engine_fnc()
            self.stop_event.wait(self.sleeptime)
        # If we get here, we aren't running any more. Wait for processes to
        # finish, then exit.
        sys.stderr.write('Engine stopping, waiting for jobs to conclude.\n')
        while len(self.current_procs) > 0:
            self.remove_finished_procs()
            if (type(self.engine_fnc) == types.FunctionType or
                inspect.ismethod(self.engine_fnc)):
                self.engine_fnc()
            time.sleep(self.sleeptime)
        sys.stderr.write('Jobs concluded, engine stopped.\n')
        self.running = False

class FLCVariantCallingEngine(ReadsFromIntervalsEngine):
    # This class extends the ReadsFromIntervalsEngine to both get the bam file
    # and call variants. However, this class is a little weird because it is
    # specifically set up to work on the Frazer lab cluster system. Here, the
    # variant calling is performed on a different server (flc) that has access
    # to the same filesystem as the server running the engine (hence the FLC
    # part of the name). This is specific to my use case but could easily
    # modified for other use cases.
    
    # TODO
    # Some way to watch for jobs that finish. Maybe I can just watch to see when
    # the output file is copied somewhere, then process it further? Or just not
    # worry after I submit? It would be nice to know even if just for the
    # purpose of monitoring progress. I would almost need another engine for
    # this, and that might be fine. I can have a monitor engine that is separate
    # of the GTFuse engine.  Something to keep track of which analysis_ids and
    # which intervals have been completed.

    def __init__(self, tumor_normal_ids, bed, java, mutect, fasta, dbsnp,
                 cosmic, external_server='flc.ucsd.edu', variant_outdir='.',
                 bam_outdir='.', threads=10, sleeptime=10, bam_fnc=None):
        """
        Initialize engine for obtaining reads for given intervals/IDs
        
        Parameters
        ----------
        tumor_normal_ids : dict
            Dict mapping whose keys are tumor analysis_ids and whose values are
            the corresponding analysis_id for the normal.
        bed : path
            Bed file with intervals.
        java : path
            Path to java installation.
        mutect : path
            Path to mutect jar.
        fasta : path
            Path to genome fasta.
        dbsnp : path
            Path to dbsnp vcf.
        cosmic : path
            Path to cosmic vcf.
        external_server : hostname
            Hostname of external server. Username assumed to be the same.
            Assumed to have access through key.
        variant_outdir : directory
            Directory to write variant calling results to.
        bam_outdir : directory
            Directory to write bam files to.
        threads : int
            Number of different processes/threads to use.
        sleeptime : int
            Number of seconds to sleep between engine updates.
        bam_fnc : function
            Function to apply to each bam file (path) after the bam file is
            made. For instance, you could make a function that uses the bam file
            path as input for variant calling, etc. This function is called as
            soon as the engine knows that the bam file is created.

        """
        self.tumor_normal_ids = tumor_normal_ids
        analysis_ids = self.get_id_list()
        # Paths to various files needed for mutect.
        self.java = java
        self.mutect = mutect
        self.fasta = fasta
        self.dbsnp = dbsnp
        self.cosmic = cosmic
        # External server to submit jobs to.
        self.external_server = external_server
        # Directory to store variant calling results.
        self.variant_outdir = variant_outdir
        # analysis_ids that we have begun calling variants for.
        self.variant_calling_started = []
        ReadsFromIntervalsEngine.__init__(self, analysis_ids, bed,
                                          bam_outdir=bam_outdir,
                                          threads=threads, sleeptime=sleeptime,
                                          engine_fnc=self.variant_calling_worker
                                         )

    def get_id_list(self):
        """Make list of ids where the tumor and normal ids as defined by
        self.tumor_normal_ids are next to each other in the list."""
        out = []
        for k in self.tumor_normal_ids.keys():
            out.append(k)
            out.append(self.tumor_normal_ids[k])
        return out

    def variant_calling_worker(self):
        import sys
        for t in ((set(self.finished) - set(self.variant_calling_started)) &
                  set(self.tumor_normal_ids.keys())):
            n = self.tumor_normal_ids[t]
            if n in self.finished:
                self.call_variants(t, n)
                self.variant_calling_started.append(t)
                self.variant_calling_started.append(n)

    def call_variants(self, tumor, normal):
        pbs = self.write_pbs_script(tumor, normal)
        self.submit_pbs_script(pbs)

    def submit_pbs_script(self, pbs):
        import subprocess
        subprocess.check_call(['ssh', self.external_server, 'qsub', pbs])

    def write_pbs_script(self, tumor, normal):
        tumor_bam = '{}.bam'.format(os.path.join(self.bam_outdir, tumor))
        normal_bam = '{}.bam'.format(os.path.join(self.bam_outdir, normal))
        # Make directory to store results if it doesn't already exist.
        analysis_dir = os.path.join(self.variant_outdir,
                                    '{}_{}'.format(tumor, normal))
        try:
            os.makedirs(analysis_dir)
        except OSError:
            pass
        bed_name = os.path.splitext(os.path.split(self.bed)[1])[0]
        out = os.path.join(analysis_dir, '{}_variants.txt'.format(bed_name))
        wig = os.path.join(analysis_dir, '{}.wig'.format(bed_name))
        pbs = os.path.join(analysis_dir, '{}_mutect.pbs'.format(bed_name))
        pbs_out = '{}_mutect.out'.format(os.path.join(analysis_dir, bed_name))
        pbs_err = '{}_mutect.err'.format(os.path.join(analysis_dir, bed_name))
        pbs_tempdir = '/scratch/{}_{}_{}'.format(tumor, normal, bed_name)
        pbs_temp_tumor = os.path.join(pbs_tempdir, '{}.bam'.format(normal))
        pbs_temp_normal = os.path.join(pbs_tempdir, '{}.bam'.format(tumor))
        pbs_tumor_sorted = os.path.join(pbs_tempdir, 
                                        '{}_sorted.bam'.format(normal))
        pbs_normal_sorted = os.path.join(pbs_tempdir, 
                                         '{}_sorted.bam'.format(tumor))
        with open(pbs, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('\n'.join(['#PBS -q high', 
                               '#PBS -N {}_{}_{}'.format(tumor, normal,
                                                         bed_name),
                               '#PBS -l nodes=1:ppn=8',
                               '#PBS -o {}'.format(pbs_out),
                               '#PBS -e {}'.format(pbs_err)]) + '\n\n')
            f.write('\n'.join(['mkdir -p {}'.format(pbs_tempdir),
                               'cd {}'.format(pbs_tempdir)]) + '\n\n')
            # Copy bam files to scratch.
            f.write('\n'.join(['rsync -avz {} {}'.format(tumor_bam, 
                                                         pbs_temp_tumor),
                               'rsync -avz {} {}'.format(normal_bam, 
                                                         pbs_temp_normal)]) 
                    + '\n\n')
            f.write('samtools sort -o {} tempt > {}\n'.format(pbs_temp_tumor,
                                                              pbs_tumor_sorted))
            f.write('samtools sort -o {} tempn > {}\n\n'.format(pbs_temp_normal,
                                                                pbs_normal_sorted))
            f.write('samtools index {}\n'.format(pbs_tumor_sorted))
            f.write('samtools index {}\n\n'.format(pbs_normal_sorted))
            # Run mutect.
            f.write('\\\n'.join(['{} -jar {}'.format(self.java, self.mutect),
                                 '\t-T MuTect',
                                 '\t-L {}'.format(self.bed),
                                 '\t-R {}'.format(self.fasta),
                                 '\t--dbsnp {}'.format(self.dbsnp),
                                 '\t--cosmic {}'.format(self.cosmic),
                                 '\t-I:normal {}'.format(pbs_normal_sorted),
                                 '\t-I:tumor {}'.format(pbs_tumor_sorted),
                                 '\t--out out.txt',
                                 '\t--coverage_file out.wig']) + '\n\n')
            # Copy results to analysis_dir.
            f.write('\n'.join(['rsync -avz out.txt {}'.format(out),
                               'rsync -avz out.wig {}'.format(wig)]) + '\n\n')
            # Remove bam files and temporary directory on scratch.
            f.write('rm -r {} {} {}\n'.format(tumor_bam, normal_bam, 
                                              pbs_tempdir))
        return pbs
