import logging
import os

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
    f.close()
    return intervals

class GTFuseBam:
    def __init__(self, analysis_id, mountpoint='.',
                 cache='./fusecache'):
        """
        Parameters
        ----------
        analysis_id : str
            TCGA analysis_id for sample of interest
        mountpoint : directory
            Directory to make temporary directory and files.
        cache : directory
            Directory to store GTFuse cache files.
        
        """
        import pwd
        username = pwd.getpwuid(os.getuid()).pw_name
        self.analysis_id = analysis_id
        self.mountpoint = mountpoint
        self.cache = cache
        self._make_tempdir()
        self.bam = ''
        self.mounted = False
        self.mount()

    def _make_tempdir(self):
        """Make a temp directory to hold the mounted bam file and allow us to 
        store other temp files."""
        self.tempdir = os.path.join(self.mountpoint,
                                    self.analysis_id)
        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)

    def mount(self):
        """
        Mount bam file for given analysis ID using GTFuse
        
        Creates a directory for this analysis ID in "mount" and mounts 
        makes mnt and cache directories then mounts the bam file using 
        GTFuse. Returns bam file path.

        """
        import glob
        import subprocess
        import time
        mnt = os.path.join(self.tempdir, 'mnt')
        # This mnt directory shouldn't already exist because the bam file
        # shouldn't already be mounted.
        os.makedirs(mnt)
        if not os.path.exists(self.cache):
            os.makedirs(self.cache)
        url = ('https://cghub.ucsc.edu/cghub/data'
               '/analysis/download/{}'.format(self.analysis_id))
        subprocess.check_call(['gtfuse', 
                               '--ssl-no-verify-ca',
                               '--cache-dir',
                               self.cache,
                               url,
                               mnt])
        files = glob.glob('{}/{}/*'.format(mnt, self.analysis_id))
        bam = [ x for x in files if os.path.splitext(x)[1] == '.bam' ][0]
        self.bam = os.path.realpath(bam)
        self.mounted = True
        
    def unmount(self, tries=100, sleeptime=10):
        """
        Unmount bam file mounted with GTFuse
        
        Parameters
        ----------
        tries : int
            Number of times to try to unmount bam file. Sometimes they don't
            unmount on the first time if they are being used or something.
        sleeptime : int
            Number of seconds to wait between tries if the bam file doesn't
            unmount on the first try.
    
        """
        import glob
        import shutil
        import subprocess
        import time
        p = os.path.sep.join(self.bam.split(os.path.sep)[0:-3] + ['*'])
        files = glob.glob(p)
        mnt = os.path.sep.join(self.bam.split(os.path.sep)[0:-2])
        # First delete other temp files that are not the mounted bam.
        for f in files:
            if f != mnt:
                os.remove(f)
        # Now clear cache.
        files = glob.glob(os.path.join(self.cache, self.analysis_id + '*'))
        for f in files:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)
        # Now unmount bam.
        t = 0
        while t < tries:
            try:
                subprocess.call('fusermount -u {}'.format(mnt),
                                shell=True)
                os.rmdir(mnt)
                os.rmdir(os.path.split(mnt)[0])
                break
            except OSError:
                time.sleep(sleeptime)
                t += 1
        self.mounted = False

class ReadsFromIntervalsBam:
    def __init__(self, gtfuse_bam, intervals, bam):
        """
        Get reads from analysis id for given intervals 
        
        Parameters
        ----------
        gtfuse_bam : GTFuseBam
            GTFuseBam object. The bam doesn't have to be mounted.
        intervals : list or path
            List of intervals of the form 1:10-200 or bed file.
        bam : str
            Path to the bam file where the reads for the intervals will be
            written.
        
        """
        self.gtfuse_bam = gtfuse_bam
        if type(intervals) == list:
            self.intervals = intervals
        else:
            self.intervals = bed_to_samtools_intervals(intervals)
        # List of intervals that we didn't obtain reads for. Hopefully this is
        # always empty.
        self.bad_intervals = []
        # Bam file with the reads from the intervals.
        self.bam = bam
        self.analysis_id = self.gtfuse_bam.analysis_id
        self.tempdir = self.gtfuse_bam.tempdir
        self.reads_from_intervals()
        
    def reads_from_intervals(self, max_intervals=10, tries=10, sleeptime=3):
        """
        Get reads from analysis id for given intervals 
        
        Parameters
        ----------
        max_intervals : int
            Maximum number of intervals to obtain with one samtools view call.
        tries : int
            Number of times to try to get a group of regions from the bam file. 
            Sometimes there is an error for some of the regions.
        sleeptime : int
            Number of seconds to wait between tries if the reads aren't obtaiend
            on the first try.
        
        """
        import copy
        import shutil
        import subprocess
        import time
        if not self.gtfuse_bam.mounted:
            self.gtfuse_bam.mount()
        temp_bams = []
        stderr = open(os.path.join(self.tempdir,
                                   '{}.err'.format(self.analysis_id)), 'w')
        for i in range(0, len(self.intervals), max_intervals):
            ints = ' '.join(self.intervals[i:i + max_intervals])
            temp_bam = os.path.join(
                self.tempdir, '{}_{}.bam'.format(self.analysis_id, i))
            c = 'samtools view -b {} {} > {}'.format(
                self.gtfuse_bam.bam, ints, temp_bam)
            t = 0
            while t < tries:
                try:
                    subprocess.check_call(c, shell=True, stderr=stderr)
                    break
                except subprocess.CalledProcessError:
                    t += 1
                    self.gtfuse_bam.unmount()
                    time.sleep(sleeptime)
                    self.gtfuse_bam.mount()
            if t == tries:
                self.bad_intervals.append(self.intervals[i:i + max_intervals])
            temp_bams.append(temp_bam)
        if len(temp_bams) > 1:
            c = 'samtools merge -f {} {}'.format(self.bam,
                                                 ' '.join(temp_bams))
            subprocess.check_call(c, shell=True, stderr=stderr)
            stderr.close()
            for b in temp_bams:
                os.remove(b)
            stderr.close()
        else:
            shutil.move(temp_bams[0], self.bam)
        os.remove(stderr.name)

class ReadsFromIntervalsEngine:
    # This class creates an "engine" that runs in the background and gets reads
    # for intervals from different bam files. The goal here was to make
    # something that could be monitored and controlled while running so that I
    # can gracefully exit a job that may take a long time. To do this, the
    # engine runs on a thread from the threading module. This thread can run in
    # the background but I can set an event to stop it. The thread creates
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

    def __init__(self, analysis_ids, bed, bam_outdir='.', bed_name=None,
                 threads=10, sleeptime=10, reads_started=[], reads_finished=[],
                 bam_fnc=None, engine_fnc=None):
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
        bed_name : str
            Name of the bed file containing the intervals. Used to name output
            files. If not provided, defaults to the file name of the bed file.
        threads : int
            Number of different processes/threads to use.
        sleeptime : int
            Number of seconds to sleep between engine updates.
        reads_started : list
            List of analysis IDs that we have already started getting reads for.
            This is used when restarting a job that a previous engine was
            working on.
        reads_finished : list
            List of analysis IDs that we have already obtained the reads for.
            This is used when restarting a job that a previous engine was
            working on.
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
        assert bam_fnc == None # Currently not implemented, may remove for now.
        self.engine_fnc = engine_fnc
        self.intervals = bed_to_samtools_intervals(bed)
        # Bed file name
        if not bed_name:
            self.bed_name = os.path.splitext(os.path.split(self.bed)[1])[0]
        else:
            self.bed_name = bed_name
        # analysis_ids that the engine has started getting reads for.
        self.reads_started = reads_started
        # analysis_ids that the engine has not started getting reads for.
        self.reads_remaining = analysis_ids
        # analysis_ids that have finished.
        self.reads_finished = reads_finished
        for a in self.reads_finished:
            self.reads_remaining.remove(a)
        self.processes = []
        # ReadsFromIntervalsBam objects.
        self.reads_from_intervals_bams = []
        # Queue used when multiprocessing.
        self.in_queue = None
        # We set this event when we want to stop the engine.
        self._stop_event = None
        # Directory to store mounted bam files.
        self.mountpoint = os.path.join(self.bam_outdir, 'tmp')
        # Directory to store GTFuse cache.
        self.fusecache = os.path.join(self.bam_outdir, 'tmp', 'fusecache')
        self.running = False
        # Process that the engine is running on.
        self.engine_thread = None
        self.start()

    def start(self):
        import threading
        # We set this event when we want to stop the engine.
        self._stop_event = threading.Event()
        t = threading.Thread(target=self._reads_from_intervals_parent)
        self.engine_thread = t
        t.start()
        self.running = True

    def stop(self):
        self._stop_event.set()

    def _reads_from_intervals_parent(self):
        import inspect
        import multiprocessing
        import Queue
        import sys
        import time
        import types

        self.in_queue = multiprocessing.JoinableQueue()
        self.out_queue = multiprocessing.Queue()
        for aid in self.analysis_ids:
            self.in_queue.put(aid)
        for i in xrange(self.threads):
            self.add_process()
            self.in_queue.put('STOP')

        while sum([p.is_alive() for p in self.processes]) > 0:
            while True:
                try:
                    bam = self.out_queue.get(timeout=self.sleeptime)
                    self.reads_finished.append(bam.analysis_id)
                except:
                    Queue.Empty:
                        break

            if (type(self.engine_fnc) == types.FunctionType or
                inspect.ismethod(self.engine_fnc)):
                self.engine_fnc()
        
        self.running = False

    def _reads_from_intervals_worker(self, in_queue, out_queue):
        import sys
        analysis_id = in_queue.get()
        while analysis_id != 'STOP':
            bam_path = os.path.join(self.bam_outdir,
                                    '{}_{}.bam'.format(self.bed_name, 
                                                       analysis_id))
            gtfuse_bam = GTFuseBam(analysis_id, mountpoint=self.mountpoint, 
                                   cache=self.fusecache)
            bam = ReadsFromIntervalsBam(gtfuse_bam, self.intervals, bam_path)
            gtfuse_bam.unmount()
            self.out_queue.put(bam)
            in_queue.task_done()
            analysis_id = in_queue.get()

        in_queue.task_done()

    def get_reads(self, analysis_id):
        import sys
        bam_path = os.path.join(self.bam_outdir,
                                '{}_{}.bam'.format(self.bed_name, analysis_id))
        logging.warning('Mounting bam for {}\n'.format(analysis_id))
        gtfuse_bam = GTFuseBam(analysis_id, mountpoint=self.mountpoint, 
                               cache=self.fusecache)
        logging.warning('Getting reads for {}\n'.format(analysis_id))
        bam = ReadsFromIntervalsBam(gtfuse_bam, self.intervals, bam_path)
        logging.warning('Unmounting bam for {}\n'.format(analysis_id))
        gtfuse_bam.unmount()
        logging.warning('Putting results in queue for {}\n'.format(analysis_id))
        self.queue.put(bam)
        logging.warning('Finished get_reads for {}\n'.format(analysis_id))

    def add_process(self):
        import multiprocessing
        p = multiprocessing.Process(target=self._reads_from_intervals_worker,
                                    args=[self.in_queue, self.out_queue]) 
        p.daemon = True
        self.processes.append(p)
        p.start()

class TumorNormalVariantCall:
    def __init__(self, tumor_id, normal_id, intervals_name, bam_outdir=',',
                 variant_outdir='.'):
        self.tumor_id = tumor_id
        self.normal_id = normal_id
        self.intervals_name = intervals_name
        self.tumor_bam = os.path.join(bam_outdir, '{}_{}.bam'.format(
            self.intervals_name, self.tumor_id))
        self.normal_bam = os.path.join(bam_outdir, '{}_{}.bam'.format(
            self.intervals_name, self.normal_id))
        self.name = '{}_{}'.format(self.tumor_id, self.normal_id)
        self.variant_dir = os.path.join(variant_outdir, self.name)
        self.variants = os.path.join(self.variant_dir, 
                                     '{}_variants.txt'.format(self.name))
        self.wig = os.path.join(self.variant_dir, '{}.wig'.format(self.name))
        self.pbs = os.path.join(self.variant_dir, 
                                '{}_mutect.pbs'.format(self.name))
        self.stdout = '{}_mutect.out'.format(os.path.join(self.variant_dir, 
                                                          self.name))
        self.stderr = '{}_mutect.err'.format(os.path.join(self.variant_dir, 
                                                          self.name))

class FLCVariantCallingEngine(ReadsFromIntervalsEngine):
    # This class extends the ReadsFromIntervalsEngine to both get the bam file
    # and call variants. However, this class is a little weird because it is
    # specifically set up to work on the Frazer lab cluster system. Here, the
    # variant calling is performed on a different server (flc) that has access
    # to the same filesystem as the server running the engine (hence the FLC
    # part of the name). This is specific to my use case but could easily
    # modified for other use cases.

    def __init__(self, tumor_normal_ids, bed, java, mutect, fasta, dbsnp,
                 cosmic, name=None, external_server='flc.ucsd.edu',
                 variant_outdir='.', bam_outdir='.', threads=10, sleeptime=10,
                 bam_fnc=None, variant_engine_fnc=None):
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
        name : str
            Name of the intervals that we are calling variants for with this
            engine. Used to name output files. Defaults to bed file name if not
            provided as argument.
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
        variant_engine_fnc : function
            Function to execute every time the engine cycles (aka every
            sleeptime number of seconds while the engine is running). 

        """
        self.tumor_normal_ids = tumor_normal_ids
        self.bed = bed
        self.analysis_ids = self._get_id_list()
        # Paths to various files needed for mutect.
        self.java = java
        self.mutect = mutect
        self.fasta = fasta
        self.dbsnp = dbsnp
        self.cosmic = cosmic
        if not name:
            self.name = os.path.splitext(os.path.split(bed)[1])[0]
        else:
            self.name = name
        # The ReadsFromIntervalsEngine has engine_fnc which I set here as a
        # function that does variant calling. However, variant_engine_fnc allows
        # for another function to run inside of engine_fnc. This allows for more
        # things to be chained together. It's turtles all the way down!
        self.variant_engine_fnc = variant_engine_fnc
        # External server to submit jobs to.
        self.external_server = external_server
        # Directory to store variant calling results.
        self.variant_outdir = variant_outdir
        self.bam_outdir = bam_outdir
        # TumorNormalVariantCall objects that we need to call variants for.
        self._init_vcs()
        # TumorNormalVariantCalls that we have begun calling variants for.
        self.variant_calling_started = []
        # Directory that holds information about this variant calling run.
        self.infodir = os.path.join(bam_outdir,
                                    '{}_variant_calling'.format(self.name))
        self.html_status = os.path.join(
            self.infodir, '{}_status.html'.format(self.name)
        )
        self.running = False
        self._setup()
        ReadsFromIntervalsEngine.__init__(
            self, self.analysis_ids, self.bed, bam_outdir=bam_outdir,
            bed_name=self.name, threads=threads, sleeptime=sleeptime,
            reads_started=self.reads_started,
            reads_finished=self.reads_finished,
            engine_fnc=self._variant_calling_worker
        )

    def _init_vcs(self):
        self.tumor_normal_variant_calls = []
        for t in self.tumor_normal_ids.keys():
            n = self.tumor_normal_ids[t]
            vc = TumorNormalVariantCall(t, n,
                                        self.name,
                                        bam_outdir=self.bam_outdir,
                                        variant_outdir=self.variant_outdir)
            self.tumor_normal_variant_calls.append(vc)
    
    def _setup(self):
        """
        Check whether a varaiant calling job with this name has been run before.
        If so, check to make sure we have the same intervals here and pick up
        where we left off. If this is the first time we've run variant calling
        for these intervals, make a directory to hold some information about
        this variant calling run and populate it.
        """
        self.reads_started = []
        self.reads_finished = []
        if os.path.exists(self.html_status):
            self._exist_setup()
        else:
            self._not_exist_setup()

    def _exist_setup(self):
        """Set up the engine given that an engine has already worked on these
        samples and intervals in the past"""
        # Update analysis ids based on which samples have already been
        # completed.
        import pandas as pd
        df = pd.read_html(self.html_status)[0]
        for vc in self.tumor_normal_variant_calls:
            t = vc.tumor_id
            n = vc.normal_id
            ind = vc.name
            if df.ix[ind, 'tumor reads'] == 'finished':
                self.reads_started.append(t)
                self.reads_finished.append(t)
            if df.ix[ind, 'normal reads'] == 'finished':
                self.reads_started.append(n)
                self.reads_finished.append(n)
            if df.ix[ind, 'variant calling'] == 'finished':
                self.variant_calling_started.append(vc)

    def _not_exist_setup(self):
        import pandas as pd
        import shutil
        os.makedirs(self.infodir)
        # Copy the bed file into the info directory for posterity.
        shutil.copy(self.bed, 
                    os.path.join(self.infodir, '{}.bed'.format(self.name)))
        self._make_html_status()

    def _make_html_status(self):
        import pandas as pd
        index = []
        tumors = [x.tumor_id for x in self.tumor_normal_variant_calls]
        normals = [x.normal_id for x in self.tumor_normal_variant_calls]
        for vc in self.tumor_normal_variant_calls:
            index.append(vc.name)
        columns = ['tumor reads', 'normal reads', 'variant calling']
        df = pd.DataFrame(index=index, columns=columns)
        df.to_html(self.html_status, na_rep='')
        self._update_html_status()

    def _update_html_status(self):
        import pandas as pd
        df = pd.read_html(self.html_status,
                          index_col=0, header=0)[0]
        pairs_finished = 0
        for vc in self.tumor_normal_variant_calls:
            t = vc.tumor_id
            n = vc.normal_id
            ind = vc.name
            if t in self.reads_finished:
                df.ix[ind, 'tumor reads'] = 'finished'
            if n in self.reads_finished:
                df.ix[ind, 'normal reads'] = 'finished'
            if os.path.exists(vc.variants):
                df.ix[ind, 'variant calling'] = 'finished'
                pairs_finished += 1
        df.to_html(self.html_status, na_rep='')
        f = open(self.html_status, 'r')
        lines = f.readlines()
        f.close()
        f = open(self.html_status, 'w')
        head = '<head>\n<title>{} status</title></head>\n'.format(self.name)
        header = ('<header>\n<h1>' + 
                  'Status of variant calling for "{}"'.format(self.name) + 
                  '</h1>\n</header>\n')
        r = ['no', 'yes'][self.running]
        para = ['<p>Currently running: {}<br>'.format(r),
                'Pairs finished: {:,}<br>'.format(pairs_finished),
                'Pairs remaining: {:,}</p>'.format(len(self.tumor_normal_ids) - 
                                                   pairs_finished)]
        lines = [head, header, '\n'.join(para) + '\n\n'] + lines
        f.write(''.join(lines))
        f.close()

    def _get_id_list(self):
        """Make list of ids where the tumor and normal ids as defined by
        self.tumor_normal_ids are next to each other in the list."""
        out = []
        for k in self.tumor_normal_ids.keys():
            out.append(k)
            out.append(self.tumor_normal_ids[k])
        return out

    def _variant_calling_worker(self):
        import inspect
        import pandas as pd
        import sys
        import time
        import types
        vc_started = set([ x.tumor_id for x in self.variant_calling_started])
        for vc in (set(self.tumor_normal_variant_calls) -
                   set(self.variant_calling_started)):
            t = vc.tumor_id
            n = vc.normal_id
            if (t in self.reads_finished) and (n in self.reads_finished):
                self._call_variants(vc)
                self.variant_calling_started.append(vc)
        if (type(self.variant_engine_fnc) == types.FunctionType or
            inspect.ismethod(self.variant_engine_fnc)):
            self.variant_engine_fnc()
        self._update_html_status()
        # If the engine is done, wait until all variant calls are done.
        if len(self.processes) == 0 and self._stop_event.is_set():
            df = pd.read_html(self.html_status,
                              index_col=0, header=0)[0]
            while set(df['variant calling']) != set(['finished']):
                time.sleep(self.sleeptime)    
                self._update_html_status()
                df = pd.read_html(self.html_status,
                                  index_col=0, header=0)[0]

    def _call_variants(self, vc):
        self._write_pbs_script(vc)
        self._submit_pbs_script(vc)

    def _submit_pbs_script(self, vc):
        import subprocess
        subprocess.check_call(['ssh', self.external_server, 'qsub', vc.pbs])

    def _write_pbs_script(self, vc):
        # Make directory to store results if it doesn't already exist.
        try:
            os.makedirs(vc.variant_dir)
        except OSError:
            pass
        pbs_tempdir = '/scratch/{}_{}_{}'.format(vc.tumor_id,
                                                 vc.normal_id,
                                                 vc.intervals_name)
        pbs_temp_tumor = os.path.join(pbs_tempdir,
                                      '{}.bam'.format(vc.normal_id))
        pbs_temp_normal = os.path.join(pbs_tempdir, 
                                       '{}.bam'.format(vc.tumor_id))
        pbs_tumor_sorted = os.path.join(pbs_tempdir, 
                                        '{}_sorted.bam'.format(vc.normal_id))
        pbs_normal_sorted = os.path.join(pbs_tempdir, 
                                         '{}_sorted.bam'.format(vc.tumor_id))
        with open(vc.pbs, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('\n'.join(['#PBS -q high', 
                               '#PBS -N {}_{}_{}'.format(vc.tumor_id, 
                                                         vc.normal_id,
                                                         vc.name),
                               '#PBS -l nodes=1:ppn=8',
                               '#PBS -o {}'.format(vc.stdout),
                               '#PBS -e {}'.format(vc.stderr)]) + '\n\n')
            f.write('\n'.join(['mkdir -p {}'.format(pbs_tempdir),
                               'cd {}'.format(pbs_tempdir)]) + '\n\n')
            # Copy bam files to scratch.
            f.write('\n'.join(['rsync -avz {} {}'.format(vc.tumor_bam, 
                                                         pbs_temp_tumor),
                               'rsync -avz {} {}'.format(vc.normal_bam, 
                                                         pbs_temp_normal)]) 
                    + '\n\n')
            f.write(' '.join(['samtools', 'sort', '-o', pbs_temp_tumor, 'tempt',
                              '>', pbs_tumor_sorted, '\n']))
            f.write(' '.join(['samtools', 'sort', '-o', pbs_temp_normal,
                              'tempn', '>', pbs_normal_sorted, '\n\n']))
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
            f.write('\n'.join(['rsync -avz out.txt {}'.format(vc.variants),
                               'rsync -avz out.wig {}'.format(vc.wig)]) 
                    + '\n\n')
            # Remove bam files and temporary directory on scratch.
            f.write('rm -r {} {} {}\n'.format(vc.tumor_bam, vc.normal_bam, 
                                              pbs_tempdir))
