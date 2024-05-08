#!/bin/env python
# USAGE: python3.8 scripts/submit_batch.py $PWD "[3792-3796]" --outdir cosmics_1stset_261120 --dry-run
# USAGE 2.0: python3.8 scripts/submit_batch.py $PWD "[4455-4456]" --outdir test_batch --mh 8 --nev 200 24 --dry-run
# USAGE with resubmission and check in runlog: python scripts/submit_batch.py  $PWD "[3010-4100]" --outdir prod-autumn22-18112022 -q cygno --nthreads 8 -r --runlog pedestals/runlog_LNGS.csv --dry-run
# reading run list from a filem: python3.8 scripts/submit_batch.py $PWD datasets/lime_April2021.txt --outdir test_batch --mh 8 --nev 200 24 --dry-run
# remove --dry-run to submit for real (otherwise only the scripts are created and commands are printed)

jobstring  = '''#!/bin/bash
ulimit -c 0 -S
ulimit -c 0 -H
set -e
cd CYGNOBASE
source scripts/activate_cygno_lngs.sh
RECOSTRING
'''

import os, sys, re, csv, subprocess

def prepare_jobpack(jobdir,logdir,workdir,cmd,maxresub,resub=False,ijob=0):
    job_file_name = jobdir+'/job{ij}_run{r}.sh'.format(ij=ijob,r=run)
    log_file_name = logdir+'/job{ij}_run{r}.log'.format(ij=ijob,r=run)
    if resub:
        ir=0
        resub_log_file_name = log_file_name.replace('.log','_resub{ir}.log'.format(ir=ir))
        while os.path.isfile(resub_log_file_name):
            resub_log_file_name = log_file_name.replace('.log','_resub{ir}.log'.format(ir=ir))
            ir+=1
        log_file_name = resub_log_file_name
        if ir>maxresub:
            maxresub=ir
    tmp_file = open(job_file_name, 'w')
    tmp_filecont = jobstring

    tmp_filecont = tmp_filecont.replace('RECOSTRING',cmd)
    tmp_filecont = tmp_filecont.replace('CYGNOBASE',workdir+'/')
    tmp_file.write(tmp_filecont)
    tmp_file.close()
    sub_cmd = 'qsub -q {queue} -l {ssd}ncpus={nt},mem={ram}mb -d {dpath} -e localhost:{logf} -o localhost:{logf} -j oe {jobf}'.format(dpath=workdir,logf=log_file_name,jobf=job_file_name,nt=nThreads,ram=RAM,ssd=ssdcache_opt,queue=options.queue)
    return sub_cmd,maxresub

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog workdir runs [options] ')
    parser.add_option(        '--dry-run',  dest='dryRun',   action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option(        '--outdir',   dest='outdir',   type="string", default=None, help='outdirectory');
    parser.add_option(        '--nthreads', dest='nthreads', type="string", default=None, help='number of threads / job. If not given, it will be decided based on the queue');
    parser.add_option('-q',   '--queue',    dest='queue',    type="string", default='cygno-custom', help='queue to be used for the jobs');
    parser.add_option('--mh'  '--max-hours',dest='maxHours', default=-1, type='float', help='Kill a subprocess if hanging for more than given number of hours.')
    parser.add_option('--nev' '--event-chunks',dest='eventChunks', default=[], nargs=2,  type='int', help='T C: Total number of events to process and events per job')
    parser.add_option('--cfg' '--config-file',dest='configFile', default="configFile_LNGS.txt",  type='string', help='the config file to be run')
    parser.add_option('-t',   '--tmp',dest='tmpdir', default="/tmp/",  type='string', help='the input directory')
    parser.add_option('-r',   '--resubmit',dest='resubmit', action='store_true', default=False, help='check the existence of the output in the outdir and resubmit the requested run range')
    parser.add_option(        '--pedestals-only',dest='pedOnly', action='store_true', default=False, help='run only on pedestal runs')
    parser.add_option(        '--runlog', dest='runlog', default=None, type='string', help='if given, check in the runlog that the runs exists before creating the job')
    (options, args) = parser.parse_args()

    if len(args)<2:
        parser.print_help()
        exit(1)

    abswpath  = os.path.abspath(args[0]) 
    if not os.path.isdir(abswpath):
        raise RuntimeError('ERROR: {p} is not a valid directory. This is the base dir where the jobs run'.format(p=abswpath))

    runr = args[1]
    runs = []
    if os.path.exists(runr):
        print("Opening the text file %s with run ranges to process "%runr)
        f = open(runr, "r")
        runs.append(exec(f.read()))
    # if a simple run range with the format [minr-maxr] is given
    else:
        p = re.compile('\[(\d+)-(\d+)\]')
        m = p.match(runr)
        if m:
            minr=m.group(1); maxr=m.group(2)
            if not str.isdigit(minr) or not str.isdigit(maxr):
                raise RuntimeError ("not a good formatted run range: {rr}. Aborting.".format(rr=runr))
            minr=int(minr); maxr=int(maxr)
            if minr>maxr:
                raise RuntimeError ("Warning: first run > end run: {rr}. Aborting.".format(rr=runr))
            runs = list(range(minr,maxr+1))
        else:
            raise RuntimeError ("not a good formatted run range string: {rr}. Aborting.".format(rr=runr))
    print ("SUBMITTING JOBS ON RUN RANGE = ",runs)

    if not options.outdir:
        raise RuntimeError('ERROR: give at least an output directory. This is where the jobs configs and logs are stored')
    else:
        absopath  = os.path.abspath(options.outdir)
        if not os.path.isdir(absopath):
            print ('making a directory ',absopath,' and running in it')
            os.system('mkdir -p {od}'.format(od=absopath))

    jobdir = absopath+'/jobs/'
    if not os.path.isdir(jobdir):
        os.system('mkdir {od}'.format(od=jobdir))
    logdir = absopath+'/logs/'
    if not os.path.isdir(logdir):
        os.system('mkdir {od}'.format(od=logdir))

    # typically for LIME images (~4MP) MEM ~1GB/core
    if options.nthreads:
        nThreads = int(options.nthreads)
        RAM = min(max(nThreads * 1200, 4000),9000)
    else:
        if options.queue=='cygno':
            nThreads = 7
            RAM = nThreads * 1250 # max is 9 GB
        elif options.queue=='cygno-custom':
            nThreads = 64
            RAM = nThreads * 1500 # max is 110 GB
        else:
            print("WARNING: queue ",options.queue," not foreseen. Using 4 threads and 1 GB/thread by default")
            nThreads = 4
            RAM = nThreads * 1000
    ssdcache_opt = 'nodes=1:disk10,' if options.queue=='cygno-custom' else ''
    # cygno-custom mounts /mnt/ssdcache, not the other queues. It seems that sometimes testing its presence works, sometimes not,
    # so when not in cygno-custom, force the usage to /tmp
    #tmpdir_opt = '' if  options.queue=='cygno-custom' else ' --tmp /tmp/'
    # IT SEEMS THAT /mnt/ssdcache is not mounted even in cygno-custom
    #tmpdir_opt = ' --tmp /tmp/'
    maxtime_opt = '' if options.maxHours < 0 else '--max-hours {hr}'.format(hr=options.maxHours)
    commands = []
    maxresub=-1
    for run in runs:
        if options.runlog:
            existing_runs=[]
            with open(options.runlog,"r") as csvfile:
                csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
                next(csvreader) # skips header
                existing_runs_list = [(int(row[0]),row[1]) for row in list(csvreader)]
            existing_runs = dict(existing_runs_list)
            if run not in existing_runs.keys():
                print ("\t=> Run %d not in %s runlog, so skipping it" % (run,options.runlog))
                continue
            if options.pedOnly and "PED" not in existing_runs[run]:
                continue
        if options.resubmit:
            recofile = '%s/reco_run%05d_3D.root' % (absopath,run)
            if os.path.exists(recofile): continue
            else: print ("Run %d has no good output file. Resubmitting it" %run )
        if len(options.eventChunks)>0:
            (totEv,evPerJob) = options.eventChunks
            print ("Preparing jobs for run {r}. The task subdivides a total of {nT} events in chunks of {nJ} events per job.".format(r=run,nT=totEv,nJ=evPerJob))
            for ij,firstEvent in enumerate(range(0,totEv,evPerJob)):
                print ("Will submit job #{ij}, processing event range: [{fev}-{lev}]".format(ij=ij,fev=firstEvent,lev=min(firstEvent+evPerJob,totEv)))
                cmd = 'python3.8 reconstruction.py {cfg} -r {r} -o reco_job{ijob} --first-event {fev} --max-entries {me} -j {nt} -t {tmpopt}  {maxtimeopt} -d {outdiropt}'.format(r=run,nt=nThreads,tmpopt=options.tmpdir,maxtimeopt=maxtime_opt,fev=firstEvent,me=evPerJob,ijob=ij,cfg=options.configFile,outdiropt=options.outdir)
                sub_cmd,maxresub = prepare_jobpack(jobdir,logdir,abswpath,cmd,maxresub,options.resubmit,ij)
                commands.append(sub_cmd)
        else:
            cmd = 'python3.8 reconstruction.py {cfg} -r {r} -j {nt} -t {tmpopt} {maxtimeopt} -d {outdiropt}'.format(r=run,nt=nThreads,tmpopt=options.tmpdir,maxtimeopt=maxtime_opt,cfg=options.configFile,outdiropt=options.outdir)
            sub_cmd,maxresub = prepare_jobpack(jobdir,logdir,abswpath,cmd,maxresub,options.resubmit)
            commands.append(sub_cmd)

    maxq=999999
    queue_properties = subprocess.run(['qstat -Qf {q}'.format(q=options.queue)], stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8')
    for row in iter(queue_properties.splitlines()):
        result = re.search(r"\s+max_user_queuable = (\d+)",row)
        if result:
            maxq = int(result.group(1))
            break

    if(len(commands)>maxq):
        raise RuntimeError("Requested to run %d jobs, byt maximum queuable per user in queue '%s' is %d" % (len(commands),options.queue,maxq))
        
    if options.dryRun:
        for c in commands:
            print (c)
    else:
        for c in commands:
            os.system(c)
        
    print ("==================== SUMMARY ==================================")
    print ("\t Requested to run over %d runs" % len(runs))
    print ("\t after filtering will submit %d jobs" % len(commands))
    if maxresub>-1:
        print ("\t this is resubmission # %d of this run range" % maxresub)
    print ("===============================================================")

    
    
    print ("DONE")


        
