from prometheus_client import start_http_server, Counter, Gauge
import os
import pathlib
import time
import math

def findProcesses(userid, executable):
    ret = []
    for p in pathlib.Path("/proc").glob("[0-9]*"):
        if p.is_dir() and p.owner() == userid:
            if (p / "exe").resolve().name == executable:
                ret.append(p)
    return ret

def gatherInfo(procpath, timeunits=0.01, pagesize=4096):
    fields = ((procpath / "stat").open()).readline().split()
    return dict(
        state = fields[2],
        time = sum(int(f) for f in fields[13:17])*timeunits,
        threads = int(fields[19]),
        vsize = int(fields[22]),
        rss = int(fields[23])*pagesize,
    )

class myGauge:
    def __init__(self, name, help):
        self._gauge = Gauge(name, help)
    def update(self, value):
        self._gauge.set(value)

class myCounter:
    def __init__(self, name, help):
        self._counter = Counter(name+"_total", help)
        self._value = 0
    def update(self, value):
        if value > self._value: # might decrease if some processes drop out
            self._counter.inc(value - self._value)
        self._value = value

def runMonitor(options):
    user = options.username if options.username else os.getlogin()
    timeunits = 1.0/os.sysconf(os.sysconf_names['SC_CLK_TCK'])
    pagesize = os.sysconf(os.sysconf_names['SC_PAGESIZE'])
    items = dict(
        time = myCounter("processing_cpu_seconds", "CPU seconds of processing"),
        threads = myGauge("processing_nthreads", "Number of threads"),
        vsize = myGauge("processing_vsize_bytes", "Virtual memory size"),
        rss = myGauge("processing_rss_bytes", "RSS memory size"),
    )
    print(f"Exporting process metrix for {options.progName} for user {user} on http://{options.address}:{options.port}/metrics, updated every {options.updateTime:.1f}s\n")
    start_http_server(options.port, addr=options.address)
    fastIters = int(math.ceil(options.scanTime/options.updateTime))
    started = False
    t0 = time.monotonic()
    while True:
        procs = findProcesses(user, options.progName)
        if not(procs):
            if started and time.monotonic() - t0 > options.timeOut:
                print(f"Exiting after {options.timeOut:.0f}s after last process is gone\n")
                break
            elif time.monotonic()- t0 > options.warmupTime:
                print(f"Exiting after {options.warmupTime:.0f}s without the first process appearing\n")
                break
        else:
            started = True
            t0 = time.monotonic()
        for _ in range(fastIters):
            allgood = True
            totals = dict((k,0) for k in items.keys())
            for p in procs:
                try:
                    info = gatherInfo(p, timeunits=timeunits, pagesize=pagesize)
                    for k in totals.keys():
                        totals[k] += info[k]
                except IOError:
                    allgood = False
                    pass # maybe process is gone
            for k,v in items.items():
                v.update(totals[k])
            time.sleep(options.updateTime)
            if not allgood:
                break # trigger a rescan anyway

def runBackgroundMonitor(options):
    if os.fork() != 0:
        return
    runMonitor(options)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="%prog")
    parser.add_option("-b", dest="background", action="store_true", default=False, help="Run in background")
    parser.add_option("-a", dest="address", type=str, default="127.0.0.1", help="Address to bind to (127.0.0.1)")
    parser.add_option("-p", dest="port", type=int, default=8084, help="Port to expose metrics on (8084)")
    parser.add_option("-P", dest="progName", type=str, default="cmsRun", help="Program to monitor for (defaults to cmsRun)")
    parser.add_option("-u", dest="username", type=str, default=None, help="Username to look for processes for (default to current user)")
    parser.add_option("-t", dest="updateTime", type=float, default=1.0, help="Frequency to update process info (in seconds)")
    parser.add_option("-T", dest="scanTime", type=float, default=15.0, help="Frequency to update scan for new processes (in seconds)")
    parser.add_option("-w", dest="warmupTime", type=float, default=300.0, help="Timeout to wait for the first process to appear (in seconds)")
    parser.add_option("-z", dest="timeOut", type=float, default=30.0, help="Timeout to exit if no more processes are found (in seconds)")
    (options,args) = parser.parse_args()
    if options.background:
        from multiprocessing import Process
        p = Process(target=runBackgroundMonitor, args=(options,))
        p.start()
        p.join()
    else:
        runMonitor(options)
