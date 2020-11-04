#!/usr/bin/env python
"""
Base class for calling obztak through the AUTOOBS named pipe mechanism.

Based on the Scheduler class by Eric Neilsen

TODO: 
- Replace ConfigParser with Yaml for consistency
"""
import os
from datetime import datetime, timedelta
import logging
import subprocess
from ConfigParser import ConfigParser
import json
import tempfile

from obztak.field import SISPI_DICT

class AutoObz(object):
    """ Automated obztak scheduler """

    def __init__(self, config_fname):
        self.configure(config_fname)

    def configure(self, config_fname):
        """Configure the auto object"""
        config = ConfigParser() 
        config.read(config_fname)

        self.stale_time_delta = timedelta(0,config.getfloat('timeouts', 'fifo'))
        
        self.output_fname = config.get('paths', 'outbox')
        self.queue_fname = config.get('paths', 'current_queue')
        self.previous_queue_fname = config.get('paths', 'previous_queue')
        self.in_progress_fname = config.get('paths', 'inprogress')
        self.fifo_fname = config.get('paths', 'fifo')  
        self.publish_script = config.get('publish','script')

        self.chunk = 10
        self.mode = None
        self.min_queue_len = 30
        self.min_queue_time = 70


    def make_script(self):
        """Create the observing script"""
        with open(self.queue_fname, 'r') as fp:
            sispi_queue = json.load(fp)

        logging.debug("Found %d exposure(s) on the SISPI/OCS queue." % len(sispi_queue))

        with open(self.in_progress_fname, 'r') as fp:
            in_progress = json.load(fp)

        for exp in in_progress:
            if exp is None: continue
            for k,v in SISPI_DICT.items():
                exp.setdefault(k,v)
            # Not great...
            exp['RA'] = 0.0
            exp['dec'] = 0.0
            exp['date'] = datetime.utcnow().strftime('%Y/%m/%d %H:%M:%S')

        tmp = tempfile.NamedTemporaryFile(suffix='.json',delete=False)
        json.dump(in_progress, tmp)
        tmp.flush()

        # Generic time for in-progress exposure
        delay = (90 + 30)*len(in_progress)
        # Sum the exposure times for the current queue
        delay += sum([q.get('expTime',0)+30 for q in sispi_queue])
        logging.debug("Total queue time: %g seconds"%delay)

        # If we don't want to add anything, return an empty list
        if len(sispi_queue) >= self.min_queue_len or delay/60. >= self.min_queue_time:
            logging.info("Queue contains %d exposures and a runtime of %d minutes; waiting..."
                         % (len(sispi_queue),delay//60))
            # Add an empty script so AUTOOBS knows the scheduler "passed"
            with open(self.output_fname, 'w') as fp:
                json.dump([], fp, indent=4)
            os.chmod(self.output_fname, 0o666)
            return

        logging.info(datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S'))
        start = datetime.utcnow() + timedelta(seconds=delay)
        utc = start.strftime('%Y-%m-%dT%H:%M:%S')
        params = dict(utc=utc,output=self.output_fname,chunk=self.chunk,
                      current=self.queue_fname,
                      previous=self.previous_queue_fname,
                      progress=tmp.name, mode=self.mode)

        # Schedule the next chunk of exposures
        cmd = "schedule_chunk -k %(chunk)i --utc %(utc)s -o %(output)s"%params
        cmd += " -c %(progress)s"%params # needs to be first
        cmd += " -c %(previous)s -c %(current)s"%params
        if self.mode: cmd += " -m %(mode)s"%params

        logging.info(cmd)

        # Generate the script
        logging.info("Calling scheduler")
        subprocess.check_call(cmd, shell=True)

    def publish(self):
        """Publish the current queue to the web"""
        logging.info("Publishing current queue...")
        script = self.publish_script
        cmd = 'ssh sispi@observer2 "%s -v $PWD/autoobs.conf"'%(script)
        try:
            raise Exception("Not publishing")
            #subprocess.call(cmd, shell=True)
        except:
            logging.warn("Publish failed.")

    def __call__(self):
        """Execute the loop to check the fifo"""
        logging.info("Scheduler starting")
        while True:
            # open block until something is sent to the fifo
            # (should by sent by obstac)
            logging.info("Waiting for AUTOOBS")
            with open(self.fifo_fname, 'r') as fp:
                time_string = fp.readline().strip()

            logging.info("Triggered by AUTOOBS")

            if len(time_string) == 0:
                continue
            
            try:
                queue_time = datetime.strptime(time_string,'%Y-%m-%d %H:%M:%S')
            except ValueError:
                logging.info("Invalid marker in FIFO: %s" % time_string)
                continue

            marker_age =  datetime.now()-queue_time
            if marker_age > self.stale_time_delta:
                logging.info("FIFO has time %s, more than %s ago; not calling scheduler"%
                            (time_string, str(self.stale_time_delta)))
                continue

            # Create the script
            new_sispi_script = self.make_script()

            # Publish the queue (even if no new script is created)
            self.publish()
