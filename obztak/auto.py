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

class AutoObz(object):
    """ Automated obztak scheduler """

    def __init__(self, config_fname):
        self.configure(config_fname)

    def configure(self, config_fname):
        config = ConfigParser() 
        config.read(config_fname)

        self.stale_time_delta = timedelta(0,config.getfloat('timeouts', 'fifo'))
        
        self.output_fname = config.get('paths', 'outbox')
        self.queue_fname = config.get('paths', 'current_queue')
        self.previous_queue_fname = config.get('paths', 'previous_queue')
        self.in_progress_fname = config.get('paths', 'inprogress')
        self.fifo_fname = config.get('paths', 'fifo')  

    def make_script(self):

        with open(self.queue_fname, 'r') as fp:
            sispi_queue = json.load(fp)

        logging.debug("Found %d exposure(s) on the SISPI/OCS queue." % len(sispi_queue))

        with open(self.in_progress_fname, 'r') as fp:
            in_progress = json.load(fp)

        # If we don't want to add anything, return an empty list
        if len(sispi_queue) >= self.min_queue_len:
            logging.info("Queue is already %d exposures long, not adding anything"
                         % len(sispi_queue))
            # Add an empty script so AUTOOBS knows the scheduler "passed"
            with open(self.output_fname, 'w') as fp:
                json.dump([], fp, indent=4)
            os.chmod(self.output_fname, 0o666)
            return

        # Generic time for in-progress exposure
        delay = (90 + 30)*len(in_progress)
        # Sum the exposure times for the current queue
        delay += sum([q['expTime']+30 for q in sispi_queue])
        logging.info("Total queue time: %g seconds"%delay)

        start = datetime.now() + timedelta(seconds=delay)
        utc = start.strftime('%Y-%m-%dT%H:%M:%S')
        params = dict(utc=utc,output=self.output_fname,complete=self.queue_fname,
                      previous=self.previous_queue_fname)
        logging.info("Calling scheduler")

        # Schedule the next chunk of exposures
        cmd = "schedule_chunk -k 15 --utc %(utc)s -o %(output)s -c %(complete)s -c %(previous)s"%params
        logging.info(cmd)

        subprocess.check_call(cmd, shell=True)

        with open(self.output_fname, 'r') as fp:
            output = json.load(fp)
            logging.info("Sending %d exposure(s) to the SISPI/OCS queue." % len(output))
        os.chmod(self.output_fname, 0o666)

    def __call__(self):
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
                queue_time = datetime.strptime(time_string, '%Y-%m-%d %H:%M:%S')
            except ValueError:
                logging.info("Invalid marker in FIFO: %s" % time_string)
                continue

            marker_age =  datetime.now()-queue_time
            if marker_age > self.stale_time_delta:
                logging.info("FIFO has time %s, more than %s ago; not calling scheduler"%
                            (time_string, str(self.stale_time_delta)))
                continue

            new_sispi_script = self.make_script()
