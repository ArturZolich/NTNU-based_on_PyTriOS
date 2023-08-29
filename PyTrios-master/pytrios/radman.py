#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Autonomous operation of hyperspectral radiometers with optional rotating measurement platform, solar power supply and remote connectivity

This script provides a class to interface with radiometers.

There should be a class for each family of sensors. Currently we have TriosManager to control 3 TriOS G1 (original) spectroradiometers and TriosG2Manager for the G2 update.
The G1 manager runs a thread for each communication port, always listening for measurement triggers and for sensor output. The G2 version monitors the sensor timer to determine when a measurement has finished and then idles while waiting for a new trigger. 

Plymouth Marine Laboratory
License: under development

"""
import os
import sys
import time
import datetime
import logging
import threading
import pytrios.pytriosg1 as ps
from numpy import log2

log = logging.getLogger('rad')
log.setLevel('INFO')


class TriosManager(object):
    """
    Trios G1 manager class
    """
    def __init__(self, port):
        # import pytrios only if used
        self.ports = [port]  # list of strings
        self.coms = ps.TMonitor(self.ports, baudrate=9600)
        self.sams = []
        self.ready = False
        self.connect_sensors()
        # track reboot cycles to prevent infinite rebooting of sensors if something unexpected happens (e.g a permanent sensor failure)
        self.reboot_counter = 0
        self.last_cold_start = datetime.datetime.now()
        self.last_connectivity_check = datetime.datetime.now()
        self.lasttrigger = None  # don't use this to get a timestamp on measurements, just used as a delay timer
        self.busy = False  # check this value to see if sensors are ready to sample

    def __del__(self):
        ps.tchannels = {}
        ps.TClose(self.coms)

    def stop(self):
        ps.tchannels = {}
        ps.TClose(self.coms)

    def connect_sensors(self):
        """(re)connect all serial ports and query all sensors"""
        self.busy = True

        ps.TClose(self.coms)
        ps.tchannels = {}
        time.sleep(1)

        log.info("Connecting: Starting listening threads")
        self.coms = ps.TMonitor(self.ports, baudrate=9600)
        time.sleep(3)

        for com in self.coms:
            # set verbosity for com channel (com messages / errors)
            # 0/1/2 = none, errors, all
            com.verbosity = 1
            # query connected instruments
            ps.TCommandSend(com, commandset=None, command='query')
        time.sleep(3)  # pause to receive query results
        self._identify_sensors()

        if len(self.sams) == 0:
            ps.TClose(self.coms)
            self.ready = False
            log.critical("no SAM modules found")
        else:
            self.ready = True

        for s in self.sams:
            # 0/1/2/3/4 = none, errors, queries(default), measurements, all
            self.tc[s].verbosity = 1
            self.tc[s].failures = 0

        self.busy = False

    def _identify_sensors(self):
        """identify SAM instruments from identified channels"""
        self.tk = list(ps.tchannels.keys())
        self.tc = ps.tchannels
        self.sams = [k for k in self.tk if ps.tchannels[k].TInfo.ModuleType in ['SAM', 'SAMIP']]  # keys
        self.chns = [self.tc[k].TInfo.TID for k in self.sams]  # channel addressing
        self.sns = [self.tc[k].TInfo.serialn for k in self.sams]  # sensor ids

        log.info("found SAM modules: {0}".format(list(zip(self.chns, self.sns))))


    def sample_all(self, trigger_id, sams_included=None, inttime=0):
        """Send a command to take a spectral sample from every sensor currently detected by the program"""
        self.lasttrigger = datetime.datetime.now()  # this is not used to timestamp measurements, only to track progress
        self.busy = True
        try:
            if sams_included is None:
                sams_included = self.sams

            for s in sams_included:
                self.tc[s].startIntSet(self.tc[s].serial, inttime, trigger=self.lasttrigger)

            # follow progress
            npending = len(sams_included)
            while npending > 0:
                # pytrios has a 12-sec timeout period for sam instruments so this will not loop forever
                # triggered measurements may not be pending but also not finished (i.e. incomplete or missing data)
                finished = [k for k in sams_included if self.tc[k].is_finished()]
                pending = [k for k in sams_included if self.tc[k].is_pending()]
                nfinished = len(finished)
                npending = len(pending)
                time.sleep(0.05)

            # account failed and successful measurement attempts
            missing = list(set(sams_included) - set(finished))

            for k in finished:
                self.tc[k].failures = 0
            for k in missing:
                self.tc[k].failures +=1

            # how long did the measurements take to arrive?
            if nfinished > 0:
                if type(self.tc[k].TSAM.lastRawSAMTime) == type(self.lasttrigger) and self.tc[k].TSAM.lastRawSAMTime is not None:
                    delays = [self.tc[k].TSAM.lastRawSAMTime - self.lasttrigger for k in sams_included]
                    delaysec = max([d.total_seconds() for d in delays])
                    log.info("\t{0} spectra received, triggered at {1} ({2} s)"
                        .format(nfinished, self.lasttrigger, delaysec))

            if len(missing) > 0:
                log.warning("Incomplete or missing result from {0}".format(",".join(missing)))

            # gather succesful results
            specs = [self.tc[s].TSAM.lastRawSAM
                    for s in sams_included if self.tc[s].is_finished()]
            sids = [self.tc[s].TInfo.serialn
                    for s in sams_included if self.tc[s].is_finished()]
            itimes = [self.tc[s].TSAM.lastIntTime
                    for s in sams_included if self.tc[s].is_finished()]

            self.busy = False
            pre_incs = [None]
            post_incs = [None]
            temp_incs = [None]
            # specs, sids, itimes may be empty lists, Last three fields for forward compatibility
            return trigger_id, specs, sids, itimes, pre_incs, post_incs, temp_incs

        except Exception as m:
            ps.TClose(self.coms)
            log.exception("Exception in TriosManager: {}".format(m))
            raise
