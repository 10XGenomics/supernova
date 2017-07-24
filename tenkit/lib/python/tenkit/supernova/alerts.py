#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# "Pseudo" Martian API for reporting errors and warnings from C++
# code prior to deployment of a formal adapter.
#

import os
import json
import martian

# hack to import from the top-level of the this package
tenkit=__import__(__name__.split('.')[0])

from tenkit.constants import ALARMS_LOCATION

def check_alert(stage, alert):
    keys = ["action", "metric", "compare", "threshold", "message"]
    for key in keys:
        if not alert.has_key(key):
           print key, " is missing in "
           print alert
           martian.throw("incorrectly formatted alert, see stdout.")
    if not ( alert["compare"] == "<" or alert["compare"] == ">" ):
        print alert
        martian.throw("invalid value for compare in alert")
    if not ( type(alert["threshold"]) == int or type(alert["threshold"]) == float ):
        martian.throw("%s: invalid type for threshold" % type(alert["threshold"]))

def load_alerts():
    alerts_file = os.path.join(ALARMS_LOCATION, "alarms-supernova.json")
    json_string = open( alerts_file, "r" ).read()
    try:
        alerts = json.loads(json_string)
    except:
        martian.throw("Incorrectly formatted alarms-supernova.json file.")
    for stage, alert_list in alerts.iteritems():
        for alert in alert_list:
            check_alert(stage, alert)
    return alerts

def write_stage_alerts(stage, path, alerts_file="alerts.list"):
    alerts = load_alerts()
    out_file = os.path.join(path, alerts_file)
    if not os.path.exists(path):
        os.makedirs(path)
    out_handle = open(out_file, "w")
    keys = ["metric", "threshold", "compare", "action", "message"]
    if not alerts.has_key(stage):
        martian.throw("No alerts found for stage %s" % stage)
    for alert in alerts[stage]:
        out_handle.write("#\n")
        out_handle.write(alert["metric"]+"\n")
        out_handle.write(str(alert["threshold"])+"\n")
        out_handle.write(alert["compare"]+"\n")
        out_handle.write(alert["action"]+"\n")
        out_handle.write(alert["message"]+"\n")
    out_handle.close()

class AlertLogger:
    def __init__(self, stage):
        alerts_all = load_alerts()
        self.alerts = alerts_all[stage]

    def issue(self, metric, value, format_string=""):
        for alert in self.alerts:
            ## find the right metric
            if alert["metric"] == metric:
                ## should we trigger?
                if (alert["compare"] == ">") ^ (value < alert["threshold"]):
                    ## optional formatting of alert message with format_string or value
                    if len(format_string) == 0:
                        format_string = str(value)
                    message = alert["message"].replace("{}", format_string)
                    ## issue an alert
                    if alert["action"] == "alarm":
                        martian.alarm(message)
                    elif alert["action"] == "exit":
                        martian.exit(message)

class SupernovaAlarms:
    SN_STAGE_ALARMS  = "martian_alerts.json"
    SN_ROLLUP_ALARMS = "alerts_rollup.txt"
    SN_EXIT  = u"exit"
    SN_ALARM = u"alarm"
    SN_ALARM_HEAD      = "The following warning(s) were issued prior to encountering an error:"
    SN_UNEXPECTED_TEXT = "An unexpected error has occurred."
    SN_ALERT_HANDLERS={ SN_ALARM    : martian.alarm, \
                        u"log_info" : martian.log_info, \
                        u"log_warn" : martian.log_warn, \
                        u"throw"    : martian.throw, \
                        SN_EXIT     : martian.exit }

    def __init__(self, base_dir,
                 alarms_file = SN_STAGE_ALARMS, rollup_file = SN_ROLLUP_ALARMS,
                 delete = True ):
        self._alarms_file=os.path.join(base_dir, alarms_file)
        self._rollup_file=os.path.join(base_dir, rollup_file)
        self._posted=False
        if delete: self.check_delete()

    def exit(self, msg=None):
        exit_handler=self.SN_ALERT_HANDLERS[self.SN_EXIT]
        if msg is None:
            full_msg = self.SN_UNEXPECTED_TEXT
        else:
            full_msg = msg
        if os.path.exists(self._rollup_file):
            issued_alarms = open(self._rollup_file, "r").read().split("\n")
            ## remove duplicate alarms that could arise if a stage that
            ## issues an alarm is restarted multiple times
            issued_alarms = list(set(issued_alarms))
            rollup = "\n".join( issued_alarms ) + "\n"
            full_msg += "\n\n"
            full_msg += self.SN_ALARM_HEAD
            full_msg += "\n\n"
            full_msg += rollup
        exit_handler(full_msg)

    def post(self):
        self._posted=True

        ## if alarm file does not exist, do nothing
        if os.path.exists(self._alarms_file):
            with open(self._alarms_file, 'r') as fp:
                alerts=json.loads(fp.read())
        else:
            return

        meta_alarm=[]
        exit_str=''

        for k,v in alerts.iteritems():
            if not k in self.SN_ALERT_HANDLERS:
                meta_alarm.append("unknown key {} in {} (BUG)".format(k,self._alarms_file))
            elif k == self.SN_EXIT:
                exit_str=';'.join(v)
            else:
                handler=self.SN_ALERT_HANDLERS[k]
                for post in v:
                    handler( post )

        for alarm in meta_alarm:
            martian.alarm(alarm)

        if len(exit_str) > 0:
            self.exit(exit_str)

    def __del__(self):
        if not self._posted:
            martian.alarm("C++ alarms were not posted, but left in file (BUG).")
        else:
            self.check_delete()

    def check_delete(self, filename = None):
        if filename == None:
            filename = self._alarms_file
        if os.path.isfile( filename ):
            os.unlink(filename)

