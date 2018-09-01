#!/usr/bin/env python
import re
import subprocess as sp
import shlex
import sys
import time
import logging

#https://stackoverflow.com/questions/13839554/how-to-change-filehandle-with-python-logging-on-the-fly-with-different-classes-a
# Create my own log location
fileh = logging.FileHandler('/fslhome/fslcollab192/slurm_status.log', 'a')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fileh.setFormatter(formatter)

logger = logging.getLogger()
for hdlr in logger.handlers[:]:  # remove all old handlers
    logger.removeHandler(hdlr)
logger.addHandler(fileh)      # set the new handler
logger.setLevel(logging.DEBUG)

STATUS_ATTEMPTS = 20
jobid = sys.argv[1]

logger.info(f"Running script with jobid={jobid}")
for i in range(STATUS_ATTEMPTS):
    logger.info(f"jobid={jobid}, status_attemp_numb={i}")
    try:
        sacct_res = sp.check_output(shlex.split("sacct -P -b -j {} -n".format(jobid)))
        res = {x.split("|")[0]: x.split("|")[1] for x in sacct_res.decode().strip().split("\n")}
        break
    except sp.CalledProcessError as e:
        logger.error("sacct process error")
        logger.error(e)
    except IndexError as e:
        pass
    # Try getting job with scontrol instead in case sacct is misconfigured
    try:
        sctrl_res = sp.check_output(shlex.split("scontrol -o show job {}".format(jobid)))
        m = re.search("JobState=(\w+)", sctrl_res.decode())
        res = {jobid: m.group(1)}
        break
    except sp.CalledProcessError as e:
        logger.error("scontrol process error")
        logger.error(e)
        if i >= STATUS_ATTEMPTS - 1:
            print("failed")
            exit(0)
        else:
            time.sleep(1)

status = res[jobid]
logger.info(res)
logger.info(f"jobid={jobid}, status={status}")

if (status == "BOOT_FAIL"):
    print("failed")
elif (status == "OUT_OF_MEMORY"):
    print("failed")
elif (status == "CANCELLED"):
    print("failed")
elif (status == "COMPLETED"):
    print("success")
elif (status == "DEADLINE"):
    print("failed")
elif (status == "FAILED"):
    print("failed")
elif (status == "NODE_FAIL"):
    print("failed")
elif (status == "PREEMPTED"):
    print("failed")
elif (status == "TIMEOUT"):
    print("failed")
# Unclear whether SUSPENDED should be treated as running or failed
elif (status == "SUSPENDED"):
    print("failed")
else:
    print("running")
