# Author: Andrei Savici & Jiao Lin
#
# Purpose: find out correlation between scattered intensities and transmited intensities
#
# Crystal at certain orientation may strongly scatter the incident neutrons.
# This script plot transmission intensity (2nd monitor) and scattering intensity
# at certain bragg angle (from masked detector pixels) vs the sample angles and
# look for correlations

# Inputs
runs300K=range(40402,40419) +range(40422,40427)+range(40430,40492)+range(40497,40534)
time_interval1=[5400,5500]
time_interval2=[8400,8600]
# this mask file was created by hand using Mantidplot GUI
mask_path = "monvsAngle_mask.xml"
IPTS_no = 9265
asciiout = "corr.dat"


from mantid.simpleapi import *
import numpy as np
# load mask
MaskWs=LoadMask('ARCS', mask_path)


def getVE(r):
    "get the only value and error"
    return r.readY(0)[0], r.readE(0)[0]

def integrateMonitor(ws, mon_no, toflimits):
    "integrate a portion of the monitor data"
    lower, upper = toflimits
    r = Integration(
        InputWorkspace=ws,RangeLower=lower,RangeUpper=upper,
        StartWorkspaceIndex=mon_no, EndWorkspaceIndex=mon_no)
    return getVE(r)

def process(runs):
    "process all runs"
    for runnumber in runs:
        print runnumber
        path='/SNS/ARCS/IPTS-%s/data/ARCS_%s_event.nxs' % (IPTS_no, runnumber)
        # monitor
        w=LoadNexusMonitors(Filename=path)
        w=NormaliseByCurrent(w)
        mon1_value, mon1_err = integrateMonitor(w, 0, time_interval1)
        mon2_value, mon2_err = integrateMonitor(w, 1, time_interval2)
        # angle
        angle=w.run()['micas70mmRot'].getStatistics().mean
        # detector
        w=Load(Filename=path)
        w=NormaliseByCurrent(w)
        all=SumSpectra(w)
        all_v, all_e=getVE(all)
        MaskDetectors(Workspace=w,MaskedWorkspace=MaskWs)
        masked = SumSpectra(w)
        masked_v, masked_e=getVE(masked)
        yield angle, mon1_value, mon1_err, mon2_value, mon2_err, all_v, all_e, masked_v, masked_e
        continue
    return

values = np.array(list(process(runs300K)))
values = values[np.argsort(values[:, 0])] # sort by angle
header = "angle\tmon1_v\tmon1_err\tmon2_v\tmon2_err\talldet_v\talldet_err\tmaskeddet_v\tmaskeddet_err\n"
np.savetxt(asciiout, values, header=header)
