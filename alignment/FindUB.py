#!/usr/bin/env python

# originally from /SNS/ARCS/IPTS-22296/shared/scripts_Arnab/FindUB_v1_MNFO.py
# modified by Jiao Lin

import sys,os
sys.path.append("/opt/mantidnightly/bin")

from mantid.simpleapi import *
from mantid import logger
import numpy as np
np.seterr("ignore")

# Inputs
instrument = 'ARCS'
IPTS = '22269'
nxdir = '/SNS/{}/IPTS-{}/nexus/'.format(instrument, IPTS)
outputdir = '/SNS/{}/IPTS-{}/shared/alignment-out/'.format(instrument, IPTS)
#runs = range(216547,216647,1)
runs = [128414]
omega_nxs_entry_name = 'Omega'
a,b,c,alpha,beta,gamm = 5.27, 5.27, 13.93, 90., 90., 120.

#maskfile = ...
#mask = LoadMask(Instrument=instrument, InputFile= maskfile)
#IDFfile = ...


# Should not need to change below
totalrun = len(runs)
print "Total number of runs %d" %totalrun
print runs
toMerge1=[]

if mtd.doesExist('mdmesh'):
    DeleteWorkspace('mdmesh')

for index, r in enumerate(runs):
    print "Run %d of %d, Processing converting MD for run : %s" %(index+1, totalrun, r)
    ows=instrument + '_'+str(r)
    toMerge1.append(ows)
    filename=nxdir+instrument+'_'+str(r)+'.nxs.h5'
    LoadEventNexus(Filename=filename, OutputWorkspace=ows)
    #LoadInstrument(Workspace= ows, Filename=IDFfile,RewriteSpectraMap=False)
    dataR=mtd[ows]
    Omega = dataR.getRun().getLogData(omega_nxs_entry_name).value.mean()+ 0
    # temp= dataR.getRun().getLogData('BL14B:SE:SampleTemp').value.mean()
    # print 'Omega = %5.2f deg., Temp = %5.2f K' %(Omega, temp)
    #MaskDetectors(Workspace=ows,MaskedWorkspace=mask)
    SetGoniometer(ows,Axis0=str(Omega)+',0,1,0,1') 
    mdmesh = ConvertToMD(
        InputWorkspace=ows, QDimensions='Q3D',dEAnalysisMode='Elastic', Q3DFrames='Q_sample', 
        LorentzCorrection='1', MinValues='-6.1,-6.1,-6.1',MaxValues='6.1,6.1,6.1', OverwriteExisting = 0
    )

data=GroupWorkspaces(toMerge1)

FindPeaksMD(InputWorkspace='mdmesh', PeakDistanceThreshold=0.5, MaxPeaks=100, DensityThresholdFactor=2000, OutputWorkspace='peaks')
FindUBUsingLatticeParameters(
    PeaksWorkspace='peaks', a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, Tolerance=0.45, FixParameters=True)
IndexPeaks(PeaksWorkspace='peaks', Tolerance=0.45, RoundHKLs=False)


#'options:  Cubic, Tetragonal, Orthorhombic, Hexagonal, Rhombohedral, Monoclinic, Triclinic'
OptimizeLatticeForCellType(PeaksWorkspace='peaks', CellType='Hexagonal', Apply=True, Tolerance=0.2)
IndexPeaks(PeaksWorkspace='peaks', Tolerance=0.2, RoundHKLs=False)

UBfile = os.path.join(outputdir, 'optub.mat')
SaveIsawUB('peaks', UBfile)

w = CreateSingleValuedWorkspace()
LoadIsawUB(w,filename)
ol=w.sample().getOrientedLattice()
print ol.getuVector(), ol.getvVector()

#define the projection directions
# proj=['1,0,0', '0, 1, 0', '0,0,1']
# LoadIsawUB(InputWorkspace='data',Filename=UBfile)
# ConvertToMD(InputWorkspace='data',OutputWorkspace='mdHKL',QDimensions='Q3D',dEAnalysisMode='Elastic', Q3DFrames='HKL',Uproj=proj[0],Vproj=proj[1],Wproj=proj[2], 
#   QConversionScales='HKL',LorentzCorrection='1', MinValues='-7.1,-7.1,-1.6',MaxValues='7.1,7.1,1.6')
# mdHKLmesh=MergeMD('mdHKL')
