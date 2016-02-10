# Author: Andrei Savici
#
# Purpose: find out correlation between scattered intensities and transmited intensities
#
# Crystal at certain orientation may strongly scatter the incident neutrons.
# This script plot transmission intensity (2nd monitor) and scattering intensity
# at certain bragg angle (from masked detector pixels) vs the sample angles and
# look for correlations

# Inputs
runs300K=range(40402,40419)+range(40422,40427)+range(40430,40492)+range(40497,40534)
time_interval1=[5400,5500]
time_interval2=[8400,8600]


from mantid.simpleapi import *

# Results
angles=[]
mon1_array=[]
mon2_array=[]
mon1e_array=[]
mon2e_array=[]
data_array=[]
error_array=[]
masked_data_array=[]
masked_error_array=[]

# this mask file was created by hand using Mantidplot GUI
MaskWs=LoadMask('ARCS','monvsAngle_mask.xml')

runs=runs300K
for runnumber in runs:
    print runnumber
    filename='/SNS/ARCS/IPTS-9265/data/ARCS_'+str(runnumber)+'_event.nxs'
    # monitor
    w=LoadNexusMonitors(Filename=filename)
    w=NormaliseByCurrent(w)
    mon1=Integration(InputWorkspace=w,RangeLower=time_interval1[0],RangeUpper=time_interval1[1],StartWorkspaceIndex=0,EndWorkspaceIndex=0)
    mon2=Integration(InputWorkspace=w,RangeLower=time_interval2[0],RangeUpper=time_interval2[1],StartWorkspaceIndex=1,EndWorkspaceIndex=1)
    mon1_value=mon1.readY(0)[0]
    mon2_value=mon2.readY(0)[0]
    mon1_err=mon1.readE(0)[0]
    mon2_err=mon2.readE(0)[0]
    angle=w.run()['micas70mmRot'].getStatistics().mean
    angles.append(angle)
    mon1_array.append(mon1_value)
    mon2_array.append(mon2_value)
    mon1e_array.append(mon1_err)
    mon2e_array.append(mon2_err)
    # detector
    w=Load(Filename=filename)
    w=NormaliseByCurrent(w)
    all=SumSpectra(w)
    data_array.append(all.readY(0)[0])
    error_array.append(all.readE(0)[0])
    MaskDetectors(Workspace=w,MaskedWorkspace=MaskWs)
    masked=SumSpectra(w)
    masked_data_array.append(masked.readY(0)[0])
    masked_error_array.append(masked.readE(0)[0])
    continue    

# save as workspace
Monitor1=CreateWorkspace(angles,mon1_array,mon1e_array)
Monitor2=CreateWorkspace(angles,mon2_array,mon2e_array)
AllData=CreateWorkspace(angles,data_array,error_array)
MaskedData=CreateWorkspace(angles,masked_data_array,masked_error_array)
Monitor1=SortXAxis(Monitor1)
Monitor2=SortXAxis(Monitor2)
AllData=SortXAxis(AllData)
MaskedData=SortXAxis(MaskedData)
