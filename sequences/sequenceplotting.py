"""Networking - Sequence Plotting
Using the translator module we can load in an xml sequence, and extract the timestamps of every event, and the analouge voltages
This code produces a plot of this data.

First the time axis of the plot is extracted

"""
import sys
sys.path.append('') # otherwise cwd isn't in sys.path 
from translator import translate, event_list, header_cluster, channel_names, analogue_cluster
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from sequencePreviewer import Previewer
try:
    from PyQt4.QtGui import QApplication
except ImportError:
    from PyQt5.QtWidgets import QApplication
import pandas as pd
import time
# Plot settings
font = {'family' : 'normal',
        'size'   : 12} 
plt.rcParams.update({'font.size':6})

# PyDex module 
t = translate()
t.load_xml('sequences\\SequenceFiles\\BECSequence_200302.xml')
t.load_xml(r"C:\Users\Jonathan\Documents\PhD\Experiment\PyDexData\41K Molasses_200310.xml")
esc = t.seq_dic['Experimental sequence cluster in'] # Sequence cluster (data)

timebase = [0] # This list will store the time stamps of each sequence
totaltime = 0 

### Range of data to Plot
steprange = range(0,30) # time step range NB this is the range [start,stop)
fastanalouges = [1,2] # fast analouges to plot
slowanalouges = [8,9,14] # slow analouges to plot

fastanalouges = []
slowanalouges = [1,3,5,7,8] # slow analouges to plot


### Figure set up ###
SequenceFig , axes = plt.subplots(len(fastanalouges)+len(slowanalouges),1,sharex=True)
axes[-1].set_xlabel('Time/ms') # Set label on bottom plot
color=iter(cm.tab20(np.linspace(0,1,10))) # An iterable of colors 
### Data frame ###
df = pd.DataFrame()


### Time Base ###
for i in steprange[0:-1]: #Because the last timestep values aren't used
    print(i)
    header  = esc['Sequence header top'][i]
    if header['Skip Step'] == '0':
        
        if header['Time unit'] == '2': #Seconds
            newtime = float(header['Time step length'])*1000
        if header['Time unit'] == '1': #milliseconds
            newtime = float(header['Time step length'])
        if header['Time unit'] == '0':
            newtime = float(header['Time step length'])/1000

        totaltime += newtime #float(header['Time step length'])  
        timebase.append(totaltime)
        timebase.append(totaltime)

df['time-ms'] = timebase

## Fast Anlouges ###
j = 0
for channel in fastanalouges:
    c = next(color)
    oldvoltage= float(esc['Fast analogue array'][channel]['Voltage'][steprange[0]])
    #Set up arrays to store
    ramps = [esc['Fast analogue array'][channel]['Ramp?'][steprange[0]]]
    voltages = [oldvoltage]

    for i in steprange[1:]: # Because the first timestep values have already been assigned 
        if esc['Sequence header top'][i]['Skip Step'] == '0':
            ramps.append(esc['Fast analogue array'][channel]['Ramp?'][i])
            print(ramps[-1])
            # print(esc['Fast analogue array'][channel]['Voltage'][i])
            newvoltage = float(esc['Fast analogue array'][channel]['Voltage'][i])
            print(i,ramps[-2],oldvoltage,newvoltage,esc['Sequence header top'][i]['Time step length'])
            if ramps[-2] =='1': #NB becuse of how ramping works this needs to be the previous step
                print('ramping')
                voltages.append(newvoltage)
            elif ramps[-2] == '0':
                voltages.append(oldvoltage)
            voltages.append(newvoltage)
            oldvoltage = newvoltage
    
    print(len(timebase),len(voltages))
    print(timebase,voltages)
    channelname = esc['Fast analogue names']['Name'][channel]

    axes[j].plot(timebase,voltages,color = c)
    axes[j].set_ylabel(channelname)
    axes[j].set_xlim(10000,10300)
    # print(esc['Fast analogue names'].keys())
    df[channelname] = voltages
    j += 1

for channel in slowanalouges:
    c = next(color)
    oldvoltage= float(esc['Slow analogue array'][channel]['Voltage'][steprange[0]])
    #Set up arrays to store
    ramps = [esc['Slow analogue array'][channel]['Ramp?'][steprange[0]]]
    voltages = [oldvoltage]

    for i in steprange[1:]: # Because the first timestep values have already been assigned 
        if esc['Sequence header top'][i]['Skip Step'] == '0':
            ramps.append(esc['Slow analogue array'][channel]['Ramp?'][i])
            print(ramps[-1])
            # print(esc['Slow analogue array'][channel]['Voltage'][i])
            newvoltage = float(esc['Slow analogue array'][channel]['Voltage'][i])
            print(i,ramps[-2],oldvoltage,newvoltage,esc['Sequence header top'][i]['Time step length'])
            if ramps[-2] =='1': #NB becuse of how ramping works this needs to be the previous step
                print('ramping')
                voltages.append(newvoltage)
            elif ramps[-2] == '0':
                voltages.append(oldvoltage)
            voltages.append(newvoltage)
            oldvoltage = newvoltage
    
    print(len(timebase),len(voltages))
    print(timebase,voltages)
    channelname = esc['Slow analogue names']['Name'][channel]

    axes[j].plot(timebase,voltages,color = c)
    axes[j].set_ylabel(channelname)
    axes[j].set_xlim(10000,10300)

    # print(esc['Slow analogue names'].keys())
    j += 1 #plot incrementer
    df["SAO "+str(channel)+" "+channelname] = voltages
# plt.plot(esc['Fast analogue array'][2]['Voltage'])
print(df.head())

df.to_csv("Ksequence.csv")
plt.show()


