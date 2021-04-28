"""Networking - Sequence Plotting
Using the translator module we can load in an xml sequence, and extract the timestamps of every event, and the analogue voltages
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
# font = {'family' : 'normal',
#         'size'   : 12} 

def sequence_to_dataframe(filepath,steprange = range(0,30),fastanalogues = [],slowanalogues = [],timeunits = "ms",alldata = True,annotate = True):
    """Function to extract info from pydex xml files into a data frame
    """
    # PyDex module 
    t = translate()

    ### Range of data to Plot ###
    t.load_xml(filepath)

    esc = t.seq_dic['Experimental sequence cluster in'] # Sequence cluster (data)
    print(len(esc['Sequence header top']))

    ### Catch input errors ###
    numsteps = len(esc['Sequence header top'])
    numfast = len(esc['Fast analogue array'])
    numslow = len(esc["Slow analogue array"])
    print("Time steps",[h['Time step name'] for h in esc['Sequence header top']])
    print("Fast analogue channel Names",[i for i in enumerate(esc['Fast analogue names']['Name'])])
    print("Slow analogue channel Names",[i for i in enumerate(esc['Slow analogue names']['Name'])])

    if max(steprange)>numsteps:
        raise ValueError("Max Step value out of range")
    if max(fastanalogues,default=0)>numfast:
        raise ValueError("Fast analogue out of range")
    if max(slowanalogues,default =0)>numslow:
        raise ValueError("Slow analogue out of range")
    if alldata:
        steprange = range(len(esc['Sequence header top']))
        fastanalogues =  range(len(esc['Fast analogue array']))
        slowanalogues =  range(len(esc['Slow analogue array']))
        fastdigitals = range(len(esc['Fast digital channels']))
        slowdigitals = range(len(esc['Slow digital channels']))
    else:
        fastdigitals= []
        slowdigitals = []
    ### Time units ###
    if timeunits=="ms":
        time_conversions = [1e-3,1,1e-3]
    elif timeunits=="s":
        time_conversions = [1e-6,1e-3,1]
    elif timeunits=="us":
        time_conversions = [1,1e3,1e6]
    else:
        raise NotImplementedError("Unknown time unit. 'ms' 's' and 'us' are avaliable ")

    ### Data frame ###
    df = pd.DataFrame()


    ### Time Base ###
    timebase = [0] # This list will store the time stamps of each sequence
    totaltime = 0 
    names = []
    for i in steprange[0:-1]: #Because the last timestep values aren't used
        print(i)
        header  = esc['Sequence header top'][i]
        if header['Skip Step'] == '0':
            
            if header['Time unit'] == '0': #Seconds
                newtime = float(header['Time step length'])*time_conversions[0]
            if header['Time unit'] == '1': #milliseconds
                newtime = float(header['Time step length'])*time_conversions[1]
            if header['Time unit'] == '2': #microseconds
                newtime = float(header['Time step length'])*time_conversions[2]

            totaltime += newtime #float(header['Time step length'])  
            timebase.append(totaltime)
            timebase.append(totaltime)
            names.append(header['Time step name']+" Start") #NB note that names is shorter than time based at this point
            names.append(header['Time step name']+" End")
    names.append("last?")
    print(len(timebase),len(names))
    df['time-'+str(timeunits)] = timebase
    df['names']  = names

    ## Fast Anlouges ###
    j = 0
    for channel in fastanalogues:
        oldvoltage= float(esc['Fast analogue array'][channel]['Voltage'][steprange[0]])
        #Set up arrays to store
        ramps = [esc['Fast analogue array'][channel]['Ramp?'][steprange[0]]]
        voltages = [oldvoltage]

        for i in steprange[1:]: # Because the first timestep values have already been assigned 
            if esc['Sequence header top'][i]['Skip Step'] == '0':
                ramps.append(esc['Fast analogue array'][channel]['Ramp?'][i])
                # print(ramps[-1])
                # print(esc['Fast analogue array'][channel]['Voltage'][i])
                newvoltage = float(esc['Fast analogue array'][channel]['Voltage'][i])
                # print(i,ramps[-2],oldvoltage,newvoltage,esc['Sequence header top'][i]['Time step length'])
                if ramps[-2] =='1': #NB becuse of how ramping works this needs to be the previous step
                    # print('ramping')
                    voltages.append(newvoltage)
                elif ramps[-2] == '0':
                    voltages.append(oldvoltage)
                voltages.append(newvoltage)
                oldvoltage = newvoltage
        
        # print(len(timebase),len(voltages))
        # print(timebase,voltages)
        channelname = esc['Fast analogue names']['Name'][channel]
        if not channelname:
            channelname ='Unnamed'
        else:
            df["FAO "+str(channel)+": "+channelname] = voltages
        j += 1
    ###Slow analogues###
    print("j at slow",j)
    for channel in slowanalogues:
        oldvoltage= float(esc['Slow analogue array'][channel]['Voltage'][steprange[0]])
        #Set up arrays to store
        ramps = [esc['Slow analogue array'][channel]['Ramp?'][steprange[0]]]
        voltages = [oldvoltage]

        for i in steprange[1:]: # Because the first timestep values have already been assigned 
            if esc['Sequence header top'][i]['Skip Step'] == '0':
                ramps.append(esc['Slow analogue array'][channel]['Ramp?'][i])
                # print(ramps[-1])
                # print(esc['Slow analogue array'][channel]['Voltage'][i])
                newvoltage = float(esc['Slow analogue array'][channel]['Voltage'][i])
                # print(i,ramps[-2],oldvoltage,newvoltage,esc['Sequence header top'][i]['Time step length'])
                if ramps[-2] =='1': #NB becuse of how ramping works this needs to be the previous step
                    # print('ramping')
                    voltages.append(newvoltage)
                elif ramps[-2] == '0':
                    voltages.append(oldvoltage)
                voltages.append(newvoltage)
                oldvoltage = newvoltage
        
        # print(len(timebase),len(voltages))
        # print(timebase,voltages)
        channelname = esc['Slow analogue names']['Name'][channel]
        if not channelname:
            channelname ='Unnamed'
        else:
            df["SAO "+str(channel)+": "+channelname] = voltages
    for channel in slowdigitals:
        values = []

        # values.append(int(esc['Slow digital channels'][steprange[0]][channel]))

        for i in steprange[0:-1]:
            if esc['Sequence header top'][i]['Skip Step'] == '0':
                values.append(int(esc['Slow digital channels'][i][channel]))
                values.append(int(esc['Slow digital channels'][i][channel]))
        print(esc["Slow digital names"]['Name'])
        channelname = esc["Slow digital names"]['Name'][channel]
        values.append(int(esc['Slow digital channels'][steprange[-1]][channel])) #append last step

        if not channelname:
            channelname ='Unnamed'
        else:
            df["SDO "+str(channel)+": "+channelname] = values
    for channel in fastdigitals:
        values = []
        for i in steprange[0:-1]:
            if esc['Sequence header top'][i]['Skip Step'] == '0':
                values.append(int(esc['Fast digital channels'][i][channel]))
                values.append(int(esc['Fast digital channels'][i][channel]))
        print(esc["Fast digital names"]['Name'])
        channelname = esc["Fast digital names"]['Name'][channel]
        values.append(int(esc['Fast digital channels'][steprange[-1]][channel])) #append last step

        if not channelname:
            channelname ='Unnamed'
        else:
            df["FDO "+str(channel)+": "+channelname] = values
    print(df.head())
    if alldata:
        df.to_csv(filepath[:-4]+'NAMEDDATA.csv')
    else:
        df.to_csv(filepath[:-4]+".csv")
    return df

def plot_sequence_from_df(dataframe):
    """Function to plot sequences from dataframes created by sequence_to_dataframe
        dataframe --dataframe object or path to csv file
    """
    #read in csv
    if type(dataframe)==str:
        df = pd.read_csv(dataframe,index_col=0)
    else:
        df = dataframe
    plt.rcParams.update({'font.size':8})
    keys = df.keys()
    for key in keys:
        if "time" in key:
            timebase = df[key]
            timeunits = key.split('-')[1]

    ### Figure set up ###
    SequenceFig , axes = plt.subplots(len(keys)-2,1,sharex=True)
    axes[-1].set_xlabel('Time/'+timeunits) # Set label on bottom plot
    color=iter(cm.tab20(np.linspace(0,1,10))) # An iterable of colors 

    i = 0 #iterator for plots
    for key in keys:
        #look for channels in csv file
        if "SAO" in key or "FAO" in key: 
            channel_name = key.split(": ")[1]
            axes[i].plot(timebase,df[key],color = next(color))
            axes[i].set_ylabel(channel_name.replace(" ","\n"),rotation = 'horizontal',ha = "right")
            i+=1 #next plot

    print(df.head())

    
    plt.show()


def sequence_plot(filepath,steprange = range(0,30),fastanalogues = [],slowanalogues = [],timeunits = "ms",annotate = True):
    """Function to plot sequences from pydex xml files
    """
    plt.rcParams.update({'font.size':8})
    # plt.style.use("tufte")
    # PyDex module 
    t = translate()



    ### Range of data to Plot ###
    t.load_xml(filepath)
    # steprange = range(0,30) # time step range NB this is the range [start,stop)
    # fastanalogues = []
    # slowanalogues = [1,3,5,7,8] # slow analogues to plot

    esc = t.seq_dic['Experimental sequence cluster in'] # Sequence cluster (data)
    print(len(esc['Sequence header top']))

    ### Catch input errors ###
    numsteps = len(esc['Sequence header top'])
    numfast = len(esc['Fast analogue array'])
    numslow = len(esc["Slow analogue array"])
    print("Time steps",[h['Time step name'] for h in esc['Sequence header top']])
    print("Fast analogue channel Names",[i for i in enumerate(esc['Fast analogue names']['Name'])])
    print("Slow analogue channel Names",[i for i in enumerate(esc['Slow analogue names']['Name'])])

    if max(steprange)>numsteps:
        raise ValueError("Max Step value out of range")
    if max(fastanalogues,default=0)>numfast:
        raise ValueError("Fast analogue out of range")
    if max(slowanalogues,default =0)>numslow:
        raise ValueError("Slow analogue out of range")



    ### Time units ###
    if timeunits=="ms":
        time_conversions = [1e-3,1,1e-3]
    elif timeunits=="s":
        time_conversions = [1e-6,1e-3,1]
    elif timeunits=="us":
        time_conversions = [1,1e3,1e6]
    else:
        raise NotImplementedError("Unknown time unit. 'ms' 's' and 'us' are avaliable ")
    ### Figure set up ###
    SequenceFig , axes = plt.subplots(len(fastanalogues)+len(slowanalogues),1,sharex=True)
    axes[-1].set_xlabel('Time/ms') # Set label on bottom plot
    color=iter(cm.tab20(np.linspace(0,1,10))) # An iterable of colors 
    ### Data frame ###
    df = pd.DataFrame()


    ### Time Base ###
    timebase = [0] # This list will store the time stamps of each sequence
    totaltime = 0 

    for i in steprange[0:-1]: #Because the last timestep values aren't used
        print(i)
        header  = esc['Sequence header top'][i]
        if header['Skip Step'] == '0':
            
            if header['Time unit'] == '0': #Seconds
                newtime = float(header['Time step length'])*time_conversions[0]
            if header['Time unit'] == '1': #milliseconds
                newtime = float(header['Time step length'])*time_conversions[1]
            if header['Time unit'] == '2': #microseconds
                newtime = float(header['Time step length'])*time_conversions[2]

            totaltime += newtime #float(header['Time step length'])  
            timebase.append(totaltime)
            timebase.append(totaltime)
            # if annotate:
                # axes[0].text(timebase[-2], 5,header['Time step name'])
    
    df['time-'+str(timeunits)] = timebase

    ## Fast Anlouges ###
    j = 0
    for channel in fastanalogues:
        c = next(color)
        oldvoltage= float(esc['Fast analogue array'][channel]['Voltage'][steprange[0]])
        #Set up arrays to store
        ramps = [esc['Fast analogue array'][channel]['Ramp?'][steprange[0]]]
        voltages = [oldvoltage]

        for i in steprange[1:]: # Because the first timestep values have already been assigned 
            if esc['Sequence header top'][i]['Skip Step'] == '0':
                ramps.append(esc['Fast analogue array'][channel]['Ramp?'][i])
                # print(ramps[-1])
                # print(esc['Fast analogue array'][channel]['Voltage'][i])
                newvoltage = float(esc['Fast analogue array'][channel]['Voltage'][i])
                # print(i,ramps[-2],oldvoltage,newvoltage,esc['Sequence header top'][i]['Time step length'])
                if ramps[-2] =='1': #NB becuse of how ramping works this needs to be the previous step
                    # print('ramping')
                    voltages.append(newvoltage)
                elif ramps[-2] == '0':
                    voltages.append(oldvoltage)
                voltages.append(newvoltage)
                oldvoltage = newvoltage
        
        # print(len(timebase),len(voltages))
        # print(timebase,voltages)
        channelname = esc['Fast analogue names']['Name'][channel]

        axes[j].plot(timebase,voltages,color = c)
        axes[j].set_ylabel(channelname.replace(" ","\n"),rotation = 'horizontal',ha = "right")
        # axes[j].set_xlim(10000,10300)
        # print(esc['Fast analogue names'].keys())
        df[channelname] = voltages
        j += 1
    ###Slow analogues###
    print("j at slow",j)
    for channel in slowanalogues:
        c = next(color)
        oldvoltage= float(esc['Slow analogue array'][channel]['Voltage'][steprange[0]])
        #Set up arrays to store
        ramps = [esc['Slow analogue array'][channel]['Ramp?'][steprange[0]]]
        voltages = [oldvoltage]

        for i in steprange[1:]: # Because the first timestep values have already been assigned 
            if esc['Sequence header top'][i]['Skip Step'] == '0':
                ramps.append(esc['Slow analogue array'][channel]['Ramp?'][i])
                # print(ramps[-1])
                # print(esc['Slow analogue array'][channel]['Voltage'][i])
                newvoltage = float(esc['Slow analogue array'][channel]['Voltage'][i])
                # print(i,ramps[-2],oldvoltage,newvoltage,esc['Sequence header top'][i]['Time step length'])
                if ramps[-2] =='1': #NB becuse of how ramping works this needs to be the previous step
                    # print('ramping')
                    voltages.append(newvoltage)
                elif ramps[-2] == '0':
                    voltages.append(oldvoltage)
                voltages.append(newvoltage)
                oldvoltage = newvoltage
        
        # print(len(timebase),len(voltages))
        # print(timebase,voltages)
        channelname = esc['Slow analogue names']['Name'][channel]

        axes[j].plot(timebase,voltages,color = c)
        axes[j].set_ylabel(channelname.replace(" ","\n"),rotation = 'horizontal',ha = "right")
        # axes[j].set_xlim(10000,10300)

        # print(esc['Slow analogue names'].keys())
        j += 1 #plot incrementer
        df["SAO "+str(channel)+" "+channelname] = voltages
    # plt.plot(esc['Fast analogue array'][2]['Voltage'])
    print(df.head())

    # df.to_csv("")
    
    plt.show()


if __name__== "__main__":
    path = r"C:\Users\jonat\Documents\PhD\Experiment\PyDex\sequences\SequenceFiles\BECSequence_200302.xml"
    dataframe = sequence_to_dataframe(path,slowanalogues= [6,8,9,14],fastanalogues=[1,2],timeunits='s')
    plot_sequence_from_df("TestBEC.csv")
    # sequence_plot(path,slowanalogues= [6,8,9,14],fastanalogues=[1,2],timeunits='s')

