import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
import pandas as pd


def plot_sequence_from_df(dataframe,timerange = None):
    """Function to plot sequences from dataframes created by sequence_to_dataframe
        dataframe --dataframe object or path to csv file
        timerange --  [start time,end time] to pick a specific time interval to plot over
    """
    #read in csv
    if type(dataframe)==str:
        df = pd.read_csv(dataframe,index_col=0)
    else:
        df = dataframe
    #reduce font size for labels
    plt.rcParams.update({'font.size':6})
    keys = df.keys()
    for key in keys:
        if "time" in key:
            timebase = df[key]
            timeunits = key.split('-')[1] #extract time units

    ### Figure set up ###
    SequenceFig , axes = plt.subplots(len(keys)-2,1,sharex=True)
    axes[-1].set_xlabel('Time/'+timeunits) # Set label on bottom plot
    color=iter(cm.tab20(np.linspace(0,1,10))) # An iterable of colors 

    i = 0 #iterator for plots
    for key in keys:
        #look for channels in csv file
        if "SAO" in key or "FAO" in key: 
            channel_name = key.split(": ")[1]
            axes[i].plot(timebase,df[key],color = next(color)) #plot data 
            axes[i].set_ylabel(channel_name.replace(" ","\n"),rotation = 'horizontal',ha = "right") #multi line label
            if timerange:
                axes[i].set_xlim(timerange[0],timerange[1])
            
            i+=1 #next plot
    print("Data:\n",df.head())    
    plt.show()

if __name__=='__main__':
    plot_sequence_from_df("BECSequence_200302.csv",[15,20.3]) #edit path to csv file as needed