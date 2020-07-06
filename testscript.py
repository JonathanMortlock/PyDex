from PyQt5.QtWidgets import *
import sys
from microscope.microscopeSettingsGui import settings_window
import numpy as np
import time
import os


os.chdir(r'C:\Users\Jonathan\Documents\PhD\Experiment\PyDex')
import master

"""Initiate an app to run the program
if running in Pylab/IPython then there may already be an app instance"""
app = QApplication.instance()
standalone = app is None # false if there is already an app instance
if standalone: # if there isn't an instance, make one
    app = QApplication(sys.argv) 

def fakeimage(boss,dexternumber):
    boss.rn.server.dxnum.emit(dexternumber)
    boss.rn.receive(image)
    

boss = master.Master(image_analysis=settings_window)
print("before show")
boss.show()
print("after show")
image_path  = r"C:\Users\Jonathan\Documents\PhD\Experiment\MicroscopeAnalysis\AtomTestImages\MIwithsparse2data.npy"
image = np.load(image_path)
    # for i in range(10):
    # time.sleep(2)
dexternumber = 0
fakeimage(boss,'0')






