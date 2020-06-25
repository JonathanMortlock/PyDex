from master import Master
from PyQt5.QtWidgets import *
import sys
from microscope.microscopeSettingsGui import settings_window
import numpy as np
import time
def run():
    """Initiate an app to run the program
    if running in Pylab/IPython then there may already be an app instance"""
    app = QApplication.instance()
    standalone = app is None # false if there is already an app instance
    if standalone: # if there isn't an instance, make one
        app = QApplication(sys.argv) 
        
    boss = Master(image_analysis=settings_window,experiment="microscope")
    print("before show")
    boss.show()
    print("after show")
    image_path  = r"C:\Users\Jonathan\Documents\PhD\Experiment\MicroscopeAnalysis\AtomTestImages\MIwithsparse2data.npy"
    image = np.load(image_path)
        # for i in range(10):
        # time.sleep(2)
    boss.rn.receive(image)
        # print("Image sent?")
    if standalone: # if an app instance was made, execute it
        sys.exit(app.exec_()) # when the window is closed, python code stops
    


if __name__ == "__main__":
    run()