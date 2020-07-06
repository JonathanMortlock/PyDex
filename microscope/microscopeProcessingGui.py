import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['image.cmap'] = 'inferno'
matplotlib.rc('font',family ='serif')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
import pyqtgraph as pg    # not as flexible as matplotlib but works a lot better with qt


sys.path.append('.')
sys.path.append('..')
import microscope.AtomImaging as ai
import imageanalysis.imageHandler as ih
from matplotlib.figure import Figure
from skimage import restoration, util, filters
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from functools import  partial
import time
import traceback, sys

# validators for user input
double_validator = QDoubleValidator() # floats
int_validator    = QIntValidator()    # integers
int_validator.setBottom(0) # don't allow -ve numbers
nat_validator    = QIntValidator()    # natural numbers 
nat_validator.setBottom(1) # > 0
def remove_slot(signal, slot, reconnect=True): #TODO This seems wrong to me, can I think of a way to replace?
    """Make sure all instances of slot are disconnected
    from signal. Prevents multiple connections to the same 
    slot. If reconnect=True, then reconnect slot to signal."""
    while True: # make sure that the slot is only connected once 
        try: signal.disconnect(slot)
        except TypeError: break
    if reconnect: signal.connect(slot)


class MplCanvas(FigureCanvasQTAgg):
    """
        How to include matplotlib in pyqt...

        This is probably slower than pyqtgraph but easier to work with and one less dependancy

    """
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111) #Could change??
        super(MplCanvas, self).__init__(fig)

class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data
    
    error
        `tuple` (exctype, value, traceback.format_exc() )
    
    result
        `object` data returned from processing, anything

    progress
        `int` indicating % progress 

    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)


class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and 
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()    

        # Add the callback to our kwargs
        # self.kwargs['progress_callback'] = self.signals.progress        

    @pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        
        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done
        



class main_window(QMainWindow):
    """
        Main GUI, inherits from the pyqt mainwindow class 
        init defines class variables
        init_UI is where the Graphical bits are defined.
        various functions to enable functionality,

        The most important thing to not here is that everything is event driven
        When an event happens data is passed by Signals to Slots. This is how you deal with asynchronus stuff and threading



    """
    event_im = pyqtSignal([np.ndarray, bool])

    def __init__(self,results_path='.', im_store_path='.', name='',
                im_handler=None, hist_handler=None, edit_ROI=False):
        super().__init__()
        self.name = name  # name is displayed in  the window title
        self.image_handler = im_handler if im_handler else ih.image_handler() 
        
        self.multirun = ''
        self.date = time.strftime("%d %b %B %Y", time.localtime()).split(" ") # day short_month long_month year
        self.log_file_name = results_path + 'log.dat' # in case init_log fails
        self.last_path = results_path # path history helps user get to the file they want
        self.init_log(results_path) # write header to the log file that collects histograms
        self.image_storage_path = im_store_path # used for loading image files
        #TODO wrap this in a "testing"
        self.image_path  = r"C:\Users\Jonathan\Documents\PhD\Experiment\MicroscopeAnalysis\AtomTestImages\MIwithsparse2data.npy"
        self.current_image = np.zeros((512,512))
        self.lattice_image = np.ones((512,512))
        remove_slot(self.event_im,partial(self.connectSelf,'current_image'))
        remove_slot(self.event_im,partial(self.connectSelf,'lattice_image'))

        # self.event_im.connect(partial(self.connectSelf,'current_image'))
        # self.event_im.connect(partial(self.connectSelf,'lattice_image'))
        # self.truthlist = np.load(self.image_path[:-8]+"truthandlattice.npy")
       

        self.threadpool = QThreadPool()

        # Processing Settings
        self.pic_width = 512
        self.pic_height = 512

        self.window = 30
        self.psf_radius = 7
        self.RLiterations = 50
        self.guess_angle = 0.0
        self.lattice_constant = 4.655
        self.lattice_vectors = np.array(((np.cos(self.guess_angle)*self.lattice_constant,np.sin(self.guess_angle)*self.lattice_constant),(-np.sin(self.guess_angle)*self.lattice_constant,np.cos(self.guess_angle)*self.lattice_constant)))
        #Variables
        self.offset = (0,0)
        self.blobs = []
        self.histogram_data = []
        self.threshold = 0
        self.atoms = []
        self.lattice = []
        self.fidelity = 0
        
        self.updated = {'current_image': False,'blobs': False,'lattice_vectors': False,'lattice': False, 'deconvolved': False,'histogram': False}
        self.init_UI(edit_ROI)
        self.set_im_show(True)
        # self.event_im.emit(self.current_image,True)
    def init_log(self, results_path='.'):
        """Create a directory for today's date as a subdirectory in the log file path
        then write the header to the log file path defined in config.dat"""
        # make subdirectory if it doesn't already exist
        results_path = os.path.join(results_path, 
                    r'%s\%s\%s'%(self.date[3],self.date[2],self.date[0]))
        try:
            os.makedirs(results_path, exist_ok=True)
        except PermissionError:  # couldn't access the path, start a log file here
            results_path = r'.\%s\%s\%s'%(self.date[3],self.date[2],self.date[0])
            os.makedirs(results_path, exist_ok=True)

        # log is saved in a dated subdirectory and the file name also has the date
        self.last_path = results_path
        self.log_file_name = os.path.join(results_path, 
                   self.name+'log'+self.date[0]+self.date[1]+self.date[3]+'.dat')  
        # # write the header to the log file commented out by JM
        # if not os.path.isfile(self.log_file_name): # don't overwrite if it already exists
        #     with open(self.log_file_name, 'w+') as f:
        #         f.write('#Single Atom Image Analyser Log File: collects histogram data\n')
        #         f.write('#include --[]\n')
        #         f.write('#'+', '.join(self.histo_handler.stats.keys())+'\n')
      
    def init_UI(self,edit_ROI):
        """
        Sets up GUI both visually and in terms of what events are connected to what functions

        """

        self.centre_widget = QWidget()
        self.tabs = QTabWidget()       # make tabs for each main display 
        self.centre_widget.layout = QVBoxLayout()
        self.centre_widget.layout.addWidget(self.tabs)
        self.centre_widget.setLayout(self.centre_widget.layout)
        self.setCentralWidget(self.centre_widget)

        double_validator = QDoubleValidator() # floats
        int_validator    = QIntValidator()    # integers
        int_validator.setBottom(0) # don't allow -ve numbers



        """
            Menubars
        """
        menubar = self.menuBar()
        
        # file menubar allows you to save/load data
        file_menu = menubar.addMenu('File')
        load_im = QAction('Load Image', self) # display a loaded image
        load_im.triggered.connect(self.load_image_from_file)
        set_lattice_image= QAction('Set lattice image ',self)
        set_lattice_image.triggered.connect(self.set_lattice_file)

        file_menu.addAction(load_im)
        file_menu.addAction(set_lattice_image)


        """ 
            Tabs
        """

        #### Settings Tab ####
        settings_tab = QWidget()
        settings_grid = QGridLayout()
        settings_tab.setLayout(settings_grid)
        
        self.auto_checkbox = QCheckBox("Auto Update?")
        self.file_path_label = QLabel(self.image_path)
        self.lattice_file_path_label = QLabel("lattice image path here")
        self.use_lattice_file_checkbox = QCheckBox("Use Lattice Image?")
        self.full_plots_checkbox = QCheckBox("Full Plots?")
        angle_label = QLabel("Lattice angle search Range")
        self.angle_min_edit = QLineEdit("0")
        self.angle_min_edit.setValidator(double_validator)
        self.angle_max_edit = QLineEdit("1.5707")
        self.angle_max_edit.setValidator(double_validator)
        self.angle_step_edit = QLineEdit("3000")
        self.angle_step_edit.setValidator(int_validator)

        self.manual_angle = QLineEdit("Angle?")
        self.manual_offset = QLineEdit("offset")
        self.manual_offset1 = QLineEdit("offset1")
        self.manual_offset.returnPressed.connect(self.manual_lattice)

        self.manual_psf_edit = QLineEdit("manual psf")
        self.manual_psf_edit.returnPressed.connect(self.manual_psf)
        # self.auto_process= self.auto_checkbox.isChecked()
        ### Grid layout
        settings_grid.addWidget(self.auto_checkbox,0,0)
        settings_grid.addWidget(QLabel("Image File path"),1,0)
        settings_grid.addWidget(self.file_path_label,1,1)
        settings_grid.addWidget(QLabel("Lattice image File path"),2,0)
        settings_grid.addWidget(self.lattice_file_path_label,2,1)
        settings_grid.addWidget(self.use_lattice_file_checkbox,3,0)
        settings_grid.addWidget(self.full_plots_checkbox,3,1)
        settings_grid.addWidget(angle_label,4,0)
        settings_grid.addWidget(self.angle_min_edit,4,1)
        settings_grid.addWidget(self.angle_max_edit,4,2)
        settings_grid.addWidget(self.angle_step_edit,4,3)
        # settings_grid.addWidget(self.loadbyid_edit,5,0)  #TODO db 
        settings_grid.addWidget(self.manual_angle,6,0)
        settings_grid.addWidget(self.manual_offset,6,1)
        settings_grid.addWidget(self.manual_offset1,6,2)
        settings_grid.addWidget(self.manual_psf_edit,7,0)

        self.tabs.addTab(settings_tab, "Settings")
        



        ###  Lattice finding tab ### 
        lattice_tab = QWidget()
        lattice_grid = QGridLayout()
        lattice_tab.setLayout(lattice_grid)
        self.tabs.addTab(lattice_tab,"Lattice Determination")

        # self.fftcanvas = MplCanvas(self)
        # ffttoolbar = NavigationToolbar(self.fftcanvas, self)
        # self.fftcanvas.axes.text(0.5,0.5,"FFT Shown here",fontsize = 20)


        self.blobcanvas = MplCanvas(self)
        blobtoolbar =  NavigationToolbar(self.blobcanvas, self)  
        self.blobcanvas.axes.text(0.5,0.5,"Blobz shown here",fontsize=20)

        self.update_blobs_button = QPushButton("Find Blobs")
        self.update_blobs_button.clicked.connect(self.find_blobs)
        
        self.find_lattice_vector_button = QPushButton("Find Lattice Vectors")
        self.find_lattice_vector_button.clicked.connect(self.find_lattice_vectors)

        self.find_offset_button = QPushButton("Find offset")
        self.find_offset_button.clicked.connect(self.find_offset)

        lattice_grid.addWidget(blobtoolbar,0,0,1,3)
        lattice_grid.addWidget(self.blobcanvas,1,0,1,3)
        # lattice_grid.addWidget(ffttoolbar,2,0)
        # lattice_grid.addWidget(self.fftcanvas,3,0)
        lattice_grid.addWidget(self.update_blobs_button,2,0)
        lattice_grid.addWidget(self.find_lattice_vector_button,2,1)
        lattice_grid.addWidget(self.find_offset_button,2,2)

        ### Deconvolved Tab ###
        self.deconvolved_canvas = MplCanvas(self, width=5, height=4, dpi=100)
        #self.find_deconvolved(self.deconvolved_canvas)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        deconvolution_toolbar = NavigationToolbar(self.deconvolved_canvas, self)
        
        deconvolved_layout = QGridLayout()
        deconvolved_layout.addWidget(deconvolution_toolbar,0,0,1,2)
        deconvolved_layout.addWidget(self.deconvolved_canvas,1,0,1,2)

        plot_tab = QWidget()
        plot_tab.setLayout(deconvolved_layout)
        self.tabs.addTab(plot_tab,"Deconvolution")

        # RL interations
        rliterationslabel = QLabel("RL iterations")
        deconvolved_layout.addWidget(rliterationslabel,2,0,1,1)
        self.rliterations_edit = QLineEdit()
        self.rliterations_edit.setText(str(self.RLiterations))
        self.rliterations_edit.returnPressed.connect(self.find_deconvolved)
        self.rliterations_edit.setValidator(int_validator)
        deconvolved_layout.addWidget(self.rliterations_edit,2,1,1,1)

        #Histogram Tab##
        histogram_tab = QWidget()
        histogram_grid = QGridLayout()
        histogram_tab.setLayout(histogram_grid)
        self.tabs.addTab(histogram_tab,"Thresholding")

        self.histogram_canvas = MplCanvas(self)
        histogram_toolbar = NavigationToolbar(self.histogram_canvas,self)

        self.bin_button = QPushButton("Find histogram")
        self.bin_button.clicked.connect(self.find_histogram)

        threshold_label = QLabel("Manual Threshold")
        self.threshold_edit = QLineEdit()
        self.threshold_edit.returnPressed.connect(self.set_threshold)
        self.threshold_edit.setValidator(double_validator)

        histogram_grid.addWidget(histogram_toolbar,0,0,1,2)
        histogram_grid.addWidget(self.histogram_canvas,1,0,1,2)
        histogram_grid.addWidget(self.bin_button,2,0,1,2)
        histogram_grid.addWidget(threshold_label,3,0)
        histogram_grid.addWidget(self.threshold_edit,3,1)

        ### Settings Tab

        image_settings_tab = QWidget()
        image_settings_grid = QGridLayout()
        image_settings_tab.setLayout(image_settings_grid)
        self.tabs.addTab(image_settings_tab, "image Settings")

        # allow user to change window name
        name_label = QLabel('Window name: ', self)
        image_settings_grid.addWidget(name_label, 0,0, 1,1)
        self.name_edit = QLineEdit(self)
        image_settings_grid.addWidget(self.name_edit, 0,1, 1,1)
        self.name_edit.setText(self.name) # default
        self.name_edit.textChanged[str].connect(self.reset_name)

        # get user to set the image size in pixels
        # get user to set the image size in pixels
        # self.pic_width_edit = QLineEdit(self)
        # self.pic_height_edit = QLineEdit(self)
        # for i, label in enumerate([['Image width: ', self.pic_width_edit, 512], 
        #         ['Image height', self.pic_height_edit, 512]]): #HARD CODED 512 for pic height #TODO
        #     textlabel = QLabel(label[0], self)
        #     image_settings_grid.addWidget(textlabel, 1,2*i, 1,1)
        #     image_settings_grid.addWidget(label[1], 1,2*i+1, 1,1)
        #     label[1].setText(str(label[2])) # default
        #     label[1].textChanged.connect(self.pic_size_text_edit)
        #     label[1].setValidator(nat_validator)

        # get image size from loading an image
        # load_im_size = QPushButton('Load size from image', self)
        # load_im_size.clicked.connect(self.load_im_size) # load image size from image
        # load_im_size.resize(load_im_size.sizeHint())
        # image_settings_grid.addWidget(load_im_size, 1,2, 1,1)

        # get user to set ROI:
        # centre of ROI x position
        roi_xc_label = QLabel('ROI x_c: ', self)
        image_settings_grid.addWidget(roi_xc_label, 2,0, 1,1)
        self.roi_x_edit = QLineEdit(self)
        image_settings_grid.addWidget(self.roi_x_edit, 2,1, 1,1)
        self.roi_x_edit.setText('1')  # default
        self.roi_x_edit.textEdited[str].connect(self.roi_text_edit)
        self.roi_x_edit.setValidator(int_validator) # only numbers
        # whether or not the user can change this window's ROI
        self.roi_x_edit.setEnabled(edit_ROI) 
        
        # centre of ROI y position
        roi_yc_label = QLabel('ROI y_c: ', self)
        image_settings_grid.addWidget(roi_yc_label, 2,2, 1,1)
        self.roi_y_edit = QLineEdit(self)
        image_settings_grid.addWidget(self.roi_y_edit, 2,3, 1,1)
        self.roi_y_edit.setText('1')  # default
        self.roi_y_edit.textEdited[str].connect(self.roi_text_edit)
        self.roi_y_edit.setValidator(int_validator) # only numbers
        self.roi_y_edit.setEnabled(edit_ROI)
        
        # ROI size
        roi_l_label = QLabel('ROI size: ', self)
        image_settings_grid.addWidget(roi_l_label, 4,0, 1,1)
        self.roi_l_edit = QLineEdit(self)
        image_settings_grid.addWidget(self.roi_l_edit, 4,1, 1,1)
        self.roi_l_edit.setText('1')  # default
        self.roi_l_edit.textEdited[str].connect(self.roi_text_edit)
        self.roi_l_edit.setValidator(nat_validator) # only numbers
        self.roi_l_edit.setEnabled(edit_ROI)

        # EMCCD bias offset
        bias_offset_label = QLabel('EMCCD bias offset: ', self)
        image_settings_grid.addWidget(bias_offset_label, 5,0, 1,1)
        self.bias_offset_edit = QLineEdit(self)
        image_settings_grid.addWidget(self.bias_offset_edit, 5,1, 1,1)
        self.bias_offset_edit.setText(str(0)) # default
        self.bias_offset_edit.editingFinished.connect(self.CCD_stat_edit)
        self.bias_offset_edit.setValidator(int_validator) # only ints

        # label to show last file analysed
        self.recent_label = QLabel('', self)
        image_settings_grid.addWidget(self.recent_label, 7,0, 1,4)
        
        #### tab for viewing images ####
        im_tab = QWidget()
        im_grid = QGridLayout()
        im_tab.setLayout(im_grid)
        self.tabs.addTab(im_tab, 'Image')
        # display the pic size widgets on this tab as well
        im_size_label = QLabel('Image Size in Pixels: ', self)
        im_grid.addWidget(im_size_label, 0,0, 1,1)
        self.pic_size_label = QLabel('', self)
        im_grid.addWidget(self.pic_size_label, 0,1, 1,1)
        self.pic_size_label.setText('(%s,%s)'%(512,512)) # default

        # toggle to continuously plot images as they come in
        self.im_show_toggle = QPushButton('Auto-display last image', self)
        self.im_show_toggle.setCheckable(True)
        self.im_show_toggle.clicked[bool].connect(self.set_im_show)
        im_grid.addWidget(self.im_show_toggle, 0,2, 1,1)
        
        im_grid_pos = 0 # starting column. 
        # centre of ROI x position
        self.xc_label = QLabel('ROI x_c: 0', self)
        im_grid.addWidget(self.xc_label, 7,im_grid_pos, 1,1)
        
        # centre of ROI y position
        self.yc_label = QLabel('ROI y_c: 0', self)
        im_grid.addWidget(self.yc_label, 7,im_grid_pos+2, 1,1)
        
        # ROI size
        self.l_label = QLabel('ROI size: 1', self)
        im_grid.addWidget(self.l_label, 7,im_grid_pos+4, 1,1)
        
        # display last image if toggle is True
        im_widget = pg.GraphicsLayoutWidget() # containing widget
        viewbox = im_widget.addViewBox() # plot area to display image
        self.im_canvas = pg.ImageItem() # the image
        # # Set colormap for pg plot
        # colormap = matplotlib.cm.get_cmap("inferno")  # cm.get_cmap("CMRmap")
        # colormap._init()
        # lut = (colormap._lut * 255).view(np.ndarray)  # Convert matplotlib colormap from 0-1 to 0 -255 for Qt
        # # Apply the colormap
        # self.im_canvas.setLookupTable(lut)

        viewbox.addItem(self.im_canvas)
        im_grid.addWidget(im_widget, 1,im_grid_pos, 6,8)
        # make an ROI that the user can drag
        self.roi = pg.ROI([0,0], [1,1], movable=False) 
        self.roi.sigRegionChangeFinished.connect(self.user_roi)
        viewbox.addItem(self.roi)
        self.roi.setZValue(10)   # make sure the ROI is drawn above the image
        # make a histogram to control the intensity scaling
        self.im_hist = pg.HistogramLUTItem()
        self.im_hist.setImageItem(self.im_canvas)
        im_widget.addItem(self.im_hist)
        
        # edits to allow the user to fix the intensity limits
        vmin_label = QLabel('Min. intensity: ', self)
        im_grid.addWidget(vmin_label, 8,im_grid_pos, 1,1)
        self.vmin_edit = QLineEdit(self)
        im_grid.addWidget(self.vmin_edit, 8,im_grid_pos+1, 1,1)
        self.vmin_edit.setText('')  # default auto from image
        self.vmin_edit.setValidator(int_validator) # only integers
        vmax_label = QLabel('Max. intensity: ', self)
        im_grid.addWidget(vmax_label, 8,im_grid_pos+2, 1,1)
        self.vmax_edit = QLineEdit(self)
        im_grid.addWidget(self.vmax_edit, 8,im_grid_pos+3, 1,1)
        self.vmax_edit.setText('')  # default auto from image
        self.vmax_edit.setValidator(int_validator) # only integers


        ### Result Tab ##
        result_tab = QWidget()
        result_layout = QGridLayout()
        result_tab.setLayout(result_layout)

        self.result_canvas = MplCanvas()
        result_toolbar = NavigationToolbar(self.result_canvas,self)

        result_layout.addWidget(result_toolbar,0,0)
        result_layout.addWidget(self.result_canvas,1,0)
        
        self.tabs.addTab(result_tab,"Result")
        ### Dashboard Tab TBI ###
        # dash_tab = QWidget()
        # dash_grid = QGridLayout()
        # dash_tab.setLayout(dash_grid)       
        # self.tabs.addTab(dash_tab,"Dashboard")
        
        # # dash_grid.addWidget(self.blobcanvas,0,0)
        # # dash_grid.addWidget(self.deconvolved_canvas,0,1)
        # # dash_grid.addWidget(self.histogram_canvas,1,0)


        """
            Status Bar
        """

        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Hi, welcome to the microscope gui")

        self.setGeometry(100, 150, 850, 500)
        self.setWindowTitle(self.name+'Microscope Imaging')
        self.setWindowIcon(QIcon('molcules.png'))        
    
        self.show()
    def set_im_show(self, toggle):
        """If the toggle is True, always update the widget with the last image."""
        remove_slot(self.event_im, self.update_im, toggle)
    def update_im(self, im, include=True):
        """Receive the image array emitted from the event signal
        display the image in the image canvas.
        event_im: [image (np.ndarray), include? (bool)]"""
        self.im_canvas.setImage(im)
        vmin, vmax = np.min(im), np.max(im)
        if self.vmin_edit.text():
            vmin = int(self.vmin_edit.text())
        if self.vmax_edit.text():
            vmax = int(self.vmax_edit.text())
        self.im_hist.setLevels(vmin, vmax)
    # def load_from_db(self):
    #     self.db_document = imgs.find_one({"_id":ObjectId(self.loadbyid_edit.text())})
    #     if type(self.db_document):
    #         print(self.db_document.keys())
    #         self.current_image = (np.frombuffer(self.db_document["raw_image"],dtype =np.float64))
    #         self.current_image.resize(self.db_document["shape"])
    #         self.truthlist =np.frombuffer(self.db_document["truthandlattice"],dtype = np.float64)
    #         self.truthlist.resize((self.db_document["latticeshape"][0],3))
    #         self.psf_radius = float(self.db_document["psf_radius"])
    #         self.status_bar.showMessage("Image loaded")
    #         self.status_bar.repaint()
    #         for key in self.updated.keys():
    #                 self.updated[key] = False
    #         self.updated['current_image'] = True
    #         if not self.use_lattice_file_checkbox.isChecked():
    #                 self.lattice_image = self.current_image
                
    #         if self.auto_checkbox.isChecked():
    #             self.find_blobs()
    #     else: 
    #         self.status_bar.showMessage("No match found :'( ")
    #         self.status_bar.repaint()
    def load_image_from_file(self):
        """Open a file dialog to get image. 
        """


        self.image_path = QFileDialog.getOpenFileName(self,"openfile")[0]
        print(self.image_path)
    
        if self.image_path.endswith("data.npy"):
            self.current_image = np.load(self.image_path)
            print(self.image_path[:-8])
            self.truthlist = np.load(self.image_path[:-8]+"truthandlattice.npy")
            self.status_bar.showMessage("{} loaded sucessfully".format(os.path.basename(self.image_path)))
            self.status_bar.repaint()
            for key in self.updated.keys():
                self.updated[key] = False
            self.updated['current_image'] = True
            self.file_path_label.setText(str(self.image_path))
            if not self.use_lattice_file_checkbox.isChecked():
                self.lattice_image = self.current_image
            
            if self.auto_checkbox.isChecked():
                self.find_blobs()
        else:
            self.status_bar.showMessage("WARNING: path doesn't end in data.npy!")
            self.status_bar.repaint()
    def CCD_stat_edit(self, emg=1, pag=4.5, Nr=8.8, acq_change=False):
        #TODO implement CCD logging   
        """Update the values used for the EMCCD bias offset and readout noise"""
        # if self.bias_offset_edit.text(): # check the label isn't empty
        #     self.image_handler.bias = int(self.bias_offset_edit.text())

        # if acq_change: # If the acquisition settings have been changed by the camera
        #     self.histo_handler.emg, self.histo_handler.pag, self.histo_handler.Nr = emg, pag, Nr
        #     self.histo_handler.dg = 2.0 if self.histo_handler.emg > 1 else 1.0 # multiplicative noise factor
        
    def set_lattice_file(self):
        #Specific to microsocpe
        self.lattice_file_path = QFileDialog.getOpenFileName(self,"Set Lattice file")[0]
        if self.lattice_file_path.endswith("data.npy"):
            self.lattice_image = np.load(self.lattice_file_path)
            self.status_bar.showMessage("{} loaded sucessfully as Lattice image".format(os.path.basename(self.lattice_file_path)))
            self.status_bar.repaint()
            self.lattice_file_path_label.setText(self.lattice_file_path)
            if self.auto_checkbox.isChecked():
                self.find_blobs()
        else:
            self.status_bar.showMessage("WARNING: path doens't end in data.npy")
            self.status_bar.repaint()
    def find_blobs(self):
        """Wrapper for ai.find_blobs 
        """
        print("Finding blobs...")
        self.status_bar.showMessage("Finding blobs...")
        self.status_bar.repaint()
        blobworker = Worker(ai.find_blobs,self.lattice_image)
        blobworker.signals.result.connect(partial(self.connectSelf,'blobs'))
        blobworker.signals.finished.connect(self.show_blobs)
        self.threadpool.start(blobworker) 
    def show_blobs(self):
        """ Plot update function for blobs
        """
        self.updated['blobs'] = True
        self.blobcanvas.axes.cla()
        self.blobcanvas.axes.imshow(self.lattice_image)
        self.status_bar.showMessage("Blobs found; Plotting Blobs")
        self.status_bar.repaint()
        print("blob plotting")

        blobarray = np.array(self.blobs)
        self.blobcanvas.axes.scatter(blobarray[:,0],blobarray[:,1],s = 100, facecolors='none', edgecolors='white')
        self.blobcanvas.draw()
        self.status_bar.showMessage("Blobs Updated")
        self.status_bar.repaint()
        if self.auto_checkbox.isChecked():
            self.find_lattice_vectors()
    def find_lattice_vectors(self):
        """
        """
        print("Finding lattice...")
        self.status_bar.showMessage("Finding Lattice angles...")
        self.status_bar.repaint()
        worker  = Worker(ai.find_lattice_vectors_new,self.blobs,float(self.angle_min_edit.text()),float(self.angle_max_edit.text()),int(self.angle_step_edit.text()))#,Plot = self.full_plots_checkbox.isChecked())
        worker.signals.result.connect(partial(self.connectSelf,'lattice_vectors'))
        worker.signals.result.connect(self.show_lattice_vectors)
        self.threadpool.start(worker)
    def show_lattice_vectors(self):
        self.updated['lattice_vectors'] = True
        self.status_bar.showMessage("Lattice vectors Found... See terminal")
        self.status_bar.repaint()
        print(self.lattice_vectors)
        if self.auto_checkbox.isChecked():
            self.find_offset()
    def find_offset(self):
        print("Finding lattice...")
        self.status_bar.showMessage("Finding Lattice Phase...")
        self.status_bar.repaint()
        worker = Worker(ai.find_offset,self.current_image,self.blobs,self.lattice_vectors)
        worker.signals.result.connect(partial(self.connectSelf,'offset'))
        worker.signals.finished.connect(self.show_offset)
        self.threadpool.start(worker)
    def manual_lattice(self):
        angle = float(self.manual_angle.text())
        self.lattice_vectors = np.array(((np.cos(angle)*self.lattice_constant,np.sin(angle)*self.lattice_constant),(-np.sin(angle)*self.lattice_constant,np.cos(angle)*self.lattice_constant)))
        self.offset = (float(self.manual_offset.text()),float(self.manual_offset1.text()))
        self.show_offset()
    def manual_psf(self):
        self.psf_radius = float(self.manual_psf_edit.text())

    def show_offset(self):
        print("Offset",self.offset)
        self.updated['lattice'] = True
        self.status_bar.showMessage("Lattice offset found... See terminal")
        self.status_bar.repaint()
        self.blobcanvas.axes.cla()
        self.lattice = ai.generate_lattice(self.current_image.shape,self.lattice_vectors,self.offset)
        self.blobcanvas.axes.imshow(self.current_image)
        self.blobcanvas.axes.scatter(self.lattice[:,0],self.lattice[:,1],color = "cyan",s = 0.5)
        self.blobcanvas.draw()
        if self.auto_checkbox.isChecked():
            self.find_deconvolved()
    # def plot_lattice(self,canvas):
    #     for point in self.lattice:
    #         for point, truth, atom in zip(unscaledlattice,truthlist,atoms):
    #         if atom == truth:
    #             col = "green"
    #         else:
    #             col = "white"
    #             wrong+=1
    #         c = plt.Circle(point,radius=0.3,color = col)
    #         ax.add_patch(c)
    #         for i in range(2):
    #             if i ==0:
    #                 j = 1
    #             else:
    #                 j = 0
    #             plt.plot([point[0]+latticevectors[0,i]/2,point[0]+latticevectors[0,i]/2+latticevectors[0,j]],[point[1]+latticevectors[1,i]/2,point[1]+latticevectors[1,i]/2+latticevectors[1,j]],color = "cyan")


    # def show_fft(self):
    #     self.fftcanvas.axes.cla()
    #     self.fftcanvas.axes.text(0.5,0.5,str(self.lattice_vectors))
  
    # Deconvolution functions
    def find_deconvolved(self):
        
        self.status_bar.showMessage("D e c o n v o l u t i o n ...")
        self.status_bar.repaint()
        self.RLiterations = int(self.rliterations_edit.text())
        print("Devonvolving {} iterations".format(self.RLiterations))
        worker = Worker(self.deconvolve,self.current_image,self.RLiterations)
        worker.signals.result.connect(partial(self.connectSelf,'deconvolved'))
        worker.signals.finished.connect(self.show_deconvolved)
        self.threadpool.start(worker) 
    def show_deconvolved(self):
        self.updated['deconvolved'] = True
        self.deconvolved_canvas.axes.cla()
        #print(self.deconvolved)
        self.deconvolved_canvas.axes.imshow(self.deconvolved)
        self.deconvolved_canvas.draw()
        self.status_bar.showMessage("Finished updating plots")
        if self.auto_checkbox.isChecked() and self.updated['lattice']== True:
            self.find_histogram()  
    def deconvolve(self,image,iterations):
        psf = ai.airy_psf(self.window,self.psf_radius)
        image = np.pad(image,int(psf.shape[0]/2),mode= "reflect")
        deconvolved = restoration.richardson_lucy(image,psf,iterations,False)
        deconvolved = util.crop(deconvolved,int(psf.shape[0]/2))
        return deconvolved        
    #Microscope histogram Functions
    def find_histogram(self):
        self.status_bar.showMessage("Binning image...")
        self.status_bar.repaint()

        worker = Worker(ai.hist_and_thresh,self.deconvolved,self.lattice, self.lattice_vectors)
        worker.signals.result.connect(self.storeHist)
        worker.signals.finished.connect(self.show_histogram)
        self.threadpool.start(worker)
    def storeHist(self,data):
        "Specific function to store histogram data"
        self.updated['threshold'] = True
        self.histogram_data = data[0]
        self.threshold = data[1]
        self.threshold_edit.setText(str(self.threshold))
        self.atoms = data[2]
        # print("atoms shape",self.atoms.shape)
    def set_threshold(self):
        self.threshold = float(self.threshold_edit.text())
        self.atoms = (np.array(self.histogram_data)> self.threshold)
        self.show_histogram()
    def show_histogram(self):
        """Update histogram canvas. Slot for hist_and_thresh finished
        """
        self.status_bar.showMessage("Histogram Found")
        self.status_bar.repaint()
        canvas = self.histogram_canvas
        canvas.axes.cla()
        canvas.axes.hist(self.histogram_data, 100)
        canvas.axes.axvline(float(self.threshold),color = "red")
        canvas.axes.set_xlabel("Counts")
        canvas.axes.set_ylabel("Occurances")
        canvas.draw()

        self.result_canvas.axes.cla()
        self.result_canvas.axes.imshow(self.current_image)
        self.result_canvas.axes.scatter(self.lattice[:,0],self.lattice[:,1],s = 2,c = self.atoms,cmap = "Reds")

        # if self.auto_checkbox.isChecked():
        #     self.find_fidelity()
    

   

    def plot_raw_image(self,canvas):
        im = np.load(self.image_path)
        canvas.axes.imshow(im)
        canvas.draw()
    



    def connectSelf(self,variable,signal):
        """ use setattr to assign a signal to a self variable
            use as signal.connect(partial(self.connectSelf),'variablename') 
        """
        setattr(self,variable,signal)
    def reset_name(self, text=''):
        """Take a new name for the window, display it in the window title."""
        self.name = self.name_edit.text()
        self.setWindowTitle(self.name+' - Microscope Processing -')
    def roi_text_edit(self, text):
        """Update the ROI position and size every time a text edit is made by
        the user to one of the line edit widgets"""
        xc, yc, l = [self.roi_x_edit.text(),
                            self.roi_y_edit.text(), self.roi_l_edit.text()]
        if any(v == '' for v in [xc, yc, l]):
            xc, yc, l = 1, 1, 1 # default 
        else:
            xc, yc, l = list(map(int, [xc, yc, l])) # crashes if the user inputs float
            
        if any(v > max(self.pic_width, self.pic_height) for v in [xc, yc, l]):
            xc, yc, l = 1, 1, 1
        
        if (xc - l//2 < 0 or yc - l//2 < 0 
            or xc + l//2 > self.pic_width 
            or yc + l//2 > self.pic_height):
            l = 2*min([xc, yc])  # can't have the boundary go off the edge
        if int(l) == 0:
            l = 1 # can't have zero width
        # self.image_handler.set_roi(dimensions=list(map(int, [xc, yc, l])))
        self.xc_label.setText('ROI x_c = '+str(xc)) 
        self.yc_label.setText('ROI y_c = '+str(yc))
        self.l_label.setText('ROI size = '+str(l))
        # update ROI on image canvas
        # note: setting the origin as top left because image is inverted
        self.roi.setPos(xc - l//2, yc - l//2)
        self.roi.setSize((l, l))

    def user_roi(self, pos):
        """Update position of ROI"""
        x0, y0 = self.roi.pos()  # lower left corner of bounding rectangle
        xw, yw = self.roi.size() # widths
        l = int(0.5*(xw+yw))  # want a square ROI
        # note: setting the origin as bottom left but the image has origin top left
        xc, yc = int(x0 + l//2), int(y0 + l//2)  # centre
        # self.image_handler.set_roi(dimensions=[xc, yc, l])
        self.xc_label.setText('ROI x_c = '+str(xc)) 
        self.yc_label.setText('ROI y_c = '+str(yc))
        self.l_label.setText('ROI size = '+str(l))
        self.roi_x_edit.setText(str(xc))
        self.roi_y_edit.setText(str(yc))
        self.roi_l_edit.setText(str(l))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = main_window()
    sys.exit(app.exec_())
