
"""Networking - Sequence Previewer
Stefan Spence 17/12/19

 - Preview DExTer sequences and facilitate loading/saving
The functions here are specific to the format of 
sequences that DExTer generates.

Note: to make it faster we could use xml.dom.minidom 
instead of python dictionaries. Since LabVIEW uses
some unicode characters, would need to parse it like
with open('filename', 'r') as f:
 dm = xml.dom.minidom.parseString(f.read().replace('\n','').replace('\t','').encode('utf-8'))
"""
import sys
import numpy as np
try:
    from PyQt4.QtCore import QThread, pyqtSignal, QEvent, QRegExp, QTimer, Qt
    from PyQt4.QtGui import (QApplication, QPushButton, QWidget, QLabel,
        QAction, QGridLayout, QMainWindow, QMessageBox, QLineEdit, QIcon, 
        QFileDialog, QDoubleValidator, QIntValidator, QComboBox, QMenu, 
        QActionGroup, QTabWidget, QVBoxLayout, QHBoxLayout, QFont, QRegExpValidator, 
        QInputDialog, QTableWidget, QTableWidgetItem, QScrollArea) 
except ImportError:
    from PyQt5.QtCore import QThread, pyqtSignal, QEvent, QRegExp, QTimer, Qt
    from PyQt5.QtGui import QFont, QIcon
    from PyQt5.QtWidgets import (QApplication, QPushButton, QWidget, 
        QTabWidget, QAction, QMainWindow, QLabel, QInputDialog, QGridLayout,
        QMessageBox, QLineEdit, QFileDialog, QComboBox, QActionGroup, QMenu,
        QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem, QScrollArea)
from translator import translate
from multirunEditor import multirun_widget
import logging
logger = logging.getLogger(__name__)


def bl(string):
    """Convert a string of a boolean to a boolean.
    This corrects for bool('0')=True."""
    try: return bool(int(string))
    except ValueError: return bool(string)

#### #### Edit sequences #### ####

# class Editor(QMainWindow):
#     """Provide a GUI for quickly editing DExTer sequences.
#     """
#     def __init__(self, num_steps=1):
#         super().__init__()
#         self.tr = translate(num_steps)
#         self.pre = Previewer(self.tr)
#         self.init_UI()

#     def make_label_edit(self, label_text, layout, position=[0,0, 1,1],
#             default_text='', validator=None):
#         """Make a QLabel with an accompanying QLineEdit and add them to the 
#         given layout with an input validator. The position argument should
#         be [row number, column number, row width, column width]."""
#         label = QLabel(label_text, self)
#         layout.addWidget(label, *position)
#         line_edit = QLineEdit(self)
#         if np.size(position) == 4:
#             position[1] += 1
#         layout.addWidget(line_edit, *position)
#         line_edit.setText(default_text) 
#         line_edit.setValidator(validator)
#         return label, line_edit
        
#     def init_UI(self):
#         """Create all of the widget objects required"""
#         self.centre_widget = QWidget()
#         self.centre_widget.layout = QGridLayout()
#         self.centre_widget.setLayout(self.centre_widget.layout)
#         self.setCentralWidget(self.centre_widget)
        
#         #### validators for user input ####
#         # reg_exp = QRegExp(r'([0-9]+(\.[0-9]+)?,?)+')
#         # comma_validator = QRegExpValidator(reg_exp) # floats and commas
#         double_validator = QDoubleValidator() # floats
#         int_validator = QIntValidator()       # integers
        
#         #### menubar at top gives options ####
#         # menubar = self.menuBar()
#         # show_windows = menubar.addMenu('Windows')
#         # menu_items = []
#         # for window_title in ['Image Analyser', 'Camera Status', 
#         #     'Image Saver', 'Monitoring']:
#         #     menu_items.append(QAction(window_title, self)) 
#         #     menu_items[-1].triggered.connect(self.show_window)
#         #     show_windows.addAction(menu_items[-1])

#         #### choose event indices ####
#         # by name
#         # by index
#         self.make_label_edit('Event index', self.centre_widget.layout, 
#             position=[1,0, 1,1], default_text='0', validator=int_validator)
        
#         #### choose channel ####
#         self.make_label_edit('Channel', self.centre_widget.layout, 
#             position=[2,0, 1,1], default_text='')

#         #### choose new value ####
#         self.make_label_edit('New value', self.centre_widget.layout, 
#             position=[3,0, 1,1], default_text='0', validator=double_validator)

#         #### preview sequence ####
#         self.preview_button = QPushButton('Preview sequence', self)
#         self.preview_button.resize(self.preview_button.sizeHint())
#         self.preview_button.clicked.connect(self.pre.show)
#         self.centre_widget.layout.addWidget(self.preview_button, 5,0, 1,1)

#         #### save to file ####
        
#         #### choose main window position and dimensions: (xpos,ypos,width,height)
#         self.setGeometry(60, 60, 900, 800)
#         self.setWindowTitle('DExTer Sequence Editor')
#         self.setWindowIcon(QIcon('docs/translatoricon.png'))

#### #### Preview sequences #### ####

class Previewer(QMainWindow):
    """Provide a display of a sequence, reminiscent
    of DExTer main view.
    """
    def __init__(self, tr=translate()):
        super().__init__()
        self.tr = tr
        self.init_UI()
        self.set_sequence()
    
    def reset_table(self, table, digital=1):
        """Set empty table items in all of the cells of the
        given table. 
        digital -- 1: Set the background colour red
                -- 0: Set the text as ''."""
        for i in range(table.rowCount()):
            for j in range(table.columnCount()):
                table.setItem(i, j, QTableWidgetItem())
                if digital:
                    table.item(i, j).setBackground(Qt.red)
                else:
                    table.item(i, j).setText('')

    def init_UI(self):
        """Create all of the widget objects required"""
        self.centre_widget = QWidget()
        self.tabs = QTabWidget()       # make tabs for each main display 
        self.centre_widget.layout = QVBoxLayout()
        self.centre_widget.layout.addWidget(self.tabs)
        self.centre_widget.setLayout(self.centre_widget.layout)
        self.setCentralWidget(self.centre_widget)
        
        num_e = len(self.tr.seq_dic['Event list array in'])
        num_s = len(self.tr.seq_dic['Experimental sequence cluster in']['Sequence header top'])
        menubar = self.menuBar()

        # save/load a sequence file
        menubar.clear() # prevents recreating menubar if init_UI() is called again 
        seq_menu = menubar.addMenu('Sequence')
        load = QAction('Load Sequence', self) 
        load.triggered.connect(self.load_seq_from_file)
        seq_menu.addAction(load)
        save = QAction('Save Sequence', self) 
        save.triggered.connect(self.save_seq_file)
        seq_menu.addAction(save)

        #### tab for previewing sequences ####
        preview_tab = QWidget()
        prv_layout = QVBoxLayout()
        preview_tab.setLayout(prv_layout)
        scroll_widget = QWidget()
        prv_layout.addWidget(scroll_widget)
        prv_vbox = QVBoxLayout()
        scroll_widget.setLayout(prv_vbox)
        self.tabs.addTab(preview_tab, "Sequence")

        # position the widgets on the layout:
        # metadata
        self.routine_name = QLabel('', self)
        self.routine_desc = QLabel('', self)
        for label, name in [[self.routine_name, 'Routine name: '], 
                [self.routine_desc, 'Routine description: ']]:
            layout = QHBoxLayout()
            title = QLabel('Routine name: ', self)
            title.setFixedWidth(200)
            layout.addWidget(title)    
            label.setStyleSheet('border: 1px solid black')
            label.setFixedWidth(400)
            layout.addWidget(label)
            prv_vbox.addLayout(layout)

        # list of event descriptions
        self.e_list = QTableWidget(4, num_e)
        self.e_list.setVerticalHeaderLabels(['Event name: ', 
            'Routine specific event? ', 'Event indices: ', 'Event path: '])
        self.e_list.setFixedHeight(250)
        self.reset_table(self.e_list, 0)
        prv_vbox.addWidget(self.e_list)
        
        # event header top 
        self.head_top = QTableWidget(14, num_s)
        self.head_top.setVerticalHeaderLabels(['Skip Step: ', 
            'Event name: ', 'Hide event steps: ', 
            'Event ID: ', 'Time step name: ', 'Populate multirun: ',
            'Time step length: ', 'Time unit: ', 'D/A trigger: ',
            'Trigger this time step? ', 'Channel: ', 'Analogue voltage (V): ',
            'GPIB event name: ', 'GPIB on/off? '])
        self.head_top.setFixedHeight(600)
        self.reset_table(self.head_top, 0)
        prv_vbox.addWidget(self.head_top)
          
        # fast digital channels
        fd_head = QLabel('Fast Digital', self) 
        prv_vbox.addWidget(fd_head)

        self.fd_chans = QTableWidget(self.tr.nfd, num_s)
        self.fd_chans.setFixedHeight(400)
        self.reset_table(self.fd_chans, 1)
        prv_vbox.addWidget(self.fd_chans)
          
        # fast analogue channels
        fa_head = QLabel('Fast Analogue', self) 
        prv_vbox.addWidget(fa_head)
        self.fa_chans = QTableWidget(self.tr.nfa, num_s*2)
        self.fa_chans.setFixedHeight(330)
        self.reset_table(self.fa_chans, 0)
        prv_vbox.addWidget(self.fa_chans)
        
        # event header middle
        self.head_mid = QTableWidget(14, num_s)
        self.head_mid.setVerticalHeaderLabels(['Skip Step: ', 
            'Event name: ', 'Hide event steps: ', 
            'Event ID: ', 'Time step name: ', 'Populate multirun: ',
            'Time step length: ', 'Time unit: ', 'D/A trigger: ',
            'Trigger this time step? ', 'Channel: ', 'Analogue voltage (V): ',
            'GPIB event name: ', 'GPIB on/off? '])
        self.head_mid.setFixedHeight(560)
        self.reset_table(self.head_mid, 0)
        prv_vbox.addWidget(self.head_mid)
        
        # slow digital channels
        sd_head = QLabel('Slow Digital', self) 
        prv_vbox.addWidget(sd_head)

        self.sd_chans = QTableWidget(self.tr.nsd, num_s)
        self.sd_chans.setFixedHeight(400)
        self.reset_table(self.sd_chans, 1)
        prv_vbox.addWidget(self.sd_chans)
        
        # slow analogue channels
        sa_head = QLabel('Slow Analogue', self) 
        prv_vbox.addWidget(sa_head)

        self.sa_chans = QTableWidget(self.tr.nsa, num_s*2)
        self.sa_chans.setFixedHeight(400)
        self.reset_table(self.sa_chans, 0)
        prv_vbox.addWidget(self.sa_chans)
        
        # place scroll bars if the contents of the window are too large
        scroll = QScrollArea(self)
        scroll.setWidget(scroll_widget)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(800)
        prv_layout.addWidget(scroll)

        #### tab for multi-run settings ####
        self.mr = multirun_widget(self.tr)
        self.tabs.addTab(self.mr, "Multirun")

        mr_menu = menubar.addMenu('Multirun')
        mrload = QAction('Load Parameters', self) 
        mrload.triggered.connect(self.mr.load_mr_params)
        mr_menu.addAction(mrload)
        mrsave = QAction('Save Parameters', self) 
        mrsave.triggered.connect(self.mr.save_mr_params)
        mr_menu.addAction(mrsave)
        
        # choose main window position and dimensions: (xpos,ypos,width,height)
        self.setGeometry(60, 60, 1000, 800)
        self.setWindowTitle('Sequence Preview')
        self.setWindowIcon(QIcon('docs/previewicon.png'))


    def reset_UI(self):
        """After loading in a new sequence, adjust the UI
        so that the tables have the right number of rows and columns. """
        num_e = len(self.tr.seq_dic['Event list array in'])
        num_s = len(self.tr.seq_dic['Experimental sequence cluster in']['Sequence header top'])
        for table, rows, cols, dig in [[self.e_list, 4, num_e, 0], [self.head_top, 14, num_s, 0],
            [self.fd_chans, self.tr.nfd, num_s, 1], [self.fa_chans, self.tr.nfa, num_s*2, 0],
            [self.head_mid, 14, num_s, 0], [self.sd_chans, self.tr.nsd, num_s, 1],
            [self.sa_chans, self.tr.nsa, num_s*2, 0]]:
            table.setRowCount(rows)
            table.setColumnCount(cols)
            self.reset_table(table, dig)
        
        self.mr.reset_sequence(self.tr)
        
    def try_browse(self, title='Select a File', file_type='all (*)', 
                open_func=QFileDialog.getOpenFileName):
        """Open a file dialog and retrieve a file name from the browser.
        title: String to display at the top of the file browser window
        default_path: directory to open first
        file_type: types of files that can be selected
        open_func: the function to use to open the file browser"""
        try:
            if 'PyQt4' in sys.modules:
                file_name = open_func(self, title, '', file_type)
            elif 'PyQt5' in sys.modules:
                file_name, _ = open_func(self, title, '', file_type)
            return file_name
        except OSError: return '' # probably user cancelled

    def load_seq_from_file(self, fname=''):
        """Choose a file name, load the sequence and then show it in the previewer."""
        if not fname: fname = self.try_browse(file_type='XML (*.xml);;all (*)')
        if fname:
            self.tr.load_xml(fname)
            self.reset_UI()
            self.set_sequence()

    def save_seq_file(self, fname=''):
        """Save the current sequence to an xml file."""
        if not fname: fname = self.try_browse(title='Choose a file name', 
                file_type='XML (*.xml);;all (*)', open_func=QFileDialog.getSaveFileName)
        if fname:
            self.tr.write_to_file(fname)

    def set_sequence(self):
        """Fill the labels with the values from the sequence"""
        seq = self.tr.seq_dic
        self.routine_name.setText(seq['Routine name in'])
        self.routine_desc.setText(seq['Routine description in'])
        ela = seq['Event list array in'] # shorthand
        esc = seq['Experimental sequence cluster in']
        self.fd_chans.setVerticalHeaderLabels(map(str.__add__, esc['Fast digital names']['Hardware ID'],
                [': '+name if name else '' for name in esc['Fast digital names']['Name']]))
        self.fa_chans.setVerticalHeaderLabels(map(str.__add__, esc['Fast analogue names']['Hardware ID'],
                [': '+name if name else '' for name in esc['Fast analogue names']['Name']]))
        self.sd_chans.setVerticalHeaderLabels(map(str.__add__, esc['Slow digital names']['Hardware ID'],
                [': '+name if name else '' for name in esc['Slow digital names']['Name']]))
        self.sa_chans.setVerticalHeaderLabels(map(str.__add__, esc['Slow analogue names']['Hardware ID'],
                [': '+name if name else '' for name in esc['Slow analogue names']['Name']]))
        for i in range(len(ela)):
            self.e_list.item(0, i).setText(ela[i]['Event name'])
            self.e_list.item(1, i).setText(str(ela[i]['Routine specific event?']))
            self.e_list.item(2, i).setText(','.join(map(str, ela[i]['Event indices'])))
            self.e_list.item(3, i).setText(ela[i]['Event path'])
        for i in range(len(esc['Sequence header top'])):
            for j, key in enumerate(['Skip Step', 'Event name', 'Hide event steps', 
                    'Event ID', 'Time step name', 'Populate multirun', 'Time step length', 
                    'Time unit', 'Digital or analogue trigger?', 'Trigger this time step?', 
                    'Channel', 'Analogue voltage (V)', 'GPIB event name', 'GPIB on/off?']):
                self.head_top.item(j, i).setText(str(esc['Sequence header top'][i][key]))
                self.head_mid.item(j, i).setText(str(esc['Sequence header middle'][i][key]))
            for j in range(self.tr.nfd):
                self.fd_chans.item(j, i).setBackground(Qt.green if bl(esc['Fast digital channels'][i][j]) else Qt.red)
            for j in range(self.tr.nfa):
                self.fa_chans.item(j, 2*i).setText(str(esc['Fast analogue array'][j]['Voltage'][i]))
                self.fa_chans.item(j, 2*i+1).setText(
                    'Ramp' if bl(esc['Fast analogue array'][j]['Ramp?'][i]) else '')
            for j in range(self.tr.nsd):
                self.sd_chans.item(j, i).setBackground(Qt.green if bl(esc['Slow digital channels'][i][j]) else Qt.red)
            for j in range(self.tr.nsa):
                self.sa_chans.item(j, 2*i).setText(str(esc['Slow analogue array'][j]['Voltage'][i]))
                self.sa_chans.item(j, 2*i+1).setText(
                    'Ramp' if bl(esc['Slow analogue array'][j]['Ramp?'][i]) else '')


    def choose_multirun_dir(self):
        """Allow the user to choose the directory where the histogram .csv
        files and the measure .dat file will be saved as part of the multi-run"""
        try:
            dir_path = QFileDialog.getExistingDirectory(self, "Select Directory", '')
            self.multirun_save_dir.setText(dir_path)
        except OSError:
            pass # user cancelled - file not found
        

####    ####    ####    #### 

def run():
    """Initiate an app to run the program
    if running in Pylab/IPython then there may already be an app instance"""
    app = QApplication.instance()
    standalone = app is None # false if there is already an app instance
    if standalone: # if there isn't an instance, make one
        app = QApplication(sys.argv) 
        
    boss = Previewer()
    boss.show()
    if standalone: # if an app instance was made, execute it
        sys.exit(app.exec_()) # when the window is closed, the python code also stops
   
if __name__ == "__main__":
    run()