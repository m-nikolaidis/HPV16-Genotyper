################################################################################
##
## BY: WANDERSON M.PIMENTA
## PROJECT MADE WITH: Qt Designer and PySide2
## V: 1.0.0
##
## This project can be used freely for all uses, as long as they maintain the
## respective credits only in the Python scripts, any information in the visual
## interface (GUI) can be modified without any implication.
##
## There are limitations on Qt licenses if you want to use your products
## commercially, I recommend reading them on the official website:
## https://doc.qt.io/qtforpython/licenses.html
##
################################################################################

import sys

import files_rc
########## New imports
import os
import re
import sys
# from typing import Sequence
import wrapper # Tool module
import simplot # Tool module
import pathlib
import logging
import pandas as pd
import pyqtgraph as pg
import plotly.express as px
from datetime import date
from Bio import SeqIO
from QLed import QLed
from PyQt5.QtGui import (QPixmap, 
	QFont, QFontDatabase, QIcon, QPalette,
	QBrush, QColor
)
from PyQt5.QtCore import ( Qt, QObject, QThread, 
	pyqtSignal, QUrl, QCoreApplication, QSize, QRect
	
)
from PyQt5.QtWidgets import (
	QApplication, QLabel, QMainWindow,
	QMenuBar, QSizePolicy, QWidget, QMenu, QAction, 
	QStatusBar, QTableWidget, QTableWidgetItem, 
	QFileDialog, QTabWidget, QGridLayout, QVBoxLayout,
	QPushButton, QScrollArea, QHBoxLayout, QButtonGroup,
	QPlainTextEdit, QMessageBox, QListWidget, QSpinBox,
	QLineEdit, QFrame, QStackedWidget, QSpacerItem, QCheckBox,
	QRadioButton, QSlider, QScrollBar, QComboBox,
	QCommandLinkButton, QAbstractScrollArea, QAbstractItemView

)

from PyQt5.QtWebEngineWidgets import QWebEngineView

from ete3 import ( 
	Tree, TreeStyle, TextFace
)

#TODO: Use QErrorMessage for errors

class Worker(QObject):
	"""Use worker classes to give a multithreading behaviour to the application
	So it does not freeze when a long task is executed
	Step 1: Create a worker class
	Step 2: Create a QThread object
	Step 3: Create a worker object
	Step 4: Move worker to the thread
	Step 5: Connect signals and slots
	Step 6: Start the thread

	"""
	def __init__(self, parent=None):
		super().__init__(parent)
		self.logW = LogWindow()
		
	finished = pyqtSignal()

	def updateWidgets(self, incident):
		if incident == "starting":
			self.logW.setCentralWidget(self.logW.pipelineWidget)
		if incident == "finished":
			self.logW.setCentralWidget(self.logW.finishedWidget)
			
		
	def runPipeline(self):
		paramsdf = mainW.paramsdf
		# TODO: Uncomment after testing
		wrapper.main(paramsdf)
		self.finished.emit()		
		return

class MainWindow(QMainWindow):
	def __init__(self):
		QMainWindow.__init__(self)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)

		startSize = QSize(1000, 720)
		self.resize(startSize)
		self.setMinimumSize(startSize)

		self.ui.stackedWidget.setMinimumWidth(20)
		GuiFunctions.addNewMenu(self, "Main Page", "homeButton", "url(:/16x16/icons/16x16/cil-home.png)", True)
		GuiFunctions.addHomeButtons(self)


		# Select home menu on init
		GuiFunctions.selectStandardMenu(self, "homeButton")

		self.ui.stackedWidget.setCurrentWidget(self.ui.homePage)

		# Remove the content margins from the window
		self.ui.horizontalLayout.setContentsMargins(0, 0, 0, 0)
		
		# TODO: Possible error in windows os from the initialization of the default parameters (forward slash)
		params = wrapper._defaultparams()
		self.paramsdf = pd.DataFrame.from_dict(params,orient='index')
		self.paramsdf.rename(columns={0:"Value"},inplace=True)

		self.show()
		# To keep the pages that are opened as extra windows
		self.openWindows = {}

	def Button(self):
		"""
		Choose what happens once i click any button
		"""
		btnWidget = self.sender()

		if btnWidget.objectName() == "homeButton":
			self.ui.stackedWidget.setCurrentWidget(self.ui.homePage)
			GuiFunctions.resetStyle(self, "homeButton")
			GuiFunctions.updatePageLabel(self, "Home")
			btnWidget.setStyleSheet(GuiFunctions.selectMenu(btnWidget.styleSheet()))

		if btnWidget.objectName() == "resultsButton":
			self.ui.stackedWidget.setCurrentWidget(self.ui.resultsPage)
			# Will customise ui.resultsPage to results
			GuiFunctions.resetStyle(self, "resultsButton")
			GuiFunctions.updatePageLabel(self, "Results")
			btnWidget.setStyleSheet(GuiFunctions.selectMenu(btnWidget.styleSheet()))
		
		if btnWidget.objectName() == "loadFastaButton":
			GuiFunctions.getFastaFile(self)
		if btnWidget.objectName() == "selectOutdirButton":
			GuiFunctions.getOutdir(self)
		if btnWidget.objectName() == "loadResultsButton":
			GuiFunctions.getResdir(self)
		if btnWidget.objectName() == "runPipelineButton":
			GuiFunctions.execPipeline(self)
		if btnWidget.objectName() == "openHelpButton":
			pass
			# GuiFunctions.
			# TODO: Implement
	
	def closeEvent(self, event):
		self.quitMsg = QMessageBox()
		reply = self.quitMsg.question(self, "Close Window", 
		"Do you want to close this window?\nThe application will be terminated",
		QMessageBox.Yes | QMessageBox.No
		)

		if reply == QMessageBox.Yes:
			event.accept()
			logging.shutdown()
		else:
			event.ignore()

class GuiFunctions(MainWindow):
	
	def toggleMenu(self, enable):
		if enable:
			pass

	# LABEL TITLE
	def labelTitle(self, text):
		self.ui.label_title_bar_top.setText(text)

	# Dynamic menus
	def addNewMenu(self, name, objName, icon, isTopMenu):
		font = QFont()
		font.setFamily(u"Segoe UI")
		button = QPushButton(str(1),self)
		button.setObjectName(objName) # TODO: This can help me with the tabs in my original version
		button.setMinimumSize(QSize(0, 70))
		button.setLayoutDirection(Qt.LeftToRight)
		button.setFont(font)
		button.setStyleSheet(Style.style_bt_standard.replace('ICON_REPLACE', icon))
		button.setText(name)
		button.setToolTip(name)
		button.clicked.connect(self.Button)

		if isTopMenu:
			self.ui.menusLayout.addWidget(button)
		else:
			self.ui.layout_menu_bottom.addWidget(button)

	## Select / deselect menus
	def selectMenu(getStyle):
		select = getStyle + ("QPushButton { border-right: 7px solid rgb(44, 49, 60); }")
		return select

	def deselectMenu(getStyle):
		deselect = getStyle.replace("QPushButton { border-right: 7px solid rgb(44, 49, 60); }", "")
		return deselect

	# Starting selections
	# TODO: Use it for starting selection in QWebEngineView?
	def selectStandardMenu(self, widget):
		for w in self.ui.frame_left_menu.findChildren(QPushButton):
			if w.objectName() == widget:
				w.setStyleSheet(GuiFunctions.selectMenu(w.styleSheet()))

	## Reset selection
	def resetStyle(self, widget):
		for w in self.ui.frame_left_menu.findChildren(QPushButton):
			if w.objectName() != widget:
				w.setStyleSheet(GuiFunctions.deselectMenu(w.styleSheet()))

	def updatePageLabel(self, text):
		# Update page label text
		newText = f" Viewing - {text.upper()}"
		self.ui.pageNameInfo.setText(newText)

	def addHomeButtons(self):
		homeButtonsObj = ["loadFastaButton", "selectOutdirButton" ,"loadResultsButton",
		"runPipelineButton","openHelpButton"		
		]
		homeButtonNames = ["Load the fasta file to analyze", "Select directory to write the output",
		"Load the results directory", "Execute the pipeline", "View help videos"
		]
		self.homeButtons = []
		indeces = range(len(homeButtonsObj))
		for idx in indeces:
			objName = homeButtonsObj[idx]
			name = homeButtonNames[idx]
			font = QFont()
			font.setFamily(u"Segoe UI")
			homeButton = QPushButton(str(1),self)
			homeButton.setObjectName(objName) 
			homeButton.setMinimumSize(QSize(0, 70))
			homeButton.setLayoutDirection(Qt.LeftToRight)
			homeButton.setFont(font)
			homeButton.setStyleSheet(Style.style_bt_standard.replace('ICON_REPLACE', "url(:/16x16/icons/16x16/cil-home.png)"))
			homeButton.setText(name)
			homeButton.setToolTip(name)
			homeButton.clicked.connect(self.Button)
			self.ui.verticalLayout_5.addWidget(homeButton)
			self.homeButtons.append(homeButton)
			if objName == "runPipelineButton":
				homeButton.setDisabled(True)
				homeButton.setStyleSheet(Style.style_bt_disabled.replace('ICON_REPLACE', "url(:/16x16/icons/16x16/cil-home.png)"))
			# For some reason the spacer does not work
			# spacer = QSpacerItem(1, 1)
			# self.ui.verticalLayout_5.addWidget(spacer)
	

	############# Home button actions
	def getFastaFile(self):
		def _enablePipeline(self):
			homeButton = self.homeButtons[3]
			homeButton.setStyleSheet(Style.style_bt_standard.replace('ICON_REPLACE', "url(:/16x16/icons/16x16/cil-home.png)"))
			homeButton.setDisabled(False)
		file_filter = "Fasta file (*.fa *.fasta);;Text files (*.txt)"
		response = QFileDialog.getOpenFileName(
			parent=self,
			caption="Select fasta file",
			directory=str(pathlib.Path.cwd()),
			filter= file_filter,
			initialFilter="Fasta file (*.fa *.fasta)"
			)
		fasta_path = pathlib.Path(response[0])
		if fasta_path != pathlib.Path("."):
			self.fasta_path = fasta_path
			self.outdir = pathlib.Path.cwd() / pathlib.Path(str(fasta_path.stem) + "_" + str(date.today()).replace("-","_"))
			self.paramsdf.loc["query"] = str(self.fasta_path)
			self.paramsdf.loc["in"] = str(self.fasta_path.parent)
			self.paramsdf.loc["out"] = str(self.outdir)
			self.analyzedSeqs = list(SeqIO.index(str(self.fasta_path),"fasta").keys())
			_enablePipeline(self)
		else:
			return
		#TODO WRITE IN MD. If no output directory is specified the results are thrown into the script directory with the Fasta input
		# filename prefix (without the extension) and the current date
	
	def getOutdir(self):
		"""
		Update the default output directory from user input
		"""
		response = QFileDialog.getExistingDirectory(
			parent=self,
			caption="Select output directory",
			)
		if response == "":
			return
		else:
			self.outdir = pathlib.Path(response)
			self.paramsdf.loc["out"] = str(self.outdir)
	
	def getResdir(self):
		response = QFileDialog.getExistingDirectory(
			parent=self,
			caption="Select the results directory",
			)
		if response == "":
			return
		else:
			self.outdir = pathlib.Path(response)
			self.paramsdf.loc["out"] = str(self.outdir)
			logfile = self.outdir / pathlib.Path("logfile.log")
			if not logfile.exists:
				raise Exception("Cannot load results. LOGFILE has been deleted")
				# TODO: Make this appear in the warning window
				# TODO: Make warning window a panel for the main page?
			# When the results are fetched from the completed dir
			# paramsdf needs to be initialized.  The logfile comes handy
			self.paramsdf = pd.DataFrame(index=[],columns=["Value"])
			fin = open(logfile,"r")
			for line in fin.readlines():
				m = re.match(r".+\tparam: (\S+),(\S+)",line)
				if m:
					index, val = m.group(1), m.group(2)
					if re.match(r"\/|\\",val):
						val = pathlib.Path(val)	
					self.paramsdf.loc[index] = val
			self.fasta_path = self.paramsdf.loc["query","Value"]
			self.analyzedSeqs = list(SeqIO.index(str(self.fasta_path),"fasta").keys())
			GuiFunctions.addNewMenu(self, "Results", "resultsButton", "url(:/16x16/icons/16x16/cil-user-follow.png)", True)

	
	def execPipeline(self, exe=True):
		"""
		By default the whole pipeline will be executed. It will create another thread to run, so the application doesn't freeze
		Step 1: Create the worker class that will execute all the necessary functions
		Step 2: Create the QThread
		Step 3: Create the worker object
		Step 4: Move the worker to the QThread
		Step 5: Connect the various signals from and to the worker
		Step 6: Start the thread
		Step 7: Clear the central widget and print the final results
		If the results are just needed to be loaded, the pipeline won't be executed but the function will be used to visualize the results 
		"""
		# TODO: If i use a dictionary for the various panels (widgets) could i display each one i want each time without destroying the previous
		self.thread = QThread() # Step 2
		self.worker = Worker() # Step 3
		self.worker.moveToThread(self.thread) # Step 4
		self.thread.started.connect(self.worker.runPipeline) # Step 5
		# Create a new window for the progress
		self.worker.finished.connect(self.thread.quit) # Step 5
		self.worker.finished.connect(self.worker.deleteLater) # Step 5
		self.thread.finished.connect(self.thread.deleteLater) # Step 5
		self.thread.start() # Step 6
		self.worker.updateWidgets("starting")
		self.thread.finished.connect(lambda: self.worker.updateWidgets("finished")) # Step 7
		self.thread.finished.connect(lambda: GuiFunctions.addNewMenu(self, "Results", "resultsButton", "url(:/16x16/icons/16x16/cil-user-follow.png)", True)) # Step 7
	
	##### Results page functions
	# def filterList(self, text):
	# 	for i in range(self.listWidget.count()):
	# 		item = self.listWidget.item(i)
	# 		item.setHidden(text not in item.text())
	def initBlastResTable():
		pass
		"""  pseudocode
		headers = df.columns
		headercount = len(df.columns)
		for i in range(headercount):
			header = headers[i]
			_blastResTableHeader = QTableWidgetItem(header)
			blastResTable.setHorizontalHeaderItem(i,_blastResTableHeader)
		"""
	def updateResultsTable():
		pass
		"""
		this going to update the tables depending on their object name
		if 
		"""

class LogWindow(QMainWindow):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.resize(300,300)
		self.setMinimumSize(300,300)
		self.setWindowTitle("Log Window")
		self.setStyleSheet(
			"""
			background-color: rgb(44, 49, 60);
			color: rgb(210, 210, 210);
			"""
		)
		self.centralWidget = QLabel("Hello world")
		self.centralWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
		self.setCentralWidget(self.centralWidget)
		
		self.pipelineWidget = QLabel("Executing pipeline... please wait")
		self.pipelineWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)

		self.finishedWidget = QLabel("Finished! \n Please don't delete the logfile\n You can safely close this window")
		self.finishedWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)

		self.show()

class Simplot_page(QMainWindow):
	####
	# NEED TO DO STYLING CHANGES
	####

			# Need it for Simplot
		# self.pushButton = QPushButton(self.frame_content_wid_1)
		# self.pushButton.setObjectName(u"pushButton")
		# self.pushButton.setMinimumSize(QSize(150, 30))
		# font8 = QFont()
		# font8.setFamily(u"Segoe UI")
		# font8.setPointSize(9)
		# self.pushButton.setFont(font8)
		# self.pushButton.setStyleSheet(u"QPushButton {\n"
		# "	border: 2px solid rgb(52, 59, 72);\n"
		# "	border-radius: 5px;	\n"
		# "	background-color: rgb(52, 59, 72);\n"
		# "}\n"
		# "QPushButton:hover {\n"
		# "	background-color: rgb(57, 65, 80);\n"
		# "	border: 2px solid rgb(61, 70, 86);\n"
		# "}\n"
		# "QPushButton:pressed {	\n"
		# "	background-color: rgb(35, 40, 49);\n"
		# "	border: 2px solid rgb(43, 50, 61);\n"
		# "}")
		# icon3 = QIcon()
		# icon3.addFile(u":/16x16/icons/16x16/cil-folder-open.png", QSize(), QIcon.Normal, QIcon.Off)
		# self.pushButton.setIcon(icon3)
	def __init__(self, sequence, window_size, step, aln_f, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Similarity plot")
		self.resize(1000,1000)
		self.setAutoFillBackground(True)
		self.setStyleSheet("background: 'white'")
		self.qseq = sequence
		self.aln_f = pathlib.Path(aln_f)

		# Simplot creation parameters
		self.wsize = window_size
		self.step = step

		# Create the pyqtgraph widget and embed it into the central Widget
		# Next create the legend widget and follow the same steps
		self.centralWidget = QWidget()
		self.centralWidget.mainLayout = QGridLayout()
		self.setCentralWidget(self.centralWidget)

		# Similarity plot graph
		self.simplotGraphWidget = pg.GraphicsLayoutWidget(parent=self)
		self.simplotGraphWidget.setWindowTitle("Similarity plot")
		
		# Initialize plot
		self.plot = self.simplotGraphWidget.addPlot()
		self.graphicsLayoutWidgetlayout = self.simplotGraphWidget.addLayout()

		# Create extra view boxes to put the legend etc
		self.viewbox = pg.ViewBox(enableMenu=False, enableMouse=False, border=pg.mkPen(color=(0,0,0),width=0.5))
		self.graphicsLayoutWidgetlayout.addItem(self.viewbox, row=0,col=1)

		self.viewbox.setMaximumWidth(200)
		self.viewbox.setMinimumWidth(200)
		self.viewbox.setBackgroundColor((191, 189, 189))
		self.viewbox_label = pg.LabelItem("<span style=\"color:black;font-size:11pt\">Legend: </span>", parent=self.viewbox)

		self.legend = pg.LegendItem(labelTextColor=(0,0,0),labelTextSize="11pt")
		self.viewbox.addItem(self.legend)
		self.legend.setParentItem(self.viewbox)
		self.legend.anchor((0,0), (0,0), offset=(0,12))

		# Prettify the widget
		self.simplotGraphWidget.setBackground("w")

		# Initialize the different pen colors
		self.pen_green = pg.mkPen(color=(40, 237, 43),  width=1.5)
		self.pen_blue = pg.mkPen(color=(66, 170, 245),  width=1.5)
		self.pen_orange = pg.mkPen(color=(237, 155, 40),  width=1.5)
		self.pen_red = pg.mkPen(color=(255, 0, 0),  width=1.5)

		# Plot the curves
		self.xaxisData, self.yaxisData, self.alnSize, self.qseqSize = simplot.calculate_similarities(self.qseq, self.aln_f,self.step,self.wsize)
		self.c1 = self.plot.plot(self.xaxisData, self.yaxisData["NC_001526.4_A1"], name="Lineage A", pen=self.pen_green)
		self.c2 = self.plot.plot(self.xaxisData, self.yaxisData["AF536180.1_B1"], name="Lineage B", pen=self.pen_blue)
		self.c3 = self.plot.plot(self.xaxisData, self.yaxisData["AF472509.1_C1"], name="Lineage C", pen=self.pen_orange)
		self.c4 = self.plot.plot(self.xaxisData, self.yaxisData["HQ644257.1_D1"], name="Lineage D", pen=self.pen_red)

		self.legend = pg.LegendItem(labelTextColor=(0,0,0),labelTextSize="11pt")
		self.viewbox.addItem(self.legend)
		self.legend.setParentItem(self.viewbox)
		self.legend.anchor((0,0), (0,0), offset=(0,12))
		self.legend.addItem(self.c1, name="Lineage A")
		self.legend.addItem(self.c2, name="Lineage B")
		self.legend.addItem(self.c3, name="Lineage C")
		self.legend.addItem(self.c4, name="Lineage D")

		# Configure the mouse functions of the plot
		# self.c1.sigClicked.connect(self.plotClicked)
		# TODO: Figure out what i need to do
		# Can leave this for now on
		# /Configure the mouse functions of the plot

		# Modify plot
		self.title = ("<span style=\"color:black;font-size:11pt;background-color:#cffaf9\">Sequence: %s Sequence size: %s Alignment size: %s Window: %s Step: %s</span>"
		%(self.qseq,self.qseqSize,self.alnSize,self.wsize,self.step)
		)

		self.plot.setTitle(self.title)
		self.plot.setLabel('left', 'Identity %', color="b")
		self.plot.setLabel('bottom', 'Alignment position', color="b")
		self.plot.axes["bottom"]["item"].setPen(pg.mkPen(color=(0,0,0)))
		self.plot.axes["bottom"]["item"].setTextPen(pg.mkPen(color=(0,0,0)))
		self.plot.axes["left"]["item"].setPen(pg.mkPen(color=(0,0,0)))
		self.plot.axes["left"]["item"].setTextPen(pg.mkPen(color=(0,0,0)))

		self.saveButton = QPushButton("Save Image")

		self.centralWidget.mainLayout.addWidget(self.simplotGraphWidget, 0, 0, 1, 6 )
		self.centralWidget.mainLayout.addWidget(self.saveButton, 1, 0, 1, 1)
		self.centralWidget.setLayout(self.centralWidget.mainLayout)
		self.setCentralWidget(self.centralWidget)

		self.saveButton.clicked.connect(self.screenCapture)

	def screenCapture(self):
		self.screenshot = screen.grabWindow( self.winId())
		file_filter = "Images (*.jpg)"
		self.response = QFileDialog.getSaveFileName(
			parent = self,
			caption = "Save similarity plot image",
			directory = f"{self.qseq} - Window size: {self.wsize} - Step: {str(self.step)}.jpg" ,
			filter = file_filter,
			initialFilter = file_filter
		)
		self.file_path = self.response[0]
		if self.file_path:
			self.screenshot.save(self.file_path, 'jpg')
		else:
			return

	def closeEvent(self,event):
		event.accept()
		simplot.rm_tmpfile(self.aln_f)

class Ui_MainWindow(QMainWindow):
	def setupUi(self, MainWindow):
		MainWindow.resize(1000, 720)
		MainWindow.setMinimumSize(QSize(1000, 720))
		MainWindow.setWindowTitle("HPV16 genotyping tool")
	   
		# Fonts        
		font1 = QFont("Segoe UI", pointSize= 12, weight= 75)
		font2 = QFont("Segoe UI", pointSize= 40)
		font3 = QFont("Segoe UI", pointSize= 14)
		font4 = QFont("Segoe UI", pointSize= 10, weight= 75)

		MainWindow.setStyleSheet(
			"""
			QMainWindow {background: transparent; }
			QToolTip {
				color: #FFFFFF;"
				background-color: rgba(27, 29, 35, 160);"
				border: 1px solid rgb(40, 40, 40);"
				border-radius: 2px;"
				}
			"""
		)
		self.centralwidget = QWidget(MainWindow)
		self.centralwidget.setStyleSheet(
				"background: transparent;\n"
				"color: rgb(210, 210, 210);"
		)

		self.horizontalLayout = QHBoxLayout(self.centralwidget)
		self.horizontalLayout.setSpacing(0)
		self.horizontalLayout.setContentsMargins(10, 10, 10, 10)
		
		self.frame_main = QFrame(self.centralwidget)
		self.horizontalLayout.addWidget(self.frame_main)

		self.verticalLayout = QVBoxLayout(self.frame_main)
		self.verticalLayout.setSpacing(0)
		self.verticalLayout.setContentsMargins(0, 0, 0, 0)
		
		self.frame_top_info = QFrame()
		self.frame_top_info.setMinimumSize(QSize(0, 45))
		self.frame_top_info.setMaximumSize(QSize(16777215, 45))
		self.frame_top_info.setStyleSheet(
			"background-color: rgb(39, 44, 54);"
		)

		self.horizontalLayout_2 = QHBoxLayout(self.frame_top_info)

		self.pageNameInfo = QLabel()
		self.pageNameInfo.setFont(font1)
		self.pageNameInfo.setStyleSheet(u"color: rgb(98, 103, 111);")

		self.verticalLayout.addWidget(self.frame_top_info)
		self.horizontalLayout_2.addWidget(self.pageNameInfo)
	
		self.frame_center = QFrame(self.frame_main)
		self.frame_center.setStyleSheet("background-color: rgb(40, 44, 52);")

		self.horizontalLayout_1 = QHBoxLayout(self.frame_center)
		self.horizontalLayout_1.setSpacing(0)
		self.horizontalLayout_1.setContentsMargins(0, 0, 0, 0)
		self.frame_left_menu = QFrame(self.frame_center)
		self.frame_left_menu.setMaximumSize(QSize(70, 16777215))

		self.frame_left_menu.setLayoutDirection(Qt.LeftToRight)
		self.frame_left_menu.setStyleSheet(u"background-color: rgb(27, 29, 35);")
		
		self.verticalLayout_3 = QVBoxLayout(self.frame_left_menu)
		self.verticalLayout_3.setSpacing(1)
		self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
		
		self.frame_menus = QFrame(self.frame_left_menu)
		
		# Will add the buttons
		self.menusLayout = QVBoxLayout(self.frame_menus)
		self.menusLayout.setSpacing(0)
		self.menusLayout.setContentsMargins(0, 0, 0, 0)

		self.verticalLayout_3.addWidget(self.frame_menus, 0, Qt.AlignTop)

		self.horizontalLayout_1.addWidget(self.frame_left_menu)

		# The right panel is where the pages are going to be displayed
		self.frameContentRightPanel = QFrame(self.frame_center)
		self.frameContentRightPanel.setStyleSheet(u"background-color: rgb(44, 49, 60);")

		self.verticalLayout_2 = QVBoxLayout(self.frameContentRightPanel)
		self.verticalLayout_2.setSpacing(0)
		self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
		# self.frame_content = QFrame(self.frameContentRightPanel)
		self.frame_content = QFrame()

		self.verticalLayout_4 = QVBoxLayout(self.frame_content)
		self.verticalLayout_4.setSpacing(0)
		self.verticalLayout_4.setContentsMargins(5, 5, 5, 5)
		self.stackedWidget = QStackedWidget(self.frame_content)
		self.stackedWidget.setStyleSheet(u"background: transparent;")
		
		# Home page widget
		self.homePage = QWidget()
		self.verticalLayout_5 = QVBoxLayout(self.homePage)
		self.label_1 = QLabel(self.homePage)

		self.label_1.setFont(font2)
		# self.label_1.setStyleSheet(u"")
		self.label_1.setAlignment(Qt.AlignCenter)

		self.verticalLayout_5.addWidget(self.label_1)

		self.label = QLabel(self.homePage)

		self.label.setFont(font3)
		self.label.setAlignment(Qt.AlignCenter)

		self.verticalLayout_5.addWidget(self.label)

		self.verticalLayout_4.addWidget(self.stackedWidget)
		self.verticalLayout_2.addWidget(self.frame_content)

		# Results page
		self.resultsPage = QWidget()
		self.resultsPage.setObjectName("resultsPage")
		self.putativeRecombinants = QCheckBox() # Will use this to "transorm" the page into a putative recombinants only
		# TODO: Implement this change when checked (just clear() and add putative recs in listWidget, and clear() the rest of the widgets)

		self.resultsPageVericalLayoutForFrames = QVBoxLayout(self.resultsPage) # This will be used to put all the frames
		self.resultsPageVericalLayoutForFrames.setObjectName("resultsPageVericalLayoutForFrames")
		
		########## QLine and list list frame
		self.listFrame = QFrame(self.resultsPage)
		self.listFrame.setObjectName("listFrame")
		self.listFrame.setStyleSheet("border-radius: 5px;")
		self.listFrame.setFrameShape(QFrame.StyledPanel)
		self.listFrame.setFrameShadow(QFrame.Raised)

		self.verticalLayout_15 = QVBoxLayout(self.listFrame)
		self.verticalLayout_15.setSpacing(0)
		self.verticalLayout_15.setObjectName(u"verticalLayout_15")
		self.verticalLayout_15.setContentsMargins(0, 0, 0, 0)
		
		self.frame_div_content_1 = QFrame(self.listFrame)
		self.frame_div_content_1.setObjectName("frame_div_content_1")
		self.frame_div_content_1.setMinimumSize(QSize(0, 110))
		self.frame_div_content_1.setMaximumSize(QSize(16777215, 110))
		self.frame_div_content_1.setStyleSheet(
			"""
			background-color: rgb(41, 45, 56);
			border-radius: 5px;
			"""
		)

		self.frame_div_content_1.setFrameShape(QFrame.NoFrame)
		self.frame_div_content_1.setFrameShadow(QFrame.Raised)
		
		self.verticalLayout_7 = QVBoxLayout(self.frame_div_content_1)
		self.verticalLayout_7.setSpacing(0)
		self.verticalLayout_7.setObjectName(u"verticalLayout_7")
		self.verticalLayout_7.setContentsMargins(0, 0, 0, 0)
		
		self.frame_title_wid_1 = QFrame(self.frame_div_content_1)
		self.frame_title_wid_1.setObjectName(u"frame_title_wid_1")
		self.frame_title_wid_1.setMaximumSize(QSize(16777215, 35))
		self.frame_title_wid_1.setStyleSheet(u"background-color: rgb(39, 44, 54);")
		self.frame_title_wid_1.setFrameShape(QFrame.StyledPanel)
		self.frame_title_wid_1.setFrameShadow(QFrame.Raised)
		
		self.verticalLayout_7.addWidget(self.frame_title_wid_1)
		
		self.frame_content_wid_1 = QFrame(self.frame_div_content_1)
		self.frame_content_wid_1.setObjectName(u"frame_content_wid_1")
		self.frame_content_wid_1.setFrameShape(QFrame.NoFrame)
		self.frame_content_wid_1.setFrameShadow(QFrame.Raised)
		
		self.horizontalLayout_9 = QHBoxLayout(self.frame_content_wid_1)
		self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
		
		self.gridLayout = QGridLayout()
		self.gridLayout.setObjectName(u"gridLayout")
		self.gridLayout.setContentsMargins(-1, -1, -1, 0)
		
		# Widgets
		self.lineEdit = QLineEdit(self.frame_content_wid_1)
		self.lineEdit.setObjectName(u"lineEdit")
		self.lineEdit.setMinimumSize(QSize(0, 30))
		self.lineEdit.setStyleSheet(
			"""
			QLineEdit {
				background-color: rgb(27, 29, 35);
				border-radius: 5px;
				border: 2px solid rgb(27, 29, 35);
				padding-left: 10px;
			}
			QLineEdit:hover {
				border: 2px solid rgb(64, 71, 88);
			}
			QLineEdit:focus {
				border: 2px solid rgb(91, 101, 124);
			}
			"""
		)
		self.listWidget = QListWidget()
		self.listWidget.resize(100,100)
		self.listWidget.setStyleSheet(
			"""
			 QScrollBar:vertical {
				border: none;
				background: rgb(52, 59, 72);
				width: 14px;
				margin: 5px 0 5px 0;
				border-radius: 0px;
			 }
			 QScrollBar::handle:vertical {	
				background: rgb(85, 170, 255);
				min-height: 25px;
				border-radius: 7px;
			 }
			 QScrollBar::add-line:vertical {
				border: none;
				background: rgb(55, 63, 77);
				height: 20px;
				border-bottom-left-radius: 7px;
				border-bottom-right-radius: 7px;
				subcontrol-position: bottom;
				subcontrol-origin: margin;
			 }
			 QScrollBar::sub-line:vertical {
				border: none;
				background: rgb(55, 63, 77);
				height: 20px;
				border-top-left-radius: 7px;
				border-top-right-radius: 7px;
				subcontrol-position: top;
				subcontrol-origin: margin;
			 }
			 QScrollBar::up-arrow:vertical, QScrollBar::down-arrow:vertical {
				background: none;
			 }
			 QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
				 background: none;
			}
			"""		
		)
		self.seqList = ["A2_refgenome","MG850175.1","MG850176.1","MG850177.1","MG850178.1","MG850179.1","MG850180.1"] # Hard coded for testing TODO: Fix after testing
		self.listWidget.addItems(self.seqList) # TODO: Implement self.seqList (in the main Window and pass it from there)
		self.listWidget.firstItem = self.listWidget.item(0)
		self.listWidget.setCurrentItem(self.listWidget.firstItem)
		self.selectedSeq = self.listWidget.currentItem().text()
		# print(f"{type(self.lineEdit.text())} or {self.lineEdit.text}")
		# self.lineEdit.textChanged.connect(lambda _:print(F"The text is {self.lineEdit.text()}"))
		# # TODO: Implement self.filterList
		self.recSeqsCheckBox = QCheckBox(self.listFrame)
		self.recSeqsCheckBox.setStyleSheet(
			"""
			QCheckBox::indicator {\n"
				border: 3px solid rgb(52, 59, 72);\n"
				width: 15px;\n"
				height: 15px;\n"
				border-radius: 10px;\n"
				background: rgb(44, 49, 60);\n"
			}
			QCheckBox::indicator:hover {
				border: 3px solid rgb(58, 66, 81);
			}
			QCheckBox::indicator:checked {
			background: 3px solid rgb(52, 59, 72);
			border: 3px solid rgb(52, 59, 72);
			background-image: url(:/16x16/icons/16x16/cil-check-alt.png);
			}
			"""
		)
		self.recSeqsCheckBox.setText("Show\nputative recombinants\nonly") # TODO: set the word position to center
		self.recSeqsCheckBox.stateChanged.connect(lambda _:print(self.recSeqsCheckBox.isChecked()))
		# TODO: Implement a function to filter the list view: clear content, then add from recseqslist

		self.horizontalLayout_9.addLayout(self.gridLayout)
		self.verticalLayout_7.addWidget(self.frame_content_wid_1)
		self.verticalLayout_15.addWidget(self.frame_div_content_1)
		self.resultsPageVericalLayoutForFrames.addWidget(self.listFrame)


		self.gridLayout.addWidget(self.lineEdit, 0, 0, 1, 1)
		self.gridLayout.addWidget(self.recSeqsCheckBox,0,1,-1,1)
		self.gridLayout.addWidget(self.listWidget, 1, 0, 2, 1)

		############# Blast results frame #############
		self.blastResFrame = QFrame(self.resultsPage)
		self.blastResFrame.setMinimumSize(QSize(0, 150))
		self.blastResFrame.setStyleSheet("""
			background-color: rgb(39, 44, 54);
			border-radius: 5px;
		""")
		self.blastResFrame.setFrameShape(QFrame.StyledPanel)
		self.blastResFrame.setFrameShadow(QFrame.Raised)

		self.blastResVerticalLayout = QVBoxLayout(self.blastResFrame)
		self.blastResVerticalLayout.setObjectName(u"blastResVerticalLayout")
		
		self.blastResGridLayout = QGridLayout()
		
		# Widgets
		# 1. Label
		self.blastLabel = QLabel(self.blastResFrame)
		self.blastLabel.setObjectName("blastLabel")
		self.blastLabel.setFont(font4)
		self.blastLabel.setText("BLAST results")
		# 2. Table
		self.blastResTable = QTableWidget(self.blastResFrame)
		# Set Rows and Columns to a first value so i can init the view
		self.blastResTable.setColumnCount(10)
		self.blastResTable.setRowCount(10)
		self.blastResTable.setStyleSheet(Style.style_table_standard)
		# 3. LED
		self._blastLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self._blastLED.setFixedSize(18,18)
		# 5. Webengine view
		self.geneIDBrowser = QWebEngineView()

		self.blastResGridLayout.addWidget(self.blastLabel, 0, 0, 1, 1)
		self.blastResGridLayout.addWidget(self._blastLED, 0, 1, 1, 1)
		self.blastResGridLayout.addWidget(self.blastResTable, 1, 0, 3, 5)
		self.blastResGridLayout.addWidget(self.geneIDBrowser, 1, 5, 3, 3)

		self.blastResVerticalLayout.addLayout(self.blastResGridLayout)
		self.resultsPageVericalLayoutForFrames.addWidget(self.blastResFrame)
		
		############# Lineage Specific SNPs Frame
		self.lineageSnpFrame = QFrame(self.resultsPage)
		self.lineageSnpFrame.setMinimumSize(QSize(0, 150))
		self.lineageSnpFrame.setStyleSheet("""
				background-color: rgb(39, 44, 54);
				border-radius: 5px;
		""")
		self.lineageSnpFrame.setFrameShape(QFrame.StyledPanel)
		self.lineageSnpFrame.setFrameShadow(QFrame.Raised)

		self.lineageSnpVerticalLayout = QVBoxLayout(self.lineageSnpFrame)
		self.lineageSnpVerticalLayout.setObjectName(u"lineageSnpVerticalLayout")
		
		self.lineageSnpGridLayout = QGridLayout()
		
		# Widgets
		# 1. Label
		self.lineageSnpLabel = QLabel(self.lineageSnpFrame)
		self.lineageSnpLabel.setFont(font4)
		self.lineageSnpLabel.setText("Lineage specific SNPs")
		# 2. Table
		self.lineageSnpTable = QTableWidget(self.lineageSnpFrame)
		# Set Rows and Columns to a first value so i can init the view
		self.lineageSnpTable.setColumnCount(10)
		self.lineageSnpTable.setRowCount(1)
		self.lineageSnpTable.setStyleSheet(Style.style_table_standard)
		# 3. LED
		self._lineageSnpLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self._lineageSnpLED.setFixedSize(18,18)
		# 5. Webengine view
		self.snpBrowser = QWebEngineView()

		self.lineageSnpGridLayout.addWidget(self.lineageSnpLabel, 0, 0, 1, 1)
		self.lineageSnpGridLayout.addWidget(self._lineageSnpLED, 0, 1, 1, 1)
		self.lineageSnpGridLayout.addWidget(self.lineageSnpTable, 1, 0, 3, 5)
		self.lineageSnpGridLayout.addWidget(self.snpBrowser, 1, 5, 3, 3)

		self.lineageSnpVerticalLayout.addLayout(self.lineageSnpGridLayout)
		self.resultsPageVericalLayoutForFrames.addWidget(self.lineageSnpFrame)

		############# Cancer SNPs Frame
		self.cancerSnpFrame = QFrame(self.resultsPage)
		self.cancerSnpFrame.setMinimumSize(QSize(0, 150))
		self.cancerSnpFrame.setStyleSheet("""
				background-color: rgb(39, 44, 54);
				border-radius: 5px;
		""")
		self.cancerSnpFrame.setFrameShape(QFrame.StyledPanel)
		self.cancerSnpFrame.setFrameShadow(QFrame.Raised)

		self.cancerSnpVerticalLayout = QVBoxLayout(self.cancerSnpFrame)
		
		self.cancerSnpGridLayout = QGridLayout()
		
		# Widgets
		# 1. Label
		self.lineageSnpLabel = QLabel(self.cancerSnpFrame)
		self.lineageSnpLabel.setFont(font4)
		self.lineageSnpLabel.setText("Increased cancer risk SNPs")
		# 2. Table
		self.cancerSnpTable = QTableWidget(self.cancerSnpFrame)
		# Set Rows and Columns to a first value so i can init the view
		self.cancerSnpTable.setColumnCount(10)
		self.cancerSnpTable.setRowCount(10)
		self.cancerSnpTable.setStyleSheet(Style.style_table_standard)

		self.cancerSnpGridLayout.addWidget(self.lineageSnpLabel, 0, 0, 1, 1)
		self.cancerSnpGridLayout.addWidget(self.cancerSnpTable, 1, 0, 3, 5)

		self.cancerSnpTable.horizontalHeader().setCascadingSectionResizes(True)
		self.cancerSnpTable.horizontalHeader().setDefaultSectionSize(150)
		# TODO: Check if i can use the default section size a fraction of
		# the total table length (table.columns)
		# when updating
		self.cancerSnpTable.horizontalHeader().setStretchLastSection(True)
		self.cancerSnpTable.verticalHeader().setStretchLastSection(True)

		self.cancerSnpVerticalLayout.addLayout(self.cancerSnpGridLayout)
		self.resultsPageVericalLayoutForFrames.addWidget(self.cancerSnpFrame)

		############# Trees Frame

		self.treeFrame = QFrame(self.resultsPage)
		self.treeFrame.setMinimumSize(QSize(0, 150))
		self.treeFrame.setStyleSheet("""
			background-color: rgb(39, 44, 54);
			border-radius: 5px;
		""")
		self.treeFrame.setFrameShape(QFrame.StyledPanel)
		self.treeFrame.setFrameShadow(QFrame.Raised)
		self.treeGridLayout = QGridLayout(self.treeFrame)
		# Widgets
		# 1. Label
		self.treeLabel = QLabel(self.treeFrame)
		self.treeLabel.setFont(font4)
		self.treeLabel.setText("Explore trees interactively")
		# 2. Table
		self.RecTreesRenderButtonGroup = QButtonGroup()
		self.RecTreesRenderButtonGroup.buttonClicked[int].connect(lambda _:print(F" Just pressed {self.RecTreesRenderButtonGroup.buttonClicked[int]} "))
		row = 1
		col = 0
		self.genes = ["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
		self.treeGridLayout.addWidget(self.treeLabel, 0, 0, 1, 1)
		for gene in self.genes:
			self.treeButton = QPushButton(gene)
			self.treeButton.setObjectName(F"treeButton{gene}")
			self.treeButton.setStyleSheet(Style.style_push_button)
			self.RecTreesRenderButtonGroup.addButton(self.treeButton)
			self.treeGridLayout.addWidget(self.treeButton, row, col, 1, 1)
			col += 1
			if col == 4:
				row += 1
				col = 0
		self.resultsPageVericalLayoutForFrames.addWidget(self.treeFrame)
			

		"""
		Implement to add colors in the interactive trees:
		Each lineage will be colored according to our schema for easier inspection
		The selected sequences will be colored a bright color that the user can change from the ETE interactive window
		"""


		# / Add the results page here
		
		self.stackedWidget.addWidget(self.homePage)
		self.stackedWidget.addWidget(self.resultsPage)

		# Necessary to display the UI correctly
		self.horizontalLayout_1.addWidget(self.frameContentRightPanel)
		self.verticalLayout.addWidget(self.frame_center)
		MainWindow.setCentralWidget(self.centralwidget)

		self.retranslateUi(MainWindow)
		# self.stackedWidget.setCurrentIndex(1)
		# QMetaObject.connectSlotsByName(MainWindow)

	def retranslateUi(self, MainWindow):
		# Maybe i can use this to initialize the QWebEngines?
		self.pageNameInfo.setText(QCoreApplication.translate("MainWindow", f" Viewing - HOME", None))
		self.label_1.setText(QCoreApplication.translate("MainWindow", "HPV Genotyping Tool ", None))
		self.label.setText(QCoreApplication.translate("MainWindow", "Please use the following menu", None))
	def retranslateResults(self, MainWindow):
		pass
		self.blastLabel.setText(QCoreApplication.translate("MainWindow", u"BLENDER INSTALLATION", None))

#       # TODO: Keeping them to use later
		# self.pushButton.setText(QCoreApplication.translate("MainWindow", u"Load results", None)) # Button text
		# Hover text to button
		# self.btn_close.setToolTip(QCoreApplication.translate("MainWindow", u"Close", None))
		# self.btn_close.setText("") # Button text




class Style():
	"""
	Class to hold all the styling i need to do
	"""
	style_bt_standard = (
	"""
	QPushButton {
		background-image: ICON_REPLACE;
		background-position: left center;
		background-repeat: no-repeat;
		border: none;
		border-left: 28px solid rgb(27, 29, 35);
		background-color: rgb(27, 29, 35);
		text-align: left;
		padding-left: 45px;
	}
	QPushButton[Active=true] {
		background-image: ICON_REPLACE;
		background-position: left center;
		background-repeat: no-repeat;
		border: none;
		border-left: 28px solid rgb(27, 29, 35);
		border-right: 5px solid rgb(44, 49, 60);
		background-color: rgb(27, 29, 35);
		text-align: left;
		padding-left: 45px;
	}
	QPushButton:hover {
		background-color: rgb(33, 37, 43);
		border-left: 28px solid rgb(33, 37, 43);
	}
	QPushButton:pressed {
		background-color: rgb(85, 170, 255);
		border-left: 28px solid rgb(85, 170, 255);
	}
	"""
	)
	style_bt_disabled = (
	"""
	QPushButton {
		background-image: ICON_REPLACE;
		background-position: left center;
		background-repeat: no-repeat;
		border: none;
		border-left: 28px solid #5B6481;
		background-color: #5B6481;
		text-align: left;
		padding-left: 45px;
	}
	QPushButton[Active=true] {
		background-image: ICON_REPLACE;
		background-position: left center;
		background-repeat: no-repeat;
		border: none;
		border-left: 28px solid #5B6481;
		border-right: 5px solid #5B6481;
		background-color: #5B6481;
		text-align: left;
		padding-left: 45px;
	}
	QPushButton:hover {
		background-color: rgb(33, 37, 43);
		border-left: 28px solid rgb(33, 37, 43);
	}
	QPushButton:pressed {
		background-color: rgb(85, 170, 255);
		border-left: 28px solid rgb(85, 170, 255);
	}
	"""
	)
	style_table_standard = ("""
			QTableWidget {
				border: none;
				border-radius: 5px
			}
			 QScrollBar:vertical {
				border: none;
				background: rgb(52, 59, 72);
				width: 14px;
				margin: 5px 0 5px 0;
				border-radius: 0px;
			 }
			 QScrollBar::handle:vertical {	
				background: rgb(85, 170, 255);
				min-height: 25px;
				border-radius: 7px;
			 }
			 QScrollBar::add-line:vertical {
				border: none;
				background: rgb(55, 63, 77);
				height: 20px;
				border-bottom-left-radius: 7px;
				border-bottom-right-radius: 7px;
				subcontrol-position: bottom;
				subcontrol-origin: margin;
			 }
			 QScrollBar::sub-line:vertical {
				border: none;
				background: rgb(55, 63, 77);
				height: 20px;
				border-top-left-radius: 7px;
				border-top-right-radius: 7px;
				subcontrol-position: top;
				subcontrol-origin: margin;
			 }
			 QScrollBar::up-arrow:vertical, QScrollBar::down-arrow:vertical {
				background: none;
			 }
			 QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
				 background: none;
			}
			# TODO: Add the horizontal bar style
			"""
	)
	style_push_button = (
		"""
		QPushButton {
			border: 2px solid rgb(52, 59, 72);
			border-radius: 5px;
			background-color: rgb(52, 59, 72);
		}
		QPushButton:hover {
			background-color: rgb(57, 65, 80);
			border: 2px solid rgb(61, 70, 86);
		}
		QPushButton:pressed {
			background-color: rgb(35, 40, 49);
			border: 2px solid rgb(43, 50, 61);
		}
		"""
	)

if __name__ == "__main__":
	app = QApplication(sys.argv)
	mainW = MainWindow()
	QFontDatabase.addApplicationFont('fonts/segoeui.ttf')
	QFontDatabase.addApplicationFont('fonts/segoeuib.ttf')
	screen = app.primaryScreen()
	sys.exit(app.exec_())

# TODO: Add all stylesheets into the Style fuction and call them from there