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
from PyQt5.QtGui import QPixmap, QFont, QFontDatabase
from PyQt5.QtCore import ( Qt, QObject, QThread, 
	pyqtSignal, QUrl, QCoreApplication, QSize
)
from PyQt5.QtWidgets import (
	QApplication, QLabel, QMainWindow,
	QMenuBar, QWidget, QMenu, QAction, 
	QStatusBar, QTableWidget, QTableWidgetItem, 
	QFileDialog, QTabWidget, QGridLayout, QVBoxLayout,
	QPushButton, QScrollArea, QHBoxLayout, QButtonGroup,
	QPlainTextEdit, QMessageBox, QListWidget, QSpinBox,
	QLineEdit, QFrame, QStackedWidget, QSpacerItem
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

		self.ui.stackedWidget.setCurrentWidget(self.ui.page_home)

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
			self.ui.stackedWidget.setCurrentWidget(self.ui.page_home)
			GuiFunctions.resetStyle(self, "homeButton")
			GuiFunctions.updatePageLabel(self, "Home")
			btnWidget.setStyleSheet(GuiFunctions.selectMenu(btnWidget.styleSheet()))

		if btnWidget.objectName() == "resultsButton":
			# self.ui.stackedWidget.setCurrentWidget(self.ui.page_widgets)
			# Will customise ui.page_widgets to results
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

class Ui_MainWindow(QMainWindow):
	def setupUi(self, MainWindow):
		MainWindow.resize(1000, 720)
		MainWindow.setMinimumSize(QSize(1000, 720))
		MainWindow.setWindowTitle("HPV16 genotyping tool")
	   
		# Fonts        
		font1 = QFont(u"Segoe UI",pointSize=12, weight=75)
		font2 = QFont(u"Segoe UI",pointSize=40)
		font3 = QFont(u"Segoe UI",pointSize=14)

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
		################### Checking this
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
		#################################

		#################### Checking this
		self.frame_center = QFrame()
		# self.frame_center = QFrame(self.frame_main)
		self.frame_center.setStyleSheet("background-color: rgb(40, 44, 52);")

		# horizontalLayout_1 will hold the menu frame
		self.horizontalLayout_1 = QHBoxLayout(self.frame_center)
		self.horizontalLayout_1.setSpacing(0)
		self.horizontalLayout_1.setContentsMargins(0, 0, 0, 0)
		self.frame_left_menu = QFrame(self.frame_center)
		self.frame_left_menu.setMaximumSize(QSize(70, 16777215))
		# TODO: Frame maximum and minimum size will fill the stupid gaps that i see
		# when deleting stuff

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
		self.page_home = QWidget()
		self.verticalLayout_5 = QVBoxLayout(self.page_home)
		self.label_1 = QLabel(self.page_home)

		self.label_1.setFont(font2)
		# self.label_1.setStyleSheet(u"")
		self.label_1.setAlignment(Qt.AlignCenter)

		self.verticalLayout_5.addWidget(self.label_1)

		self.label = QLabel(self.page_home)

		self.label.setFont(font3)
		self.label.setAlignment(Qt.AlignCenter)

		self.verticalLayout_5.addWidget(self.label)

		self.verticalLayout_4.addWidget(self.stackedWidget)
		self.verticalLayout_2.addWidget(self.frame_content)
		############################ / Checking this

		# Add the results page here
		# / Add the results page here
		
		self.stackedWidget.addWidget(self.page_home)

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


if __name__ == "__main__":
	app = QApplication(sys.argv)
	mainW = MainWindow()
	QFontDatabase.addApplicationFont('fonts/segoeui.ttf')
	QFontDatabase.addApplicationFont('fonts/segoeuib.ttf')
	screen = app.primaryScreen()
	sys.exit(app.exec_())
