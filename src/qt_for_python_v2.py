import os
import re
import sys
from typing import Sequence
import wrapper # Tool module
import simplot # Tool module
import pathlib
import pandas as pd
import pyqtgraph as pg
import plotly.express as px
from datetime import date
from QLed import QLed
from pyqtgraph import PlotWidget, plot
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import ( 
	Qt, QObject, QThread, 
	pyqtSignal, QUrl
)
from PyQt5.QtWidgets import (
	QApplication, QLabel, QMainWindow,
	QMenuBar, QWidget, QMenu, QAction, 
	QStatusBar, QTableWidget, QTableWidgetItem, 
	QFileDialog, QTabWidget, QGridLayout, QVBoxLayout,
	QPushButton, QScrollArea, QHBoxLayout, QButtonGroup,
	QPlainTextEdit, QMessageBox, QListWidget, QSpinBox
)

from PyQt5.QtWebEngineWidgets import QWebEngineView

from ete3 import ( 
	Tree, TreeStyle, TextFace, NodeStyle, CircleFace
)

#TODO: Better use QErrorMessage for errors

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
	finished = pyqtSignal()

	def runPipeline(self):
		paramsdf = mainW.paramsdf
		# TODO: Uncomment after testing
		wrapper.main(paramsdf)
		self.finished.emit()		
		return 

class Main_page(QMainWindow):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("HPV16 genotyping tool")
		self.resize(400, 200)
		self.centralWidget = QLabel("Hello, World")
		self.centralWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
		self.setCentralWidget(self.centralWidget)
		
		self.pipelineWidget = QLabel("Executing pipeline... please wait")
		self.pipelineWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)

		self.finishedWidget = QLabel("Finished!")
		self.finishedWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
		
		self.table = QTableWidget()
		
		self._createActions()
		self._createMenus()
		self._connectActions()
		self._createStatusBar()

		# Init exec parameters
		self.system = sys.platform
		if "win" in self.system:
			self.paramsdf = pd.read_excel(pathlib.Path(os.getcwd()) / pathlib.Path('../resources/params_win.xlsx'),index_col=0, engine="openpyxl")
		if "linux" in self.system:
			self.paramsdf = pd.read_excel(pathlib.Path(os.getcwd()) / pathlib.Path('../resources/params_linux.xlsx'),index_col=0, engine="openpyxl")

		self.paramsdf.loc["system","Value"] = self.system
		self.show()
		# To keep the pages that are opened as extra windows
		self.openWindows = {}

	def _createMenus(self):
		menuBar = QMenuBar(self)
		self.setMenuBar(menuBar)
		# File
		fileMenu = QMenu("&File", self)
		menuBar.addMenu(fileMenu)
		fileMenu.addAction(self.loadFasta)
		fileMenu.addAction(self.loadOutdir)
		fileMenu.addAction(self.loadResults)
		fileMenu.addAction(self.exit)
		# Execute
		executeMenu = QMenu("&Execute", self)
		menuBar.addMenu(executeMenu)
		executeMenu.addAction(self.pipelineAnalysis)
		executeMenu.addAction(self.individualAnalysis)
		# Options
		optionsMenu = QMenu("&Options", self)
		menuBar.addMenu(optionsMenu)
		optionsMenu.addAction(self.cleanTmp)
		# Help
		helpMenu = QMenu("&Help", self)
		menuBar.addMenu(helpMenu)
		helpMenu.addAction(self.documentation)
		helpMenu.addAction(self.helpVideos)

	def _createActions(self):
		# File actions
		self.loadFasta = QAction("&Load Fasta file",self)
		self.loadOutdir = QAction("&Select output directory",self)
		self.loadResults = QAction("Load &results directory",self)
		self.exit = QAction("&Exit",self)
		# Execute actions
		self.pipelineAnalysis = QAction("Run &pipeline")
		self.individualAnalysis = QAction("Run &individual analysis")
		self.pipelineAnalysis.setEnabled(False)
		self.individualAnalysis.setEnabled(False)
		# Options actions
		self.cleanTmp = QAction("Clean tmp directories", self)
		self.cleanTmp.isCheckable()
		self.cleanTmp.setCheckable(True)
		self.cleanTmp.isChecked()
		
		# Help actions
		self.documentation = QAction("Open &documentation",self)
		self.helpVideos = QAction("Open help &videos",self)

		# Adding shortcuts
		self.loadFasta.setShortcut("Ctrl+O")

		# Adding help tips
		cleanTmpTip = "Delete temporary directories after executing the pipeline"
		self.cleanTmp.setStatusTip(cleanTmpTip)
		self.cleanTmp.setToolTip(cleanTmpTip)

		helpTip = "Display the documentation file"
		self.documentation.setStatusTip(helpTip)
		self.documentation.setToolTip(helpTip)

	def _connectActions(self):
		# File connected actions		
		self.loadFasta.triggered.connect(lambda _: self.getFastaFile())
		self.loadOutdir.triggered.connect(lambda _: self.getOutdir())
		self.loadResults.triggered.connect(lambda _: self.getResdir())
		self.exit.triggered.connect(sys.exit)

		# Execute connected actions
		self.pipelineAnalysis.triggered.connect(lambda _: self.execPipeline())

		# Options connected actions
		# TODO: Need to read about this
		# Help connected actions
		# self.documentation.triggered.connect()

	def _createStatusBar(self):
		self.statusBar = QStatusBar()
		self.setStatusBar(self.statusBar)
		self.statusBar.showMessage("Please provide Fasta file to begin...", 13000)

	def getFastaFile(self):
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
			self.pipelineAnalysis.setEnabled(True)
			self.individualAnalysis.setEnabled(True)
		else:
			return
		
		#TODO WRITE IN MD. If no output directory is specified the results are thrown into the script directory with the Fasta input
		# filename and the current date

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
			caption="Select output directory",
			)
		if response == "":
			return
		else:
			self.resdir = pathlib.Path(response)
			self.paramsdf.loc["out"] = str(self.resdir)

	def execPipeline(self):
		
		# Step 2: Create a QThread object
		self.thread = QThread()
		# Step 3: Create a worker object
		self.worker = Worker()
		# Step 4: Move worker to the thread
		self.worker.moveToThread(self.thread)
		# Step 5: Connect signals and slots

		self.thread.started.connect(self.worker.runPipeline)
		self.worker.finished.connect(self.thread.quit)
		self.worker.finished.connect(self.worker.deleteLater)
		self.thread.finished.connect(self.thread.deleteLater)
		# Step 6: Start the thread
		self.thread.start()
		self.setCentralWidget(self.pipelineWidget)

		# Step 7: Print finished message and show the results
		self.thread.finished.connect(lambda: self.setCentralWidget(self.finishedWidget))
		self.thread.finished.connect(lambda: self.showResults())		
		
	def showError(self, error):
		pass
	
	# def show_results(self):
	# 	resW = Results_page()
	# 	resW.show()

	def closeEvent(self, event):
		reply = QMessageBox.question(self, 'Window Close', 'Do you want to close this window and terminate the application?',
				QMessageBox.Yes | QMessageBox.No)
		if reply == QMessageBox.Yes:
			event.accept()
		else:
			event.ignore()

	def showResults(self):
		self.page = Results_page(self)
		self.openWindows["Results"] = self.page
		self.page.show()
		
class Progress_widget(QWidget):
	"""
	Create a progress page and embed it in the main window
	Show the progress of the pipeline
	"""
	# TODO: Finish this
	pass


class Error_page(QWidget):
	"""
	Open a window and display the last line of the log file
	The last line contains the error
	"""
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Error message")
		self.resize(1000, 200)
		self.mainLayout = QGridLayout()
		self.logFile = "logfile.log"

		self.ErrorMsg = QLabel("An error has occured please check the log below ...")
		self.ErrorMsg.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
		self.FileRead = QPlainTextEdit()
		self.text= list(open(self.logFile).readlines())
		self.FileRead.setPlainText(self.text[-1])
		self.mainLayout.addWidget(self.ErrorMsg)
		self.mainLayout.addWidget(self.FileRead)
		self.setLayout(self.mainLayout)

# TODO: Need to fix in results page
"""
1. Once i click in a recombinant sequence it automatically fixes the window and cannot resize
2. Add box to search threw the list widget
"""
class Results_page(QMainWindow):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Results")

		self.centralWidget = QWidget()

		self.centralWidget.mainLayout = QGridLayout()
		self.genes = ["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
		self.gene_name_regex = re.compile(r"^\S+_(\S+)_aln.+$")
		# Accept a list of pathlib Paths for the phylogenetic trees
		
		trees_dir = "tmp/trees"
		trees = os.listdir(trees_dir)
		trees = [pathlib.Path(trees_dir) / t for t in trees]
		
		self.trees_l = trees
		self.trees = {}
		for t in self.trees_l:
			m = re.match(self.gene_name_regex, t.name)
			gene = m.group(1)
			self.trees[gene] = Tree(str(t))

		self.outdir = mainW.outdir

		# Add tabs
		self.tabs = QTabWidget()
		self.tabs.addTab(self.BlastResTab(),"BLAST results")
		self.tabs.addTab(self.SnpResTab(),"SNP results")
		self.tabs.addTab(self.TreesResTab(),"Phylogenetic trees")
		self.tabs.addTab(self.RecSeqsTab(),"Recombinant sequences")
		
		self.centralWidget.mainLayout.addWidget(self.tabs, 0, 0)
		self.centralWidget.setLayout(self.centralWidget.mainLayout)
		self.setCentralWidget(self.centralWidget)

	def BlastResTab(self):
		# BLAST results table tab
		self.BlastResDataFile = self.outdir / "Filt_mock_GeneID_Blastn_res_GeneID.xlsx"
		self.geneIdFile = self.outdir / "Graphics/GeneSim/All_sequences.html"
		
		self.BlastTab = QWidget()
		self.BlastTabLayout = QVBoxLayout()
		self.BlastResTable = QTableWidget()
		self.displayBlastResData(str(self.BlastResDataFile))
		self.geneIdPlot = QWebEngineView()
		self.geneIdPlotUrl = QUrl.fromLocalFile(str(self.geneIdFile)) # Needs complete path as string
		self.geneIdPlot.load(self.geneIdPlotUrl)
		self.BlastTabLayout.addWidget(self.BlastResTable)
		self.BlastTabLayout.addWidget(self.geneIdPlot)
		self.BlastTab.setLayout(self.BlastTabLayout)
		return self.BlastTab
	
	def SnpResTab(self):
		self.SnpResDataFile = self.outdir / "Filt_mock_Probe_Blastn_res_SNP_nucl.xlsx"
		self.SnpPlotPath = self.outdir / "Graphics/SNP_detection/All_sequences.html"
		# SNP results tab
		self.SnpTab = QWidget()
		self.SnpTabLayout = QVBoxLayout()

		self.SnpResTable = QTableWidget()
		self.displaySNPResData(str(self.SnpResDataFile))
		self.SnpPlot = QWebEngineView()
		self.SnpPlotUrl = QUrl.fromLocalFile(str(self.SnpPlotPath)) # Needs complete path as string
		self.SnpPlot.load(self.SnpPlotUrl)
		self.SnpTabLayout.addWidget(self.SnpResTable)
		self.SnpTabLayout.addWidget(self.SnpPlot)
		self.SnpTab.setLayout(self.SnpTabLayout)
		return self.SnpTab
	
	def TreesResTab(self):
		# Trees results tab
		self.TreeTab = QWidget()
		self.TreeTabLayout = QVBoxLayout()

		self.TreesScrollArea = QScrollArea(widgetResizable=True)
		self.TreesScrollArea.setStyleSheet(
			"""
			QWidget{ background-color: white } 
			QScrollBar{ background-color: none } 
			"""
		)
		self.TreesImgContentWidget = QWidget()
		self.TreesScrollArea.setWidget(self.TreesImgContentWidget)
		self.TreesScrollArea.layout = QVBoxLayout(self.TreesImgContentWidget)
		
		# Sort images correctly
		highlight_dir = self.outdir / pathlib.Path("Graphics/Trees_images")
		self.tree_files = os.listdir(highlight_dir)
		self.tree_files = [pathlib.Path(f) for f in self.tree_files]
		self.tree_files_d = {}
		for file in self.tree_files:
			m = re.match(self.gene_name_regex, file.name)
			gene = m.group(1)
			self.tree_files_d[gene] = file.name
		for gene in self.genes:
			file = self.tree_files_d[gene]	
			pixmap = QPixmap(os.path.join(highlight_dir, file))
			if not pixmap.isNull():
				label = QLabel(pixmap=pixmap)
				self.TreesScrollArea.layout.addWidget(label)
		
		self.InteractiveTrees = QWidget()
		self.TreesHbox = QHBoxLayout()
		self.InteractiveTreesLabel = QLabel("Explore trees interactively")
		self.TreesHbox.addWidget(self.InteractiveTreesLabel)
		self.TreesRenderButtonGroup = QButtonGroup()
		self.TreesRenderButtonGroup.buttonClicked[int].connect(self.eteInteractive)
		
		# Gene buttons
		for gene in self.genes:
			self.button = QPushButton(gene)
			self.TreesRenderButtonGroup.addButton(self.button)
			self.TreesHbox.addWidget(self.button)
		self.InteractiveTrees.setLayout(self.TreesHbox)
		
		self.TreeTabLayout.addWidget(self.TreesScrollArea)
		self.TreeTabLayout.addWidget(self.InteractiveTrees)
		self.TreeTab.setLayout(self.TreeTabLayout)

		return self.TreeTab

	def RecSeqsTab(self):
		
		self.RecTab = QWidget()
		self.RecTab.mainLayout = QGridLayout()
		# Create the scrollable list view
		self.RecTab.listview = QListWidget()
		self.RecTab.listview.resize(100,100)
		self.RecTab.listview.addItems(filenames)
		# Add action in double click signal
		self.RecTab.listview.itemDoubleClicked.connect(self.displaySeqBlastResData)
		self.RecTab.listview.itemDoubleClicked.connect(self.displayGeneIDGraph)
		self.RecTab.listview.itemDoubleClicked.connect(self.displaySeqSNPResData)
		self.RecTab.listview.itemDoubleClicked.connect(self.displaySNPGraph)
		self.RecTab.listview.itemDoubleClicked.connect(self.updateRecLED)
		
		# # HardCoded stuff for the testing
		self.genes = ["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
		self.gene_name_regex = re.compile(r"^\S+_(\S+)_aln.+$")
		self.trees = {}
		self.trees_l = trees
		for t in self.trees_l:
			m = re.match(self.gene_name_regex, t.name)
			gene = m.group(1)
			self.trees[gene] = Tree(str(t))
		self.recombinationStatus = {
		"MG850175.1":{"Genes":1,"SNP":1},
		"MG850176.1":{"Genes":0,"SNP":1}
		}
		self.geneIDcolor_dict = {
		"A":"#00ff00",
		"B":"#0080ff",
		"C":"#ff8000",
		"D":"#ff0000"
		}
		self.snp_color_dict = {
		"Lin_A":"#00ff00",
		"Lin_B":"#0080ff",
		"Lin_C":"#ff8000",
		"Lin_D":"#ff0000",
		"Lin_BCD":"#ebeb34",
		"Other":"#ffffff"
		}

		self.RecTabLayout = QHBoxLayout()
		
		# Some tmp variables to test
		self.RecTab.RecSeqsList = ["MG850175","MG850176","MG850177","MG850178","MG850179","MG850180","MG850181","MG850182","MG850183","MG850184","MG850185","MG850186","MG850187","MG850188","MG850189","MG850190","MG850191","MG850192","MG850193","MG850194","MG850195","MG850196","MG850197","MG850198","MG850199","MG850200","MG850201","MG850202","MG850203","MG850204","MG850205","MG850206","MG850207","MG850208","MG850209","MG850210","MG850211","MG850212","MG850213","MG850214","MG850215","MG850216","MG850217","MG850218","MG850219","MG850220","MG850221","MG850222","MG850223","MG850224","MG850225","MG850226","MG850227","MG850228","MG850229","MG850230","MG850231","MG850232","MG850233","MG850234","MG850235","MG850236","MG850237","MG850238","MG850239","MG850240","MG850241","MG850242","MG850243","MG850244","MG850245","MG850246","MG850247","MG850248","MG850249","MG850250","MG850251","MG850252","MG850253","MG850254","MG850255","MG850256","MG850257","MG850258","MG850259","MG850260","MG850261","MG850262","MG850263","MG850264","MG850265","MG850266","MG850267","MG850268","MG850269","MG850270","MG850271","MG850272","MG850273","MG850274","MG850275","MG850276","MG850277","MG850278","MG850279","MG850280","MG850281","MG850282","MG850283","MG850284","MG850285","MG850286","MG850287","MG850288","MG850289"]
		self.RecTab.RecSeqsList = ["MG850175","MG850176","MG850177"]
		
		# Create scroll area to add sequence names as buttons to show the different widgets
		# Scroll area
		# Set up the buttons
		#############################################
		
		# Label of each box
		self.RecTab.blastResLabelBox = QWidget()
		self.RecTab.blastResLabelBoxLayout = QHBoxLayout()
		self.RecTab.blastResLabel = QLabel("BLAST results")
		# Led widget for Recombination results
		self.RecTab._blastLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self.RecTab._blastLED.setFixedSize(18,18)
		self.RecTab.blastResLabelBoxLayout.addWidget(self.RecTab.blastResLabel)
		self.RecTab.blastResLabelBoxLayout.addWidget(self.RecTab._blastLED)
		self.RecTab.blastResLabelBox.setLayout(self.RecTab.blastResLabelBoxLayout)
		# TODO: Check stylsheet method

		self.RecTab.snpResLabelBox = QWidget()
		self.RecTab.snpResLabelBoxLayout = QHBoxLayout()
		self.RecTab.snpResLabel = QLabel("SNP results")
		# Led widget for Recombination results
		self.RecTab._snpLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self.RecTab._snpLED.setFixedSize(18,18)
		
		self.RecTab.snpResLabelBoxLayout.addWidget(self.RecTab.snpResLabel)
		self.RecTab.snpResLabelBoxLayout.addWidget(self.RecTab._snpLED)
		self.RecTab.snpResLabelBox.setLayout(self.RecTab.snpResLabelBoxLayout)
		# TODO: Check stylsheet method

		# Add Blast results table widget
		self.RecTab.blastResTable = QTableWidget()
		# Gene results plotly graph
		self.RecTab.geneIDBrowser = QWebEngineView()
		# Add SNP results table widget
		self.RecTab.SnpResTable = QTableWidget()
		# SNP results plotly graph
		self.RecTab.snpBrowser = QWebEngineView()

		# Widget with Simplot Button
		self.RecTab.simplotButtonAreaWidget = QWidget()
		# self.simplotButtonAreaWidget.setFixedSize(150,150)
		self.RecTab.simplotButtonAreaWidget.setMaximumSize(250,200)
		self.RecTab.simplotButtonAreaWidgetSimplotLabel = QLabel("Similarity plot parameters")
		self.RecTab.simplotButtonAreaWidgetWindowLabel = QLabel("Window size")
		self.RecTab.simplotButtonAreaWidgetStepLabel = QLabel("Step")
		self.RecTab.simplotButtonAreaWidgetLayout = QGridLayout()
		self.RecTab.simplotButton = QPushButton("Calculate Simplot")
		self.RecTab.windowSizeSpinBox = QSpinBox(maximum=10000,value=300)
		self.RecTab.stepSizeSpinBox = QSpinBox(maximum=10000, value=100)
		self.RecTab.simplotButtonAreaWidgetLayout.addWidget(self.RecTab.simplotButtonAreaWidgetSimplotLabel,0,0)
		self.RecTab.simplotButtonAreaWidgetLayout.addWidget(self.RecTab.simplotButtonAreaWidgetWindowLabel,1,1)
		self.RecTab.simplotButtonAreaWidgetLayout.addWidget(self.RecTab.simplotButtonAreaWidgetStepLabel,2,1)
		self.RecTab.simplotButtonAreaWidgetLayout.addWidget(self.RecTab.simplotButton,3,0,1,2)
		self.RecTab.simplotButtonAreaWidgetLayout.addWidget(self.RecTab.windowSizeSpinBox,1,0)
		self.RecTab.simplotButtonAreaWidgetLayout.addWidget(self.RecTab.stepSizeSpinBox,2,0)
		self.RecTab.simplotButtonAreaWidget.setLayout(self.RecTab.simplotButtonAreaWidgetLayout)

		self.RecTab.simplotButton.clicked.connect(self.createSimplot)
		#TODO: Implement createSimplot

		# Trees with highlighted sequence in QScrollArea 
		self.RecTab.RecSeqsTreesWidget = QWidget()
		self.RecTab.RecSeqsTreesWidgetLayout = QVBoxLayout()
		self.RecTab.RecSeqsTreesScollArea = QScrollArea(widgetResizable=True)
		self.RecTab.RecSeqsTreesScollArea.setStyleSheet(
			"""
			QWidget{ background-color: white } 
			QScrollBar{ background-color: none } 
			"""
		)
		self.RecTab.RecTreesImgContentWidget = QWidget()
		self.RecTab.RecSeqsTreesScollArea.setWidget(self.RecTab.RecTreesImgContentWidget)
		self.RecTab.RecSeqsTreesScollArea.layout = QVBoxLayout(self.RecTab.RecTreesImgContentWidget)

		# Sort images correctly
		# highlight_dir = self.outdir / pathlib.Path("Graphics\\Trees_images")
		highlight_dir = pathlib.Path("mock_2021_08_24/Graphics/Trees_images")
		self.tree_files = os.listdir(highlight_dir)
		self.tree_files = [pathlib.Path(f) for f in self.tree_files]
		self.tree_files_d = {}
		for file in self.tree_files:
			m = re.match(self.gene_name_regex, file.name)
			gene = m.group(1)
			self.tree_files_d[gene] = file.name
		for gene in self.genes:
			file = self.tree_files_d[gene]	
			pixmap = QPixmap(os.path.join(highlight_dir, file))
			if not pixmap.isNull():
				label = QLabel(pixmap=pixmap)
				self.RecTab.RecSeqsTreesScollArea.layout.addWidget(label)
		
		self.RecTab.RecInteractiveTrees = QWidget()
		self.RecTab.RecTreesHbox = QHBoxLayout()
		self.RecTab.RecInteractiveTreesLabel = QLabel("Explore trees interactively")
		self.RecTab.RecTreesHbox.addWidget(self.RecTab.RecInteractiveTreesLabel)
		self.RecTab.RecTreesRenderButtonGroup = QButtonGroup()
		self.RecTab.RecTreesRenderButtonGroup.buttonClicked[int].connect(self.eteInteractive)
		
		# Gene buttons
		for gene in self.genes:
			self.RecTab.button = QPushButton(gene)
			self.RecTab.RecTreesRenderButtonGroup.addButton(self.RecTab.button)
			self.RecTab.RecTreesHbox.addWidget(self.RecTab.button)
		self.RecTab.RecInteractiveTrees.setLayout(self.RecTab.RecTreesHbox)
		
		self.RecTab.RecSeqsTreesWidgetLayout.addWidget(self.RecTab.RecSeqsTreesScollArea)
		self.RecTab.RecSeqsTreesWidgetLayout.addWidget(self.RecTab.RecInteractiveTrees)
		self.RecTab.RecSeqsTreesWidget.setLayout(self.RecTab.RecSeqsTreesWidgetLayout)

		self.RecTab.mainLayout.addWidget(self.RecTab.listview,1,0) # Row, column
		self.RecTab.mainLayout.addWidget(self.RecTab.simplotButtonAreaWidget,2,0) # Row, column
		self.RecTab.mainLayout.addWidget(self.RecTab.RecSeqsTreesWidget,3,0,3,1)
		self.RecTab.mainLayout.addWidget(self.RecTab.blastResLabelBox,0,1)
		self.RecTab.mainLayout.addWidget(self.RecTab.blastResTable,1,1)
		self.RecTab.mainLayout.addWidget(self.RecTab.geneIDBrowser,2,1)
		self.RecTab.mainLayout.addWidget(self.RecTab.snpResLabelBox,3,1)
		self.RecTab.mainLayout.addWidget(self.RecTab.SnpResTable,4,1)
		self.RecTab.mainLayout.addWidget(self.RecTab.snpBrowser,5,1)
		
		self.RecTab.setLayout(self.RecTab.mainLayout)

		return self.RecTab
		
	def displayBlastResData(self, excel_path):
		df = pd.read_excel(excel_path,engine="openpyxl")
		if df.size == 0:
			return
			
		df.fillna('', inplace=True)
		self.BlastResTable.setRowCount(df.shape[0])
		self.BlastResTable.setColumnCount(df.shape[1])
		self.BlastResTable.setHorizontalHeaderLabels(df.columns)
		
		# returns pandas array object
		for row in df.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.BlastResTable.setItem(row[0], col_index, tableItem)
		
		self.BlastResTable.setColumnWidth(2, 300)

	def displaySNPResData(self, excel_path):
		df = pd.read_excel(excel_path,engine="openpyxl")
		if df.size == 0:
			return
			
		df.fillna('', inplace=True)
		self.SnpResTable.setRowCount(df.shape[0])
		self.SnpResTable.setColumnCount(df.shape[1])
		self.SnpResTable.setHorizontalHeaderLabels(df.columns)
		
		# returns pandas array object
		for row in df.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.SnpResTable.setItem(row[0], col_index, tableItem)
		
		self.SnpResTable.setColumnWidth(2, 300)

	def eteInteractive(self, id):
		"""
		Open tree interactive window with the press of a button
		for the user to manipulate
		"""
		for button in self.TreesRenderButtonGroup.buttons():
			if button is self.TreesRenderButtonGroup.button(id):
				button_id = button.text()
				# Add faces etc
				t = self.trees[button_id]
				t.ladderize(direction=1)
				ts = TreeStyle()
				ts.title.add_face(TextFace(button_id),column=1)
				ts.show_branch_support = True
				t.ladderize(direction=1)
				return t.show(tree_style=ts)

	def displaySeqBlastResData(self, excel_path):
		self.RecTab.recseq = self.RecTab.listview.currentItem().text()
		self.RecTab.blastResTable.clear()
		if dfBlast.size == 0:
			return
		dfBlast.fillna('', inplace=True)
		tmpdf = dfBlast[dfBlast["qaccver"] == self.RecTab.recseq]
		self.RecTab.blastResTable.setRowCount(tmpdf.shape[0])
		self.RecTab.blastResTable.setColumnCount(tmpdf.shape[1])
		self.RecTab.blastResTable.setHorizontalHeaderLabels(tmpdf.columns)
		
		# returns pandas array object
		row_num = 0
		for row in tmpdf.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.RecTab.blastResTable.setItem(row_num, col_index, tableItem)
			row_num += 1
		
		self.RecTab.blastResTable.setColumnWidth(2, 300)

	def displayGeneIDGraph(self):
		self.RecTab.geneIDBrowser = QWebEngineView()
		tmpdf = recgID_for_graphdf[recgID_for_graphdf["Sequences"] == self.RecTab.recseq]
		tmpdf = tmpdf.sort_values("qstart")
		fig = px.scatter(tmpdf, x="Gene", y="Sequences", size="pident", color="Lineage",text="Database",
		hover_data=["Database","qstart","qend","pident"], title="Gene identification recombinant sequences", 
		color_discrete_map = self.geneIDcolor_dict
		)
		self.RecTab.geneIDBrowser.setHtml(fig.to_html(include_plotlyjs='cdn'))
		self.RecTab.geneIDBrowser.setFixedHeight(300)
		self.RecTab.mainLayout.addWidget(self.RecTab.geneIDBrowser,2,1)
		self.RecTab.setLayout(self.RecTab.mainLayout)

	def displaySeqSNPResData(self, excel_path):
		self.RecTab.SnpResTable.clear()
		if dfSNP.size == 0:
			return
		dfSNP.fillna('', inplace=True)
		tmpdf = dfSNP[dfSNP["Index"] == self.RecTab.recseq]
		self.RecTab.SnpResTable.setRowCount(tmpdf.shape[0])
		self.RecTab.SnpResTable.setColumnCount(tmpdf.shape[1])
		self.RecTab.SnpResTable.setHorizontalHeaderLabels(tmpdf.columns)
		
		# returns pandas array object
		row_num = 0
		for row in tmpdf.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.RecTab.SnpResTable.setItem(row_num, col_index, tableItem)
			row_num += 1
		
		self.RecTab.SnpResTable.setColumnWidth(2, 300)
		self.RecTab.SnpResTable.setFixedHeight(80)

	def displaySNPGraph(self):
		self.RecTab.snpBrowser = QWebEngineView()
		tmpdf = recSNP_for_graphdf[recSNP_for_graphdf["Sequences"] == self.RecTab.recseq]
		fig = px.scatter(tmpdf, x="SNP reference pos", y="Sequences", color="Lineage",
			hover_data=["SNP reference pos","Lineage","Nucleotide"], title="SNP identification recombinant sequences", 
			color_discrete_map=self.snp_color_dict
		)
		fig.update_layout(xaxis_type = 'linear')
		self.RecTab.snpBrowser.setHtml(fig.to_html(include_plotlyjs='cdn'))
		self.RecTab.snpBrowser.setFixedHeight(400)
		self.RecTab.mainLayout.addWidget(self.RecTab.snpBrowser,5,1)
		self.RecTab.setLayout(self.RecTab.mainLayout)

	def updateRecLED(self):
		self.RecTab._blastLED.setValue(self.recombinationStatus[self.RecTab.recseq]["Genes"])
		self.RecTab._snpLED.setValue(self.recombinationStatus[self.RecTab.recseq]["SNP"])

	def createSimplot(self):
		test = "A2_refgenome"
		# simplotW = Simplot_page(self.RecTab.recseq)
		self.simplotW = Simplot_page(test, 300, 20) # TODO: Need to make these variable
		print(self.simplotW)
		mainW.openWindows["Simplot"] = self.simplotW
		self.simplotW.show()
		print(mainW.openWindows)

class Simplot_page(QMainWindow):
	def __init__(self, sequence, window_size, step, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Similarity plot")
		self.resize(1000,1000)
		self.setAutoFillBackground(True)
		self.setStyleSheet("background: 'white'")
		# Hard coded for testing
		self.qseq = sequence
		print(self.qseq)
		self.aln_f = pathlib.Path("tmp/A2_refgenome_aln.fa")
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
		self.viewbox_label = pg.LabelItem("<span style=\"color:black;font-size:10pt\">Legend: </span>", parent=self.viewbox)


		self.legend = pg.LegendItem(labelTextColor=(0,0,0),labelTextSize="10pt")
		self.viewbox.addItem(self.legend)
		self.legend.setParentItem(self.viewbox)
		self.legend.anchor((0,0), (0,0), offset=(0,12))

		# Prettify the widget
		self.simplotGraphWidget.setBackground("w")

		# Initialize the different pen colors
		self.pen_green = pg.mkPen(color=(40, 237, 43),  width=1)
		self.pen_blue = pg.mkPen(color=(66, 170, 245),  width=1)
		self.pen_orange = pg.mkPen(color=(237, 155, 40),  width=1)
		self.pen_red = pg.mkPen(color=(255, 0, 0),  width=1)

		# Plot the curves
		self.xaxisData, self.yaxisData, self.alnSize, self.qseqSize = simplot.calculate_similarities(self.qseq, self.aln_f,self.step,self.wsize)
		self.c1 = self.plot.plot(self.xaxisData, self.yaxisData["NC_001526.4_A1"], name="Lineage A", pen=self.pen_green)
		self.c2 = self.plot.plot(self.xaxisData, self.yaxisData["AF536180.1_B1"], name="Lineage B", pen=self.pen_blue)
		self.c3 = self.plot.plot(self.xaxisData, self.yaxisData["AF472509.1_C1"], name="Lineage C", pen=self.pen_orange)
		self.c4 = self.plot.plot(self.xaxisData, self.yaxisData["HQ644257.1_D1"], name="Lineage D", pen=self.pen_red)

		self.legend = pg.LegendItem(labelTextColor=(0,0,0),labelTextSize="10pt")
		self.viewbox.addItem(self.legend)
		self.legend.setParentItem(self.viewbox)
		self.legend.anchor((0,0), (0,0), offset=(0,12))
		self.legend.addItem(self.c1, name="Lineage A")
		self.legend.addItem(self.c2, name="Lineage B")
		self.legend.addItem(self.c3, name="Lineage C")
		self.legend.addItem(self.c4, name="Lineage D")

		# # Configure the mouse functions of the plot
		# self.c1.sigClicked.connect(self.plotClicked)
		# # TODO: Figure out what i need to do
		# Modify plot
		
		self.title = ("<span style=\"color:black;font-size:10pt;background-color:#cffaf9\">Sequence: %s Sequence size: %s Alignment size: %s Window: %s Step: %s</span>" %(self.qseq,self.qseqSize,self.alnSize,self.wsize,self.step))
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
			directory = f"{self.qseq} - Window size: {self.wsize} - Step: {str(self.step)}.jpg " ,
			filter = file_filter,
			initialFilter = file_filter
		)
		self.file_path = self.response[0]
		if self.file_path:
			self.screenshot.save(self.file_path, 'jpg')
		else:
			return

	def closeEvent(self, event):
		reply = QMessageBox.question(self, 'Window Close', 'Are you sure you want to close the window?',
				QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
		if reply == QMessageBox.Yes:
			event.accept()
			# simplot.rm_tmpfile(self.aln_f)
			# TODO: Uncomment after testing
		else:
			event.ignore()


if __name__ == "__main__":
	app = QApplication(sys.argv)
	# There is a problem with the thread in the main page
	mainW = Main_page()
	# TODO: Main Window should be able to be closed
	screen = app.primaryScreen()

	# For testing
	filenames = ["MG850175.1","MG850176.1","MG850177.1","MG850178","MG850179","MG850180","MG850181","MG850182","MG850183","MG850184","MG850185","MG850186","MG850187","MG850188","MG850189","MG850190","MG850191","MG850192","MG850193","MG850194","MG850195","MG850196","MG850197","MG850198","MG850199","MG850200","MG850201","MG850202","MG850203","MG850204","MG850205","MG850206","MG850207","MG850208","MG850209","MG850210","MG850211","MG850212","MG850213","MG850214","MG850215","MG850216","MG850217","MG850218","MG850219","MG850220","MG850221","MG850222","MG850223","MG850224","MG850225","MG850226","MG850227","MG850228","MG850229","MG850230","MG850231","MG850232","MG850233","MG850234","MG850235","MG850236","MG850237","MG850238","MG850239","MG850240","MG850241","MG850242","MG850243","MG850244","MG850245","MG850246","MG850247","MG850248","MG850249","MG850250","MG850251","MG850252","MG850253","MG850254","MG850255","MG850256","MG850257","MG850258","MG850259","MG850260","MG850261","MG850262","MG850263","MG850264","MG850265","MG850266","MG850267","MG850268","MG850269","MG850270","MG850271","MG850272","MG850273","MG850274","MG850275","MG850276","MG850277","MG850278","MG850279","MG850280","MG850281","MG850282","MG850283","MG850284","MG850285","MG850286","MG850287","MG850288","MG850289"]
	dfBlast = pd.read_excel("tmp/TestExcel.xlsx", engine="openpyxl")
	dfSNP = pd.read_excel("tmp/TestExcelSNP.xlsx", engine="openpyxl")
	recgID_for_graphdf = pd.read_excel("tmp/RecOnlyGeneID_to_usedf.xlsx", engine="openpyxl", index_col=0)
	recSNP_for_graphdf = pd.read_excel("tmp/RecOnlySNP_to_usedf.xlsx", engine="openpyxl", index_col=0)
	
	trees_dir = "tmp/trees"
	trees = os.listdir(trees_dir)
	trees = [pathlib.Path(trees_dir) / t for t in trees]
	sys.exit(app.exec_())
