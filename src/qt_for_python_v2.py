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
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import ( 
	QSortFilterProxyModel, Qt, QObject, QThread, 
	pyqtSignal, QUrl
)
from PyQt5.QtWidgets import (
	QApplication, QLabel, QMainWindow,
	QMenuBar, QWidget, QMenu, QAction, 
	QStatusBar, QTableWidget, QTableWidgetItem, 
	QFileDialog, QTabWidget, QGridLayout, QVBoxLayout,
	QPushButton, QScrollArea, QHBoxLayout, QButtonGroup,
	QPlainTextEdit, QMessageBox, QListWidget, QSpinBox,
	QLineEdit
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
		self.resize(1000, 720)

		
		self.centralWidget = QLabel("Hello, World")
		self.centralWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
		self.setCentralWidget(self.centralWidget)
		
		self.pipelineWidget = QLabel("Executing pipeline... please wait")
		self.pipelineWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)

		self.finishedWidget = QLabel("Finished! \n Please don't delete the logfile")
		self.finishedWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
		
		self._createActions()
		self._createMenus()
		self._connectActions()
		self._createStatusBar()

		# Init exec parameters
		# if "win" in self.system:
			# self.paramsdf = pd.read_excel(pathlib.Path(os.getcwd()) / pathlib.Path('../resources/params_win.xlsx'),index_col=0, engine="openpyxl")
		# if "linux" in self.system:
			# self.paramsdf = pd.read_excel(pathlib.Path(os.getcwd()) / pathlib.Path('../resources/params_linux.xlsx'),index_col=0, engine="openpyxl")
	
		params = wrapper._defaultparams()
		# TODO: Possible error in windows from the initialization of the default parameters (forward slash)
		self.paramsdf = pd.DataFrame.from_dict(params,orient='index')
		self.paramsdf.rename(columns={0:"Value"},inplace=True)

		self.show()
		# To keep the pages that are opened as extra windows
		self.openWindows = {}


	def _createMenus(self):
		menuBar = QMenuBar(self)
		self.setMenuBar(menuBar)
		# File menu
		fileMenu = QMenu("&File", self)
		menuBar.addMenu(fileMenu)
		fileMenu.addAction(self.loadFasta)
		fileMenu.addAction(self.loadOutdir)
		fileMenu.addSeparator()
		fileMenu.addAction(self.loadResultsAction)
		fileMenu.addSeparator()
		fileMenu.addAction(self.exit)
		# Execute menu
		executeMenu = QMenu("&Execute", self)
		menuBar.addMenu(executeMenu)
		executeMenu.addAction(self.pipelineAnalysis)
		executeMenu.addAction(self.individualAnalysis)
		# Options menu
		optionsMenu = QMenu("&Options", self)
		menuBar.addMenu(optionsMenu)
		optionsMenu.addAction(self.cleanTmp)
		# Help menu
		helpMenu = QMenu("&Help", self)
		menuBar.addMenu(helpMenu)
		helpMenu.addAction(self.documentation)
		helpMenu.addAction(self.helpVideos)

	def _createActions(self):
		# File actions
		self.loadFasta = QAction("&Load Fasta file",self)
		self.loadOutdir = QAction("&Select output directory",self)
		self.loadResultsAction = QAction("Load &results directory",self)
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
		self.pipelineAnalysis.setShortcut("Ctrl+E")

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
		self.loadResultsAction.triggered.connect(lambda _: self.getResdir())
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
			self.analyzedSeqs = list(SeqIO.index(str(self.fasta_path),"fasta").keys())
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
			self.showResults()

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
		self.worker.finished.connect(self.thread.quit) # Step 5
		self.worker.finished.connect(self.worker.deleteLater) # Step 5
		self.thread.finished.connect(self.thread.deleteLater) # Step 5
		self.thread.start() # Step 6
		self.setCentralWidget(self.pipelineWidget) # Step 7
		self.thread.finished.connect(lambda: self.setCentralWidget(self.finishedWidget)) # Step 7
		self.thread.finished.connect(lambda: self.showResults()) # Step 7
		# TODO: Uncomment after testing
		
	def showError(self, error):
		pass
	
	def showResults(self):
		if "Results" not in self.openWindows:
			self.page = Results_page(self)
			self.openWindows["Results"] = self.page
		self.openWindows["Results"].show()
		# This memoization actually yields a bug
		# When i close the results and open another folder i get the same page
		# If i add del or something i always get 
		#QQuickWidget: Failed to make context current
		#QQuickWidget::resizeEvent() no OpenGL context
		# The error is produced when i initialize the Results QMainWindow
		# after closing it and not terminating the app
		# Maybe i can use the OpenGL widget?
	
	def closeEvent(self, event):
		self.quitMsg = QMessageBox()
		reply = self.quitMsg.question(self, "Window Close", 
		"Do you want to close this window?\nThe application will be terminated",
		QMessageBox.Yes | QMessageBox.No
		)

		if reply == QMessageBox.Yes:
			event.accept()
			logging.shutdown()
		else:
			event.ignore()
		
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

class Results_page(QMainWindow):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Results")
		self.centralWidget = QWidget()
		self.centralWidget.mainLayout = QGridLayout()

		# Initialize parameters that are needed
		self.outdir = mainW.outdir
		self.fasta_f = mainW.fasta_path
		self.system = mainW.paramsdf.loc["system","Value"]
	
		self.genes = ["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
		self.gene_name_regex = re.compile(r"^\S+_(\S+)_aln.+$")

		self.trees_dir = self.outdir/"Phylogenetic_Trees"
		trees = os.listdir(self.trees_dir)
		trees = [pathlib.Path(self.trees_dir) / t for t in trees]
		
		self.trees_l = trees
		self.trees = {}
		for t in self.trees_l:
			m = re.match(self.gene_name_regex, t.name)
			gene = m.group(1)
			self.trees[gene] = Tree(str(t))
		
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

		# TODO: Filt_mock is the filtered query_file name, should make this variable
		self.dfBlast = pd.read_excel(self.outdir / "Filt_mock_GeneID_Blastn_res_GeneID.xlsx", engine="openpyxl")
		self.dfSNP = pd.read_excel(self.outdir / "SNPs_results.xlsx", engine="openpyxl")
		self.dfCancerSNP = pd.read_excel(self.outdir / "Filt_mock_cProbe_Blastn_res_cSNP_nucl.xlsx", engine="openpyxl")
		self.recgID_for_graphdf = pd.read_excel(self.outdir / "RecOnlyGeneID_to_usedf.xlsx", engine="openpyxl", index_col=0)
		self.recSNP_for_graphdf = pd.read_excel(self.outdir / "RecOnlySNP_to_usedf.xlsx", engine="openpyxl", index_col=0)
		self.recSeqsList = ["A2_refgenome","MG850175.1"] #TODO: Make this variable
		
		self.seqList = mainW.analyzedSeqs
		self.totalSequenceRecDict = SeqIO.index(str(self.fasta_f),"fasta")
		self.binaries = wrapper._init_binaries(self.system)
		
		self.muscle_bin= self.binaries[2]
		self.simplotDBFasta = mainW.paramsdf.loc["SimplotRef_database","Value"]
		
		# Add tabs
		self.tabs = QTabWidget()
		self.tabList = []
		self.tabs.addTab(self.AnalyzedSeqsTab(),"Analyzed sequences")
		self.tabs.addTab(self.RecSeqsTab(),"Potential recombinant sequences")
		
		# Add the tabs to a list so i can easily call the selected sequence
		# of each one
		self.tabList.append(self.SeqsTab)
		self.tabList.append(self.RecTab)
		self.tabIdx = self.tabs.currentIndex()
		self.tabs.currentChanged.connect(self.changeTabIdx)

		self.centralWidget.mainLayout.addWidget(self.tabs, 0, 0)
		self.centralWidget.setLayout(self.centralWidget.mainLayout)
		self.setCentralWidget(self.centralWidget)

	# TODO: I think much of the tabs code is redundant
	# Maybe i can find a way to generalise it? with the tabIdx function
	
	def AnalyzedSeqsTab(self):
		self.SeqsTab = QWidget()
		self.SeqsTab.mainLayout = QGridLayout()

		# Create the scrollable list and a search field
		self.SeqsTab.seqListWidget = QWidget()
		self.SeqsTab.seqListWidget.mainLayout = QVBoxLayout()

		self.SeqsTab.listWidget = QListWidget()
		self.SeqsTab.listWidget.resize(100,100)
		self.SeqsTab.listWidget.addItems(self.seqList)
	
		self.SeqsTab.listWidget.firstItem = self.SeqsTab.listWidget.item(0)
		self.SeqsTab.listWidget.setCurrentItem(self.SeqsTab.listWidget.firstItem)
		self.SeqsTab.selectedSeq = self.SeqsTab.listWidget.currentItem().text()

		# The search field
		self.SeqsTab.listSearchField = QLineEdit()
		self.SeqsTab.listSearchField.setStyleSheet("font-size: 11")
		self.SeqsTab.listSearchField.textChanged.connect(self.filterList)
		
		# Update the view when double clicking the list items
		self.SeqsTab.listWidget.itemDoubleClicked.connect(self.changeTabSelectedSeq)
		self.SeqsTab.listWidget.itemDoubleClicked.connect(self.displaySeqBlastResData)
		self.SeqsTab.listWidget.itemDoubleClicked.connect(self.displayGeneIDGraph)
		self.SeqsTab.listWidget.itemDoubleClicked.connect(self.displaySeqSNPResData)
		self.SeqsTab.listWidget.itemDoubleClicked.connect(self.displaySeqCancerSNPResData)
		self.SeqsTab.listWidget.itemDoubleClicked.connect(self.displaySNPGraph)
		
		self.SeqsTab.seqListWidget.mainLayout.addWidget(self.SeqsTab.listSearchField)
		self.SeqsTab.seqListWidget.mainLayout.addWidget(self.SeqsTab.listWidget)
		self.SeqsTab.seqListWidget.setLayout(self.SeqsTab.seqListWidget.mainLayout)
		
		# Blast results (Gene identification)
		self.SeqsTab.blastResWidget = QWidget()
		self.SeqsTab.blastResWidgetLayout = QGridLayout()
		self.SeqsTab.blastResLabel = QLabel("BLAST results")

		# Blast results table widget
		self.SeqsTab.blastResTable = QTableWidget()

		self.SeqsTab.blastResWidgetLayout.addWidget(self.SeqsTab.blastResLabel,0,0,1,1)
		self.SeqsTab.blastResWidgetLayout.addWidget(self.SeqsTab.blastResTable,1,0,2,2)
		self.SeqsTab.blastResWidget.setLayout(self.SeqsTab.blastResWidgetLayout)

		# SNP results
		self.SeqsTab.snpResWidget = QWidget()
		self.SeqsTab.snpResWidgetLayout = QGridLayout()
		self.SeqsTab.snpResLabel = QLabel("Lineage specific SNPs")
		
		# Add SNP results table widget
		self.SeqsTab.SnpResTable = QTableWidget()		
		self.SeqsTab.snpResWidgetLayout.addWidget(self.SeqsTab.snpResLabel,0,0,1,1)
		self.SeqsTab.snpResWidgetLayout.addWidget(self.SeqsTab.SnpResTable,1,0,2,2)
		self.SeqsTab.snpResWidget.setLayout(self.SeqsTab.snpResWidgetLayout)

		# Cancer SNP results
		self.SeqsTab.cancersnpResWidget = QWidget()
		self.SeqsTab.cancersnpResWidgetLayout = QGridLayout()
		self.SeqsTab.cancersnpResLabel = QLabel("Lineage specific SNPs")
		
		# Add cancer SNP results table widget
		self.SeqsTab.cancerSnpResTable = QTableWidget()		
		self.SeqsTab.cancersnpResWidgetLayout.addWidget(self.SeqsTab.cancersnpResLabel,0,0,1,1)
		self.SeqsTab.cancersnpResWidgetLayout.addWidget(self.SeqsTab.cancerSnpResTable,1,0,2,2)
		self.SeqsTab.cancersnpResWidget.setLayout(self.SeqsTab.cancersnpResWidgetLayout)

		# Gene results plotly graph
		self.SeqsTab.geneIDBrowser = QWebEngineView()
		# SNP results plotly graph
		self.SeqsTab.snpBrowser = QWebEngineView()
		
		# cancer SNP results table widget
		self.SeqsTab.cancerSnpResTable = QTableWidget()

		# Trees with highlighted sequence in QScrollArea 
		self.SeqsTab.SeqsTreesWidget = QWidget()
		self.SeqsTab.SeqsTreesWidgetLayout = QVBoxLayout()
		self.SeqsTab.SeqsTreesScollArea = QScrollArea(widgetResizable=True)
		self.SeqsTab.SeqsTreesScollArea.setStyleSheet(
			"""
			QWidget{ background-color: white } 
			QScrollBar{ background-color: none } 
			"""
		)
		self.SeqsTab.RecTreesImgContentWidget = QWidget()
		self.SeqsTab.SeqsTreesScollArea.setWidget(self.SeqsTab.RecTreesImgContentWidget)
		self.SeqsTab.SeqsTreesScollArea.layout = QVBoxLayout(self.SeqsTab.RecTreesImgContentWidget)

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
				self.SeqsTab.SeqsTreesScollArea.layout.addWidget(label)
		
		self.SeqsTab.RecInteractiveTrees = QWidget()
		self.SeqsTab.RecTreesHbox = QHBoxLayout()
		self.SeqsTab.RecInteractiveTreesLabel = QLabel("Explore trees interactively")
		self.SeqsTab.RecTreesHbox.addWidget(self.SeqsTab.RecInteractiveTreesLabel)
		self.SeqsTab.RecTreesRenderButtonGroup = QButtonGroup()
		self.SeqsTab.RecTreesRenderButtonGroup.buttonClicked[int].connect(self.eteInteractive)
		
		# Gene buttons
		for gene in self.genes:
			self.SeqsTab.button = QPushButton(gene)
			self.SeqsTab.RecTreesRenderButtonGroup.addButton(self.SeqsTab.button)
			self.SeqsTab.RecTreesHbox.addWidget(self.SeqsTab.button)
		self.SeqsTab.RecInteractiveTrees.setLayout(self.SeqsTab.RecTreesHbox)
		
		self.SeqsTab.SeqsTreesWidgetLayout.addWidget(self.SeqsTab.SeqsTreesScollArea)
		self.SeqsTab.SeqsTreesWidgetLayout.addWidget(self.SeqsTab.RecInteractiveTrees)
		self.SeqsTab.SeqsTreesWidget.setLayout(self.SeqsTab.SeqsTreesWidgetLayout)

		# QGridLayout.addItem(item, row, column[, rowSpan=1[, columnSpan=1[, alignment=Qt.Alignment()]]])
		self.SeqsTab.mainLayout.addWidget(self.SeqsTab.seqListWidget,0,0) # span 2,1
		self.SeqsTab.mainLayout.addWidget(self.SeqsTab.blastResWidget,0,1) # span 2,2
		self.SeqsTab.mainLayout.addWidget(self.SeqsTab.cancerSnpResTable,1,0)

		self.SeqsTab.mainLayout.addWidget(self.SeqsTab.geneIDBrowser,2,1)
		

		self.SeqsTab.mainLayout.addWidget(self.SeqsTab.snpResWidget,3,1,1,-1)
		self.SeqsTab.mainLayout.addWidget(self.SeqsTab.snpBrowser,5,1,1,-1)

		self.SeqsTab.mainLayout.addWidget(self.SeqsTab.SeqsTreesWidget,3,0,3,1) # span 3, 1

		self.SeqsTab.setLayout(self.SeqsTab.mainLayout)
		
		return self.SeqsTab

	def RecSeqsTab(self):
		
		self.RecTab = QWidget()
		self.RecTab.mainLayout = QGridLayout()

		# Create the scrollable list and a search field
		self.RecTab.seqListWidget = QWidget()
		self.RecTab.seqListWidget.mainLayout = QVBoxLayout()

		self.RecTab.listWidget = QListWidget()
		self.RecTab.listWidget.resize(100,100)
		self.RecTab.listWidget.addItems(self.recSeqsList)
	
		self.RecTab.listWidget.firstItem = self.RecTab.listWidget.item(0)
		self.RecTab.listWidget.setCurrentItem(self.RecTab.listWidget.firstItem)
		self.RecTab.selectedSeq = self.RecTab.listWidget.currentItem().text()

		# The search field
		self.RecTab.listSearchField = QLineEdit()
		self.RecTab.listSearchField.setStyleSheet("font-size: 11")
		self.RecTab.listSearchField.textChanged.connect(self.filterList)
		
		# Update the view when double clicking the list items
		self.RecTab.listWidget.itemDoubleClicked.connect(self.changeTabSelectedSeq)
		self.RecTab.listWidget.itemDoubleClicked.connect(self.displaySeqBlastResData)
		self.RecTab.listWidget.itemDoubleClicked.connect(self.displayGeneIDGraph)
		self.RecTab.listWidget.itemDoubleClicked.connect(self.displaySeqSNPResData)
		self.RecTab.listWidget.itemDoubleClicked.connect(self.displaySeqCancerSNPResData)
		self.RecTab.listWidget.itemDoubleClicked.connect(self.displaySNPGraph)
		# self.RecTab.listWidget.itemDoubleClicked.connect(self.updateRecLED)
		# TODO: Uncomment after testing
		# Need to implement recombination status
		
		self.RecTab.seqListWidget.mainLayout.addWidget(self.RecTab.listSearchField)
		self.RecTab.seqListWidget.mainLayout.addWidget(self.RecTab.listWidget)
		self.RecTab.seqListWidget.setLayout(self.RecTab.seqListWidget.mainLayout)
		
		# Blast results (Gene identification)
		self.RecTab.blastResWidget = QWidget()
		self.RecTab.blastResWidgetLayout = QGridLayout()
		self.RecTab.blastResLabel = QLabel("BLAST results")

		# Led widget for Recombination results
		self.RecTab._blastLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self.RecTab._blastLED.setFixedSize(18,18)
		
		# Blast results table widget
		self.RecTab.blastResTable = QTableWidget()
		self.RecTab.blastResWidgetLayout.addWidget(self.RecTab.blastResLabel,0,0,1,1)
		self.RecTab.blastResWidgetLayout.addWidget(self.RecTab._blastLED,0,1,1,1)
		self.RecTab.blastResWidgetLayout.addWidget(self.RecTab.blastResTable,1,0,2,2)
		self.RecTab.blastResWidget.setLayout(self.RecTab.blastResWidgetLayout)

		# SNP results
		self.RecTab.snpResWidget = QWidget()
		self.RecTab.snpResWidgetLayout = QGridLayout()
		self.RecTab.snpResLabel = QLabel("Lineage specific SNPs")
		# Led widget for Recombination results
		self.RecTab._snpLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self.RecTab._snpLED.setFixedSize(18,18)

		# Add SNP results table widget
		self.RecTab.SnpResTable = QTableWidget()		
		self.RecTab.snpResWidgetLayout.addWidget(self.RecTab.snpResLabel,0,0,1,1)
		self.RecTab.snpResWidgetLayout.addWidget(self.RecTab._snpLED,0,1,1,1)
		self.RecTab.snpResWidgetLayout.addWidget(self.RecTab.SnpResTable,1,0,2,2)
		self.RecTab.snpResWidget.setLayout(self.RecTab.snpResWidgetLayout)
	
		# Cancer SNP results
		self.RecTab.cancersnpResWidget = QWidget()
		self.RecTab.cancersnpResWidgetLayout = QGridLayout()
		self.RecTab.cancersnpResLabel = QLabel("Lineage specific SNPs")
		
		# Add cancer SNP results table widget
		self.RecTab.cancerSnpResTable = QTableWidget()		
		self.RecTab.cancersnpResWidgetLayout.addWidget(self.RecTab.cancersnpResLabel,0,0,1,1)
		self.RecTab.cancersnpResWidgetLayout.addWidget(self.RecTab.cancerSnpResTable,1,0,2,2)
		self.RecTab.cancersnpResWidget.setLayout(self.RecTab.cancersnpResWidgetLayout)

		# Gene results plotly graph
		self.RecTab.geneIDBrowser = QWebEngineView()
		# SNP results plotly graph
		self.RecTab.snpBrowser = QWebEngineView()


		# Widget with Simplot Button
		self.RecTab.simplotButtonAreaWidget = QWidget()
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

		# QGridLayout.addItem(item, row, column[, rowSpan=1[, columnSpan=1[, alignment=Qt.Alignment()]]])
		self.RecTab.mainLayout.addWidget(self.RecTab.seqListWidget,0,0) # span 2,1
		self.RecTab.mainLayout.addWidget(self.RecTab.blastResWidget,0,1) # span 2,2
		self.RecTab.mainLayout.addWidget(self.RecTab.simplotButtonAreaWidget,0,2) # span 2,1
		self.RecTab.mainLayout.addWidget(self.RecTab.cancerSnpResTable,1,0)

		self.RecTab.mainLayout.addWidget(self.RecTab.geneIDBrowser,2,1)

		self.RecTab.mainLayout.addWidget(self.RecTab.snpResWidget,3,1,1,-1)
		self.RecTab.mainLayout.addWidget(self.RecTab.snpBrowser,5,1,1,-1)

		self.RecTab.mainLayout.addWidget(self.RecTab.RecSeqsTreesWidget,3,0,3,1) # span 3, 1

		self.RecTab.setLayout(self.RecTab.mainLayout)
		
		return self.RecTab
	
	def changeTabIdx(self):
		self.tabIdx = self.tabs.currentIndex()
	
	def filterList(self, text):
		for i in range(self.RecTab.listWidget.count()):
			item = self.RecTab.listWidget.item(i)
			item.setHidden(text not in item.text())

	def changeTabSelectedSeq(self):
		self.tabList[self.tabIdx].selectedSeq = self.tabList[self.tabIdx].listWidget.currentItem().text()
		
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

	def displaySeqBlastResData(self):
		self.tabList[self.tabIdx].blastResTable.clear()
		if self.dfBlast.size == 0:
			return
		self.dfBlast.fillna('', inplace=True)
		tmpdf = self.dfBlast[self.dfBlast["qaccver"] == self.tabList[self.tabIdx].selectedSeq]
		
		self.tabList[self.tabIdx].blastResTable.setRowCount(tmpdf.shape[0])
		self.tabList[self.tabIdx].blastResTable.setColumnCount(tmpdf.shape[1])
		self.tabList[self.tabIdx].blastResTable.setHorizontalHeaderLabels(tmpdf.columns)
		
		# returns pandas array object
		row_num = 0
		for row in tmpdf.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.tabList[self.tabIdx].blastResTable.setItem(row_num, col_index, tableItem)
			row_num += 1
		self.tabList[self.tabIdx].blastResTable.resizeColumnsToContents()

	def displayGeneIDGraph(self):
		self.tabList[self.tabIdx].geneIDBrowser = QWebEngineView()
		tmpdf = self.recgID_for_graphdf[self.recgID_for_graphdf["Sequences"] == self.tabList[self.tabIdx].selectedSeq]
		tmpdf = tmpdf.sort_values("qstart")
		fig = px.scatter(tmpdf, x="Gene", y="Sequences", size="pident", color="Lineage",text="Database",
		hover_data=["Database","qstart","qend","pident"], title="Gene identification recombinant sequences", 
		color_discrete_map = self.geneIDcolor_dict
		)
		self.tabList[self.tabIdx].geneIDBrowser.setHtml(fig.to_html(include_plotlyjs='cdn'))
		self.tabList[self.tabIdx].mainLayout.addWidget(self.tabList[self.tabIdx].geneIDBrowser,2,1)
		self.tabList[self.tabIdx].setLayout(self.tabList[self.tabIdx].mainLayout)

	def displaySeqSNPResData(self):
		self.tabList[self.tabIdx].SnpResTable.clear()
		if self.dfSNP.size == 0:
			return
		self.dfSNP.fillna('', inplace=True)
		tmpdf = self.dfSNP[self.dfSNP["Sequences"] == self.tabList[self.tabIdx].selectedSeq]
		self.tabList[self.tabIdx].SnpResTable.setRowCount(tmpdf.shape[0])
		self.tabList[self.tabIdx].SnpResTable.setColumnCount(tmpdf.shape[1])
		self.tabList[self.tabIdx].SnpResTable.setHorizontalHeaderLabels(tmpdf.columns)
		
		# returns pandas array object
		row_num = 0
		for row in tmpdf.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.tabList[self.tabIdx].SnpResTable.setItem(row_num, col_index, tableItem)
			row_num += 1
		
		self.tabList[self.tabIdx].SnpResTable.resizeColumnsToContents()
		self.tabList[self.tabIdx].SnpResTable.setFixedHeight(80)

	def displaySNPGraph(self):
		self.tabList[self.tabIdx].snpBrowser = QWebEngineView()
		tmpdf = self.recSNP_for_graphdf[self.recSNP_for_graphdf["Sequences"] == self.RecTab.selectedSeq]
		fig = px.scatter(tmpdf, x="SNP reference pos", y="Sequences", color="Lineage",
			hover_data=["SNP reference pos","Lineage","Nucleotide"], title="SNP identification recombinant sequences", 
			color_discrete_map=self.snp_color_dict
		)
		fig.update_layout(xaxis_type = 'linear')
		self.tabList[self.tabIdx].snpBrowser.setHtml(fig.to_html(include_plotlyjs='cdn'))
		self.tabList[self.tabIdx].mainLayout.addWidget(self.tabList[self.tabIdx].snpBrowser,5,1)
		self.tabList[self.tabIdx].setLayout(self.tabList[self.tabIdx].mainLayout)

	def displaySeqCancerSNPResData(self):
		self.tabList[self.tabIdx].cancerSnpResTable.clear()
		if self.dfCancerSNP.size == 0:
			return
		self.dfCancerSNP.fillna('', inplace=True)
		tmpdf = self.dfCancerSNP[self.dfCancerSNP["Sequences"] == self.tabList[self.tabIdx].selectedSeq]
		# tmpdf = tmpdf.T
		self.tabList[self.tabIdx].cancerSnpResTable.setRowCount(tmpdf.shape[0])
		self.tabList[self.tabIdx].cancerSnpResTable.setColumnCount(tmpdf.shape[1])
		self.tabList[self.tabIdx].cancerSnpResTable.setHorizontalHeaderLabels(tmpdf.columns)
		
		# returns pandas array object
		row_num = 0
		for row in tmpdf.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.tabList[self.tabIdx].cancerSnpResTable.setItem(row_num, col_index, tableItem)
			row_num += 1
		self.tabList[self.tabIdx].cancerSnpResTable.resizeColumnsToContents()
		self.tabList[self.tabIdx].cancerSnpResTable.setFixedHeight(80)

	def updateRecLED(self):
		"""
		Self explanatory
		This function applies only to the Putative recombinant sequences tab
		"""
		self.RecTab._blastLED.setValue(self.recombinationStatus[self.RecTab.selectedSeq]["Genes"])
		self.RecTab._snpLED.setValue(self.recombinationStatus[self.RecTab.selectedSeq]["SNP"])

	def createSimplot(self):
		"""
		Self explanatory
		This function applies only to the Putative recombinant sequences tab
		"""
		self.tmpOut = self.outdir / "tmp_dir"
		self.aln_f = simplot.align(self.muscle_bin, self.simplotDBFasta, self.totalSequenceRecDict,
		 self.RecTab.selectedSeq, self.tmpOut
		)
		self.window_size = self.RecTab.windowSizeSpinBox.value()
		self.step = self.RecTab.stepSizeSpinBox.value()
		self.simplotW = Simplot_page(self.RecTab.selectedSeq, self.window_size, 
		self.step, self.aln_f
		)
		mainW.openWindows["Simplot"] = self.simplotW
		mainW.openWindows["Simplot"].show()

	def closeEvent(self, event):
		self.totalSequenceRecDict.close()
		del mainW.openWindows["Results"]
	
class Simplot_page(QMainWindow):
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

if __name__ == "__main__":
	app = QApplication(sys.argv)
	mainW = Main_page()
	screen = app.primaryScreen()
	sys.exit(app.exec_())

# TODO: Main Window should be able to be closed and not terminate the app?
# TODO: Write in MD. Never remove the logfile. It is necessary for the applicationi
# to continue the process
# TODO: Each  time the tabs are changed the QWebEngines become smaller and smaller
# TODO: Need to check in windows os
