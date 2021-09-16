################################################################################
## Basic GUI idea and core were incepted and created by: WANDERSON M.PIMENTA
## USING WITH: Qt Designer and PySide2
## V: 1.0.0
##
## There are limitations on Qt licenses if you want to use your products
## commercially, I recommend reading them on the official website:
## https://doc.qt.io/qtforpython/licenses.html
##
##
## All the PySide2 API calls were rewritten in PyQt5 library
## Further GUI changes (in this script) were made by NIKOLAIDIS MARIOS
## All the CLI functions and analysis were created by NIKOLAIDIS MARIOS
################################################################################

import os
import re
import sys
import blast # Tool module
import wrapper # Tool module
import simplot # Tool module
import pathlib
import logging
import files_rc
import numpy as np
import pandas as pd
import pyqtgraph as pg
import plotly.express as px
from datetime import date
from Bio import SeqIO
from QLed import QLed
from PyQt5.QtGui import (
	QFont, QFontDatabase, QIcon,
	QDesktopServices
)
from PyQt5.QtCore import ( 
	Qt, QObject, QThread, pyqtSignal, 
	QUrl, QCoreApplication, QSize,
	QMetaObject	
)
from PyQt5.QtWidgets import (
	QApplication, QLabel, QListWidgetItem, 
	QMainWindow, QProgressBar, QSpacerItem, QWidget, 
	QTableWidget, QTableWidgetItem, 
	QFileDialog, QGridLayout, QVBoxLayout,
	QPushButton, QHBoxLayout, QButtonGroup,
	QMessageBox, QListWidget, QSpinBox,
	QLineEdit, QFrame, QStackedWidget, QCheckBox,
	QTextBrowser, QErrorMessage
)
from PyQt5.QtWebEngineWidgets import QWebEngineView
from ete3 import ( 
	Tree, TreeStyle, TextFace, NodeStyle
)

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
		# self.logW = LogWindow()
		GuiFunctions.addNewMenu(mainW, "Progress", "progressButton", "url(:/16x16/icons/16x16/cil-magnifying-glass.png)", True)
		print(mainW.buttons)
		mainW.buttons["Progress"].click()
		GuiFunctions.makeHomeUnavailable(self)
		
	finished = pyqtSignal()

	# def updateWidgets(self, incident):
	# 	if incident == "starting":
	# 		self.logW.setCentralWidget(self.logW.pipelineWidget)
	# 	if incident == "finished":
	# 		self.logW.setCentralWidget(self.logW.finishedWidget)

	def runPipeline(self) -> None:
		
		paramsdf = mainW.paramsdf
		# query_f_path, hpv16error =	wrapper.main(paramsdf)
		# if hpv16error == True:
		# 	GuiFunctions.showError(self, "No HPV16 sequences were identified. Programme execution has been stopped.")
		# 	self.logW.close()
		# mainW.paramsdf.loc["query"] = query_f_path
		self.finished.emit()		
		return

class MainWindow(QMainWindow):
	def __init__(self):
		QMainWindow.__init__(self)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		self.ui.retranslateUi(self)
		startSize = QSize(1000, 720)
		self.resize(startSize)
		self.setMinimumSize(startSize)

		self.buttons = {}		
		
		self.ui.stackedWidget.setMinimumWidth(20)
		GuiFunctions.addNewMenu(self, "Main Page", "homeButton", "url(:/16x16/icons/16x16/cil-home.png)", True)
		GuiFunctions.addHomeButtons(self)

		# Select home menu on init
		GuiFunctions.selectStandardMenu(self, "homeButton")

		self.ui.stackedWidget.setCurrentWidget(self.ui.homePage)

		# Remove the content margins from the window
		self.ui.horizontalLayout.setContentsMargins(0, 0, 0, 0)
		
		params = wrapper._defaultparams()
		self.paramsdf = pd.DataFrame.from_dict(params,orient='index')
		self.paramsdf.rename(columns={0:"Value"},inplace=True)
		self.binaries = wrapper._init_binaries(sys.platform) # Init some necessary params
		self.openWindows = {}
		self.show()

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
		if btnWidget.objectName() == "progressButton":
			self.ui.stackedWidget.setCurrentWidget(self.ui.progressPage)
			GuiFunctions.resetStyle(self, "progressButton")
			GuiFunctions.updatePageLabel(self, "Progress")
			btnWidget.setStyleSheet(GuiFunctions.selectMenu(btnWidget.styleSheet()))
		if btnWidget.objectName() == "resultsButton":
			self.ui.stackedWidget.setCurrentWidget(self.ui.resultsPage)
			GuiFunctions.resetStyle(self, "resultsButton")
			GuiFunctions.updatePageLabel(self, "Results")
			btnWidget.setStyleSheet(GuiFunctions.selectMenu(btnWidget.styleSheet()))
		if btnWidget.objectName() == "loadFastaButton":
			GuiFunctions.getFastaFile(self)
		if btnWidget.objectName() == "selectOutdirButton":
			GuiFunctions.getOutdir(self)
		if btnWidget.objectName() == "loadResultsButton":
			GuiFunctions.loadResDir(self)
		if btnWidget.objectName() == "runPipelineButton":
			GuiFunctions.execPipeline(self)
		if btnWidget.objectName() == "openHelpButton":
			pass
			# GuiFunctions.showHelp(self)
	
	def closeEvent(self, event):
		self.quitMsg = QMessageBox()
		reply = self.quitMsg.question(self, "Close Window", 
		"Do you want to close this window?\nThe application will be terminated",
		QMessageBox.Yes | QMessageBox.No
		)
		if reply == QMessageBox.Yes:
			event.accept()
			logging.shutdown()
			for w in mainW.openWindows:
				mainW.openWindows[w].close()
		else:
			event.ignore()


class GuiFunctions(MainWindow):
	
	def labelTitle(self, text):
		self.ui.label_title_bar_top.setText(text)

	# Dynamic menus
	def addNewMenu(self, name, objName, icon, isTopMenu):
		font = QFont()
		font.setFamily(u"Segoe UI")
		button = QPushButton(str(1),self)
		button.setObjectName(objName)
		button.setMinimumSize(QSize(0, 70))
		button.setLayoutDirection(Qt.LeftToRight)
		button.setFont(font)
		button.setStyleSheet(Style.style_bt_standard.replace('ICON_REPLACE', icon))
		button.setText(name)
		button.setToolTip(name)
		button.clicked.connect(self.Button)
		self.buttons[name] = button
		if isTopMenu:
			self.ui.menusLayout.addWidget(button)
		else:
			self.ui.layout_menu_bottom.addWidget(button)

	def selectMenu(getStyle):
		select = getStyle + ("QPushButton { border-right: 7px solid rgb(44, 49, 60); }")
		return select

	def deselectMenu(getStyle):
		deselect = getStyle.replace("QPushButton { border-right: 7px solid rgb(44, 49, 60); }", "")
		return deselect

	# Starting selections
	def selectStandardMenu(self, widget):
		for w in self.ui.frame_left_menu.findChildren(QPushButton):
			if w.objectName() == widget:
				w.setStyleSheet(GuiFunctions.selectMenu(w.styleSheet()))
	
	# Reset selection
	def resetStyle(self, widget):
		for w in self.ui.frame_left_menu.findChildren(QPushButton):
			if w.objectName() != widget:
				w.setStyleSheet(GuiFunctions.deselectMenu(w.styleSheet()))

	def updatePageLabel(self, text):
		newText = f" Viewing - {text.upper()}"
		self.ui.pageNameInfo.setText(newText)

	def addHomeButtons(self):
		self.homeButtonsObj = ["loadFastaButton", "selectOutdirButton" ,"loadResultsButton",
		"openHelpButton"		
		]
		self.homeButtonNames = ["Load the fasta file to analyze", "Select directory to write the output",
		"Load results", "View help videos  (Opens a separate window)"
		]
		# self.homeButtons = []
		indeces = range(len(self.homeButtonsObj))
		menuSpacer = QSpacerItem(20, 20)
		for idx in indeces:
			objName = self.homeButtonsObj[idx]
			name = self.homeButtonNames[idx]
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
			if objName == "loadResultsButton":
				self.ui.homeMenuGridLayout.addWidget(self.ui.orLabel)
			if objName == "openHelpButton":
				self.ui.homeMenuGridLayout.addItem(menuSpacer)
			self.ui.homeMenuGridLayout.addWidget(homeButton)
			# self.homeButtons.append(homeButton)
	
	def enablePipeline(self):
		# Delete all the existing widgets
		for i in reversed(range(self.ui.homeMenuGridLayout.count())): 
			if self.ui.homeMenuGridLayout.itemAt(i).widget() != None:
				self.ui.homeMenuGridLayout.itemAt(i).widget().setParent(None)

		self.ui.homeMenuGridLayout.addWidget(self.ui.toolLabel)
		self.ui.homeMenuGridLayout.addWidget(self.ui.menuLabel)

		self.homeButtonsObj = ["runPipelineButton", "selectOutdirButton" ,"loadResultsButton",
		"openHelpButton"		
		]
		self.homeButtonNames = ["Execute the pipeline", "Select directory to write the output",
		"Load results", "View help videos  (Opens a separate window)"
		]
		# self.homeButtons = []
		indeces = range(len(self.homeButtonsObj))
		menuSpacer = QSpacerItem(20, 20)
		
		for idx in indeces:
			objName = self.homeButtonsObj[idx]
			name = self.homeButtonNames[idx]
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
			if objName == "runPipelineButton":
				homeButton.setStyleSheet(Style.style_bt_highlight.replace('ICON_REPLACE', "url(:/16x16/icons/16x16/cil-home.png)"))
			if objName == "loadResultsButton":
				self.ui.homeMenuGridLayout.addWidget(self.ui.orLabel)
			if objName == "openHelpButton":
				self.ui.homeMenuGridLayout.addItem(menuSpacer)
			self.ui.homeMenuGridLayout.addWidget(homeButton)
			# self.homeButtons.append(homeButton)

	def makeHomeUnavailable(self):
		pass

	def showError(self,text: str, windowTitle: str = "Error") -> None:
		errorMsg = QErrorMessage(parent = mainW)
		errorMsg.setWindowTitle(windowTitle)
		errorMsg.setWindowModality(Qt.WindowModal)
		errorMsg.showMessage(text)
		errorMsg.setFixedSize(QSize(500,200))
		return 

	############# Home button actions
	def getFastaFile(self: QMainWindow) -> None:
		# Self is mainW
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
			mainW.fasta_path = fasta_path
			mainW.paramsdf.loc["out"] = mainW.paramsdf.loc["out"] / pathlib.Path(str(fasta_path.stem) + "_" + str(date.today()).replace("-","_"))
			mainW.paramsdf.loc["query"] = str(mainW.fasta_path)
			mainW.paramsdf.loc["in"] = str(mainW.fasta_path.parent)
			mainW.outdir = mainW.paramsdf.loc["out", "Value"]
			GuiFunctions.enablePipeline(self)
			return 
		else:
			return
	
	def getOutdir(self) -> None:
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
			mainW.outdir = pathlib.Path(response)
			mainW.paramsdf.loc["out"] = str(mainW.outdir)
			return
	
	def loadResDir(self) -> None:
		response = QFileDialog.getExistingDirectory(
			parent=self,
			caption="Select the results directory",
			)
		if response == "":
			return
		else:
			mainW.outdir = pathlib.Path(response)
			mainW.paramsdf.loc["out"] = str(mainW.outdir)
			logfile = mainW.outdir / pathlib.Path(".logfile.log")
			if logfile.exists() == False:
				GuiFunctions.showError(self, F"Cannot load results from specified directory.\nThe file .logfile.log has been deleted.\nPlease re-run the analysis")
				return 
			mainW.paramsdf = pd.DataFrame(index=[],columns=["Value"])
			fin = open(logfile,"r")
			for line in fin.readlines():
				m = re.match(r".+\tparam: (\S+),(\S+)",line)
				if m:
					index, val = m.group(1), m.group(2)
					if re.match(r"\/|\\",val):
						val = pathlib.Path(val)	
					mainW.paramsdf.loc[index] = val
			mainW.fasta_path = mainW.paramsdf.loc["query","Value"]
			if mainW.outdir != mainW.fasta_path.parent:
				mainW.paramsdf.loc["query","Value"] = mainW.outdir / mainW.fasta_path.name
				mainW.fasta_path = mainW.paramsdf.loc["query","Value"]
			mainW.totalSequenceRecDict = SeqIO.index(str(mainW.fasta_path),"fasta")
			mainW.analyzedSeqs = list(mainW.totalSequenceRecDict.keys())
			GuiFunctions.addNewMenu(self, "Results", "resultsButton", "url(:/16x16/icons/16x16/cil-magnifying-glass.png)", True)
			GuiFunctions.initVariables(self)
			mainW.showMaximized()
			mainW.buttons["Results"].click()
			return

	def execPipeline(self) -> None:
		"""
		Self is GuiFunctions, child of MainWindow class
		By default the whole pipeline will be executed. It will create another thread to run, so the application doesn't freeze
		Step 1: Create the worker class that will execute all the necessary functions
		Step 2: Create the QThread
		Step 3: Create the worker object
		Step 4: Move the worker to the QThread
		Step 5: Connect the various signals from and to the worker
		Step 6: Start the thread
		Step 7: Clear the central widget and load the results
		"""
		self.thread = QThread() # Step 2
		self.worker = Worker() # Step 3
		self.worker.moveToThread(self.thread) # Step 4
		
		self.thread.started.connect(self.worker.runPipeline) # Step 5
		# Create a new window for the progress
		self.worker.finished.connect(self.thread.quit) # Step 5
		self.worker.finished.connect(self.worker.deleteLater) # Step 5
		self.thread.finished.connect(self.thread.deleteLater) # Step 5
		self.thread.start() # Step 6
		# self.worker.updateWidgets("starting")
		# self.thread.finished.connect(lambda: self.worker.updateWidgets("finished")) # Step 7
		
		# self.thread.finished.connect(lambda: GuiFunctions.initVariables(self))
		# self.thread.finished.connect(lambda: GuiFunctions.addNewMenu(self, "Results", "resultsButton", "url(:/16x16/icons/16x16/cil-magnifying-glass.png)", True)) # Step 7
		# self.thread.finished.connect(lambda: mainW.showMaximized())
		# self.thread.finished.connect(lambda: mainW.buttons["Results"].click())

	##### Results page functions
	def initVariables(self: QMainWindow) -> None:
		# self is mainW
		self.dfBlast = pd.read_excel(self.outdir / "GeneIdentification_results.xlsx", engine="openpyxl")
		self.dfLineageSnp = pd.read_excel(self.outdir / "LineageSpecificSNPs.xlsx", engine="openpyxl", index_col=0)
		self.dfCancerSnp = pd.read_excel(self.outdir / "cancerSNP_results.xlsx", engine="openpyxl")
		self.dfGeneIdentification = pd.read_excel(self.outdir / "GeneIdentification_results.xlsx", engine="openpyxl", index_col=0)
		self.fasta_path = self.paramsdf.loc["query", "Value"]
		self.totalSequenceRecDict = SeqIO.index(str(self.fasta_path), "fasta")
		self.analyzedSeqs = list(self.totalSequenceRecDict.keys())
		self.muscle_bin= self.binaries[2]
		self.simplotDBFasta = self.paramsdf.loc["SimplotRef_database","Value"]
		self.geneIDColorDict = {
		"A":"#00ff00",
		"B":"#0080ff",
		"C":"#ff8000",
		"D":"#ff0000"
		}
		self.lineageSpecificSnpColorDict = {
		"Lin_A":"#00ff00",
		"Lin_B":"#0080ff",
		"Lin_C":"#ff8000",
		"Lin_D":"#ff0000",
		"Lin_BCD":"#ebeb34",
		"Other":"#ffffff"
		}
		
		self.genes = ["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
		self.gene_name_regex = re.compile(r"^\S+_([E|L]\d+)_aln.+$")
		self.trees_dir = self.outdir/"Phylogenetic_Trees"
		self.trees_l = os.listdir(self.trees_dir)
		self.trees_l = [pathlib.Path(self.trees_dir) / t for t in self.trees_l]
		if len(self.trees_l) == 0:
			pass
			# Throwing two errors makes the app crash
			# Should be different threads
			# GuiFunctions.showError(self, "No tree (newick) files have been found", "Information")
		else:
			self.trees = {}
			for t in self.trees_l:
				self.trees[t.name.replace("_aln_NJ_tree.nwk","")] = Tree(str(t))
			# TODO: Maybe optimize for memory?
		mainW.ui.translateResultsUI(self)
		self.putRecSeqs = GuiFunctions.loadRecList(self)
		GuiFunctions.showError(self, F"{len(self.putRecSeqs)} putative recombinants have been found", "Information")
		return

	def updateListWidget(self) -> None:
		"""
		Update the listWidget to display all the analyzed sequences
		or the putative recombinant sequences
		"""
		status = mainW.ui.recSeqsCheckBox.isChecked()
		mainW.ui.listWidget.clear()
		if status:
			mainW.ui.listWidget.addItems(mainW.putRecSeqs)
		else:
			mainW.ui.listWidget.addItems(mainW.analyzedSeqs)
		return

	def filterList(self) -> None:
		for i in range(mainW.ui.listWidget.count()):
			item = mainW.ui.listWidget.item(i)
			item.setHidden(self not in item.text())
		return 

	def updateResultsTables(self) -> None:
		"""
		Update the tables depending on their object name
		Self is listWidgetItem
		"""
		mainW.selectedSeq = self.text()
		tables = [mainW.ui.blastResTable,mainW.ui.lineageSnpTable,mainW.ui.lineageSumsTable,mainW.ui.cancerSnpTable]
		indeces = range(len(tables))
		for idx in indeces:
			table = tables[idx]
			tableName = table.objectName()
			if tableName == "blastResTable":
				if mainW.dfBlast.size == 0: return
				tmpdf = mainW.dfBlast[mainW.dfBlast["Query sequence"] == mainW.selectedSeq]
				tmpdf = tmpdf.drop("Query sequence",axis = 1)
				tmpdf.sort_values("Query start", inplace = True)
				tmpdf["Perc. identity"] = tmpdf["Perc. identity"].apply(lambda x: "{:0.2f}".format(x))
				tmpdf["E-value"] = tmpdf["E-value"].apply(lambda x: F"{x:.2e}")
			if tableName == "lineageSnpTable":
				if mainW.dfLineageSnp.size == 0: return
				tmpdf = mainW.dfLineageSnp[mainW.dfLineageSnp.index == mainW.selectedSeq].tail(1).T.head(67) # 67 == num_snps
				tmpdf["Lineage specific SNP"] = tmpdf.index
				tmpdf.rename({tmpdf.columns[0]:"Lineage"}, inplace = True, axis = 1)
				tmpdf = tmpdf[["Lineage specific SNP","Lineage"]]
			if tableName == "lineageSumsTable":
				tmpdf = mainW.dfLineageSnp[mainW.dfLineageSnp.index == mainW.selectedSeq].tail(1).T
				numRows = tmpdf.shape[0] - 67 # Num snps
				tmpdf = tmpdf.tail(numRows)
				tmpdf["Lineage"] = tmpdf.index
				tmpdf.index = range(len(tmpdf))
				tmpdf.rename({tmpdf.columns[0]:"Proportion"}, inplace=True, axis=1)
				tmpdf = tmpdf[["Lineage","Proportion"]]
			if tableName == "cancerSnpTable":
				if mainW.dfCancerSnp.size == 0: return
				tmpdf = mainW.dfCancerSnp[mainW.dfCancerSnp["Query sequence"] == mainW.selectedSeq]
				tmpdf.sort_values(by=["Query position"],inplace=True,ignore_index=True)
				tmpdf.drop("Unnamed: 0",axis=1, inplace=True)
				tmpdf = tmpdf[["SNP","Query nucleotide", "Query position", "E-value"]]
			if tableName == "lineageSumsTable":
				tmp = tmpdf["Proportion"].values
				zeroes = 0
				for i in tmp:
					if i == 0: zeroes += 1
				table.setRowCount(tmpdf.shape[0] - zeroes)
			else:
				table.setRowCount(tmpdf.shape[0])
			
			table.setColumnCount(tmpdf.shape[1])
			table.setHorizontalHeaderLabels(tmpdf.columns)
			rowNum = 0
			for row in tmpdf.iterrows():
				values = row[1]
				if tableName == "lineageSumsTable" and values["Proportion"] == 0:
					continue
				for colIdx, value in enumerate(values):
					tableItem = QTableWidgetItem(str(value))
					table.setItem(rowNum, colIdx, tableItem)
				rowNum += 1
			table.resizeColumnsToContents()
			
			table.horizontalHeader().setCascadingSectionResizes(True)
			if tableName != "lineageSnpTable":
				table.horizontalHeader().setDefaultSectionSize(int(1100 / tmpdf.shape[1]))
			if tableName == "lineageSumsTable":
				table.horizontalHeader().setDefaultSectionSize(120)
			if tableName == "cancerSnpTable":
				table.horizontalHeader().setDefaultSectionSize(int(900 / tmpdf.shape[1]))
		return

	def displayGraphs(self: QListWidgetItem) -> None:
		mainW.ui.graphicsList = []
		browsers = [mainW.ui.geneIDBrowser, mainW.ui.snpBrowser]
		data = {"Gene":"", "Lineage" : "", "Sublineage" : "",
			"Query start": 0, "Query end": 0, "Subject sequence": "", 
			"Subject start": 0, "Subject end": 0, "Aln length": 0, "Perc. identity": 0, "E-value": 1.0}
		for browser in browsers:
			browserName = browser.objectName()
			if browserName == "geneIDBrowser":
				tmpdf = mainW.dfGeneIdentification[mainW.dfGeneIdentification.index == mainW.selectedSeq]
				tmpdf = tmpdf.sort_values("Query start")
				idx = tmpdf.tail(1).index
				for gene in mainW.ui.genes:
					if gene not in tmpdf["Gene"].values:
						geneDF = pd.DataFrame(index=idx, data = data)
						geneDF["Gene"] = gene
						tmpdf = pd.concat([tmpdf, geneDF])		
				fig = px.scatter(tmpdf, x="Gene", y=tmpdf.index, size="Perc. identity", color="Lineage", # text="Sublineage",
				hover_data=["Sublineage","Query start","Query end","Perc. identity"], 
				color_discrete_map = mainW.geneIDColorDict, category_orders = {
					"Gene":["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
					}
				)
				fig.update_layout(height=200, legend_y=1.5)
			if browserName == "lineageSnpBrowser":
				tmpdf = mainW.dfLineageSnp[mainW.dfLineageSnp.index == mainW.selectedSeq]
				tmpdf = tmpdf.T.head(67) # num of lineage specific snps
				tmpdf.columns=["Nucleotide","Lineage"]
				tmpdf["SNP reference pos"] = tmpdf.index
				tmpdf["SNP reference pos"] = tmpdf["SNP reference pos"].apply(lambda x: int(x.split("_")[2]) )
				tmpdf.index = range(len(tmpdf))
				tmpdf["Query sequence"] = mainW.selectedSeq
				fig = px.scatter(tmpdf, x="SNP reference pos", y="Query sequence", color="Lineage",
				hover_data=["SNP reference pos","Lineage","Nucleotide"], 
				color_discrete_map=mainW.lineageSpecificSnpColorDict
				)
				fig.update_layout(xaxis_type = 'linear', height=200, legend_y=1.5)
			mainW.ui.graphicsList.append(fig)
			browser.setHtml(fig.to_html(include_plotlyjs='cdn'))
		return

	def displayCancerSnpInfo(self: QTableWidgetItem) -> None:
		cSnp = mainW.ui.cancerSnpTable.item(self.row(),0).text()
		m = re.match(r"(\S+_NuclPos_\d+).+", cSnp)
		cSnp = m.group(1) + ".md"
		fin = pathlib.Path(__file__).parent / pathlib.Path("resources") / pathlib.Path("cSNPinfo") / cSnp 
		mainW.ui.cancerSnpTextBrowser.setSource(QUrl(str(fin)))
		mainW.ui.cancerSnpTextBrowser.setAlignment(Qt.AlignJustify)
		return 

	def openPubmedUrl(self: QUrl) -> None:
		if "#citations" in self.toString():
			return
		else:
			QDesktopServices.openUrl(self)
		return 

	def eteInteractive(self) -> None:
		buttonId = mainW.ui.TreesRenderButtonGroup.button(self).text()
		if hasattr(mainW, "selectedSeq") == False:
			GuiFunctions.showError(self, F"Please select a sequence first")
			return
		if mainW.selectedSeq + "_" + buttonId not in mainW.trees:
			GuiFunctions.showError(self, F"Sequence {mainW.selectedSeq}\ndoes not have {buttonId} gene")
			return
		t = mainW.trees[mainW.selectedSeq + "_" + buttonId]
		t.ladderize(direction=1)
		ts = TreeStyle()
		ts.title.add_face(TextFace(buttonId + " Gene", fsize = 13), column = 1)

		# Styling for certain clades
		selectedSeqstyle = NodeStyle()
		selectedSeqstyle["bgcolor"] = "Gray"
		defaultStyle = NodeStyle()
		defaultStyle["bgcolor"] = "White"
		linAStyle = NodeStyle()
		linAStyle["bgcolor"] = "Green"
		linBStyle = NodeStyle()
		linBStyle["bgcolor"] = "SteelBlue"
		linCStyle = NodeStyle()
		linCStyle["bgcolor"] = "Orange"
		linDStyle = NodeStyle()
		linDStyle["bgcolor"] = "FireBrick"
		for leaf in t.iter_leaves():
			if leaf.name == mainW.selectedSeq:
				leaf.img_style = selectedSeqstyle
			if leaf.name != mainW.selectedSeq:
				leaf.img_style = defaultStyle
			if re.match(r"^A\d+_\S\d",leaf.name):
				leaf.img_style=linAStyle			
			if re.match(r"^B\d+_\S\d",leaf.name):
				leaf.img_style=linBStyle			
			if re.match(r"^C\d+_\S\d",leaf.name):
				leaf.img_style=linCStyle			
			if re.match(r"^D\d+_\S\d",leaf.name):
				leaf.img_style=linDStyle
		ts.show_branch_support = True
		t.ladderize(direction = 1)
		t.show(tree_style = ts, child_app = True)
		return

	def updateLed(self) -> None:
		if mainW.selectedSeq in mainW.putRecSeqs:
			mainW.ui._blastLED.setValue(1)
			mainW.ui._lineageSnpLED.setValue(1)
		else:
			mainW.ui._blastLED.setValue(0)
			mainW.ui._lineageSnpLED.setValue(0)
		return

	def createSimplot(self) -> None:
		mainW.tmpOut = mainW.outdir / ".tmp"
		mainW.aln_f = simplot.align(
			mainW.muscle_bin, mainW.simplotDBFasta, 
			mainW.totalSequenceRecDict, mainW.selectedSeq, 
			mainW.tmpOut
		)
		mainW.window_size = mainW.ui.windowSizeSpinBox.value()
		mainW.step = mainW.ui.stepSizeSpinBox.value()
		mainW.simplotW = Simplot_page(
			mainW.selectedSeq, mainW.window_size, 
			mainW.step, mainW.aln_f
		)
		mainW.simplotW.show()
		return 

	def saveGraphics(self) -> None:
		for i in range(len(mainW.ui.graphicsList)):
			fig = mainW.ui.graphicsList[i]
			name = pathlib.Path(mainW.outdir) / pathlib.Path("Graphics")
			if i == 0:
				name = name / pathlib.Path("GeneIdentification")
			if i == 1:
				name = name / pathlib.Path("LineageSpecificSNPs")
			name = name / (mainW.selectedSeq + ".html")
			fig.write_html(str(name),include_plotlyjs='cdn')
		return

	def loadRecList(self) -> list:
		f = mainW.outdir / "putativeRecombinants.txt"
		fin = open(str(f),"r")
		recList = [line.rstrip() for line in fin.readlines()]
		return recList

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

		self.finishedWidget = QLabel("Finished! \n Please don't delete the (hidden) logfile.log \n You can safely close this window")
		self.finishedWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
		
		self.widgets = {
			"pipelineWidget" : self.pipelineWidget,
			"finishedWidget" : self.finishedWidget,
			"progressWidget" : None
		}
		mainW.openWindows["log"] = self
		self.show()

class Simplot_page(QMainWindow):
	def __init__(self, sequence, window_size, step, aln_f, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Similarity plot")
		self.resize(1000,1000)
		self.setMinimumSize(QSize(1000,500))
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
		self.pen_green = pg.mkPen(color=(40, 237, 43),  width=2.5)
		self.pen_blue = pg.mkPen(color=(66, 170, 245),  width=2.5)
		self.pen_orange = pg.mkPen(color=(237, 155, 40),  width=2.5)
		self.pen_red = pg.mkPen(color=(255, 0, 0),  width=2.5)

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
		mainW.openWindows["simplot"] = self

	def screenCapture(self) -> None:
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
		return 

	def closeEvent(self,event):
		event.accept()
		simplot.rm_tmpfile(self.aln_f)

class Ui_MainWindow(QMainWindow):
	def setupUi(self, MainWindow):
		MainWindow.resize(1000, 720)
		MainWindow.setMinimumSize(QSize(1000, 720))
		MainWindow.setWindowTitle("HPV-16 genotyper")
	   
		# Fonts        
		self.font1 = QFont("Segoe UI", pointSize= 12, weight= 75)
		self.font2 = QFont("Segoe UI", pointSize= 40)
		self.font3 = QFont("Segoe UI", pointSize= 14)
		self.font4 = QFont("Segoe UI", pointSize= 10, weight= 75)
		self.font5 = QFont("Segoe UI", pointSize= 25)
		
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
		self.pageNameInfo.setFont(self.font1)
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
		self.homeMenuGridLayout = QGridLayout(self.homePage)
		self.toolLabel = QLabel(self.homePage)

		self.toolLabel.setFont(self.font2)
		# self.toolLabel.setStyleSheet(u"")
		self.toolLabel.setAlignment(Qt.AlignCenter)

		self.homeMenuGridLayout.addWidget(self.toolLabel)

		self.menuLabel = QLabel(self.homePage)
		self.menuLabel.setFont(self.font3)
		self.menuLabel.setAlignment(Qt.AlignCenter)

		self.orLabel = QLabel(self.homePage)
		self.orLabel.setFont(self.font3)
		self.orLabel.setAlignment(Qt.AlignCenter) # This label will be added later

		self.homeMenuGridLayout.addWidget(self.menuLabel)

		self.verticalLayout_4.addWidget(self.stackedWidget)
		self.verticalLayout_2.addWidget(self.frame_content)

		# Progress page
		self.progressPage = QWidget()
		self.progressPage.setObjectName("progressPage")
		
		self.progressPageVboxLayout = QVBoxLayout(self.progressPage)
		
		self.progressPageGridLayout = QGridLayout()

		self.progressLabel = QLabel(self.progressPage)
		self.progressLabel.setText("Executing the pipeline")
		self.progressLabel.setFont(self.font5)
		self.toolLabel.setAlignment(Qt.AlignCenter)
		spacer = QSpacerItem(20, 20)
		self.progressPageGridLayout.addWidget(self.progressLabel)
		self.progressPageGridLayout.addItem(spacer)
		
		self.statusLines = ["HPV-16 filtering", "Gene Identification", 
			"Lineage specific SNPs scan", "cancer specific SNPs scan",
			"Gene alignment", "Trees calculation"
		]
		
		pos = 1
		for status in self.statusLines:
			print(status)
			label = QLabel()
			label.setText(status)
			label.setFont(self.font3)
			self.progressPageGridLayout.addWidget(label, pos, 0, 2, 1)
			self.progressPageGridLayout.addItem(spacer)
			pos += 1
		
		self.progressBars = []
		pos = 1
		for status in self.statusLines:
			progressBar = QProgressBar()
			progressBar.setMinimum = 0
			progressBar.setMaximum = 100
			progressBar.setValue(100)
			progressBar.setStyleSheet(Style.style_progress_bar)
			self.progressBars.append(progressBar)
			self.progressPageGridLayout.addWidget(progressBar, pos, 1, 2, 1)
			self.progressPageGridLayout.addItem(spacer)
			pos += 1

		self.progressPageVboxLayout.addLayout(self.progressPageGridLayout)

		# Results page
		self.resultsPage = QWidget()
		self.resultsPage.setObjectName("resultsPage")
		self.putativeRecombinants = QCheckBox() # Will use this to "transform" the list into a putative recombinants only

		# This will be used to put all the frames
		self.resultsPageGridLayoutForFrames = QGridLayout(self.resultsPage)
		
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
		# 1. Label
		self.listFrameLabel = QLabel(self.listFrame)
		self.listFrameLabel.setFont(self.font4)
		self.listFrameLabel.setText("Analyzed sequences")
		# 2. Line edit
		self.lineEdit = QLineEdit(self.frame_content_wid_1)
		self.lineEdit.setObjectName(u"lineEdit")
		self.lineEdit.setMinimumSize(QSize(0, 30))
		self.lineEdit.setStyleSheet(Style.style_line_edit)
		self.listWidget = QListWidget()
		self.listWidget.resize(100,100)
		self.listWidget.setStyleSheet(Style.style_list_widget)
		self.lineEdit.textChanged.connect(GuiFunctions.filterList)

		self.listWidget.itemDoubleClicked.connect(GuiFunctions.updateResultsTables)
		self.listWidget.itemDoubleClicked.connect(GuiFunctions.displayGraphs)
		self.listWidget.itemDoubleClicked.connect(GuiFunctions.updateLed)
		#3. Check box
		self.recSeqsCheckBox = QCheckBox(self.listFrame)
		self.recSeqsCheckBox.setStyleSheet(Style.style_checkbox)
		self.recSeqsCheckBox.setText("Show\nputative recombinants\nonly")
		self.recSeqsCheckBox.stateChanged.connect(GuiFunctions.updateListWidget)

		self.horizontalLayout_9.addLayout(self.gridLayout)
		self.verticalLayout_7.addWidget(self.frame_content_wid_1)
		self.verticalLayout_15.addWidget(self.frame_div_content_1)

		self.gridLayout.addWidget(self.listFrameLabel, 0, 0, 1, 1)
		self.gridLayout.addWidget(self.lineEdit, 1, 0, 1, 1)
		self.gridLayout.addWidget(self.recSeqsCheckBox,1,1,-1,1)
		self.gridLayout.addWidget(self.listWidget, 2, 0, 2, 1)

		############# Blast results frame #############
		self.blastResFrame = QFrame(self.resultsPage)
		self.blastResFrame.setMinimumSize(QSize(0, 150))
		self.blastResFrame.setStyleSheet(Style.style_results_frame)
		self.blastResFrame.setFrameShape(QFrame.StyledPanel)
		self.blastResFrame.setFrameShadow(QFrame.Raised)

		self.blastResVerticalLayout = QVBoxLayout(self.blastResFrame)
		self.blastResVerticalLayout.setObjectName(u"blastResVerticalLayout")
		
		self.blastResGridLayout = QGridLayout()
		
		# Widgets
		# 1. Label
		self.blastLabel = QLabel(self.blastResFrame)
		self.blastLabel.setObjectName("blastLabel")
		self.blastLabel.setFont(self.font4)
		self.blastLabel.setText("Gene identification")
		# 2. Table
		self.blastResTable = QTableWidget(self.blastResFrame)
		self.blastResTable.setObjectName("blastResTable")
		# Set Rows and Columns to a first value so i can init the view
		self.blastResTable.setColumnCount(10)
		self.blastResTable.setRowCount(10)
		self.blastResTable.setStyleSheet(Style.style_table_standard)
		self.blastResTable.horizontalHeader().setDefaultSectionSize(150)
		# 3. LED
		self._blastLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self._blastLED.setFixedSize(18,18)
		# 5. Webengine view
		self.geneIDBrowser = QWebEngineView()
		self.geneIDBrowser.setObjectName("geneIDBrowser")

		self.blastResGridLayout.addWidget(self.blastLabel, 0, 0, 1, 1)
		self.blastResGridLayout.addWidget(self.blastResTable, 1, 0, 3, 5)
		self.blastResGridLayout.addWidget(self.geneIDBrowser, 0, 5, -1, 3)
		self.blastResGridLayout.addWidget(self._blastLED, 0, 5, Qt.AlignLeft)

		self.blastResVerticalLayout.addLayout(self.blastResGridLayout)

		############# Lineage Specific SNPs Frame
		self.lineageSnpFrame = QFrame(self.resultsPage)
		self.lineageSnpFrame.setMinimumSize(QSize(0, 150))
		self.lineageSnpFrame.setStyleSheet(Style.style_results_frame)
		self.lineageSnpFrame.setFrameShape(QFrame.StyledPanel)
		self.lineageSnpFrame.setFrameShadow(QFrame.Raised)

		self.lineageSnpVerticalLayout = QVBoxLayout(self.lineageSnpFrame)
		self.lineageSnpVerticalLayout.setObjectName(u"lineageSnpVerticalLayout")
		
		self.lineageSnpGridLayout = QGridLayout()
		
		# Widgets
		# 1. Label
		self.lineageSnpLabel = QLabel(self.lineageSnpFrame)
		self.lineageSnpLabel.setFont(self.font4)
		self.lineageSnpLabel.setText("Lineage specific SNPs")
		# 2. Basic Table
		self.lineageSnpTable = QTableWidget(self.lineageSnpFrame)
		self.lineageSnpTable.setObjectName("lineageSnpTable")
		# Set Rows and Columns to a first value so i can init the view
		self.lineageSnpTable.setColumnCount(2)
		self.lineageSnpTable.setRowCount(4)
		self.lineageSnpTable.setMaximumWidth(285)
		self.lineageSnpTable.setStyleSheet(Style.style_table_standard)
		# 3. LED
		self._lineageSnpLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self._lineageSnpLED.setFixedSize(18,18)
		# 4. Dominant lineage table
		self.lineageSumsTable = QTableWidget(self.lineageSnpFrame)
		self.lineageSumsTable.setObjectName("lineageSumsTable")
		self.lineageSumsTable.setColumnCount(2)
		self.lineageSumsTable.setRowCount(4)
		self.lineageSumsTable.setMaximumWidth(270)
		self.lineageSumsTable.setStyleSheet(Style.style_table_standard)
		# 5. Dominant lineage label
		self.lineageSumsLabel = QLabel(self.lineageSnpFrame)
		self.lineageSumsLabel.setFont(self.font4)
		self.lineageSumsLabel.setText("Proportion of lineages")
		# 6. Webengine view
		self.snpBrowser = QWebEngineView()
		self.snpBrowser.setObjectName("lineageSnpBrowser")

		self.lineageSnpGridLayout.addWidget(self.lineageSnpLabel, 0, 0)
		self.lineageSnpGridLayout.addWidget(self.lineageSnpTable, 1, 0, 1, 1)
		self.lineageSnpGridLayout.addWidget(self.lineageSumsLabel, 0, 1, 1, 1)
		self.lineageSnpGridLayout.addWidget(self.lineageSumsTable, 1, 1, 1, 1)
		self.lineageSnpGridLayout.addWidget(self.snpBrowser, 0, 2,-1, 6)
		self.lineageSnpGridLayout.addWidget(self._lineageSnpLED, 0, 2, Qt.AlignLeft)

		self.lineageSnpVerticalLayout.addLayout(self.lineageSnpGridLayout)
		
		############# Cancer SNPs Frame
		self.cancerSnpFrame = QFrame(self.resultsPage)
		self.cancerSnpFrame.setMinimumSize(QSize(0, 150))
		self.cancerSnpFrame.setStyleSheet(Style.style_results_frame)
		self.cancerSnpFrame.setFrameShape(QFrame.StyledPanel)
		self.cancerSnpFrame.setFrameShadow(QFrame.Raised)

		self.cancerSnpVerticalLayout = QVBoxLayout(self.cancerSnpFrame)
		self.cancerSnpGridLayout = QGridLayout()
		
		# Widgets
		# 1. Label
		self.cancnerSnpLabel = QLabel(self.cancerSnpFrame)
		self.cancnerSnpLabel.setFont(self.font4)
		self.cancnerSnpLabel.setText("Increased cancer risk SNPs")
		# 2. Table
		self.cancerSnpTable = QTableWidget(self.cancerSnpFrame)
		self.cancerSnpTable.setObjectName("cancerSnpTable")
		# Set Rows and Columns to a first value so i can init the view
		self.cancerSnpTable.setColumnCount(5)
		self.cancerSnpTable.setRowCount(5)
		self.cancerSnpTable.setStyleSheet(Style.style_table_standard)
		self.cancerSnpTable.itemClicked.connect(GuiFunctions.displayCancerSnpInfo)
		# Text browser
		self.cancerSnpTextBrowser = QTextBrowser()
		self.cancerSnpTextBrowser.setMinimumSize(QSize(0,100))
		self.cancerSnpTextBrowser.setMaximumWidth(800)
		self.cancerSnpTextBrowser.setStyleSheet(Style.style_text_browser)
		self.cancerSnpTextBrowser.setOpenLinks(False)
		self.cancerSnpTextBrowser.anchorClicked.connect(GuiFunctions.openPubmedUrl)

		self.cancerSnpGridLayout.addWidget(self.cancnerSnpLabel, 0, 0, 1, 1)
		self.cancerSnpGridLayout.addWidget(self.cancerSnpTable, 1, 0, 2, 1)
		self.cancerSnpGridLayout.addWidget(self.cancerSnpTextBrowser, 1, 1,-1,-1)

		self.cancerSnpTable.horizontalHeader().setCascadingSectionResizes(True)
		self.cancerSnpTable.horizontalHeader().setDefaultSectionSize(100)
	
		self.cancerSnpTable.horizontalHeader().setStretchLastSection(True)
		self.cancerSnpTable.verticalHeader().setStretchLastSection(True)

		self.cancerSnpVerticalLayout.addLayout(self.cancerSnpGridLayout)
		
		############# Trees Frame
		self.treeFrame = QFrame(self.resultsPage)
		self.treeFrame.setMinimumSize(QSize(0, 150))
		self.treeFrame.setStyleSheet(Style.style_results_frame)
		self.treeFrame.setFrameShape(QFrame.StyledPanel)
		self.treeFrame.setFrameShadow(QFrame.Raised)
		self.treeGridLayout = QGridLayout(self.treeFrame)
		# Widgets
		# 1. Label
		self.treeLabel = QLabel(self.treeFrame)
		self.treeLabel.setFont(self.font4)
		self.treeLabel.setText("Explore trees interactively")
		# 2. Table
		self.TreesRenderButtonGroup = QButtonGroup()
		self.TreesRenderButtonGroup.buttonClicked[int].connect(GuiFunctions.eteInteractive)
		row = 1
		col = 0
		self.genes = ["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
		self.treeGridLayout.addWidget(self.treeLabel, 0, 0)
		for gene in self.genes:
			self.treeButton = QPushButton(gene)
			self.treeButton.setObjectName(F"treeButton{gene}")
			self.treeButton.setStyleSheet(Style.style_push_button)
			self.TreesRenderButtonGroup.addButton(self.treeButton)
			self.treeGridLayout.addWidget(self.treeButton, row, col)
			col += 1
			if col == 4:
				row += 1
				col = 0

		############# Simplot Frame
		self.SimplotFrame = QFrame(self.resultsPage)
		self.SimplotFrame.setMinimumSize(QSize(0, 150))
		self.SimplotFrame.setStyleSheet(
			"""
			background-color: rgb(39, 44, 54);
			border-radius: 5px;
			"""
		)
		self.SimplotFrame.setFrameShape(QFrame.StyledPanel)
		self.SimplotFrame.setFrameShadow(QFrame.Raised)

		self.SimplotVerticalLayout = QVBoxLayout(self.SimplotFrame)
		
		self.SimplotGridLayout = QGridLayout()
		
		# Widgets
		# 1. Label
		self.SimplotLabel = QLabel(self.SimplotFrame)
		self.SimplotLabel.setFont(self.font4)
		self.SimplotLabel.setText(" Create similarity plot ")
		
		# Simplot Widget
		self.simplotButtonAreaWidget = QWidget()
		self.simplotButtonAreaWidget.setMaximumSize(250,200)
		self.simplotButtonAreaWidgetSimplotLabel = QLabel("Similarity plot parameters")
		self.simplotButtonAreaWidgetWindowLabel = QLabel("Window size")
		self.simplotButtonAreaWidgetStepLabel = QLabel("Step")
		self.simplotButtonAreaWidgetLayout = QGridLayout()
		self.simplotButton = QPushButton("Calculate Simplot")
		self.simplotButton.setStyleSheet(Style.style_push_button)
		self.windowSizeSpinBox = QSpinBox(maximum=100000,value=300)
		self.windowSizeSpinBox.setStyleSheet(Style.style_spinbox)
		self.stepSizeSpinBox = QSpinBox(maximum=100000, value=100)
		self.stepSizeSpinBox.setStyleSheet(Style.style_spinbox)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButtonAreaWidgetSimplotLabel,0,0,1,2)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButtonAreaWidgetWindowLabel,1,1)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButtonAreaWidgetStepLabel,2,1)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButton,3,0,1,2)
		self.simplotButtonAreaWidgetLayout.addWidget(self.windowSizeSpinBox,1,0)
		self.simplotButtonAreaWidgetLayout.addWidget(self.stepSizeSpinBox,2,0)
		self.simplotButtonAreaWidget.setLayout(self.simplotButtonAreaWidgetLayout)
		self.simplotButton.clicked.connect(GuiFunctions.createSimplot)

		self.SimplotGridLayout.addWidget(self.SimplotLabel, 0, 0, 1, 1)
		self.SimplotGridLayout.addWidget(self.simplotButtonAreaWidget, 1, 0, 3, 5)

		self.SimplotVerticalLayout.addLayout(self.SimplotGridLayout)
		
		# Save graphics and list frame
		self.saveGraphicsFrame = QFrame(self.resultsPage)
		self.saveGraphicsFrame.setMinimumSize(QSize(0,150))
		self.saveGraphicsFrame.setStyleSheet(Style.style_results_frame)
		self.saveGraphicsFrame.setFrameShape(QFrame.StyledPanel)
		self.saveGraphicsFrame.setFrameShadow(QFrame.Raised)

		self.saveGraphicsVBoxLayout = QVBoxLayout(self.saveGraphicsFrame)
		self.saveGraphicsFrameButton = QPushButton(self.saveGraphicsFrame)
		self.saveGraphicsFrameButton.setMinimumSize(QSize(150,30))
		self.saveGraphicsFrameButton.setFont(QFont("Segoe UI",9))
		self.saveGraphicsFrameButton.setStyleSheet(Style.style_push_button)
		icon = QIcon()
		icon.addFile("url(:/16x16/icons/16x16/cil-magnifying-glass.png)", QSize(), QIcon.Normal, QIcon.Off)
		self.saveGraphicsFrameButton.setIcon(icon)
		self.saveGraphicsFrameButton.setText("Save graphics")
		self.saveGraphicsFrameButton.clicked.connect(GuiFunctions.saveGraphics)
		self.saveGraphicsVBoxLayout.addWidget(self.saveGraphicsFrameButton)
		# Add frames
		self.resultsPageGridLayoutForFrames.addWidget(self.listFrame,0,0,1,-1)
		self.resultsPageGridLayoutForFrames.addWidget(self.cancerSnpFrame,1,0,2,-1)
		self.resultsPageGridLayoutForFrames.addWidget(self.lineageSnpFrame,3,0,2,-1)
		self.resultsPageGridLayoutForFrames.addWidget(self.blastResFrame,5,0,2,-1)
		self.resultsPageGridLayoutForFrames.addWidget(self.treeFrame,7,0,1,4)
		self.resultsPageGridLayoutForFrames.addWidget(self.SimplotFrame,7,4,1,1)
		self.resultsPageGridLayoutForFrames.addWidget(self.saveGraphicsFrame,7,5,1,1)
		
		# / Add the results page here
		
		self.stackedWidget.addWidget(self.homePage)
		self.stackedWidget.addWidget(self.progressPage)
		self.stackedWidget.addWidget(self.resultsPage)

		# Necessary to display the UI correctly
		self.horizontalLayout_1.addWidget(self.frameContentRightPanel)
		self.verticalLayout.addWidget(self.frame_center)
		MainWindow.setCentralWidget(self.centralwidget)

		self.stackedWidget.setCurrentIndex(1)
		QMetaObject.connectSlotsByName(MainWindow)


	def retranslateUi(self, MainWindow):
		self.pageNameInfo.setText(QCoreApplication.translate("MainWindow", f" Viewing - HOME", None))
		self.toolLabel.setText(QCoreApplication.translate("MainWindow", "HPV-16 genotyper ", None))
		self.menuLabel.setText(QCoreApplication.translate("MainWindow", "Please use the following menu", None))
		self.orLabel.setText(QCoreApplication.translate("MainWindow", "OR", None))

	def translateResultsUI(self, mainWindow):
		self.listWidget.addItems(mainW.analyzedSeqs)

class Style():
	"""
	Class to hold all the stylings
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
	style_bt_highlight = (
		"""
		QPushButton {
			background-image: ICON_REPLACE;
			background-position: left center;
			background-repeat: no-repeat;
			border: none;
			border-left: 28px solid #3c57b0;
			background-color: #3c57b0;
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
			background-color: #486899;
			border-left: 28px solid #486899;
		}
		QPushButton:pressed {
			background-color: rgb(85, 170, 255);
			border-left: 28px solid rgb(85, 170, 255);
		}
		"""
	)
	style_table_standard = ("""
			QTableWidget {
				border: 2px;
				border-radius: 5px;
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
	style_spinbox = (
		"""
		QSpinBox {
			border: 2px solid rgb(52, 59, 72);
			border-radius: 5px;
			background-color: rgb(52, 59, 72);
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
		"""
		#TODO: Implement this
	)
	style_text_browser = (
		"""
		QTextBrowser {
			border: 2px solid rgb(52, 59, 72);
			border-radius: 5px;
			background-color: rgb(52, 59, 72);
		}
		"""
	)
	style_list_widget = (
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
	style_checkbox = (
		"""
		QCheckBox {
			text-align: center;
		}
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
	style_line_edit = (
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
	style_results_frame = (
		"""
		background-color: rgb(39, 44, 54);
		border-radius: 5px;
		"""
	)
	style_progress_bar = (
		"""
		QProgressBar {
			text-align: center;
		}
		QProgressBar::chunk {
			background-color: rgb(85, 170, 255);
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
