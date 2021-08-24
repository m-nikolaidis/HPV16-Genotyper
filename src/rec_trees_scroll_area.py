import re
import os
import pathlib
import sys
import pandas as pd
import plotly.express as px
from QLed import QLed
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtGui import QPixmap, QStandardItemModel, QStandardItem
from PyQt5.QtCore import Qt, QSortFilterProxyModel, QStringListModel, QModelIndex, QItemSelectionModel
from PyQt5.QtWidgets import (
	QWidget, QApplication, QListWidget, QListView, QHBoxLayout, QPushButton,
	QTableWidget, QTableWidgetItem, QGridLayout, QLabel, QSpinBox,
	QVBoxLayout, QScrollArea, QButtonGroup, QLineEdit, QAbstractItemView, 
	QCheckBox
)
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace

recseqs = ["MG850175.1","MG850176.1","MG850177.1","MG850178","MG850179","MG850180","MG850181","MG850182","MG850183","MG850184","MG850185","MG850186","MG850187","MG850188","MG850189","MG850190","MG850191","MG850192","MG850193","MG850194","MG850195","MG850196","MG850197","MG850198","MG850199","MG850200","MG850201","MG850202","MG850203","MG850204","MG850205","MG850206","MG850207","MG850208","MG850209","MG850210","MG850211","MG850212","MG850213","MG850214","MG850215","MG850216","MG850217","MG850218","MG850219","MG850220","MG850221","MG850222","MG850223","MG850224","MG850225","MG850226","MG850227","MG850228","MG850229","MG850230","MG850231","MG850232","MG850233","MG850234","MG850235","MG850236","MG850237","MG850238","MG850239","MG850240","MG850241","MG850242","MG850243","MG850244","MG850245","MG850246","MG850247","MG850248","MG850249","MG850250","MG850251","MG850252","MG850253","MG850254","MG850255","MG850256","MG850257","MG850258","MG850259","MG850260","MG850261","MG850262","MG850263","MG850264","MG850265","MG850266","MG850267"]
dfBlast = pd.read_excel("tmp/TestExcel.xlsx", engine="openpyxl")
dfSNP = pd.read_excel("tmp/TestExcelSNP.xlsx", engine="openpyxl")
recgID_for_graphdf = pd.read_excel("tmp/RecOnlyGeneID_to_usedf.xlsx", engine="openpyxl", index_col=0)
recSNP_for_graphdf = pd.read_excel("tmp/RecOnlySNP_to_usedf.xlsx", engine="openpyxl", index_col=0)

trees_dir = "tmp/trees"
trees = os.listdir(trees_dir)
trees = [pathlib.Path(trees_dir) / t for t in trees]
class RecTab(QWidget):
	""" A Fast test gui show how to create buttons in a ScrollArea"""
	def __init__(self):
		super(RecTab, self).__init__()
		self.mainLayout = QGridLayout()
		
		# Create the scrollable list view and Simplot parameter next to it
		self.recSeqsLabel = QLabel("Recombinant sequences")
		self.recSeqsLabel.setStyleSheet("font: bold 12px;")
		self.listViewBox = QWidget()
		self.listViewBoxLayout = QGridLayout()
		self.listViewListofStringModel = QStringListModel(recseqs)
		# Filter box
		self.model = QStandardItemModel(len(recseqs),1)
		for num, recseq in enumerate(recseqs):
			item = QStandardItem(recseq)
			self.model.setItem(num,0, item)
		self.firstRecSeqItem = self.model.item(0,0)
		self.firstRecSeqItemIdx = self.model.indexFromItem(self.firstRecSeqItem)

		self.filterProxyModel = QSortFilterProxyModel()
		self.filterProxyModel.setSourceModel(self.model)
		self.filterProxyModel.setFilterCaseSensitivity(Qt.CaseInsensitive)
		self.filterProxyModel.setFilterKeyColumn(0)

		self.listViewSearchBox = QLineEdit()
		self.listViewSearchBox.setStyleSheet("font: 12px;")
		self.listViewSearchBox.textChanged.connect(self.filterProxyModel.setFilterRegExp)       

		# Create listView and Update graphs when double clicking a recombinant sequence
		self.listView = QListView()
		self.listView.resize(100,100)
		self.listView.setModel(self.filterProxyModel)
		self.listView.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.listView.doubleClicked.connect(self.displaySeqBlastResData)
		self.listView.doubleClicked.connect(self.displayGeneIDGraph)
		self.listView.doubleClicked.connect(self.displaySeqSNPResData)
		self.listView.doubleClicked.connect(self.displaySNPGraph)
		self.listView.doubleClicked.connect(self.updateLED)

		# Widget with Simplot
		self.simplotButtonAreaWidget = QWidget()
		self.simplotButtonAreaWidget.setMaximumSize(250,200)
		self.simplotButtonAreaWidgetSimplotLabel = QLabel("Similarity plot parameters")
		self.simplotButtonAreaWidgetSimplotLabel.setStyleSheet("font: bold 12px;")
		self.simplotIgnoreGapsCheckBox = QCheckBox("Ignore gaps")
		self.simplotButtonAreaWidgetWindowLabel = QLabel("Window size")
		self.simplotButtonAreaWidgetStepLabel = QLabel("Step")
		self.simplotButtonAreaWidgetLayout = QGridLayout()
		self.simplotButton = QPushButton("Create Simplot")
		self.simplotWindowSizeSpinBox = QSpinBox(maximum=10000,value=300)
		self.simplotStepSizeSpinBox = QSpinBox(maximum=10000, value=100)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButtonAreaWidgetSimplotLabel,0,0)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotIgnoreGapsCheckBox,1,0)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButtonAreaWidgetWindowLabel,2,1)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButtonAreaWidgetStepLabel,3,1)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotButton,4,0,1,2)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotWindowSizeSpinBox,2,0)
		self.simplotButtonAreaWidgetLayout.addWidget(self.simplotStepSizeSpinBox,3,0)
		self.simplotButtonAreaWidget.setLayout(self.simplotButtonAreaWidgetLayout)

		
		# Add widgets to listViewBox widget
		self.listViewBoxLayout.addWidget(self.recSeqsLabel, 0,1)
		self.listViewBoxLayout.addWidget(self.listViewSearchBox, 1,1)
		self.listViewBoxLayout.addWidget(self.listView, 2,1)
		self.listViewBoxLayout.addWidget(self.simplotButtonAreaWidget, 2,2)
		self.listViewBox.setLayout(self.listViewBoxLayout)

		# HardCoded stuff for testing
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

		# Label of each box
		self.blastResLabelBox = QWidget()
		self.blastResLabelBoxLayout = QHBoxLayout()
		self.blastResLabel = QLabel("BLAST results")
		self.blastResLabel.setStyleSheet("font: bold 12px;")
		# Led widget for Recombination results
		self._blastLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self._blastLED.setFixedSize(18,18)
		self.blastResLabelBoxLayout.addWidget(self.blastResLabel)
		self.blastResLabelBoxLayout.addWidget(self._blastLED)
		self.blastResLabelBox.setLayout(self.blastResLabelBoxLayout)
		# TODO: Check stylsheet method

		self.snpResLabelBox = QWidget()
		self.snpResLabelBoxLayout = QHBoxLayout()
		self.snpResLabel = QLabel("SNP results")
		self.snpResLabel.setStyleSheet("font: bold 12px;")
		# Led widget for Recombination results
		self._snpLED=QLed(self, onColour=QLed.Green,offColour=QLed.Red, shape=QLed.Circle)
		self._snpLED.setFixedSize(18,18)
		
		self.snpResLabelBoxLayout.addWidget(self.snpResLabel)
		self.snpResLabelBoxLayout.addWidget(self._snpLED)
		self.snpResLabelBox.setLayout(self.snpResLabelBoxLayout)
		# TODO: Check stylsheet method


		# Add Blast results table widget
		self.blastResTable = QTableWidget()
		# Gene results plotly graph
		self.geneIDBrowser = QWebEngineView()
		# Add SNP results table widget
		self.SnpResTable = QTableWidget()
		# SNP results plotly graph
		self.snpBrowser = QWebEngineView()




		# Trees with highlighted sequence in QScrollArea 
		self.RecSeqsTreesWidget = QWidget()
		self.RecSeqsTreesWidgetLayout = QVBoxLayout()
		
		self.RecSeqsTreesScollArea = QScrollArea(widgetResizable=True)
		self.RecSeqsTreesScollArea.setStyleSheet(
			"""
			QWidget{ background-color: white } 
			QScrollBar{ background-color: none } 
			"""
		)
		self.RecTreesImgContentWidget = QWidget()
		self.RecSeqsTreesScollArea.setWidget(self.RecTreesImgContentWidget)
		self.RecSeqsTreesScollArea.layout = QVBoxLayout(self.RecTreesImgContentWidget)

		# Sort images correctly
		# highlight_dir = self.outdir / pathlib.Path("Graphics\\Trees_images")
		highlight_dir = pathlib.Path("tmp\\Images")
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
				self.RecSeqsTreesScollArea.layout.addWidget(label)
		
		self.RecInteractiveTrees = QWidget()
		self.RecTreesHbox = QHBoxLayout()
		self.RecInteractiveTreesLabel = QLabel("Explore trees interactively")
		self.RecTreesHbox.addWidget(self.RecInteractiveTreesLabel)
		self.RecTreesRenderButtonGroup = QButtonGroup()
		self.RecTreesRenderButtonGroup.buttonClicked[int].connect(self.eteInteractive)
		
		# Gene buttons
		for gene in self.genes:
			self.button = QPushButton(gene)
			self.RecTreesRenderButtonGroup.addButton(self.button)
			self.RecTreesHbox.addWidget(self.button)
		self.RecInteractiveTrees.setLayout(self.RecTreesHbox)
		
		self.RecSeqsTreesWidgetLayout.addWidget(self.RecSeqsTreesScollArea)
		self.RecSeqsTreesWidgetLayout.addWidget(self.RecInteractiveTrees)
		self.RecSeqsTreesWidget.setLayout(self.RecSeqsTreesWidgetLayout)

		self.mainLayout.addWidget(self.listViewBox,1,0) # Row, column
		self.mainLayout.addWidget(self.RecSeqsTreesWidget,2,0,4,1)
		self.mainLayout.addWidget(self.blastResLabelBox,0,1)
		self.mainLayout.addWidget(self.blastResTable,1,1)
		self.mainLayout.addWidget(self.geneIDBrowser,2,1)
		self.mainLayout.addWidget(self.snpResLabelBox,3,1)
		self.mainLayout.addWidget(self.SnpResTable,4,1)
		self.mainLayout.addWidget(self.snpBrowser,5,1)
		
		self.setLayout(self.mainLayout)

		self.initView()

		self.show()
	
	def initView(self):
		self.listView.setCurrentIndex(self.firstRecSeqItemIdx)
		print(self.listViewListofStringModel.data(self.firstRecSeqItemIdx, 0))
		selectionModel = QItemSelectionModel()
		selectionModel.setSourceModel(self.model)
		

	def displaySeqBlastResData(self, excel_path):
		self.recseqIdx = self.listView.currentIndex()
		self.recseq = self.listViewListofStringModel.data(self.recseqIdx, 0)
		self.blastResTable.clear()
		if dfBlast.size == 0:
			return
		dfBlast.fillna('', inplace=True)
		tmpdf = dfBlast[dfBlast["qaccver"] == self.recseq]
		tmpdf["pident"] = round(tmpdf["pident"] * 100,3)
		tmpdf = tmpdf.sort_values("qstart")
		self.blastResTable.setRowCount(tmpdf.shape[0])
		self.blastResTable.setColumnCount(tmpdf.shape[1])
		self.blastResTable.setHorizontalHeaderLabels(tmpdf.columns)
		
		# returns pandas array object
		row_num = 0
		for row in tmpdf.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.blastResTable.setItem(row_num, col_index, tableItem)
			row_num += 1
		
		self.blastResTable.resizeColumnsToContents()


	def displayGeneIDGraph(self):
		self.geneIDBrowser = QWebEngineView()
		tmpdf = recgID_for_graphdf[recgID_for_graphdf["Sequences"] == self.recseq]
		tmpdf = tmpdf.sort_values("qstart")
		fig = px.scatter(tmpdf, x="Gene", y="Sequences", size="pident", color="Lineage",text="Database",
		hover_data=["Database","qstart","qend","pident"], title="Gene identification recombinant sequences", 
		color_discrete_map = self.geneIDcolor_dict
		)
		self.geneIDBrowser.setHtml(fig.to_html(include_plotlyjs='cdn'))
		self.geneIDBrowser.setFixedHeight(300)
		self.mainLayout.addWidget(self.geneIDBrowser,2,1)
		self.setLayout(self.mainLayout)

	def displaySeqSNPResData(self, excel_path):
		self.SnpResTable.clear()
		if dfSNP.size == 0:
			return
		dfSNP.fillna('', inplace=True)
		tmpdf = dfSNP[dfSNP["Index"] == self.recseq]
		self.SnpResTable.setRowCount(tmpdf.shape[0])
		self.SnpResTable.setColumnCount(tmpdf.shape[1])
		self.SnpResTable.setHorizontalHeaderLabels(tmpdf.columns)
		
		# returns pandas array object
		row_num = 0
		for row in tmpdf.iterrows():
			values = row[1]
			for col_index, value in enumerate(values):
				tableItem = QTableWidgetItem(str(value))
				self.SnpResTable.setItem(row_num, col_index, tableItem)
			row_num += 1
		
		self.SnpResTable.resizeColumnsToContents()
		self.SnpResTable.setFixedHeight(80)

	def displaySNPGraph(self):
		self.snpBrowser = QWebEngineView()
		tmpdf = recSNP_for_graphdf[recSNP_for_graphdf["Sequences"] == self.recseq]
		fig = px.scatter(tmpdf, x="SNP reference pos", y="Sequences", color="Lineage",
			hover_data=["SNP reference pos","Lineage","Nucleotide"], title="SNP identification recombinant sequences", 
			color_discrete_map=self.snp_color_dict
		)
		fig.update_layout(xaxis_type = 'linear')
		self.snpBrowser.setHtml(fig.to_html(include_plotlyjs='cdn'))
		self.snpBrowser.setFixedHeight(400)
		self.mainLayout.addWidget(self.snpBrowser,5,1)
		self.setLayout(self.mainLayout)


	def updateLED(self):
		self._blastLED.setValue(self.recombinationStatus[self.recseq]["Genes"])
		self._snpLED.setValue(self.recombinationStatus[self.recseq]["SNP"])


	def eteInteractive(self, id):
		"""
		Open tree interactive window with the press of a button
		"""
		color = "#FF0000"
		for button in self.RecTreesRenderButtonGroup.buttons():
			if button is self.RecTreesRenderButtonGroup.button(id):
				button_id = button.text()
				t = self.trees[button_id]
				# Color recombinant
				rec_nodestyle = NodeStyle()
				rec_nodestyle["bgcolor"] = color
				for leaf in t.traverse():
					if leaf.name == self.recseq:
						leaf.set_style(rec_nodestyle)
				# Add faces etc
				ts = TreeStyle()
				ts.title.add_face(TextFace(button_id),column=1)
				ts.show_branch_support = True
				t.ladderize(direction=1)
				return t.show(tree_style=ts)

if __name__ == '__main__':
	app = QApplication(sys.argv)
	tg = RecTab()
	sys.exit(app.exec_())