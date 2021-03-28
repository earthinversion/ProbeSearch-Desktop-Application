################################################################################
##
# Utpal Kumar
# @Earth Inversion
################################################################################

import sys
import os

from PyQt5.QtGui import QColor
from PyQt5 import uic
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import yaml
import sys
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from io import StringIO
# from support import getProbePattern
from xml.sax import saxutils as su
from Bio import SeqIO
import pandas as pd
import os
import yaml
import regex

homeDir = str(Path.home())  # home directory


# loation of the app cache
appDir = os.path.join(homeDir, ".genomSequencingApp")
os.makedirs(appDir, exist_ok=True)
userDataFileName = os.path.join(
    appDir, "userData.yaml")  # userdata is kept here
if os.path.exists(userDataFileName):
    with open(userDataFileName) as file:
        userData = yaml.load(file, Loader=yaml.FullLoader)
else:
    userData = {}
    userData['loadDir'] = homeDir

prob_primer_seqs_file = "probe_primer_seq.yml"
if os.path.exists(prob_primer_seqs_file):
    with open(prob_primer_seqs_file) as file:
        prob_primer_seqs_dict = yaml.load(file, Loader=yaml.FullLoader)
else:
    prob_primer_seqs_dict = {}
currentDir = str(Path().absolute())

# GLOBALS
counter = 0
jumper = 10


class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        self.ui = uic.loadUi("mainwindow.ui", self)

        self.userData = userData
        self.all_loaded_files = []
        self.IconSize = QSize(50, 50)
        self.statusBar.setStyleSheet("color: red")
        if "loadedFiles" in self.userData and self.userData['loadedFiles']:
            item2Add = []
            for loadFile in self.userData['loadedFiles']:
                if os.path.exists(loadFile):
                    item2Add.append(loadFile)
                    self.all_loaded_files.append(loadFile)
                    textToAdd = self.getBaseFile()
                    self.label.setText(textToAdd)
            # add files to the list of loaded files
            self.listWidget.addItems(item2Add)
        else:
            self.label.setText('No files loaded')

        self.exec_cmd_lst = ['Quick Analyze']
        for key in prob_primer_seqs_dict.keys():
            self.exec_cmd_lst.append("Find "+key)

        # add loaded files to the list widget
        self.listWidget.setWordWrap(True)

        #
        self.plainTextEdit.setReadOnly(True)

        self.progressBar.setStyle(QStyleFactory.create("Windows"))
        self.progressBar.setTextVisible(True)
        self.lastruncmd = None
        self.label_last_run.setText(f"{self.lastruncmd}")
        # self.progressBar.hide()

        # add loaded files to the list widget
        self.listWidget_3.addItems(self.exec_cmd_lst)
        self.pxm = QIcon(os.path.join(
            'icons', 'arrow-alt-circle-right.svg')).pixmap(self.IconSize)
        self.label_5.setPixmap(self.pxm)

        self.listWidget.itemClicked.connect(self.load_file_sel_chn)
        self.listWidget_3.itemClicked.connect(self.exec_cmd_sel_clk)

        self.pushButton.clicked.connect(self.openFileDialog)
        self.actionLoad_Files.triggered.connect(self.openFileDialog)
        self.actionSave_Defaults.triggered.connect(self.saveDefaults)

        # load toolbar
        loadAct = QAction(
            QIcon(os.path.join('icons', 'plus.svg')), 'Load', self)
        loadAct.triggered.connect(self.openFileDialog)

        self.toolbar = self.addToolBar('Load')
        self.toolbar.addAction(loadAct)

        # save defaults toolbar
        saveDefltsAct = QAction(
            QIcon(os.path.join('icons', 'bookmark.svg')), 'Save Defaults', self)
        saveDefltsAct.triggered.connect(self.saveDefaults)

        self.toolbar = self.addToolBar('Save Defaults')
        self.toolbar.addAction(saveDefltsAct)

        self.pushButton_exec_3.clicked.connect(self.analyzeFasta)
        self.pushButton_dldCSV.clicked.connect(self.download_csv)
        self.pushButton_dldHTML.clicked.connect(self.download_html)
        self.pushButton_lw_clear.clicked.connect(self.clear_file_list)

        # EXECUTE commands
        self.selCmd = None
        self.selFile = None

    def clear_file_list(self):
        self.listWidget.clear()
        self.all_loaded_files = []

    def getBaseFileGen(self, loadedFile):
        loadedFile_basename = os.path.basename(loadedFile)
        loadedFile_dirname = os.path.split(os.path.dirname(
            loadedFile))[-1]

        return f".../{loadedFile_dirname}/{loadedFile_basename}"

    def load_file_sel_chn(self, item):
        self.label_4.setText(self.getBaseFileGen(item.text()))
        self.selFile = item.text()

    def exec_cmd_sel_clk(self, item):
        self.label_6.setText(item.text())
        self.selCmd = item.text()

    def openFileDialog(self):

        fname = QFileDialog.getOpenFileName(
            self, 'Open File', self.userData['loadDir'], 'Fasta files (*.fasta)')
        if len(fname[0]):
            loadedFile = fname[0]
            if loadedFile:
                # print("loaded file: ", loadedFile)
                if loadedFile not in self.all_loaded_files:
                    if os.path.exists(loadedFile):
                        self.all_loaded_files.append(loadedFile)
                        self.listWidget.addItem(loadedFile)
                        self.userData['loadedFiles'] = self.all_loaded_files
                        self.userData['loadDir'] = os.path.dirname(loadedFile)

            textToAdd = self.getBaseFile()
            self.label.setText(textToAdd)

    def getBaseFile(self):
        totalFiles = len(self.userData['loadedFiles'])
        loadedFile = self.userData['loadedFiles'][0]

        loadedFile_basename = os.path.basename(loadedFile)
        loadedFile_dirname = os.path.split(os.path.dirname(
            loadedFile))[-1]
        return f"{loadedFile_dirname} --> {loadedFile_basename} (total: {totalFiles} files)"

    def saveDefaults(self):
        with open(userDataFileName, 'w') as file:
            yaml.dump(self.userData, file)
        self.statusBar.setStyleSheet("color: green")
        self.statusBar.showMessage("Defaults saved!", 5000)

    def analyzeFasta(self):
        if self.selCmd and self.selFile:

            self.lengthSeqs = sum(
                1 for _ in SeqIO.parse(self.selFile, 'fasta'))

            if self.selCmd == self.exec_cmd_lst[0]:
                # self.progressBar.show()
                self.progressBar.setValue(0)
                self.plainTextEdit.clear()
                try:
                    self.lastruncmd = self.selCmd
                    self.sequences = SeqIO.parse(self.selFile, 'fasta')
                    self.worker0 = WorkerThread(
                        self.sequences, self.lengthSeqs)
                    self.worker0.start()
                    self.worker0.finished.connect(self.evt_worker_finished)
                    self.worker0.update_progress.connect(
                        self.evt_update_progress)
                    self.worker0.update_plainTextEdit.connect(
                        self.evt_update_plaintextedit)

                except Exception as err:
                    self.plainTextEdit.insertPlainText(err)
                # self.progressBar.hide()
            else:
                # self.progressBar.show()
                try:
                    allowedMismatch = int(self.spinBox_nmismatch.value())
                except:
                    allowedMismatch = 5
                    self.statusBar.showMessage(
                        "Number of allowed mismatch set to 5", 5000)
                self.totalSequences = None
                self.mismatchCountDict = None
                self.progressBar.setValue(0)
                self.plainTextEdit.clear()
                try:
                    self.sequences = SeqIO.parse(self.selFile, 'fasta')
                    probeName = self.selCmd.split()[1]
                    self.lastruncmd = self.selCmd

                    self.worker1 = WorkerThread_find(
                        self.sequences, self.lengthSeqs, probeName, allowedMismatch=allowedMismatch+1)
                    self.worker1.start()
                    self.worker1.finished.connect(self.evt_worker_finished)
                    self.worker1.update_progress.connect(
                        self.evt_update_progress)
                    self.worker1.update_plainTextEdit.connect(
                        self.evt_update_plaintextedit)
                    self.worker1.update_mismatch.connect(
                        self.evt_update_mismatch)
                    self.worker1.update_total_records.connect(
                        self.evt_update_total_rec)

                except Exception as err:
                    self.plainTextEdit.insertPlainText(err)
        else:
            self.statusBar.setStyleSheet("color: red")
            if not self.selCmd:
                self.statusBar.showMessage("Select a command to execute", 5000)
            elif not self.selFile:
                self.statusBar.showMessage("Select a file", 5000)

    def evt_worker_finished(self):
        # print('Worker finished')
        self.progressBar.setValue(int(100))
        self.statusBar.setStyleSheet("color: green")
        self.label_last_run.setText(f"{self.lastruncmd}")
        self.statusBar.showMessage("Finished computation", 10000)

    def evt_update_progress(self, val):
        self.progressBar.setValue(val)

    def evt_update_mismatch(self, dictval):
        self.mismatchCountDict = dictval

    def evt_update_total_rec(self, val):
        self.totalSequences = val

    def evt_update_plaintextedit(self, val):
        self.plainTextEdit.appendPlainText(
            val)
        # self.plainTextEdit.insertPlainText(
        #     val)

    def download_csv(self):
        app_csv_res = os.path.join(appDir, f'analysis_results.csv')
        try:
            textinbox = self.plainTextEdit.toPlainText()
            if not textinbox:
                self.statusBar.setStyleSheet("color: red")
                self.statusBar.showMessage("Nothing to download", 5000)
                return
            if os.path.exists(app_csv_res):
                text = pd.read_csv(app_csv_res)
            else:
                text = self.plainTextEdit.toPlainText()
            cmd_to_run = "_".join(self.selCmd.split()).lower()
            if len(text):
                name = QFileDialog.getSaveFileName(
                    self, 'Save CSV File', os.path.join(os.path.dirname(self.selFile), f'analysis_results_{cmd_to_run}.csv'), "CSV Files (*.csv);;Text Files (*.txt);;All Files (*)")
                fileName = name[0]
                if fileName:
                    if os.path.exists(app_csv_res):
                        text.to_csv(fileName, index=False)
                    else:
                        file = open(fileName, 'w')
                        file.write(text)
                        file.close()
                    self.statusBar.setStyleSheet("color: green")
                    self.statusBar.showMessage(
                        "File downloaded successfully", 5000)
            else:
                self.statusBar.setStyleSheet("color: red")
                self.statusBar.showMessage("Nothing to download", 5000)
        except Exception as err:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error")
            msg.setInformativeText(str(err))
            msg.setWindowTitle("Invalid Input")
            msg.exec_()

    def download_html(self):
        try:
            app_csv_res = os.path.join(appDir, f'analysis_results.csv')
            textinbox = self.plainTextEdit.toPlainText()
            if not textinbox:
                self.statusBar.setStyleSheet("color: red")
                self.statusBar.showMessage("Nothing to download", 5000)
                return
            if os.path.exists(app_csv_res):
                df_HTML = pd.read_csv(app_csv_res)
            else:
                text = self.plainTextEdit.toPlainText()
                df_HTML = pd.read_csv(StringIO(text), error_bad_lines=False)
            cmd_to_run = "_".join(self.selCmd.split()).lower()

            html = df_HTML.to_html()

            if len(df_HTML):
                name = QFileDialog.getSaveFileName(
                    self, 'Save HTML File', os.path.join(os.path.dirname(self.selFile), f'analysis_results_{cmd_to_run}.html'), "HTML Files (*.html);;Text Files (*.txt);;All Files (*)")
                fileName = name[0]
                if fileName:
                    file = open(fileName, 'w')
                    file.write(html)
                    file.write("\n")
                    file.write(
                        f"Total Number of sequences: {self.totalSequences}<br>")
                    spacing = "&nbsp;"*4
                    for key, val in self.mismatchCountDict.items():
                        file.write(
                            f"{spacing}Number of sequences with {key} mismatch: {val}<br>")
                    file.close()
                    self.statusBar.setStyleSheet("color: green")
                    self.statusBar.showMessage(
                        "File downloaded successfully", 5000)
            else:
                self.statusBar.setStyleSheet("color: red")
                self.statusBar.showMessage("Nothing to download", 5000)
        except Exception as err:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error")
            msg.setInformativeText(str(err))
            msg.setWindowTitle("Invalid Input")
            msg.exec_()


#########################################################################################################################################################


class WorkerThread_find(QThread):
    update_progress = pyqtSignal(int)
    update_plainTextEdit = pyqtSignal(str)
    update_mismatch = pyqtSignal(dict)
    update_total_records = pyqtSignal(int)

    def __init__(self, sequences, lengthSeqs, probeName, allowedMismatch=5):
        super().__init__()
        self.sequences = sequences
        self.lengthSeqs = lengthSeqs
        self.probeName = probeName
        # get the sequence from the yml file
        seq_to_search = prob_primer_seqs_dict[self.probeName]
        self.seq_to_search = seq_to_search
        if allowedMismatch < len(self.seq_to_search):
            self.allowedMismatch = allowedMismatch
        else:
            self.allowedMismatch = len(self.seq_to_search)
        self.records = list(self.sequences)
        self.misMatchCounts = {}

        self.update_progress.emit(0)

    def run(self):
        global appDir
        self.textToSend = f"genome,seq,Mismatch,misMatchLocs\n"

        self.update_plainTextEdit.emit(self.textToSend)

        results_list = []
        count = 0
        for recIdx, record in enumerate(self.records):
            count += 1
            recseq = str(record.seq)
            recid = record.id
            lenrecSeq = len(recseq)
            ratioN = recseq.count('N')/lenrecSeq * 100
            if ratioN > 1:  # ignore the record if the percent of N is greater than 1
                # print(f"N ratio is high! {ratioN}%")
                continue

            match = regex.search(
                f"({self.seq_to_search})"+"{s<"+f"{self.allowedMismatch}"+"}", recseq)

            if match:
                match_span = f"{match.span()[0]+1}-{match.span()[1]}"
                match_group = match.group()
                if "N" in match_group:
                    continue

                dict_ = {}
                dict_.update({"genome": record.id+f"/{match_span}"})
                dict_.update({"seq": match_group})
                # levenshtein_distance_val = levenshtein_distance(
                #     self.seq_to_search, match_group)
                mistMatchesLocs = self.get_mismatch_idx(
                    self.seq_to_search, match_group, numMismatch=match.fuzzy_counts[0])
                dict_.update({"Mismatch": match.fuzzy_counts[0]})
                dict_.update({"MismatchLocation": mistMatchesLocs})
                results_list.append(dict_)
                self.textToSend = f"{record.id}/{match_span},{match_group},{match.fuzzy_counts[0]},{mistMatchesLocs}"
                self.update_plainTextEdit.emit(self.textToSend)

            else:
                continue
            self.update_progress.emit(
                int(round(recIdx/self.lengthSeqs*100)))
        self.totalSequences = count  # total number of sequences
        self.update_plainTextEdit.emit(
            "---")
        self.update_plainTextEdit.emit(
            f"Total number of records analyzed: {self.totalSequences}")

        dff = pd.DataFrame(results_list)
        try:

            for ii in range(self.allowedMismatch):
                totmismatch = (dff['Mismatch'] == ii).sum()
                self.misMatchCounts[f'{ii}'] = totmismatch
                self.update_plainTextEdit.emit(
                    f"Number of sequences with {ii} mismatch: {totmismatch}")
            self.update_mismatch.emit(self.misMatchCounts)
            self.update_total_records.emit(self.totalSequences)
        except:
            pass
        filename = os.path.join(appDir, f'analysis_results.csv')
        dff.to_csv(filename, index=False)

    # def getPattern(self):
    #     # Create pattern to search
    #     pattern = getProbePattern(self.probeName, self.seq_to_search)
    #     return pattern

    def get_mismatch_idx(self, w1, w2, numMismatch=999):
        idxs = []
        for (i, (a, b)) in enumerate(zip(w1, w2)):
            if a != b:
                idxs.append(i)
            if len(idxs) >= numMismatch:
                break
        # print(idxs)
        idxs = [str(int(idval+1)) for idval in idxs]
        idxs = ";".join(idxs)
        return idxs


class WorkerThread(QThread):
    update_progress = pyqtSignal(int)
    update_plainTextEdit = pyqtSignal(str)

    def __init__(self, sequences, lengthSeqs):
        super().__init__()
        self.sequences = sequences
        self.lengthSeqs = lengthSeqs
        # self.dict_list = []

    def run(self):
        textToSend = f"recordID,sequence,length of seq,percent N,Total P"

        self.update_plainTextEdit.emit(textToSend)
        dict_list = []
        for idx, record in enumerate(self.sequences):
            recseq = str(record.seq)
            lenrecSeq = len(recseq)
            ratioN = recseq.count('N')/lenrecSeq * 100
            numP = recseq.count('P')

            textToSend = f"{record.id},{recseq[0:10]}...{recseq[-10:-1]},{lenrecSeq},{ratioN:.1f},{numP}"

            self.update_plainTextEdit.emit(textToSend)
            self.update_progress.emit(int(round(idx/self.lengthSeqs*100)))
            dict_ = {}
            dict_.update({"genome": record.id})
            dict_.update({"recSeq": f"{recseq[0:10]}...{recseq[-10:-1]}"})
            dict_.update({"lenSeq": lenrecSeq})
            dict_.update({"N ratio": f"{ratioN:.1f}"})
            dict_.update({"Num Ps": numP})
            dict_list.append(dict_)
        dff = pd.DataFrame(dict_list)
        filename = os.path.join(appDir, f'analysis_results.csv')
        dff.to_csv(filename, index=False)

# ==> SPLASHSCREEN WINDOW


class SplashScreen(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        # self.ui = Ui_SplashScreen()
        # self.ui.setupUi(self)
        self.ui = uic.loadUi("splash_screen.ui", self)
        # ==> SET INITIAL PROGRESS BAR TO (0) ZERO
        self.progressBarValue(0)

        # ==> REMOVE STANDARD TITLE BAR
        self.setWindowFlags(Qt.FramelessWindowHint)  # Remove title bar
        # Set background to transparent
        self.setAttribute(Qt.WA_TranslucentBackground)

        # ==> APPLY DROP SHADOW EFFECT
        self.shadow = QGraphicsDropShadowEffect(self)
        self.shadow.setBlurRadius(20)
        self.shadow.setXOffset(0)
        self.shadow.setYOffset(0)
        self.shadow.setColor(QColor(0, 0, 0, 120))
        self.ui.circularBg.setGraphicsEffect(self.shadow)

        # QTIMER ==> START
        self.timer = QTimer()
        self.timer.timeout.connect(self.progress)
        # TIMER IN MILLISECONDS
        self.timer.start(15)

        # SHOW ==> MAIN WINDOW
        ########################################################################
        self.show()
        ## ==> END ##

    # DEF TO LOANDING
    ########################################################################
    def progress(self):
        global counter
        global jumper
        value = counter

        # HTML TEXT PERCENTAGE
        htmlText = """<p><span style=" font-size:68pt;">{VALUE}</span><span style=" font-size:58pt; vertical-align:super;">%</span></p>"""

        # REPLACE VALUE
        newHtml = htmlText.replace("{VALUE}", str(jumper))

        if(value > jumper):
            # APPLY NEW PERCENTAGE TEXT
            self.ui.labelPercentage.setText(newHtml)
            jumper += 10

        # SET VALUE TO PROGRESS BAR
        # fix max value error if > than 100
        if value >= 100:
            value = 1.000
        self.progressBarValue(value)
        if counter == 10:
            self.main = MainWindow()

        # CLOSE SPLASH SCREE AND OPEN APP
        if counter > 100:
            # STOP TIMER
            self.timer.stop()

            # SHOW MAIN WINDOW
            # self.main = MainWindow()
            self.main.show()

            # CLOSE SPLASH SCREEN
            self.close()

        # INCREASE COUNTER
        counter += 1.0

    # DEF PROGRESS BAR VALUE
    ########################################################################
    def progressBarValue(self, value):

        # PROGRESSBAR STYLESHEET BASE
        styleSheet = """
        QFrame{
        	border-radius: 150px;
        	background-color: qconicalgradient(cx:0.5, cy:0.5, angle:90, stop:{STOP_1} rgba(255, 0, 127, 0), stop:{STOP_2} rgba(85, 170, 255, 255));
        }
        """

        # GET PROGRESS BAR VALUE, CONVERT TO FLOAT AND INVERT VALUES
        # stop works of 1.000 to 0.000
        progress = (100 - value) / 100.0

        # GET NEW VALUES
        stop_1 = str(progress - 0.001)
        stop_2 = str(progress)

        # SET VALUES TO NEW STYLESHEET
        newStylesheet = styleSheet.replace(
            "{STOP_1}", stop_1).replace("{STOP_2}", stop_2)

        # APPLY STYLESHEET WITH NEW VALUES
        self.ui.circularProgress.setStyleSheet(newStylesheet)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SplashScreen()
    sys.exit(app.exec_())
