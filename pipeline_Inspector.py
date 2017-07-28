#!/usr/bin/python3
import os, re, sys, stat
import json, pickle
import logging
import logging.handlers
import subprocess as sp

import pandas as pd
from GS_variants import GS_variant_effect_prediction
from glob import glob
from operator import methodcaller


class Logger(object):

    def __init__(self):
        self.Save_file_path = '/home/jsgene/JK1/NGS/exec_history/'  # log files that indicate execution time
        self.JSON_path = '/data1/home/jsgene/JSONS/'
        self.PICKLE_path = '/data1/home/jsgene/PICKLES'
        self.Desktop = '/home/jsgene/JK1/NGS/'

    def GS_TIME(self):
        now = time.localtime()
        begin_time = '%04d-%02d-%02d %02d:%02d:%02d' % (now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
        return begin_time

    def GS_exec_log(self, filename, saving_data):
        # Save history of operation
        with open('%s%s' % (self.Save_file_path, filename), 'a') as f:
            f.write(saving_data)

    def GS_progress_log(self, log_filename, saving_data):

        normal_logger = logging.getLogger()
        normal_formatter = logging.Formatter('[%(process)d %(processName)s %(thread)d %(threadName)s|%(filename)s:%(lineno)s] %(asctime)s > %(message)s')

        # create file handler which logs even debug messages
        normal_fileHandler = logging.FileHandler('%s%s.log' % (self.Save_file_path, log_filename))
        normal_streamHandler = logging.StreamHandler()

        # create formatter and add it to the handlers
        normal_streamHandler.setFormatter(normal_formatter)
        normal_fileHandler.setFormatter(normal_formatter)

        normal_logger.addHandler(normal_fileHandler)
        normal_logger.addHandler(normal_streamHandler)
        normal_logger.setLevel(logging.DEBUG)
        normal_logger.debug(saving_data)


    def GS_error_log(self, log_filename, function_name, message, exc_info=True):
        try:
            # execution of the function
            function_name
        except StandardError as e:
            stderr_logger = logging.getLogger()
            stderr_formatter = logging.Formatter('[%(process)d %(processName)s %(thread)d %(threadName)s|%(filename)s:%(lineno)s] %(asctime)s > %(message)s')

            stderr_fileHandler = logging.FileHandler('%s%s.log' % (self.Save_file_path, log_filename))
            stderr_streamHandler = logging.StreamHandler()

            # create formatter and add it to the handlers
            stderr_streamHandler.setFormatter(stderr_formatter)
            stderr_fileHandler.setFormatter(stderr_formatter)

            stderr_logger.addHandler(stderr_fileHandler)
            stderr_logger.addHandler(stderr_streamHandler)

            stderr_logger.setLevel(logging.ERROR)
            stderr_logger.error(e, exc_info)
            stderr_logger.error(message, exc_info)
        except EnvironmentError as err:
            stderr_logger = logging.getLogger()
            stderr_formatter = logging.Formatter('[%(process)d %(processName)s %(thread)d %(threadName)s|%(filename)s:%(lineno)s] %(asctime)s > %(message)s')
            stderr_fileHandler = logging.FileHandler('%s%s.log' % (self.Save_file_path, log_filename))
            stderr_streamHandler = logging.StreamHandler()

            # create formatter and add it to the handlers
            stderr_streamHandler.setFormatter(stderr_formatter)
            stderr_fileHandler.setFormatter(stderr_formatter)

            stderr_logger.addHandler(stderr_fileHandler)
            stderr_logger.addHandler(stderr_streamHandler)

            stderr_logger.setLevel(logging.ERROR)
            stderr_logger.error(err, exc_info)
            stderr_logger.error(message, exc_info)

class Inspector(object):
    def __init__(self, group, type="GS"):

        if type == "GS":
            self.groupName = "GS_" + group
            self.basedir = '/EQL7/pipeline/'
            self.bam_dir = '%s%s/' % (self.basedir, self.groupName)
            self.fastq_link = '/EQL2/%s/WXS/fastq/link/' %(self.groupName)
        elif type == "RSq":
            self.groupName = "SGI" + group + "_rsq2"
            self.basedir = '/EQL8/pipeline/'
            self.tdf_dir = '%s%sexpr/' % (self.basedir, self.groupName)
            self.bam_dir = '%s%smut/' % (self.basedir, self.groupName)
            self.fastq_link = '/EQL2/SGI_%s/RNASeq/fastq/link/' %(group)

        self.Pair = '_pair_filter_vep.dat'
        self.Single = '_single_filter_vep.dat'
        self.type = type

    def file_existance(self, fileLists, Dir):

        List_Done = []
        List_Not_Done = []

        for aa in range(len(fileLists)):
            #ReferenceDirs = glob(Dir+"*"+fileLists[aa]+"*" + self.Pair)
            ReferenceDirs = glob(Dir + "*" + fileLists[aa] + "*")
            if ReferenceDirs == []:
                if not fileLists[aa].split('_')[-1] == "B":
                    List_Not_Done.append(fileLists[aa])
            else:
                print(ReferenceDirs)
                List_Done.append(fileLists[aa])

        print('\n\n\n\n\n', List_Done, " exists!!!")
        print( '\n\n\n\n\n', List_Not_Done, " not exists")

    def inspection_for_Copynumber(self):
        CN = GS_Copynumber()
        CN.initialize_setting(self.groupName)
        CN.Select_TODO_list()

    def Sample_naming_inspection(self, matchingList):

        standards = []
        needToRename = []

        if self.type == "GS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})(_T)_(GS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(GS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+")
            standards = [standard1, standard2, standard3, standard4, standard5, standard6]
        elif self.type == "RSq":
            standard1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9]|_T0[1-9]_P0)_(RSq)+")
            standards = [standard1]
        elif self.type == "WXS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})_T_(SS|TS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(SS|TS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(SS|TS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(SS|TS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(SS|TS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(SS|TS)+")
            standards = [standard1, standard2, standard3, standard4, standard5, standard6]

        for ref in standards:
            needToRename = list(filter(lambda x: ref.match(x) == None, matchingList)) + needToRename

        return needToRename


    def Renaming_Files(self):
        sampleNames = glob("%s/*" % (self.bam_dir))
        unmatchedList = []

        ## At first, upper dir names are changed
        for sampleGroup in sampleNames:
            if sampleGroup.split(".")[-1] == "html":
                ######## Html #####################
                FirstFilter = sampleGroup.split("/")  # e.g : ['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1021_02_RSq.html']
                SecondFilter = FirstFilter[-1]  # e.g : ['IRCR_BT16_1021_02_RSq.html']
                HTMLFilter = self.Sample_naming_inspection([SecondFilter])  # e.g :
                HTMLFilter = list(set(HTMLFilter))
                CorrectedHTML = self.correcting_naming_mistakes(HTMLFilter)

                ####### Link ######################
                if len(HTMLFilter) == 1:
                    OLD = "/".join(FirstFilter[:-1] + HTMLFilter)
                    NEW = "/".join(FirstFilter[:-1] + CorrectedHTML)

                    DIRFilter = HTMLFilter[0][:-5]
                    print(DIRFilter)
                    CorrectedName = CorrectedHTML[0][:-5]
                    #os.path.join(os.getcwd(), self.fastq_link)
                    OldLink1 = '%s%s.1.fq.gz' % (self.fastq_link, DIRFilter)
                    OldLink2 = '%s%s.2.fq.gz' % (self.fastq_link, DIRFilter)
                    NewLink1 = '%s%s.1.fq.gz' % (self.fastq_link, CorrectedName)
                    NewLink2 = '%s%s.2.fq.gz' % (self.fastq_link, CorrectedName)
                    print("old %s vs new %s" % (OldLink1, NewLink1))

                    self.Rename(OLD, NEW)
                    self.Rename(OldLink1, NewLink1)
                    self.Rename(OldLink2, NewLink2)
            else:
                continue

        ## Then, sub dir names and file names are changed
        for sampleGroup in sampleNames:

            matchingFiles = glob("%s/*" % (sampleGroup))  # e.g : ['/EQL8/pipeline/SGI20170718/IRCR_BT16_1021_02_RSq/IRCR_BT16_1021_02_RSq_splice.bam', '/EQL8/pipeline/SGI20170718/IRCR_BT16_1141_RSq/IRCR_BT16_1141_RSq_splice.dedup.bai']

            FirstFilter = list(map(methodcaller("split", "/"), matchingFiles)) # e.g : [['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1021_02_RSq', 'IRCR_BT16_1021_02_RSq_splice.bam'], ['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1141_RSq', 'IRCR_BT16_1141_RSq_splice.dedup.bai']]
            SecondFilter = list(map(lambda x: x[-1], FirstFilter)) # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam', 'IRCR_BT16_1141_RSq_splice.dedup.bai']
            ThirdFilter = self.Sample_naming_inspection(SecondFilter) # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam']
            ThirdFilter = list(set(ThirdFilter))
            CorretedFileNameList = self.correcting_naming_mistakes(ThirdFilter)
            OLDFILES = list(map(lambda x: sampleGroup + '/' + x, ThirdFilter))
            NEWFILES = list(map(lambda x: sampleGroup + '/' + x, CorretedFileNameList))

            if len(OLDFILES) == len(NEWFILES):
                for i in range(len(OLDFILES)):
                    self.Rename(OLDFILES[i], NEWFILES[i])

            FourthFilter = sampleGroup.split("/")
            FifthFilter = self.Sample_naming_inspection(FourthFilter[-1:])
            NEWDIRList = self.correcting_naming_mistakes(FifthFilter)
            FourthFilter[-1] = NEWDIRList[-1] 
            self.Rename(sampleGroup, '/'.join(FourthFilter))



    def correcting_naming_mistakes(self, inputList):
        outputList = []
        mistake1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_0[1-9]|_0[1-9]_P0)_([A-Z]{2,3})+")
        mistake2 = re.compile("^S([0-9]{2})([0-9]{5,})_(GS)+")
        mistake3 = re.compile("^(|IRCR_)S([0-9]{2})_([0-9]{5,})_(GS)+")

        for NAME in inputList:
            if mistake1.match(NAME) != None:
                newNAME = NAME.replace("_0", "_T0")
            elif mistake2.match(NAME) != None:
                newNAME = NAME.replace("_GS", "_T_GS")
            elif mistake3.match(NAME) != None:
                if "IRCR_S" in NAME:
                    newNAME = NAME.replace("IRCR_S", "S")

                Delimiter = "_GS"
                newNAMEList = newNAME.split(Delimiter)
                newNAMEList[0] = ''.join(newNAMEList[0].split("_"))
                newNAME = Delimiter.join(newNAMEList)
                if mistake2.match(newNAME) != None:
                    newNAME = newNAME.replace("_GS", "_T_GS")
            else:
                newNAME = NAME
            outputList.append(newNAME)
        return outputList

    def ReLink(self,toUnlink, newLink):
        STDOUT = sp.check_output("ls -l %s" %(toUnlink), shell=True)
        TargetPoint = str(STDOUT).split("->")[-1][1:-3]
        sp.check_call("unlink %s" %(toUnlink))
        sp.check_call("ln -s %s %s"%(TargetPoint, newLink), shell=True)


    def Rename(self, A, B):
        os.system('mv %s %s' % (A, B))
        print("%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n %s " %(A, B))



    def read_DAT(self, isPAIR=True):
        VEP = GS_variant_effect_prediction()
        VEP.initialize_setting(self.groupName)
        VEP.initialize_Single_Pair(isPair=isPAIR)
        LIST = VEP.Select_TODO_list()
        INDEX = ['Sample_Name', 'CHR', 'POSITION', 'REF', 'ALT', 'n_nRef', 'n_nAlt', 't_nRef', 't_nAlt', 'Annot_type', 'canonical', 'SC']
        instant_msg = []

        for dirs in LIST:
            ID_mutect = dirs.split('/')[-1] + '.mutect'
            ID_indels = dirs.split('/')[-1] + '.indels'
            for ID in [ID_mutect, ID_indels]:
                if isPAIR:
                    file = dirs + '/' + ID + self.Pair
                else:
                    file = dirs + '/' + ID + self.Single
                try:
                    with open(file, 'r') as f:
                        for line in f:
                            (SName, CHR, POS, REF, ALT, n_nRef, n_nAlt, t_nRef, t_nAlt, AnnotGene, AnnotCH_dna,
                             AnnotCH_prot, Annot_type, canonical, refseq, SC) = line.split('\t')
                            Sample = [SName, CHR, POS, REF, ALT, n_nRef, n_nAlt, t_nRef, t_nAlt, Annot_type, canonical, SC]
                            Sample_info = pd.Series(data=Sample, index=INDEX)
                            for element in Sample:
                                if element == '':
                                    print(line, Sample_info)
                                    instant_msg.append(file)
                                    with open("/home/jsgene/JK1/NGS/exec_history/%s_missing_in_DAT.txt" % ID, 'a') as wf:
                                        wf.write(line)
                                else:
                                    print(dirs, 'is OK')
                                    continue
                except IOError as err:
                    logging.basicConfig(filename="/home/jsgene/JK1/NGS/exec_history/%s_no_DAT_file.log" % ID, level=logging.DEBUG)
                    logging.debug(err)
                    continue
            print("something wrong with %s" %instant_msg)

#    def sequence_Quality(self):

class File_manager(object):
    import os, sys
    import subprocess
    def from_EQL_to_Storage(self, type):
        departureHost = "jsgene@119.86.100.106"
        departure = "/EQL8/pipeline/SGI20150709_rsq2expr/"
        destinationHost = "smcbi@119.4.212.148"
        destination = "/media/backup"
        if type == "GS":
            destination = destination + "/WXS_pipeline/"
        elif type == "RSq":
            destination = destination + "/RNASeq_pipeline/"

        cmd1 = 'sshpass -p smcbi '
        cmd2 = 'rsync -avXz --progress %s %s:%s' % (departure, destinationHost, destination)

        #cmd = cmd1 + cmd2
        cmd3 = "SSHPASS=smcbi"
        cmd4 = "sshpass -e sftp -oBatchMode=no -b - %s <<EOF" %(destinationHost)
        cmd5 = "put %s" %(departure)
        cmd6 = "bye"
        cmd7 = "EOF"
        print(cmd)
        os.system(cmd)

    def move_whole_dir(self, FROM, TO):
        sp.check_output(["cp" "-r", FROM, TO], stderr=sp.STDOUT, shell=True)






if __name__=="__main__":
    ins = Inspector("20170711", "GS")
    ins.Renaming_Files()
    #ins.read_DAT()
    #mgr = File_manager()
    #mgr.from_EQL_to_Storage("RSq")
