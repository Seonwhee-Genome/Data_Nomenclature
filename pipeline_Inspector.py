#!/usr/bin/python3
import os, re, sys, stat
import json, pickle
import logging
import logging.handlers
import subprocess as sp
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
        elif type == "RSq":
            self.groupName = "SGI" + group + "_rsq2"
            self.basedir = '/EQL8/pipeline/'
            self.tdf_dir = '%s%sexpr/' % (self.basedir, self.groupName)
            self.bam_dir = '%s%smut/' % (self.basedir, self.groupName)

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
        print('matchingTarget', matchingList)

        if self.type == "GS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})(|_T)_(GS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(GS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+")
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
        sampleNames = glob("%s/*" % (self.tdf_dir))
        unmatchedList = []
        for sampleGroup in sampleNames:
            if sampleGroup.split(".")[-1] == "html":
                continue
            else:
                matchingFiles = glob("%s/*" % (sampleGroup))  # e.g : ['/EQL8/pipeline/SGI20170718/IRCR_BT16_1021_02_RSq/IRCR_BT16_1021_02_RSq_splice.bam', '/EQL8/pipeline/SGI20170718/IRCR_BT16_1141_RSq/IRCR_BT16_1141_RSq_splice.dedup.bai']

                FirstFilter = list(map(methodcaller("split", "/"), matchingFiles)) # e.g : [['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1021_02_RSq', 'IRCR_BT16_1021_02_RSq_splice.bam'], ['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1141_RSq', 'IRCR_BT16_1141_RSq_splice.dedup.bai']]
                SecondFilter = list(map(lambda x: x[-1], FirstFilter)) # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam', 'IRCR_BT16_1141_RSq_splice.dedup.bai']
                ThirdFilter = self.Sample_naming_inspection(SecondFilter) # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam']

                for matchingFile in matchingFiles:
                    if True in list(map(lambda x: x in matchingFile, ThirdFilter)):
                        unmatchedList.append(matchingFile)


        unmatchedList = list(set(unmatchedList))
        print(unmatchedList)
        self.correcting_naming_mistakes(unmatchedList)

        #for NAME in unmatchedList:

            #print(NAME)
            # os.system('mv %s %s' %(path+NAME, path+new_name))
    def correcting_naming_mistakes(self, inputList):
        mistake1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_0[1-9]|_0[1-9]_P0)_([A-Z]{2,3})+")
        for PATH, NAME in inputList:
            if mistake1.match(NAME) != None:
                print("%s is mistake1" %(NAME))
                newNAME = NAME.replace("_0", "_T0")
                print("and now changed into %s" %(newNAME))
                os.system('mv %s %s' % (PATH + NAME, PATH + newNAME))
            else:
                print("%s is not mistake1" %(NAME))





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

if __name__=="__main__":
    ins = Inspector("20170718", "RSq")
    ins.Renaming_Files()