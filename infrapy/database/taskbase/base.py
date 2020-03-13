'''
Created on Oct 30, 2014

@author: omarcillo
'''
import os.path as os_path
import sys
import configparser


class Base(object):
    '''
    Base class for Infrapy
    '''
    general='GeneralParams'

    def __init__(self, conf_file,task):
        '''
        Constructor
        '''
        self.general_PARAM={}
        self.task_PARAM={}
        self.db_connected=False
        self.data_retrieved=False
        self.data_processed=False
        self.result_submitted=False
        if conf_file==[]:
            print('Error: configuration file is not defined')
            sys.exit(1)
        if os_path.isfile(conf_file) == False:
            print('Error: configuration file does not exist in directory')
            sys.exit(2)
        Config = configparser.ConfigParser()
        Config.read(conf_file)
        self.db_PARAM=dict(Config.items('database'))
        self.general_PARAM=dict(Config.items(self.general))
        self.task_PARAM=dict(Config.items(task))
        print('\n GENERAL PARAMETERS')
        for il in self.general_PARAM:
            print(il,self.general_PARAM[il])
        print('\n DB CONNECTION')
        for il in self.db_PARAM:
            print(il,self.db_PARAM[il])
        print('\n TASK PARAMETERS')
        for il in self.task_PARAM:
            print(il,self.task_PARAM[il])


    @classmethod
    def database_connecting(self):
        '''
        Constructor
        '''
        print('Hello')
        return self.db_connected

    @classmethod
    def data_retrieving(self):
        '''
        Constructor
        '''
        return self.data_retrieved

    @classmethod
    def data_processing(self):
        '''
        Constructor
        '''
        return self.data_processed

    @classmethod
    def results_submitting(self):
        '''
        Constructor
        '''
        return self.result_submitted

    @classmethod
    def database_disconnecting(self):
        '''
        Constructor
        '''
        return self.db_connected
