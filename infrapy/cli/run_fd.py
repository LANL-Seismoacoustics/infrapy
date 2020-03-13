#!/usr/bin/env python
import sys

from infrapy.database.taskbase.fdet import FDet


def run(config_file, array):
    pdetect = FDet(config_file, ARRAY_NAME=array)
    pdetect.database_connecting()
    pdetect.data_processing()

if __name__ == '__main__':

    try:
        config_file = sys.argv[1]
    except Exception as ex1:
        print('A configuration file is required to start data processing')
        sys.exit()
    print('Running fk with configuration file:', config_file)

    try:
        array = sys.argv[2]
        print('Array override:', array)
    except Exception as ex1:
        array = None

    run(config_file, array)
