#!/usr/bin/env python
import sys

from infrapy.database.taskbase.assoc import AssocInfraPy_LANL


def run(config_file):
    pdetect = AssocInfraPy_LANL(config_file)
    pdetect.database_connecting()
    pdetect.data_processing()


if __name__ == '__main__':

    try:
        config_file = sys.argv[1]
    except Exception as ex1:
        print('A configuration file is required to start data processing')
        sys.exit()

    print('run FK, with configuration file:', config_file)
    run(config_file)
