#!/usr/bin/evn python
#
# Helper code to read command line arguments, read a dictionary from
# a .yaml or .json file, and work with dictionaries.
# ----------------------------------------------------------------------
# Author: J. P. Johnston <jpj13@mit.edu>
#
# Some code is based on morpho
# morpho is located at https://github.com/project8/morpho/tree/v1.1.5
# The authors of morpho are:
#   J. A. Formaggio <josephf@mit.edu>
#   T. E. Weiss <tweiss@mit.edu>
#   M. Guigue <mathieu.guigue@pnnl.gov>
#   J. N. Kofron <jared.kofron@gmail.com>
# ----------------------------------------------------------------------

from yaml import load as yload
from argparse import ArgumentParser

import logging
logger = logging.getLogger('shapegenerator')
logger.setLevel(logging.DEBUG)
base_format = '%(asctime)s[%(levelname)-8s] %(name)s(%(lineno)d) -> %(message)s'
logging.basicConfig(format=base_format, datefmt='%m/%d/%Y %H:%M:%S')

import os

# ----------------------------------------------------------------------
# Methods to parse command line arguments and read data from files
# Code is from morpho.py
def parse_args():
    p = ArgumentParser(description='''
        Preprocessing of model parameters for a stan/morpho model
    ''')
    p.add_argument('-c','--config',
                   metavar='<configuration file>',
                   help='Full path to the configuration file',
                   required=True)
    p.add_argument('param',nargs='*',
                   default=False,
                   help='Manualy change of a parameter and its value')
    return p.parse_args()
# ----------------------------------------------------------------------
def read_param(yaml_data, node, default):
        data = yaml_data
        xpath = node.split('.')
        try:
            for path in xpath:
                data = data[path]
        except Exception as exc:
            if default == 'required':
                err = """FATAL: Configuration parameter {0} required but not\
                provided in config file!
                """.format(node)
                logger.debug(err)
                raise exc
            else:
                data = default
        return data
# ----------------------------------------------------------------------
def change_and_format(b):
    if b == 'True':
        return True
    elif b == 'False':
        return False
    else:
        try:
            a = float(b)
            return a
        except:
            return b
# ----------------------------------------------------------------------
def merge(a, b, path=None):
    '''
    merges b into a
    '''
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                a[key] = change_and_format( b[key])
        else:
            a[key] = change_and_format( b[key])
    return a
# ----------------------------------------------------------------------
def update_from_arguments(the_dict,args):
    '''
    Update the dictionary extracted from the config file
    '''
    logger.debug('Update dict parameters')
    new_dict = the_dict
    for a_arg in args:
        result = a_arg.split('=')
        print(result)
        xpath = result[0].split('.')
        print(xpath)
        to_update_dict = {xpath[-1]:result[1]}
        print(to_update_dict)
        for path in reversed(xpath[:-1]):
            print("path=%s"%path)
            to_update_dict = {path:to_update_dict}
            print(to_update_dict)
        new_dict = merge(new_dict,to_update_dict)
    return new_dict
# ----------------------------------------------------------------------
# Inputs: path- file path or directory path
#         file_path- True if the given path includes a filename at the
#                    end, false if the path is to a directory
#
# Creates the directory if it does not exist
def create_path(path,file_path=True):
    if(file_path):
        path = "/".join(path.split("/")[:-1])
    if not os.path.exists(path):
        os.makedirs(path)
    return
# ----------------------------------------------------------------------

