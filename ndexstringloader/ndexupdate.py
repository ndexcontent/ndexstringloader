#! /usr/bin/env python

import argparse
import sys
import logging
from logging import config
from ndexutil.config import NDExUtilConfig
import ndexstringloader

import csv
import pandas as pd
import json

from datetime import datetime

import gzip
import shutil

import os

import ndexutil.tsv.tsv2nicecx2 as t2n

from ndexutil.tsv.streamtsvloader import StreamTSVLoader



import configparser
import urllib

import requests

import ndex2


logger = logging.getLogger(__name__)

TSV2NICECXMODULE = 'ndexutil.tsv.tsv2nicecx2'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(relativeCreated)dms " \
             "%(filename)s::%(funcName)s():%(lineno)d %(message)s"

def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_fm)

    parser.add_argument('--profile', help='Profile in configuration '
                                          'file to use to load '
                                          'NDEx credentials which means'
                                          'configuration under [XXX] will be'
                                          'used '
                                          '(default '
                                          'ndexstringloader)',
                        default='ndexstringloader')
    parser.add_argument('--logconf', default=None,
                        help='Path to python logging configuration file in '
                             'this format: https://docs.python.org/3/library/logging.config.html#logging-config-fileformat'
                             'Setting this overrides -v parameter which uses '
                             ' default logger. (default None)')

    parser.add_argument('--conf', help='Configuration file to load '
                                       '(default ~/' +
                                       NDExUtilConfig.CONFIG_FILE)
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module and'
                             'in ' + TSV2NICECXMODULE + '. Messages are '
                             'output at these python logging levels '
                             '-v = ERROR, -vv = WARNING, -vvv = INFO, '
                             '-vvvv = DEBUG, -vvvvv = NOTSET (default no '
                             'logging)')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' +
                                 ndexstringloader.__version__))

    parser.add_argument('--networkid', help='UUID of network on the server to update', required=True)

    parser.add_argument('--cx', help='CX file to use for updating netwoprk on server', required=True)

    return parser.parse_args(args)


def _setup_logging(args):
    """
    Sets up logging based on parsed command line arguments.
    If args.logconf is set use that configuration otherwise look
    at args.verbose and set logging for this module and the one
    in ndexutil specified by TSV2NICECXMODULE constant
    :param args: parsed command line arguments from argparse
    :raises AttributeError: If args is None or args.logconf is None
    :return: None
    """

    if args.logconf is None:
        level = (50 - (10 * args.verbose))
        logging.basicConfig(format=LOG_FORMAT,
                            level=level)
        logging.getLogger(TSV2NICECXMODULE).setLevel(level)
        logger.setLevel(level)
        return

    # logconf was set use that file
    logging.config.fileConfig(args.logconf,
                              disable_existing_loggers=False)


class NDExNdexstringloaderLoader(object):
    """
    Class to load content
    """
    def __init__(self, args):
        """
        :param args:
        """
        self._conf_file = args.conf
        self._profile = args.profile

        self._network_id = args.networkid
        self._update_cx_file_name = args.cx

    def _parse_config(self):
        """
        Parses config
        :return:
        """
        ncon = NDExUtilConfig(conf_file=self._conf_file)
        con = ncon.get_config()
        self._user = con.get(self._profile, NDExUtilConfig.USER)
        self._pass = con.get(self._profile, NDExUtilConfig.PASSWORD)
        self._server = con.get(self._profile, NDExUtilConfig.SERVER)


    def _update_network_on_server(self):

        print('{} - updating network {} from cx file {} on server {} for user {}...'
              .format(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                  self._network_id, self._update_cx_file_name, self._server, self._user))

        with open(self._update_cx_file_name, 'br') as network_out:

            my_client = ndex2.client.Ndex2(host=self._server, username=self._user, password=self._pass)
            my_client.update_cx_network(network_out, self._network_id)

        print('{} - updated network {} from cx file {} on server {} for user {}'
              .format(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                  self._network_id, self._update_cx_file_name, self._server, self._user))

        return

    def run(self):
        """
        Runs content loading for NDEx STRING Content Loader
        :param theargs:
        :return:
        """
        self._parse_config()
        self._update_network_on_server()





def main(args):
    """
    Main entry point for program
    :param args:
    :return:
    """
    desc = """
    Version {version}

    Loads NDEx STRING Content Loader data into NDEx (http://ndexbio.org).

    To connect to NDEx server a configuration file must be passed
    into --conf parameter. If --conf is unset the configuration
    the path ~/{confname} is examined.

    The configuration file should be formatted as follows:

    [<value in --profile (default ncipid)>]

    {user} = <NDEx username>
    {password} = <NDEx password>
    {server} = <NDEx server(omit http) ie public.ndexbio.org>


    """.format(confname=NDExUtilConfig.CONFIG_FILE,
               user=NDExUtilConfig.USER,
               password=NDExUtilConfig.PASSWORD,
               server=NDExUtilConfig.SERVER,
               version=ndexstringloader.__version__)
    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = ndexstringloader.__version__

    try:
        _setup_logging(theargs)
        loader = NDExNdexstringloaderLoader(theargs)
        loader.run()
        return 0

    except Exception as e:
        print('\n   {} {}\n'.format(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")), e))
        logger.exception('Caught exception')
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
