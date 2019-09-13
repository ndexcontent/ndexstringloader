#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ndexstringloader` package."""

import os
import tempfile
import shutil

import unittest
from unittest.mock import MagicMock
from unittest import mock
from unittest.mock import Mock
from unittest.mock import patch

from ndexutil.config import NDExUtilConfig
from ndexstringloader.ndexloadstring import NDExSTRINGLoader
import ndexstringloader
from ndexstringloader import ndexloadstring

import ndex2

import json
import pandas as pd
import ndexutil.tsv.tsv2nicecx2 as t2n

import time

from unittest.mock import patch, mock_open
import requests
import requests_mock


from requests.exceptions import Timeout

class Param(object):
    """
    Dummy object
    """
    pass


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class TestNdexstringloader(unittest.TestCase):
    """Tests for `ndexstringloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

        self._args = {
            'conf': None,
            'profile': None,
            'loadplan': None,
            'stringversion': '11.0',
            'args': None,
            'datadir': tempfile.mkdtemp(),
            'cutoffscore': 0.7,
            'iconurl': 'https://home.ndexbio.org/img/STRING-logo.png'
        }

        self._args = dotdict(self._args)
        self._network_name = 'Network for Junit Testing STRING Loader - delete it'


    def tearDown(self):
        """Remove temp directory created by setUp"""
        if os.path.exists(self._args['datadir']):
            shutil.rmtree(self._args['datadir'])

    #@unittest.skip("skip it  now - uncomment later")
    def test_0010_parse_config(self):

        temp_dir = self._args['datadir']
        try:
            p = Param()
            self._args['profile'] = 'test_conf_section'
            self._args['conf'] = os.path.join(temp_dir, 'temp.conf')

            with open(self._args['conf'], 'w') as f:
                f.write('[' + self._args['profile'] + ']' + '\n')
                f.write(NDExUtilConfig.USER + ' = aaa\n')
                f.write(NDExUtilConfig.PASSWORD + ' = bbb\n')
                f.write(NDExUtilConfig.SERVER + ' = dev.ndexbio.org\n')
                f.flush()

            loader = NDExSTRINGLoader(self._args)
            loader._parse_config()
            self.assertEqual('aaa', loader._user)
            self.assertEqual('bbb', loader._pass)
            self.assertEqual('dev.ndexbio.org', loader._server)
        finally:
            shutil.rmtree(temp_dir)


    #@unittest.skip("skip it  now - uncomment later")
    def test_0020_remove_duplicate_edges(self):

        # some duplicate records in the same format as in STRING 9606.protein.links.full.v11.0.txt
        duplicate_records = [
            '9606.ENSP00000238651 9606.ENSP00000364486 0 0 0 0 0 0 45 0 0 800 0 0 0 800',
            '9606.ENSP00000268876 9606.ENSP00000216181 0 0 0 0 0 0 73 0 381 0 0 422 203 700',
            '9606.ENSP00000242462 9606.ENSP00000276480 0 0 0 0 0 0 0 0 0 0 0 0 401 400',
            '9606.ENSP00000364486 9606.ENSP00000238651 0 0 0 0 0 0 45 0 0 800 0 0 0 800',
            '9606.ENSP00000276480 9606.ENSP00000242462 0 0 0 0 0 0 0 0 0 0 0 0 401 400',
            '9606.ENSP00000216181 9606.ENSP00000268876 0 0 0 0 0 0 73 0 381 0 0 422 203 700'
        ]
        ensembl_ids = {
            '9606.ENSP00000216181': {
                'display_name': 'MYH9',
                'alias': 'ncbigene:4627|ensembl:ENSP00000216181',
                'represents': 'uniprot:P3557'
            },
            '9606.ENSP00000238651': {
                'display_name': 'ACOT2',
                'alias': 'ncbigene:10965|ensembl:ENSP00000238651',
                'represents': 'uniprot:P49753'
            },
            '9606.ENSP00000242462': {
                'display_name': 'NEUROG3',
                'alias': 'ncbigene:50674|ensembl:ENSP00000242462',
                'represents': 'uniprot:Q9Y4Z2'
            },
            '9606.ENSP00000268876': {
                'display_name': 'UNC45B',
                'alias': 'ncbigene:146862|ensembl:ENSP00000268876',
                'represents': 'uniprot:Q8IWX7'
            },
            '9606.ENSP00000276480': {
                'display_name': 'ST18',
                'alias': 'ncbigene:9705|ensembl:ENSP00000276480',
                'represents': 'uniprot:O60284'
            },
            '9606.ENSP00000364486': {
                'display_name': 'FBP2',
                'alias': 'ncbigene:8789|ensembl:ENSP00000364486',
                'represents': 'uniprot:O00757'
            }
        }

        temp_dir = self._args['datadir']

        try:
            string_loader = NDExSTRINGLoader(self._args)
            string_loader.__setattr__('ensembl_ids', ensembl_ids)

            file_with_duplicates = os.path.join(temp_dir, string_loader._full_file_name)

            # create file with duplicate records
            with open(file_with_duplicates, 'w') as o_f:
                o_f.write('header line' + '\n') # the first line is header; don't care what its content in this test
                for line in duplicate_records:
                    o_f.write(line + '\n')
                o_f.flush()

            # validate that the file with duplicate records was written fine
            with open(file_with_duplicates, 'r') as i_f:
                next(i_f)  # skip header
                index = 0
                for line in i_f:
                    self.assertEqual(line.rstrip(), duplicate_records[index])
                    index += 1

            # generate tsv file without duplicates
            string_loader.create_output_tsv_file()


            # records that should be in the new file after calling create_output_tsv_file
            unique_records = [
               'ACOT2\tuniprot:P49753\tncbigene:10965|ensembl:ENSP00000238651\tFBP2\tuniprot:O00757\tncbigene:8789|ensembl:ENSP00000364486\t0\t0\t0\t0\t0\t0\t45\t0\t0\t800\t0\t0\t0\t800',
               'UNC45B\tuniprot:Q8IWX7\tncbigene:146862|ensembl:ENSP00000268876\tMYH9\tuniprot:P3557\tncbigene:4627|ensembl:ENSP00000216181\t0\t0\t0\t0\t0\t0\t73\t0\t381\t0\t0\t422\t203\t700',
               'CDC16\tuniprot:Q13042\tncbigene:8881|ensembl:ENSP00000353549\tANAPC5\tuniprot:Q9UJX4\tncbigene:51433|ensembl:ENSP00000261819\t0\t0\t0\t0\t0\t102\t90\t987\t260\t900\t0\t754\t622\t999'
            ]

            # open the newly-generated file and validate that all records are unique
            with open(string_loader._output_tsv_file_name, 'r') as i_f:
                index = 0
                next(i_f) # skip header
                for line in i_f:
                    self.assertEqual(line.rstrip(), unique_records[index])
                    index += 1

        finally:
            shutil.rmtree(temp_dir)


    #@unittest.skip("skip it  now - uncomment later")
    def test_0030_exception_on_duplicate_edge_with_different_scores(self):

        # some duplicate records in the same format as in STRING 9606.protein.links.full.v11.0.txt
        duplicate_records = [
            '9606.ENSP00000238651 9606.ENSP00000364486 0 0 0 0 0 0 45 0 0 800 0 0 0 800',
            '9606.ENSP00000238651 9606.ENSP00000364486 0 0 0 0 0 0 45 0 0 800 0 0 0 801'
        ]
        ensembl_ids = {
            '9606.ENSP00000238651': {
                'display_name': 'ACOT2',
                'alias': 'ncbigene:10965|ensembl:ENSP00000238651',
                'represents': 'uniprot:P49753'
            },
            '9606.ENSP00000364486': {
                'display_name': 'FBP2',
                'alias': 'ncbigene:8789|ensembl:ENSP00000364486',
                'represents': 'uniprot:O00757'
            }
        }


        for i in range(0, 2):
            temp_dir = self._args['datadir']

            try:
                string_loader = NDExSTRINGLoader(self._args)
                string_loader.__setattr__('ensembl_ids', ensembl_ids)

                file_with_duplicates = os.path.join(temp_dir, string_loader._full_file_name)

                # create file with duplicate records
                with open(file_with_duplicates, 'w') as o_f:
                    o_f.write(
                        'header line' + '\n')  # the first line is header; don't care what its content in this test
                    for line in duplicate_records:
                        o_f.write(line + '\n')
                    o_f.flush()

                # validate that the file with duplicate records was written fine
                with open(file_with_duplicates, 'r') as i_f:
                    next(i_f)  # skip header
                    index = 0
                    for line in i_f:
                        self.assertEqual(line.rstrip(), duplicate_records[index])
                        index += 1

                with self.assertRaises(ValueError):
                    string_loader.create_output_tsv_file()

            finally:
                shutil.rmtree(temp_dir)
                self._args['datadir'] = tempfile.mkdtemp()

                # re-init duplicates and re-rerun the teast
                duplicate_records = [
                    '9606.ENSP00000238651 9606.ENSP00000364486 0 0 0 0 0 0 45 0 0 800 0 0 0 801',
                    '9606.ENSP00000364486 9606.ENSP00000238651 0 0 0 0 0 0 45 0 0 800 0 0 0 800'
                ]


    #@unittest.skip("skip it  now - uncomment later")
    def test_0040_init_network_atributes(self):
        net_attributes = {}

        cutoffscore = str(self._args['cutoffscore'])

        net_attributes['name'] = 'STRING - Human Protein Links - High Confidence (Score >= ' + cutoffscore + ')'

        net_attributes['description'] = '<br>This network contains high confidence (score >= ' \
                    + cutoffscore + ') human protein links with combined scores. ' \
                    + 'Edge color was mapped to the combined score value using a gradient from light grey ' \
                    + '(low Score) to black (high Score).'

        net_attributes['rights'] = 'Attribution 4.0 International (CC BY 4.0)'

        net_attributes['rightsHolder'] = 'STRING CONSORTIUM'

        net_attributes['version'] = self._args['stringversion']

        net_attributes['organism'] = 'Homo sapiens (human)'

        net_attributes['networkType'] = ['interactome', 'ppi']

        net_attributes['reference'] = '<p>Szklarczyk D, Morris JH, Cook H, Kuhn M, Wyder S, ' \
                    + 'Simonovic M, Santos A, Doncheva NT, Roth A, Bork P, Jensen LJ, von Mering C.<br><b> ' \
                    + 'The STRING database in 2017: quality-controlled protein-protein association networks, ' \
                    + 'made broadly accessible.</b><br>Nucleic Acids Res. 2017 Jan; ' \
                    + '45:D362-68.<br> <a target="_blank" href="https://doi.org/10.1093/nar/gkw937">' \
                    + 'DOI:10.1093/nar/gkw937</a></p>'

        net_attributes['prov:wasDerivedFrom'] = \
            'https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz'

        net_attributes['prov:wasGeneratedBy'] = \
            '<a href="https://github.com/ndexcontent/ndexstringloader" target="_blank">ndexstringloader ' \
            + str(ndexstringloader.__version__) + '</a>'

        net_attributes['__iconurl'] = self._args['iconurl']


        loader = NDExSTRINGLoader(self._args)

        # get network attributes from STRING loader object
        network_attributes = loader._init_network_attributes()

        self.assertDictEqual(net_attributes, network_attributes, 'unexpected network properties')


    #@unittest.skip("skip it  now - uncomment later")
    def test_0050_check_if_data_dir_exists(self):

        self._args['datadir'] = '__temp_dir_for_testing__'
        absolute_path = os.path.abspath(self._args['datadir'])

        if os.path.exists(absolute_path):
            os.rmdir(absolute_path)

        loader = NDExSTRINGLoader(self._args)

        # _check_if_data_dir_exists will create dir if it doesn't exist
        loader._check_if_data_dir_exists()
        self.assertTrue(os.path.exists(absolute_path))

        os.rmdir(absolute_path)
        self.assertFalse(os.path.exists(absolute_path))


    def test_0060_get_package_dir(self):
        actual_package_dir = ndexloadstring.get_package_dir()
        expected_package_dir = os.path.dirname(ndexstringloader.__file__)
        self.assertEqual(actual_package_dir, expected_package_dir)


    def test_0070_get_load_plan(self):
        actual_load_plan = ndexloadstring.get_load_plan()
        expected_load_plan = os.path.join(ndexloadstring.get_package_dir(), ndexloadstring.STRING_LOAD_PLAN)
        self.assertEqual(actual_load_plan, expected_load_plan)


    def test_0080_get_style(self):
        actual_style = ndexloadstring.get_style()
        expected_style = os.path.join(ndexloadstring.get_package_dir(), ndexloadstring.STYLE)
        self.assertEqual(actual_style, expected_style)


    def test_0090_parse_args(self):
        desc = """
        Version {version}

        Loads NDEx STRING Content Loader data into NDEx (http://ndexbio.org).

        To connect to NDEx server a configuration file must be passed
        into --conf parameter. If --conf is unset, the configuration
        ~/{confname} is examined.

        The configuration file should be formatted as follows:

        [<value in --profile (default dev)>]

        {user} = <NDEx username>
        {password} = <NDEx password>
        {server} = <NDEx server(omit http) ie public.ndexbio.org>
        """.format(confname=NDExUtilConfig.CONFIG_FILE,
                   user=NDExUtilConfig.USER,
                   password=NDExUtilConfig.PASSWORD,
                   server=NDExUtilConfig.SERVER,
                   version=ndexstringloader.__version__)
        temp_dir = self._args['datadir']

        args = []
        args.append('--cutoffscore')
        args.append('0.75')
        args.append(temp_dir)

        expected_args = {}
        expected_args['conf'] = None
        expected_args['cutoffscore'] = 0.75
        expected_args['datadir'] = temp_dir
        expected_args['iconurl'] = 'https://home.ndexbio.org/img/STRING-logo.png'
        expected_args['loadplan'] = os.path.join(ndexloadstring.get_package_dir(), ndexloadstring.STRING_LOAD_PLAN)
        expected_args['logconf'] = None
        expected_args['profile'] = 'ndexstringloader'
        expected_args['skipdownload'] = False
        expected_args['stringversion'] = '11.0'
        expected_args['style'] = os.path.join(ndexloadstring.get_package_dir(), ndexloadstring.STYLE)
        expected_args['verbose'] = 0
        expected_args['template'] = None
        expected_args['update'] = None

        the_args = ndexloadstring._parse_arguments(desc, args)

        self.assertDictEqual(the_args.__dict__, expected_args)


    def test_0100_setup_logging(self):

        verbose_level = 0
        args = {
                'logconf': None,
                'verbose': verbose_level
               }

        args = dotdict(args)
        ndexloadstring._setup_logging(args)

        logger_level_set = ndexloadstring.logger.getEffectiveLevel()

        self.assertEqual((50 - (10 * verbose_level)), logger_level_set)


    def test_0110_load_style_template(self):

        self._args.style = ndexloadstring.get_style()

        loader = NDExSTRINGLoader(self._args)

        loader._load_style_template()

        style_template_actual = loader.__getattribute__('_template')

        style_template_expected = \
            ndex2.create_nice_cx_from_file(os.path.abspath(os.path.join(ndexloadstring.get_package_dir(), 'style.cx')))

        self.assertDictEqual(style_template_actual.__dict__, style_template_expected.__dict__)


    def test_0120_get_headers_headers_of_links_file(self):
        header = [
            'protein1',
            'protein2',
            'neighborhood',
            'neighborhood_transferred',
            'fusion',
            'cooccurence',
            'homology',
            'coexpression',
            'coexpression_transferred',
            'experiments',
            'experiments_transferred',
            'database',
            'database_transferred',
            'textmining',
            'textmining_transferred',
            'combined_score'
        ]

        header_str = ' '.join(header)

        temp_dir = self._args['datadir']
        tempfile = os.path.join(temp_dir, '__temp_link_file__.txt')

        with open(tempfile, 'w') as f:
            f.write(header_str + '\n')
            f.flush()

        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_full_file_name', tempfile)

        header_actual = loader._get_headers_headers_of_links_file()

        self.assertEqual(header, header_actual)


    def test_0130_init_ensembl_ids(self):
        header = [
            'protein1',
            'protein2',
            'neighborhood',
            'neighborhood_transferred',
            'fusion',
            'cooccurence',
            'homology',
            'coexpression',
            'coexpression_transferred',
            'experiments',
            'experiments_transferred',
            'database',
            'database_transferred',
            'textmining',
            'textmining_transferred',
            'combined_score'
        ]
        content = [
            '9606.ENSP00000000233 9606.ENSP00000272298 0 0 0 332 0 0 62 0 181 0 0 0 125 490',
            '9606.ENSP00000000233 9606.ENSP00000253401 0 0 0 0 0 0 0 0 186 0 0 0 56 198',
            '9606.ENSP00000000233 9606.ENSP00000401445 0 0 0 0 0 0 0 0 160 0 0 0 0 159',
            '9606.ENSP00000000233 9606.ENSP00000418915 0 0 0 0 0 0 61 0 158 0 0 542 0 606',
            '9606.ENSP00000000233 9606.ENSP00000327801 0 0 0 0 0 69 61 0 78 0 0 0 89 167',
            '9606.ENSP00000000233 9606.ENSP00000466298 0 0 0 0 0 141 0 0 131 0 0 0 98 267',
            '9606.ENSP00000000233 9606.ENSP00000232564 0 0 0 0 0 0 62 0 171 0 0 0 56 201',
            '9606.ENSP00000000233 9606.ENSP00000393379 0 0 0 0 0 0 61 0 131 0 0 0 43 150',
            '9606.ENSP00000000233 9606.ENSP00000371253 0 0 0 0 0 0 61 0 0 0 0 0 224 240',
            '9606.ENSP00000000233 9606.ENSP00000373713 0 0 0 0 0 0 63 0 63 0 0 0 237 271'
        ]

        ensembl_ids_expected = {
            '9606.ENSP00000000233': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000272298': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000253401': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000401445': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000418915': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000327801': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000466298': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000232564': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000393379': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000371253': { 'display_name': None, 'alias': None, 'represents': None },
            '9606.ENSP00000373713': { 'display_name': None, 'alias': None, 'represents': None }
        }

        header_str = ' '.join(header)

        temp_dir = self._args['datadir']
        tempfile = os.path.join(temp_dir, '__temp_link_file__.txt')

        with open(tempfile, 'w') as f:
            f.write(header_str + '\n')
            for c in content:
                f.write(c + '\n')
            f.flush()

        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_full_file_name', tempfile)

        loader._init_ensembl_ids()

        ensembl_ids_actual = loader.__getattribute__(('ensembl_ids'))

        self.assertEqual(ensembl_ids_expected, ensembl_ids_actual)


    def test_0140_populate_display_names(self):
        links_header = [
            'protein1',
            'protein2',
            'neighborhood',
            'neighborhood_transferred',
            'fusion',
            'cooccurence',
            'homology',
            'coexpression',
            'coexpression_transferred',
            'experiments',
            'experiments_transferred',
            'database',
            'database_transferred',
            'textmining',
            'textmining_transferred',
            'combined_score'
        ]
        links_content = [
            '9606.ENSP00000000233 9606.ENSP00000272298 0 0 0 332 0 0 62 0 181 0 0 0 125 490',
            '9606.ENSP00000000233 9606.ENSP00000253401 0 0 0 0 0 0 0 0 186 0 0 0 56 198',
            '9606.ENSP00000000233 9606.ENSP00000401445 0 0 0 0 0 0 0 0 160 0 0 0 0 159',
            '9606.ENSP00000000233 9606.ENSP00000418915 0 0 0 0 0 0 61 0 158 0 0 542 0 606',
            '9606.ENSP00000000233 9606.ENSP00000327801 0 0 0 0 0 69 61 0 78 0 0 0 89 167',
            '9606.ENSP00000000233 9606.ENSP00000466298 0 0 0 0 0 141 0 0 131 0 0 0 98 267',
            '9606.ENSP00000000233 9606.ENSP00000232564 0 0 0 0 0 0 62 0 171 0 0 0 56 201',
            '9606.ENSP00000000233 9606.ENSP00000393379 0 0 0 0 0 0 61 0 131 0 0 0 43 150',
            '9606.ENSP00000000233 9606.ENSP00000371253 0 0 0 0 0 0 61 0 0 0 0 0 224 240',
            '9606.ENSP00000000233 9606.ENSP00000373713 0 0 0 0 0 0 63 0 63 0 0 0 237 271'
        ]
        links_header_str = ' '.join(links_header)

        ensembl_ids_expected = {
            '9606.ENSP00000000233': { 'display_name': 'ARF5', 'alias': None, 'represents': None },
            '9606.ENSP00000272298': { 'display_name': 'CALM2', 'alias': None, 'represents': None },
            '9606.ENSP00000253401': { 'display_name': 'ARHGEF9', 'alias': None, 'represents': None },
            '9606.ENSP00000401445': { 'display_name': 'ERN1', 'alias': None, 'represents': None },
            '9606.ENSP00000418915': { 'display_name': 'CDKN2A', 'alias': None, 'represents': None },
            '9606.ENSP00000327801': { 'display_name': 'P4HB', 'alias': None, 'represents': None },
            '9606.ENSP00000466298': { 'display_name': 'STX10', 'alias': None, 'represents': None },
            '9606.ENSP00000232564': { 'display_name': 'GNB4', 'alias': None, 'represents': None },
            '9606.ENSP00000393379': { 'display_name': 'KIF5C', 'alias': None, 'represents': None },
            '9606.ENSP00000371253': { 'display_name': 'GART', 'alias': None, 'represents': None },
            '9606.ENSP00000373713': { 'display_name': 'SACM1L', 'alias': None, 'represents': None }
        }

        #  names header is '# NCBI taxid / display name / STRING'
        names_header = [
            '# NCBI taxid',
            'display name',
            'STRING'
        ]
        names_header_str = ' / '.join(names_header)

        # the last entry '9606\tARF55\t9606.ENSP00000000233' has the same id as the first
        # '9606\tARF5\t9606.ENSP00000000233'; this should never happen in reality
        names_content = [
            '9606\tARF5\t9606.ENSP00000000233',
            '9606\tCALM2\t9606.ENSP00000272298',
            '9606\tARHGEF9\t9606.ENSP00000253401',
            '9606\tERN1\t9606.ENSP00000401445',
            '9606\tCDKN2A\t9606.ENSP00000418915',
            '9606\tP4HB\t9606.ENSP00000327801',
            '9606\tSTX10\t9606.ENSP00000466298',
            '9606\tGNB4\t9606.ENSP00000232564',
            '9606\tKIF5C\t9606.ENSP00000393379',
            '9606\tGART\t9606.ENSP00000371253',
            '9606\tSACM1L\t9606.ENSP00000373713',
            '9606\tARF55\t9606.ENSP00000000233'
        ]

        temp_dir = self._args['datadir']
        temp_links_file = os.path.join(temp_dir, '__temp_link_file__.txt')
        temp_names_file = os.path.join(temp_dir, '__temp_name_file__.txt')

        with open(temp_links_file, 'w') as f:
            f.write(links_header_str + '\n')
            for l in links_content:
                f.write(l + '\n')
            f.flush()

        with open(temp_names_file, 'w') as f:
            f.write(names_header_str + '\n')
            for n in names_content:
                f.write(n + '\n')
            f.flush()


        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_full_file_name', temp_links_file)
        loader.__setattr__('_names_file', temp_names_file)

        loader._init_ensembl_ids()

        loader._populate_display_names()
        ensembl_ids_actual = loader.__getattribute__(('ensembl_ids'))
        self.assertEqual(ensembl_ids_expected, ensembl_ids_actual)

        duplicate_names = loader.__getattribute__('duplicate_display_names')
        self.assertEqual(duplicate_names, {'9606.ENSP00000000233': ['ARF5', 'ARF55']})


    def test_0150_populate_aliases(self):
        links_header = [
            'protein1',
            'protein2',
            'neighborhood',
            'neighborhood_transferred',
            'fusion',
            'cooccurence',
            'homology',
            'coexpression',
            'coexpression_transferred',
            'experiments',
            'experiments_transferred',
            'database',
            'database_transferred',
            'textmining',
            'textmining_transferred',
            'combined_score'
        ]
        links_content = [
            '9606.ENSP00000000233 9606.ENSP00000253401 0 0 0 0 0 0 0 0 186 0 0 0 56 198',
            '9606.ENSP00000000233 9606.ENSP00000401445 0 0 0 0 0 0 0 0 160 0 0 0 0 159',
            '9606.ENSP00000000233 9606.ENSP00000418915 0 0 0 0 0 0 61 0 158 0 0 542 0 606',
            '9606.ENSP00000000233 9606.ENSP00000327801 0 0 0 0 0 69 61 0 78 0 0 0 89 167',
            '9606.ENSP00000000233 9606.ENSP00000466298 0 0 0 0 0 141 0 0 131 0 0 0 98 267',
            '9606.ENSP00000000233 9606.ENSP00000232564 0 0 0 0 0 0 62 0 171 0 0 0 56 201',
            '9606.ENSP00000000233 9606.ENSP00000393379 0 0 0 0 0 0 61 0 131 0 0 0 43 150',
            '9606.ENSP00000000233 9606.ENSP00000371253 0 0 0 0 0 0 61 0 0 0 0 0 224 240',
            '9606.ENSP00000000233 9606.ENSP00000373713 0 0 0 0 0 0 63 0 63 0 0 0 237 271',
            '9606.ENSP00000000233 9606.ENSP00000479069 0 0 0 0 0 0 0 0 70 0 0 0 215 238'
        ]
        links_header_str = ' '.join(links_header)

        ensembl_ids_expected = {
            '9606.ENSP00000000233': { 'display_name': None, 'alias': 'ncbigene:381|ensembl:ENSP00000000233', 'represents': None },
            '9606.ENSP00000253401': { 'display_name': None, 'alias': 'ncbigene:23229|ensembl:ENSP00000253401', 'represents': None },
            '9606.ENSP00000401445': { 'display_name': None, 'alias': 'ncbigene:2081|ensembl:ENSP00000401445', 'represents': None },
            '9606.ENSP00000418915': { 'display_name': None, 'alias': 'ncbigene:1029|ensembl:ENSP00000418915', 'represents': None },
            '9606.ENSP00000327801': { 'display_name': None, 'alias': 'ncbigene:5034|ensembl:ENSP00000327801', 'represents': None },
            '9606.ENSP00000466298': { 'display_name': None, 'alias': 'ncbigene:8677|ensembl:ENSP00000466298', 'represents': None },
            '9606.ENSP00000232564': { 'display_name': None, 'alias': 'ncbigene:59345|ensembl:ENSP00000232564', 'represents': None },
            '9606.ENSP00000393379': { 'display_name': None, 'alias': 'ncbigene:3800|ensembl:ENSP00000393379', 'represents': None },
            '9606.ENSP00000371253': { 'display_name': None, 'alias': 'ncbigene:2618|ensembl:ENSP00000371253', 'represents': None },
            '9606.ENSP00000373713': { 'display_name': None, 'alias': 'ncbigene:22908|ensembl:ENSP00000373713', 'represents': None },
            '9606.ENSP00000479069': { 'display_name': None, \
                                      'alias': 'ncbigene:101930165|ncbigene:105369241|ncbigene:728929|ensembl:ENSP00000479069', \
                                      'represents': None }
        }

        #  entrez header is '# NCBI taxid / entrez / STRING'
        entrez_header = [
            '# NCBI taxid',
            'entrez',
            'STRING'
        ]
        entrez_header_str = ' / '.join(entrez_header)

        entrez_content = [
            '9606\t381\t9606.ENSP00000000233',
            '9606\t23229\t9606.ENSP00000253401',
            '9606\t2081\t9606.ENSP00000401445',
            '9606\t1029\t9606.ENSP00000418915',
            '9606\t5034\t9606.ENSP00000327801',
            '9606\t8677\t9606.ENSP00000466298',
            '9606\t59345\t9606.ENSP00000232564',
            '9606\t3800\t9606.ENSP00000393379',
            '9606\t2618\t9606.ENSP00000371253',
            '9606\t22908\t9606.ENSP00000373713',
            '9606\t101930165|105369241|728929\t9606.ENSP00000479069'
        ]

        temp_dir = self._args['datadir']
        temp_links_file = os.path.join(temp_dir, '__temp_link_file__.txt')
        temp_entrez_file = os.path.join(temp_dir, '__temp_entrez_file__.txt')

        with open(temp_links_file, 'w') as f:
            f.write(links_header_str + '\n')
            for l in links_content:
                f.write(l + '\n')
            f.flush()

        with open(temp_entrez_file, 'w') as f:
            f.write(entrez_header_str + '\n')
            for e in entrez_content:
                f.write(e + '\n')
            f.flush()


        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_full_file_name', temp_links_file)
        loader.__setattr__('_entrez_file', temp_entrez_file)

        loader._init_ensembl_ids()

        loader._populate_aliases()
        self.maxDiff = None
        eids = loader.__getattribute__('ensembl_ids')
        self.assertEqual(ensembl_ids_expected, loader.__getattribute__('ensembl_ids'))


    def test_0160_populate_represents(self):
        links_header = [
            'protein1',
            'protein2',
            'neighborhood',
            'neighborhood_transferred',
            'fusion',
            'cooccurence',
            'homology',
            'coexpression',
            'coexpression_transferred',
            'experiments',
            'experiments_transferred',
            'database',
            'database_transferred',
            'textmining',
            'textmining_transferred',
            'combined_score'
        ]
        links_content = [
            '9606.ENSP00000000233 9606.ENSP00000253401 0 0 0 0 0 0 0 0 186 0 0 0 56 198',
            '9606.ENSP00000000233 9606.ENSP00000401445 0 0 0 0 0 0 0 0 160 0 0 0 0 159',
            '9606.ENSP00000000233 9606.ENSP00000418915 0 0 0 0 0 0 61 0 158 0 0 542 0 606',
            '9606.ENSP00000000233 9606.ENSP00000327801 0 0 0 0 0 69 61 0 78 0 0 0 89 167',
            '9606.ENSP00000000233 9606.ENSP00000466298 0 0 0 0 0 141 0 0 131 0 0 0 98 267',
            '9606.ENSP00000000233 9606.ENSP00000232564 0 0 0 0 0 0 62 0 171 0 0 0 56 201',
            '9606.ENSP00000000233 9606.ENSP00000393379 0 0 0 0 0 0 61 0 131 0 0 0 43 150',
            '9606.ENSP00000000233 9606.ENSP00000371253 0 0 0 0 0 0 61 0 0 0 0 0 224 240',
            '9606.ENSP00000000233 9606.ENSP00000373713 0 0 0 0 0 0 63 0 63 0 0 0 237 271'
        ]
        links_header_str = ' '.join(links_header)

        ensembl_ids_expected = {
            '9606.ENSP00000000233': { 'display_name': None, 'alias': None, 'represents': 'uniprot:P84085' },
            '9606.ENSP00000253401': { 'display_name': None, 'alias': None, 'represents': 'uniprot:O43307' },
            '9606.ENSP00000401445': { 'display_name': None, 'alias': None, 'represents': 'uniprot:O75460' },
            '9606.ENSP00000418915': { 'display_name': None, 'alias': None, 'represents': 'uniprot:P42771' },
            '9606.ENSP00000327801': { 'display_name': None, 'alias': None, 'represents': 'uniprot:P07237' },
            '9606.ENSP00000466298': { 'display_name': None, 'alias': None, 'represents': 'uniprot:O60499' },
            '9606.ENSP00000232564': { 'display_name': None, 'alias': None, 'represents': 'uniprot:Q9HAV0' },
            '9606.ENSP00000393379': { 'display_name': None, 'alias': None, 'represents': 'uniprot:O60282' },
            '9606.ENSP00000371253': { 'display_name': None, 'alias': None, 'represents': 'uniprot:P22102' },
            '9606.ENSP00000373713': { 'display_name': None, 'alias': None, 'represents': 'uniprot:Q9NTJ5' }
        }


        # uniprot file doesn't have header
        uniprot_content = [
            '9606\tP84085|ARF5_HUMAN\t9606.ENSP00000000233\t100.0\t374.0',
            '9606\tO43307|ARHG9_HUMAN\t9606.ENSP00000253401\t100.0\t1078.0',
            '9606\tO75460|ERN1_HUMAN\t9606.ENSP00000401445\t100.0\t2028.0',
            '9606\tP42771|CDN2A_HUMAN\t9606.ENSP00000418915\t98.701\t298.0',
            '9606\tP07237|PDIA1_HUMAN\t9606.ENSP00000327801\t100.0\t1037.0',
            '9606\tO60499|STX10_HUMAN\t9606.ENSP00000466298\t100.0\t507.0',
            '9606\tQ9HAV0|GBB4_HUMAN\t9606.ENSP00000232564\t100.0\t703.0',
            '9606\tO60282|KIF5C_HUMAN\t9606.ENSP00000393379\t100.0\t1965.0',
            '9606\tP22102|PUR2_HUMAN\t9606.ENSP00000371253\t100.0\t2064.0',
            '9606\tQ9NTJ5|SAC1_HUMAN\t9606.ENSP00000373713\t100.0\t1230.0'
        ]

        temp_dir = self._args['datadir']
        temp_links_file = os.path.join(temp_dir, '__temp_link_file__.txt')
        temp_uniprot_file = os.path.join(temp_dir, '__temp_uniprot_file__.txt')

        with open(temp_links_file, 'w') as f:
            f.write(links_header_str + '\n')
            for l in links_content:
                f.write(l + '\n')
            f.flush()

        with open(temp_uniprot_file, 'w') as f:
            for u in uniprot_content:
                f.write(u + '\n')
            f.flush()


        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_full_file_name', temp_links_file)
        loader.__setattr__('_uniprot_file', temp_uniprot_file)

        loader._init_ensembl_ids()

        loader._populate_represents()
        self.maxDiff = None
        self.assertEqual(ensembl_ids_expected, loader.__getattribute__('ensembl_ids'))



        os.remove(temp_uniprot_file)

        # test for unlikely situation when one key maps to different Unirpot IDs:
        # the first two entries have the same id 9606.ENSP00000000233
        # mapping to different Uniprot IDs - this should never happen in reality though
        uniprot_content_duplicate_ids = [
            '9606\tP84085|ARF5_HUMAN\t9606.ENSP00000000233\t100.0\t374.0',
            '9606\tO43307|ARHG9_HUMAN\t9606.ENSP00000000233\t100.0\t1078.0',
            '9606\tO75460|ERN1_HUMAN\t9606.ENSP00000401445\t100.0\t2028.0'
        ]
        with open(temp_uniprot_file, 'w') as f:
            for u in uniprot_content_duplicate_ids:
                f.write(u + '\n')
            f.flush()

        loader._init_ensembl_ids()
        loader._populate_represents()
        duplicate_uniprot_ids = loader.__getattribute__('duplicate_uniprot_ids')
        self.assertEqual({'9606.ENSP00000000233': ['uniprot:P84085', 'uniprot:O43307']}, duplicate_uniprot_ids)


    def test_0170_get_name_rep_alias(self):

        ensembl_ids = {
            '9606.ENSP00000216181': {
                'display_name': None,
                'represents': None,
                'alias': None
            },
            '9606.ENSP00000238651': {
                'display_name': 'ACOT2',
                'represents': None,
                'alias': 'ncbigene:10965|ensembl:ENSP00000238651'
            },
            '9606.ENSP00000242462': {
                'display_name': 'NEUROG3',
                'represents': 'uniprot:Q9Y4Z2',
                'alias': None
            },
            '9606.ENSP00000268876': {
                'display_name': None,
                'represents': 'uniprot:Q8IWX7',
                'alias': 'ncbigene:146862|ensembl:ENSP00000268876'
            },
            '9606.ENSP00000276480': {
                'display_name': 'ST18',
                'represents': 'uniprot:O60284',
                'alias': 'ncbigene:9705|ensembl:ENSP00000276480'
            },
            '9606.ENSP00000364486': {
                'display_name': 'FBP2',
                'represents': 'uniprot:O00757',
                'alias': 'ncbigene:8789|ensembl:ENSP00000364486',
            }
        }

        represents_expected = {
            '9606.ENSP00000216181': 'ensembl:ENSP00000216181\tensembl:ENSP00000216181\tensembl:ENSP00000216181',
            '9606.ENSP00000238651': 'ACOT2\thgnc:ACOT2\tncbigene:10965|ensembl:ENSP00000238651',
            '9606.ENSP00000242462': 'NEUROG3\tuniprot:Q9Y4Z2\tuniprot:Q9Y4Z2',
            '9606.ENSP00000268876': 'ensembl:ENSP00000268876\tuniprot:Q8IWX7\tncbigene:146862|ensembl:ENSP00000268876',
            '9606.ENSP00000276480': 'ST18\tuniprot:O60284\tncbigene:9705|ensembl:ENSP00000276480',
            '9606.ENSP00000364486': 'FBP2\tuniprot:O00757\tncbigene:8789|ensembl:ENSP00000364486'
        }

        string_loader = NDExSTRINGLoader(self._args)
        string_loader.__setattr__('ensembl_ids', ensembl_ids)

        for key, value in ensembl_ids.items():
            name_rep_alias = string_loader._get_name_rep_alias(key)
            self.assertEqual(name_rep_alias, represents_expected[key])

    def test_0180_create_NDEx_connection(self):
        loader = NDExSTRINGLoader(self._args)

        user_name = 'aaa'
        password = 'aaa'
        server = 'dev.ndexbio.org'

        loader.__setattr__('_pass', password)
        loader.__setattr__('_server', server)
        ndex_client = loader.create_ndex_connection()
        self.assertIsNone(ndex_client)

        loader.__setattr__('_user', 'aaa')
        ndex_client = loader.create_ndex_connection()
        self.assertIsNotNone(ndex_client)

        user_obj = ndex_client.get_user_by_username(user_name)
        self.assertTrue(len(user_obj) == 15)


    def test_0190_generate_CX_file(self):

        tsv_header = [
            'name1',
            'represents1',
            'alias1',
            'name2',
            'represents2',
            'alias2'
            'neighborhood',
            'neighborhood_transferred',
            'fusion',
            'cooccurence',
            'homology',
            'coexpression',
            'coexpression_transferred',
            'experiments',
            'experiments_transferred',
            'database',
            'database_transferred',
            'textmining',
            'textmining_transferred',
            'combined_score'
        ]

        tsv_header_str = '\t'.join(tsv_header) + '\n'

        tsv_body = [
            'VCL\tuniprot:P18206\tncbigene:7414|ensembl:ENSP00000211998\tTLN1\tuniprot:Q9Y490\tncbigene:7094|ensembl:ENSP00000316029\t0\t0\t0\t0\t0\t106\t82\t870\t809\t900\t0\t701\t538\t999',
            'VCL\tuniprot:P18206\tncbigene:7414|ensembl:ENSP00000211998\tPXN\tuniprot:P49023\tncbigene:5829|ensembl:ENSP00000267257\t0\t0\t0\t0\t0\t76\t63\t888\t377\t900\t0\t957\t534\t999',
            'VCL\tuniprot:P18206\tncbigene:7414|ensembl:ENSP00000211998\tACTN1\tuniprot:P12814\tncbigene:87|ensembl:ENSP00000377941\t0\t0\t0\t0\t0\t242\t81\t870\t809\t900\t0\t556\t504\t999'
        ]

        temp_dir = self._args['datadir']
        temp_links_tsv_file = os.path.join(temp_dir, '__protein_links_tmp__.tsv')
        temp_cx_network = os.path.join(temp_dir, '__networks__.cx')


        with open(temp_links_tsv_file, 'w') as f:
            f.write(tsv_header_str)
            for t in tsv_body:
                f.write(t + '\n')
            f.flush()


        self._args.style = ndexloadstring.get_style()
        loader = NDExSTRINGLoader(self._args)

        loader.__setattr__('_output_tsv_file_name', temp_links_tsv_file)
        loader.__setattr__('_cx_network', temp_cx_network)
        loader.__setattr__('_load_plan', ndexloadstring.get_load_plan())
        loader._load_style_template()


        network_attributes = loader._init_network_attributes()
        network_attributes['name'] = self._network_name

        loader._generate_CX_file(network_attributes)

        network_in_nice_cx_1 = ndex2.create_nice_cx_from_file(temp_cx_network)

        nice_cx_path = ndexloadstring.get_package_dir() + '/../tests/test_network.cx'
        #with open(nice_cx_path, 'w') as f:
        #    f.write(json.dumps(network_in_nice_cx_1.to_cx(), indent=4))

        network_in_nice_cx_2 = ndex2.create_nice_cx_from_file(nice_cx_path)

        dict_1 = network_in_nice_cx_1.__dict__
        dict_2 = network_in_nice_cx_2.__dict__

        self.maxDiff = None
        self.assertDictEqual(dict_1, dict_2)


    def test_0200_load_network_to_server_cx_network_is_none(self):
        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_user', 'u')
        loader.__setattr__('_pass', 'p')
        loader.__setattr__('_server', 's')
        try:
            res = loader._load_network_to_server('ha')
            self.assertEqual(2, res)
            self.fail('Expected exception')
        except FileNotFoundError as fe:
            self.assertTrue('No such file' in str(fe))

    def test_0210_load_network_to_server_cx_network_no_id(self):
        cxfile = os.path.join(self._args.datadir, 'hello.cx')
        with open(cxfile, 'w') as f:
            f.write('hello')
        loader = NDExSTRINGLoader(self._args)
        mockndex = MagicMock()
        mockndex.save_cx_stream_as_new_network = MagicMock()
        loader.set_ndex_connection(mockndex)
        loader.__setattr__('_user', 'u')
        loader.__setattr__('_pass', 'p')
        loader.__setattr__('_server', 's')
        loader.__setattr__('_cx_network', cxfile)
        res = loader._load_network_to_server('ha')
        self.assertEqual(0, res)
        mockndex.save_cx_stream_as_new_network.assert_called()

    def test_0220_load_network_to_server_cx_network_raise_except(self):
        cxfile = os.path.join(self._args.datadir, 'hello.cx')
        with open(cxfile, 'w') as f:
            f.write('hello')
        loader = NDExSTRINGLoader(self._args)
        mockndex = MagicMock()
        mockndex.save_cx_stream_as_new_network = MagicMock(side_effect=Exception('hi'))
        loader.set_ndex_connection(mockndex)
        loader.__setattr__('_user', 'u')
        loader.__setattr__('_pass', 'p')
        loader.__setattr__('_server', 's')
        loader.__setattr__('_cx_network', cxfile)
        res = loader._load_network_to_server('ha')
        self.assertEqual(2, res)
        mockndex.save_cx_stream_as_new_network.assert_called()

    def test_0230_update_network_on_server_cx_network_withid(self):
        cxfile = os.path.join(self._args.datadir, 'hello.cx')
        with open(cxfile, 'w') as f:
            f.write('hello')
        loader = NDExSTRINGLoader(self._args)
        mockndex = MagicMock()
        mockndex.update_cx_network = MagicMock()
        mockndex.update_cx_network.return_value = 0
        loader.set_ndex_connection(mockndex)
        loader.__setattr__('_user', 'u')
        loader.__setattr__('_pass', 'p')
        loader.__setattr__('_server', 's')
        loader.__setattr__('_cx_network', cxfile)
        res = loader._update_network_on_server('haha', network_id='hehe')
        self.assertEqual(0, res)
        mockndex.update_cx_network.assert_called()


    def test_0240_get_network_uuid(self):

        loader = NDExSTRINGLoader(self._args)
        mockndex = MagicMock()

        network_summaries_for_mock = [
            {'name':'Network 1', 'externalId':'111-111-111'},
            {'name':'Network 2', 'externalId':'222-222-222'},
            {'name':'Network 3', 'externalId':'333-333-333'},
            {'name':'Network 4', 'externalId':'444-444-444'}
        ]

        mockndex.get_network_summaries_for_user = MagicMock(return_value=network_summaries_for_mock)

        loader.set_ndex_connection(mockndex)

        loader.__setattr__('_user', 'u')
        loader.__setattr__('_pass', 'p')
        loader.__setattr__('_server', 's')

        # test scenario where Network Name is found in network summaries
        for summary in network_summaries_for_mock:
            network_name = summary.get('name')
            network_uuid = summary.get('externalId')
            network_uuid_from_server = loader.get_network_uuid(network_name, network_summaries_for_mock)

            #mockndex.get_network_summaries_for_user.assert_called_with(loader.__getattribute__('_user'))
            self.assertEqual(network_uuid,network_uuid_from_server)

        # test scenario where Network Name is not found in network summaries
        network_name = 'Non Existant Name'
        network_uuid_from_server = loader.get_network_uuid(network_name, network_summaries_for_mock)
        #mockndex.get_network_summaries_for_user.assert_called_with(loader.__getattribute__('_user'))
        self.assertIsNone(network_uuid_from_server)

        # test scenario where Network Name is not found in network summaries
        #mockndex.get_network_summaries_for_user.side_effect = Exception('Server is Down')
        #network_uuid_from_server = loader.get_network_uuid(network_name, network_summaries_for_mock)
        #self.assertEqual(network_uuid_from_server, 2)


    def test_0250_download(self):

        entrez_url = \
            'https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz'

        local_file_name = 'entrez.tsv'
        local_downloaded_file_name_unzipped = self._args['datadir'] + '/' + local_file_name
        local_downloaded_file_name_zipped = local_downloaded_file_name_unzipped + '.gz'

        loader = NDExSTRINGLoader(self._args)

        # expect to raise Timeout exception
        with patch('ndexstringloader.ndexloadstring.requests') as mock_requests:
            mock_requests.get.side_effect = Timeout
            with self.assertRaises(Timeout):
                loader._download(entrez_url, local_downloaded_file_name_zipped)
            mock_requests.get.assert_called_once_with(entrez_url)

        with patch('ndexstringloader.ndexloadstring.requests.get') as mock_get:
            not_found_code = 404
            mock_get.return_value.status_code = not_found_code
            assert loader._download(entrez_url, local_downloaded_file_name_zipped) == not_found_code

            mock_get.assert_called_once_with(entrez_url)

        with patch('ndexstringloader.ndexloadstring.requests.get') as mock_get:
            mock_get.return_value.status_code = 200
            mock_get.return_value.content = 'hello'.encode()  # this is the 'content' to be written to file

            assert loader._download(entrez_url, local_downloaded_file_name_zipped) == 0
            mock_get.assert_called_once_with(entrez_url)


    def test_0260_download_STRING_files(self):

        loader = NDExSTRINGLoader(self._args)

        _protein_links_url = loader.__getattribute__('_protein_links_url')
        _names_file_url = loader.__getattribute__('_names_file_url')
        _entrez_ids_file_url = loader.__getattribute__('_entrez_ids_file_url')
        _uniprot_file = loader.__getattribute__('_uniprot_ids_file_url')

        not_found_code = 404

        with requests_mock.mock() as m:
            # test scenario when all fiels are downloaded fine
            m.get(_protein_links_url, text='protein data returned by the get', status_code=200)
            m.get(_names_file_url, text='name data returned by the get', status_code=200)
            m.get(_entrez_ids_file_url, text='entrez data returned by the get', status_code=200)
            m.get(_uniprot_file, text='uniprot data returned by the get', status_code=200)
            assert loader._download_STRING_files() == 0

            # test scenario when failed to download _protein_links_url
            m.get(_protein_links_url, text='data returned by the get', status_code=not_found_code)
            assert loader._download_STRING_files() == not_found_code

            # test scenario when failed to download _names_file_url
            m.get(_protein_links_url, text='data returned by the get', status_code=200)
            m.get(_names_file_url, text='name data returned by the get', status_code=not_found_code)
            assert loader._download_STRING_files() == not_found_code

            # test scenario when failed to download _entrez_ids_file_url
            m.get(_protein_links_url, text='protein data returned by the get', status_code=200)
            m.get(_names_file_url, text='name data returned by the get', status_code=200)
            m.get(_entrez_ids_file_url, text='entrez data returned by the get', status_code=not_found_code)
            assert loader._download_STRING_files() == not_found_code

            # test scenario when failed to download _uniprot_file
            m.get(_protein_links_url, text='protein data returned by the get', status_code=200)
            m.get(_names_file_url, text='name data returned by the get', status_code=200)
            m.get(_entrez_ids_file_url, text='entrez data returned by the get', status_code=200)
            m.get(_uniprot_file, text='uniprot data returned by the get', status_code=not_found_code)
            assert loader._download_STRING_files() == not_found_code


    #@unittest.skip("needs to be finished")
    def test_0270_unpack_STRING_files(self):

        loader = NDExSTRINGLoader(self._args)
        vals = {
            ('full.gz'): ndexloadstring.ERROR_CODE,
            ('entrez.gz'): ndexloadstring.ERROR_CODE,
            ('names.gz'): ndexloadstring.ERROR_CODE,
            ('uniprot.gz'): ndexloadstring.ERROR_CODE
        }

        def side_effect(args):
            ret_value = vals[args]
            return ret_value

        loader._unzip = MagicMock(return_value = 0)

        loader.__setattr__('_full_file_name', 'full')
        loader.__setattr__('_entrez_file', 'entrez')
        loader.__setattr__('_names_file', 'names')
        loader.__setattr__('_uniprot_file', 'uniprot')

        ret_value = loader._unpack_STRING_files()
        self.assertEqual(ret_value, 0)

        loader._unzip = MagicMock(side_effect=side_effect)
        ret_value = loader._unpack_STRING_files()
        self.assertEqual(ret_value, ndexloadstring.ERROR_CODE)

        # now, I want simulate a scenario where unzipping full.gz succeeds, but unzipping entrez.gz fails
        vals[('full.gz')] = ndexloadstring.SUCCESS_CODE
        ret_value = loader._unpack_STRING_files()
        self.assertEqual(ret_value, ndexloadstring.ERROR_CODE)

        vals[('entrez.gz')] = ndexloadstring.SUCCESS_CODE
        ret_value = loader._unpack_STRING_files()
        self.assertEqual(ret_value, ndexloadstring.ERROR_CODE)

        vals[('names.gz')] = ndexloadstring.SUCCESS_CODE
        ret_value = loader._unpack_STRING_files()
        self.assertEqual(ret_value, ndexloadstring.ERROR_CODE)

        vals[('uniprot.gz')] = ndexloadstring.SUCCESS_CODE
        ret_value = loader._unpack_STRING_files()
        self.assertEqual(ret_value, ndexloadstring.SUCCESS_CODE)


    @mock.patch('ndexstringloader.ndexloadstring.gzip.open')
    @mock.patch('ndexstringloader.ndexloadstring.open')
    @mock.patch('ndexstringloader.ndexloadstring.shutil.copyfileobj')
    @mock.patch('ndexstringloader.ndexloadstring.os.remove')
    def test_0280_unzip(self, mock_remove, mock_copyfileobj, mock_open, mock_gzopen):

        loader = NDExSTRINGLoader(self._args)

        full_file_name = loader.__getattribute__('_full_file_name')
        full_file_name_gz = full_file_name + '.gz'

        ret_value = loader._unzip(full_file_name_gz)
        self.assertEqual(ret_value, 0)

        mock_gzopen.assert_called_once()
        mock_gzopen.assert_called_with(full_file_name_gz, 'rb')

        mock_open.assert_called_once()
        mock_open.assert_called_with(full_file_name, 'wb')

        mock_remove.assert_called_once()
        mock_remove.assert_called_with(full_file_name_gz)

        # hmm, not sure how to test args passed to  shutil.copyfileobj(f_in, f_out)  here ...
        # will just make sure it was called ...
        #mock_copyfileobj.assert_called_with(mock_gzopen, mock_open)
        mock_copyfileobj.assert_called_once()


    def test_0290_is_valid_update_UUID(self):

        loader = NDExSTRINGLoader(self._args)

        loader.__setattr__('_update_UUID', None)
        self.assertTrue(loader._is_valid_update_UUID())

        loader.__setattr__('_update_UUID', 'a62c9252-ce13-11e9-8bd8-525400c25d22')
        self.assertTrue(loader._is_valid_update_UUID())

        loader.__setattr__('_update_UUID', 'a62c9252-ce13-11e9-8bd8-525400c25d2')
        self.assertFalse(loader._is_valid_update_UUID())



    @unittest.skip("this test actually gets human.entrez_2_string.2018.tsv.gz from STRING server; we skip it")
    def test_1000_download_and_unzip(self):

        entrez_url = \
            'https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz'

        local_file_name = 'entrez.tsv'
        local_downloaded_file_name_unzipped = self._args['datadir'] + '/' + local_file_name
        local_downloaded_file_name_zipped = local_downloaded_file_name_unzipped + '.gz'

        loader = NDExSTRINGLoader(self._args)

        loader._download(entrez_url, local_downloaded_file_name_zipped)
        self.assertTrue(os.path.exists(local_downloaded_file_name_zipped))

        loader._unzip(local_downloaded_file_name_zipped)
        self.assertTrue(os.path.exists(local_downloaded_file_name_unzipped))

    @unittest.skip("this test actually downloads files from server and unpacks them;  we skip it")
    def test_1010_download_and_unzip_STRING_files(self):

        loader = NDExSTRINGLoader(self._args)

        loader._download_STRING_files()

        full_file = loader.__getattribute__('_full_file_name') + '.gz'
        names_file = loader.__getattribute__('_names_file') + '.gz'
        entrez_file = loader.__getattribute__('_entrez_file') + '.gz'
        uniprot_file = loader.__getattribute__('_uniprot_file') + '.gz'

        self.assertTrue(os.path.exists(full_file))
        self.assertTrue(os.path.exists(names_file))
        self.assertTrue(os.path.exists(entrez_file))
        self.assertTrue(os.path.exists(uniprot_file))

        loader._unpack_STRING_files()

        full_file = loader.__getattribute__('_full_file_name')
        names_file = loader.__getattribute__('_names_file')
        entrez_file = loader.__getattribute__('_entrez_file')
        uniprot_file = loader.__getattribute__('_uniprot_file')

        self.assertTrue(os.path.exists(full_file))
        self.assertTrue(os.path.exists(names_file))
        self.assertTrue(os.path.exists(entrez_file))
        self.assertTrue(os.path.exists(uniprot_file))

    @unittest.skip("this test actually creates a small network and uploads it to server; we skip it")
    def test_1020_load_to_NDEx(self):
        user_name = 'aaa'
        password = 'aaa'
        server = 'dev.ndexbio.org'

        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_pass', password)
        loader.__setattr__('_server', server)


        self._args.style = ndexloadstring.get_style()

        # no user name was set - load_to_NDEx() is expected to fail to create connection with NDEx server
        ret_code = loader.load_to_NDEx()
        self.assertEqual(ret_code, 2)

        loader.__setattr__('_user', user_name)

        # _cutoffscore is not used in this test for filtering rows of initial file, it is used for generating network name
        loader.__setattr__('_cutoffscore', '0.999')

        tsv_header = [
            'name1',
            'represents1',
            'alias1',
            'name2',
            'represents2',
            'alias2'
            'neighborhood',
            'neighborhood_transferred',
            'fusion',
            'cooccurence',
            'homology',
            'coexpression',
            'coexpression_transferred',
            'experiments',
            'experiments_transferred',
            'database',
            'database_transferred',
            'textmining',
            'textmining_transferred',
            'combined_score'
        ]

        tsv_header_str = '\t'.join(tsv_header) + '\n'

        tsv_body = [
            'VCL\tuniprot:P18206\tncbigene:7414|ensembl:ENSP00000211998\tTLN1\tuniprot:Q9Y490\tncbigene:7094|ensembl:ENSP00000316029\t0\t0\t0\t0\t0\t106\t82\t870\t809\t900\t0\t701\t538\t999',
            'VCL\tuniprot:P18206\tncbigene:7414|ensembl:ENSP00000211998\tPXN\tuniprot:P49023\tncbigene:5829|ensembl:ENSP00000267257\t0\t0\t0\t0\t0\t76\t63\t888\t377\t900\t0\t957\t534\t999',
            'VCL\tuniprot:P18206\tncbigene:7414|ensembl:ENSP00000211998\tACTN1\tuniprot:P12814\tncbigene:87|ensembl:ENSP00000377941\t0\t0\t0\t0\t0\t242\t81\t870\t809\t900\t0\t556\t504\t999'
        ]

        temp_dir = self._args['datadir']
        temp_links_tsv_file = os.path.join(temp_dir, '__protein_links_tmp__.tsv')
        temp_cx_network = os.path.join(temp_dir, '__networks__.cx')


        with open(temp_links_tsv_file, 'w') as f:
            f.write(tsv_header_str)
            for t in tsv_body:
                f.write(t + '\n')
            f.flush()

        self._args.style = ndexloadstring.get_style()

        loader.__setattr__('_output_tsv_file_name', temp_links_tsv_file)
        loader.__setattr__('_cx_network', temp_cx_network)
        loader.__setattr__('_load_plan', ndexloadstring.get_load_plan())

        loader._load_style_template()

        ret_code = loader.load_to_NDEx()
        self.assertEqual(ret_code, 0)

    @unittest.skip("this test actually uses test_network.cx to upload and update it on server; we skip it")
    def test_0240_load_or_update_network_on_server(self):
        user_name = 'aaa'
        password = 'aaa'
        server = 'dev.ndexbio.org'

        loader = NDExSTRINGLoader(self._args)
        loader.__setattr__('_user', user_name)
        loader.__setattr__('_pass', password)
        loader.__setattr__('_server', server)


        nice_cx_path = ndexloadstring.get_package_dir() + '/../tests/test_network.cx'
        loader.__setattr__('_cx_network', nice_cx_path)

        # NDex client connection
        ndex_client = loader.create_ndex_connection()

        # delete all networks with the given name on NDEx server
        network_UUID = loader.get_network_uuid(self._network_name)
        while network_UUID and network_UUID != 2:
            ndex_client.delete_network(network_UUID)
            network_UUID = loader.get_network_uuid(self._network_name)


        # upload networks to the server 4 times
        loader._load_or_update_network_on_server('test network')
        loader._load_or_update_network_on_server('test network')
        loader._load_or_update_network_on_server('test network')
        loader._load_or_update_network_on_server('test network')

        network_UUID = loader.get_network_uuid(self._network_name)


        count = 0
        # now delete all the uploaded networks - there should be 4 of them
        while network_UUID and network_UUID != 2:
            ndex_client.delete_network(network_UUID)
            network_UUID = loader.get_network_uuid(self._network_name)
            count += 1

        self.assertEqual(count, 4)


        # now upload network, and update it (overwrite it) on the server
        time.sleep(1)
        loader._load_or_update_network_on_server('test network')
        time.sleep(1)
        network_UUID = loader.get_network_uuid(self._network_name)


        time.sleep(1)
        loader._load_or_update_network_on_server('test network', network_UUID)
        for i in range (0, 4):
            # over-write/update the existing network 4 times
            time.sleep(1)
            loader._load_or_update_network_on_server('test network', network_UUID)

        count = 0
        while network_UUID and network_UUID != 2:
            ndex_client.delete_network(network_UUID)
            network_UUID = loader.get_network_uuid(self._network_name)
            count += 1

        # there should only be one network this time
        self.assertEqual(count, 1)


        # try to get network UUID for a non-existant user; we expect to receive 2 from get_network_uuid()
        loader.__setattr__('_user', '_no_exists_')
        ret_code = loader.get_network_uuid('test network')
        self.assertEqual(ret_code, 2)

        # try to create network for non-existant user; we expect to receive 2 from _load_or_update_network_on_server()
        ndex_client.__setattr__('username', '_no_exists_')
        ret_code = loader._load_or_update_network_on_server('test network')
        self.assertEqual(ret_code, 2)
