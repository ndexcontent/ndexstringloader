#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ndexstringloader` package."""

import os
import tempfile
import shutil

import unittest
from ndexutil.config import NDExUtilConfig
from ndexstringloader.ndexloadstring import NDExSTRINGLoader
import ndexstringloader
from ndexstringloader import ndexloadstring

import ndex2


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
        args.append('--datadir')
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
            '9606.ENSP00000000233 9606.ENSP00000373713 0 0 0 0 0 0 63 0 63 0 0 0 237 271'
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
            '9606.ENSP00000373713': { 'display_name': None, 'alias': 'ncbigene:22908|ensembl:ENSP00000373713', 'represents': None }
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
            '9606\t22908\t9606.ENSP00000373713'
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
        self.assertEqual(ensembl_ids_expected, loader.__getattribute__('ensembl_ids'))



    #@unittest.skip("skip it  now - uncomment later")
    def test_0900_download_and_unzip(self):

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


    #@unittest.skip("skip it  now - uncomment later")
    def test_0910_download_and_unzip_STRING_files(self):

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


