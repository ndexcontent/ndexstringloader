==========================
NDEx STRING Content Loader
==========================


.. image:: https://img.shields.io/pypi/v/ndexstringloader.svg
        :target: https://pypi.python.org/pypi/ndexstringloader

.. image:: https://img.shields.io/travis/vrynkov/ndexstringloader.svg
        :target: https://travis-ci.org/ndexcontent/ndexstringloader

.. image:: https://coveralls.io/repos/github/ndexcontent/ndexstringloader/badge.svg?branch=master
        :target: https://coveralls.io/github/ndexcontent/ndexstringloader?branch=master

.. image:: https://readthedocs.org/projects/ndexstringloader/badge/?version=latest
        :target: https://ndexstringloader.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Python application for loading STRING data into `NDEx <http://ndexbio.org>`_.

This tool downloads and unpacks the `STRING <https://string-db.org/>`_ files below

    `9606.protein.links.full.v11.0.txt.gz <https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz>`_

    `human.name_2_string.tsv.gz <https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz>`_

    `human.entrez_2_string.2018.tsv.gz <https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz>`_

    `human.uniprot_2_string.2018.tsv.gz <https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz>`_

generates a new tsv file, transforms it to CX, and uploads it to NDEx server.



**1\)** Below is an example of a record
from `9606.protein.links.full.v11.0.txt.gz <https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz>`_

.. code-block::

9606.ENSP00000261819 9606.ENSP00000353549 0 0 0 0 0 102 90 987 260 900 0 754 622 999


To generate a STRING network, the loader reads rows from that file one by one and compares the value of the last
column combined_score with the value cutoffscore argument.  The line is not added to the network generated in case
combined_score is less than cutoffscore.


**2\)** If combined_score is no than less cutoffscore, the loader process two first columns

column 1 - protein1 (9606.ENSP00000261819)
column 2 - protein2 (9606.ENSP00000353549)

When processing first column protein1, the script

replaces Ensembl Id with a display name, for example 9606.ENSP00000261819 becomes ANAPC5. Mapping of
display names to Enseml Ids is found in
`human.name_2_string.tsv.gz <https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz>`_

uses `human.uniprot_2_string.2018.tsv.gz <https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz>`_
to create represents value.  For example, represents for 9606.ENSP00000261819 is uniprot:Q9UJX4

uses `human.entrez_2_string.2018.tsv.gz <https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz>`_
to create list of aliases for the current protein.  Thus, list of aliases for 9606.ENSP00000261819 is
ncbigene:51433|ensembl:ENSP00000261819

**3\)** The second column protein2 is processed the same way as column 1.

**4\)**  In the generated tsv file, protein1 and protein2 values from the original file are replaced with

protein_display_name_1 represents_1 alias_1 protein_display_name_2 represents_2 alias_2


Dependencies
------------

* ndex2
* ndexutil

Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block::

   git clone https://github.com/vrynkov/ndexstringloader
   cd ndexstringloader
   make dist
   pip install dist/ndexloadstring*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image 

.. code-block::

   make

Output:

.. code-block::

   clean                remove all build, test, coverage and Python artifacts
   clean-build          remove build artifacts
   clean-pyc            remove Python file artifacts
   clean-test           remove test and coverage artifacts
   lint                 check style with flake8
   test                 run tests quickly with the default Python
   test-all             run tests on every Python version with tox
   coverage             check code coverage quickly with the default Python
   docs                 generate Sphinx HTML documentation, including API docs
   servedocs            compile the docs watching for changes
   testrelease          package and upload a TEST release
   release              package and upload a release
   dist                 builds source and wheel package
   install              install the package to the active Python's site-packages
   dockerbuild          build docker image and store in local repository
   dockerpush           push image to dockerhub


Configuration
-------------

The **ndexloadstring.py** requires a configuration file to be created.
The default path for this configuration is :code:`~/.ndexutils.conf` but can be overridden with
:code:`--conf` flag.

**Configuration file**

Networks listed in **[network_ids]** section need to be visible to the **user**

.. code-block::

    [dev]
    user = joe123 
    password = somepassword123 
    server = dev.ndexbio.org
    hi_confidence = 311b0e5f-6570-11e9-8c69-525400c25d22
    ProteinLinksFile = https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz
    NamesFile = https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
    EntrezIdsFile = https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz
    UniprotIdsFile = https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz
    full_file_name = 9606.protein.links.full.v11.0.txt
    entrez_file = human.entrez_2_string.2018.tsv
    names_file = human.name_2_string.tsv
    uniprot_file = human.uniprot_2_string.2018.tsv
    output_tsv_file_name = 9606.protein.links.full.v11.0.tsv.txt
    output_hi_conf_tsv_file_name = 9606.protein.links.full.v11.0.hi_conf.tsv.txt

    [prod]
    user = joe123 _prod
    password = somepassword123_prod 
    server = prod.ndexbio.org
    hi_confidence = 311b0e5f-6570-11e9-8c69-525400c25d22
    ProteinLinksFile = https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz
    NamesFile = https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
    EntrezIdsFile = https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz
    UniprotIdsFile = https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz
    full_file_name = 9606.protein.links.full.v11.0.txt
    entrez_file = human.entrez_2_string.2018.tsv
    names_file = human.name_2_string.tsv
    uniprot_file = human.uniprot_2_string.2018.tsv
    output_tsv_file_name = 9606.protein.links.full.v11.0.tsv.txt
    output_hi_conf_tsv_file_name = 9606.protein.links.full.v11.0.hi_conf.tsv.txt


Needed files
------------

Load plan is required for running this script.  **string_plan.json**  found at **ndexstringloader/ndexstringloader** can be used for this purpose.


Usage
-----

For information invoke :code:`ndexloadstring.py -h`

**Example usage**

Here is how this command can be run for **dev** and **prod** targets:

.. code-block::

   ndexloadstring.py --profile dev

   ndexloadstring.py --profile prod


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

**TODO:** Add information about example usage


.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` vrynkov/ndexstringloader:0.1.0 ndexloadstring.py --conf conf # TODO Add other needed arguments here


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
