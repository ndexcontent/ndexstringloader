==========================
NDEx STRING Content Loader
==========================


.. image:: https://img.shields.io/pypi/v/ndexstringloader.svg
        :target: https://pypi.python.org/pypi/ndexstringloader

.. image:: https://img.shields.io/travis/vrynkov/ndexstringloader.svg
        :target: https://travis-ci.org/vrynkov/ndexstringloader

.. image:: https://readthedocs.org/projects/ndexstringloader/badge/?version=latest
        :target: https://ndexstringloader.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Python Boilerplate contains all the boilerplate you need to create a Python NDEx Content Loader package.


* Free software: BSD license
* Documentation: https://ndexstringloader.readthedocs.io.


Tools
-----

* **ndexloadstring.py** -- Loads STRING into NDEx_

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
    
    [prod]
    user = joe123 _prod
    password = somepassword123_prod 
    server = prod.ndexbio.org 
    
    [network_ids]
    style = 9c23e193-5d73-11e9-8c69-525400c25d22
    hi_confidence_style = 9c100b71-5d73-11e9-8c69-525400c25d22 
    full = a57b23c5-65fe-11e9-8c69-525400c25d22 
    hi_confidence = 311b0e5f-6570-11e9-8c69-525400c25d22 

    [source]
    ProteinLinksFile = https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz
    NamesFile = https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
    EntrezIdsFile = https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz
    UniprotIdsFile = https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz

    [input]
    full_file_name = 9606.protein.links.full.v11.0.txt
    entrez_file = human.entrez_2_string.2018.tsv
    names_file = human.name_2_string.tsv
    uniprot_file = human.uniprot_2_string.2018.tsv

    [output]
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

   ndexloadstring.py --loadplan loadplan.json  --profile dev 

   ndexloadstring.py --loadplan loadplan.json  --profile prod 


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
