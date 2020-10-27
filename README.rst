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

generates a new tsv file, transforms it to CX, and uploads it to NDEx server. Duplicate edges
(edges that have the same Source and Target nodes and the same value of :code:`combined_score`)
are included to the generated tsv and CX files only once. Name of the newly generated network includes
the value of :code:`cutoffscore` argument, for example,
:code:`STRING - Human Protein Links - High Confidence (Score >= 0.7)`. In case user didn't specify :code:`--update UUID`
argument, then the network with this name gets over-written in case if already exists on NDEx server;
otherwise, a new network is created.
Specifying :code:`--update UUID` command line argument will over-write network with this UUID if it is found.
If not, then user is asked if (s)he wants to create a new network. When network is updated, only edges and nodes are
changed; network attributes are not modified.


**1\)** Below is an example of a record
from `9606.protein.links.full.v11.0.txt.gz <https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz>`_

.. code-block::

   9606.ENSP00000261819 9606.ENSP00000353549 0 0 0 0 0 102 90 987 260 900 0 754 622 999


To generate a STRING network, the loader reads rows from that file one by one and compares the value of the last
column :code:`combined_score` with the value :code:`cutoffscore` argument.  The row is not added to the network generated in case
:code:`combined_score` is less than the commad-line argument :code:`cutoffscore`.


**2\)** If :code:`combined_score` is no than less :code:`cutoffscore`, the loader process two first columns

.. code-block::

   column 1 - protein1 (9606.ENSP00000261819)
   column 2 - protein2 (9606.ENSP00000353549)

When processing first column :code:`protein1`, the script

replaces :code:`Ensembl Id` with a :code:`display name`, for example :code:`9606.ENSP00000261819` becomes :code:`ANAPC5`. Mapping of
display names to :code:`Enseml Ids` is found in
`human.name_2_string.tsv.gz <https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz>`_

uses `human.uniprot_2_string.2018.tsv.gz <https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz>`_
to create :code:`represents` value.  For example, :code:`represents` for :code:`9606.ENSP00000261819` is :code:`uniprot:Q9UJX4`

uses `human.entrez_2_string.2018.tsv.gz <https://stringdb-static.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz>`_
to create list of aliases for the current protein.  Thus, list of aliases for :code:`9606.ENSP00000261819` is
:code:`ncbigene:51433|ensembl:ENSP00000261819`

**3\)** The second column :code:`protein2` is processed the same way as :code:`column 1`.

**4\)**  In the generated tsv file :code:`9606.protein.links.tsv`, :code:`protein1` and :code:`protein2` values from the original file are replaced with

.. code-block::

   protein_display_name_1 represents_1 alias_1 protein_display_name_2 represents_2 alias_2

So, the original

.. code-block::

   9606.ENSP00000261819 9606.ENSP00000353549 0 0 0 0 0 102 90 987 260 900 0 754 622 999

becomes

.. code-block::

   ANAPC5 uniprot:Q9UJX4 ncbigene:51433|ensembl:ENSP00000261819 CDC16 uniprot:Q13042  ncbigene:8881|ensembl:ENSP00000353549 0 0 0 0 0 102 90 987 260 900 0 754 622 999


**5\)**  The generated tsv file :code:`9606.protein.links.tsv` is then transformed to CX :code:`9606.protein.links.cx`.
The default style defined in :code:`style.cx` distributed with this loader is applied to the
generated network in case neither :code:`--style` nor :code:`--template` is specified.
User can specify style template file with either :code:`--style` argument or
style template network UUID :code:`--template UUID_of_style_template_network`.
Specifying both :code:`--template` and :code:`--style` is not allowed.

**6\)**  :code:`9606.protein.links.cx` is then uploaded to NDEx server either replacing
an existing network (in case :code:`--update UUID` is specified or network with this name already exists),
or creating a new network.


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

   git clone https://github.com/ndexcontent/ndexstringloader
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

    [ndexstringloader]
    user = joe123
    password = somepassword123
    server = dev.ndexbio.org


Needed files
------------

Load plan is required for running this script.  **string_plan.json**  found at **ndexstringloader/ndexstringloader** can be used for this purpose.


Usage
-----

For information invoke :code:`ndexloadstring.py -h`

**Example usage**

Here is how this command can be run for **dev** and **prod** targets:

.. code-block::

   ndexloadstring.py --profile dev tmpdir/

   ndexloadstring.py --profile prod tmpdir/


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
