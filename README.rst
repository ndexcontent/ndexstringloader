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

* `ndex2 <https://pypi.org/project/ndex2>`_
* `ndexutil <https://pypi.org/project/ndexutil>`_

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

The **ndexloadstring.py** requires a configuration file in the following format be created.
The default path for this configuration is :code:`~/.ndexutils.conf` but can be overridden with
:code:`--conf` flag.

**Format of configuration file**

.. code-block::

    [<value in --profile (default ndexstringloader)>]

    user = <NDEx username>
    password = <NDEx password>
    server = <NDEx server(omit http) ie public.ndexbio.org>
    style = <NDEx UUID of network to use for styling networks created>


The NDEx UUID needed for **style** can be obtained by uploading the :code:`style.cx` file found under
the :code:`data/` directory of this repository. NOTE: The network needs to be uploaded to the same
server as defined in **style** :code:`public.ndexbio.org` is NDEx_ production. Also the network needs
to be visible to the **user**

**Example configuration file**

.. code-block::

    [ndexstringloader_dev]

    user = joe123
    password = somepassword123
    server = dev.ndexbio.org
    style = 86f63bf8-1b48-11e9-a05d-525400c25d22


Needed files
------------

**TODO:** Add description of needed files


Usage
-----

For information invoke :code:`ndexloadstring.py -h`

**Example usage**

**TODO:** Add information about example usage

.. code-block::

   ndexloadstring.py # TODO Add other needed arguments here


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
