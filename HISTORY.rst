=======
History
=======

1.0.0 (11-11-2020)
------------------

* New default behavior: **force-directed-cl** layout is now applied on
  networks via py4cytoscape library and a running instance of Cytoscape.
  Alternate Cytoscape layouts and the networkx "spring" layout can be
  run by setting appropriate value via the new **--layout** flag

0.3.0 (10-28-2020)
------------------

* Added ``--skipupload`` that lets caller skip upload of network to NDEx

* Spring layout applied by default for all networks that have less then 2,000,000
  edges. This can be overridden with new flag ``--layoutedgecutoff``

0.2.4 (2019-12-01)
------------------
* Fixed defect UD-462 Verify new network attributes are correctly set in ndexstringloader (https://ndexbio.atlassian.net/browse/UD-462).

0.2.3 (2019-09-13)
------------------
* If user loads the entire STRING network (i.e., runs the script with --cutoffscore 0), the name of the resulting netwpork should be "STRING - Human Protein Links", not "STRING - Human Protein Links - High Confidence".

0.2.2 (2019-09-12)
------------------
* Added new featured specified by UD-577 Quick improvement for new String loader (added optional --update argument that allows to specify the UUID of a target network to update; added optional --template argument that allows to specify the UUID of a target network to use as style template, the update operation now only changes nodes and edges, but leaves network properties untouched).

0.2.1 (2019-08-23)
------------------
* Improved README file.
* Added new JUnit tests (JUnit test coverage is 87%).

0.2.0 (2019-07-26)
------------------
* Removed duplicate edges. Every pair of connected nodes in STRING networks had the same edge duplicated (one edge going from A to B, and another going from B to A).  Since edges in STRING are not directed, we can safely remove half of them.

* Added new arguments to command line:
   optional --cutoffscore (default is 0.7) - used to filter on combined_score column. To include edges with combined_score of 800 or higher, --cutoffscore 0.8 should be specified

   required --datadir specifies a working directory where STRING files will be downloaded to and processed style.cx file that contains style is supplied with the STRING loader and used by default. It can be overwritten with --style argument.

0.1.0 (2019-03-13)
------------------
* First release on PyPI.
