=======
History
=======


0.2.0 (2019-07-26)
------------------
* Removed duplicate edges. Every pair of connected nodes in STRING networks had the same edge duplicated (one edge going from A to B, and another going from B to A).  Since edges in STRING are not directed, we can safely remove half of them.

* Added new arguments to command line:
   optional --cutoffscore (default is 0.7) - used to filter on combined_score column. To include edges with combined_score of 800 or higher, --cutoffscore 0.8 should be specified

   required --datadir specifies a working directory where STRING files will be downloaded to and processed style.cx file that contains style is supplied with the STRING loader and used by default. It can be overwritten with --style argument.


* First release on PyPI.

0.1.0 (2019-03-13)
------------------

* First release on PyPI.
