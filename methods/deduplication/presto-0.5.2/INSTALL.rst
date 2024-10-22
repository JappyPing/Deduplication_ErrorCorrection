Installation
================================================================================

The latest stable release of pRESTO may be downloaded from
`PyPI <https://pypi.python.org/pypi/presto>`__ or
`Bitbucket <https://bitbucket.org/kleinstein/presto/downloads>`__.
The simplest way to install the latest stable release of pRESTO is via pip::

    > pip3 install presto

Requirements
--------------------------------------------------------------------------------

-  `Python 3.4.0 <http://python.org>`__
-  `setuptools 2.0 <http://bitbucket.org/pypa/setuptools>`__
-  `NumPy 1.8 <http://numpy.org>`__
-  `SciPy 0.14 <http://scipy.org>`__
-  `pandas 0.15 <http://pandas.pydata.org>`__
-  `Biopython 1.65 <http://biopython.org>`__
-  `MUSCLE v3.8 <http://www.drive5.com/muscle>`__
-  `USEARCH v7.0 <http://www.drive5.com/usearch>`__

Linux
--------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the
   full SciPy stack using the
   `instructions <http://scipy.org/install.html>`__, then install
   Biopython according to its
   `instructions <http://biopython.org/DIST/docs/install/Installation.html>`__.

2. Extract the pRESTO bundle and run::

   > python3 setup.py install --user

Mac OS X
--------------------------------------------------------------------------------

1. Install Xcode. Available from the Apple store or
   `developer downloads <http://developer.apple.com/downloads>`__.

2. Older versions Mac OS X will require you to install XQuartz 2.7.5. Available
   from the `XQuartz project <http://xquartz.macosforge.org/landing>`__.

3. Install Homebrew following the installation and post-installation
   `instructions <http://brew.sh>`__.

4. Install Python 3.4.0+ and set the path to the python3 executable::

   > brew install python3
   > echo 'export PATH=/usr/local/bin:$PATH' >> ~/.profile

5. Exit and reopen the terminal application so the PATH setting takes effect.

6. You may, or may not, need to install gfortran (required for SciPy). Try
   without first, as this can take an hour to install and is not needed on
   newer releases. If you do need gfortran to install SciPy, you can install it
   using Homebrew::

   > brew install gfortran

   If the above fails run this instead::

   > brew install --env=std gfortran

7. Install NumPy, SciPy, pandas and Biopyton using the Python package
   manager::

   > pip3 install numpy scipy pandas biopython

8. Extract the pRESTO bundle, open a terminal window, change directories
   to the pRESTO folder, and run::

   > python3 setup.py install

Windows
--------------------------------------------------------------------------------

1. Install Python 3.4.0+ from `Python <http://python.org/downloads>`__,
   selecting both the options 'pip' and 'Add python.exe to Path'.

2. Install NumPy, SciPy, pandas and Biopython using the packages
   available from the
   `Unofficial Windows binary <http://www.lfd.uci.edu/~gohlke/pythonlibs>`__
   collection.

3. Unzip the pRESTO bundle, open a Command Prompt, change directories to
   the pRESTO folder, and run::

   > python setup.py install

4. For a default installation of Python 3.4, the pRESTO scripts will be
   installed into ``C:\Python34\Scripts`` and should be directly
   executable from the Command Prompt. If this is not the case, then
   follow step 5 below.

5. Add both the ``C:\Python34`` and ``C:\Python34\Scripts`` directories
   to your ``%Path%``. On Windows 7 the ``%Path%`` setting is located
   under Control Panel -> System and Security -> System -> Advanced
   System Settings -> Environment variables -> System variables -> Path.
