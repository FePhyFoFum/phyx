Instructions for installing dependencies
========================================

Phyx requires the dependencies NLopt and Armadillo. Armadillo itself has some optional dependencies which may provide faster execution. 

### NLopt
For the time being, please use the provided source code, regardless of system. Type:

	tar -xvzf nlopt-2.4.2.tar.gz
	cd nlopt-2.4.2
	./configure --with-cxx --without-octave --without-matlab
	make
	sudo make install

LINUX
---------------
Install using existing repositories. Type the following command into a terminal prompt:

	sudo apt-get install liblapack-dev
	sudo apt-get install libatlas-cpp-0.6-dev
	sudo apt-get install libarmadillo-dev

MAC
---------------
Make sure you have compilers (e.g. from XCode, homebrew, or MacPorts) installed.
Installation will be from source code (provided in the deps directory).

### LAPACK (optional) - this can take a while

	tar -xvzf lapack-3.4.2.tgz
	cd lapack-3.4.2
	cp make.inc.example make.inc
	make blaslib
	make

### ATLAS (optional)

	tar -xvjf atlas3.10.2.tar.bz2
	cd ATLAS
	mkdir build
	cd build
	../configure
	make build
	make check
	make time
	sudo make install

### Armadillo

	tar -xvzf armadillo-7.400.2.tar.gz
	cd armadillo-7.400.2
	cmake .
	make
	sudo make install

WINDOWS
---------------
We do not have frequent access to Windows computers, but will help anyway we can.
