all: 
	cd src; make;
clean:
	cd src; make clean;

check: 
	cd tests; python tests.py