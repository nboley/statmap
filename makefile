all: 
	cd src; make;
	cd utilities; make;
clean:
	cd src; make clean;
	cd utilities; make clean;

check: 
	cd tests; python tests.py