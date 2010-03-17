all: 
	cd src; make;

clean:
	cd src; make clean;

check: 
	cd unit_tests; python unit_tests.py