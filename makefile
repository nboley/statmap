
all: 
	cd src; make;
	cd utilities; make;
clean:
	cd src; make clean;
	cd utilities; make clean;

check: 
	cd tests; rm -rf smo_* *.fa *.map; python tests.py && python simulate_chipseq.py;
