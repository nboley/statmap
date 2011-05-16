python large_scale_sim.py data/Dnase_2L.wig data/chr2L_short.fa 21146691 chr2l_shrt fl_dist.txt new_test
mkdir test_run; mv data/new_test* test_run; cd test_run;
../../../bin/statmap -g ../data/chr2L_short.fa -r new_test_IP.fastq -c new_test_CONTROL.fastq -o mapped -f ../fl_dist.txt -n 10 -a i
