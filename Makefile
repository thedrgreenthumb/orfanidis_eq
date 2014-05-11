all:
	g++ -o eq test_orfanidis_eq.cpp -I. 

clean:
	rm -f *.tstdat
	rm -f eq	
