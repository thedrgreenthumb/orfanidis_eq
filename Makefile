all:
	g++ -o eq test_orfanidis_eq.cpp -I. 

clean:
	rm -f orfanidis_eq.tstdat
	rm -f butterworth_*
	rm -f chebyshev1_*
	rm -f chebyshev2_*
	rm -f eq	
