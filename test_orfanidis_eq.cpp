/*
 * Copyright (c) 2018 Fedor Uporov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "orfanidis_eq.h"

using namespace std;
using namespace OrfanidisEq;

static const int classConversionsRangeDb = eqGainRangeDb;

static const double defaultMaxGain = 1.0;
static const double defaultMinGain = 0.0;
static const double defaultUnitImpAmp = 1.0;

/*
 * Provide auxiliary functions for EQ testing.
 */
class TestOrfanidisEq {
public:
	TestOrfanidisEq(){}
	~TestOrfanidisEq(){}

	bool setFreqGrid(FrequencyGrid& fg, unsigned int numberOfBands)
	{
		switch(numberOfBands) {
			case 5:
				fg.set5Bands(bandsGridCenterFreqHz);
				break;

			case 10:
				fg.set10Bands(bandsGridCenterFreqHz);
				break;

			case 20:
				fg.set20Bands(bandsGridCenterFreqHz);
				break;

			case 30:
				fg.set30Bands(bandsGridCenterFreqHz);
				break;

			default:
				cout <<
				    "Can not configue freq grid, provide 5,10,20 or 30."
				    << "Provided : " << numberOfBands << endl;

				return true;
		}

		return false;
	}

	void setUnitImpulse(vector<eq_double_t>& vect)
	{
		vect[0] = defaultUnitImpAmp;
	}

	void processEq(Eq& equalizer, vector<eq_double_t> &in,
	    vector<eq_double_t> &out)
	{
		for (unsigned int i = 0; i < in.size(); i++)
			equalizer.SBSProcess(&in[i], &out[i]);
	}

	void saveCSFile(string fileName, vector<eq_double_t>& testVector)
	{
		/* Write to file. */
		ofstream dataFile;

		dataFile.open(fileName.c_str());
		if (!dataFile.is_open())
			return;

		for (unsigned int j = 0; j < testVector.size() - 1; j++)
			dataFile << testVector[j] << ", ";

		dataFile << testVector[testVector.size() - 1] << endl;

		dataFile.close();
	}
};

void classConversionsTest()
{
	Conversions convs(classConversionsRangeDb);

	cout << "1. conversions class: fast conversions test : " << endl;
	cout << "db -> lin -> fast_lin" << endl;

	eq_double_t testRangeDb = classConversionsRangeDb + 3;
	eq_double_t refDb = 0;

	for (refDb = -testRangeDb; refDb < testRangeDb; refDb += 6.4)
		cout << refDb << " " << convs.db2Lin(refDb) << " " <<
		    convs.fastDb2Lin(refDb) << endl;

	cout << "lin -> db -> fast_db" << endl;

	eq_double_t refLin = 0;
	for (refLin = 0.00316 /* -50 db */; refLin < 316 /* 50 db */;
	    refLin*=convs.db2Lin(6))
		cout << refLin << " " << convs.lin2Db(refLin) << " " <<
		    convs.fastLin2Db(refLin) << endl;

	cout << "----" << endl;
}

void classFreqGridTest()
{
	cout << "2. freq_grid class: print freqs for 1/3 octave eq:" << endl;

	FrequencyGrid fg;
	fg.set30Bands();

	cout << "band freq : rounded band freq" << endl;

	for (unsigned int i = 0; i < fg.getNumberOfBands(); i++)
	    cout << fg.getFreq(i) << " : " << fg.getRoundedFreq(i) << endl;
}

void printInputConfigurationData(int sampleRate, int testVectorLength,
    int numberOfBands, const vector<double>& gains)
{
	cout << "SampleRate = " << sampleRate << " " << "Hz." << endl;
	cout << "Input vector size = " << testVectorLength << " " << "samples" << endl;
	cout << "Number of bands = " << numberOfBands << endl;

	cout << "Gains data : ";
	for (int i = 0; i < numberOfBands; i++)
		cout << gains[i] << ", ";

	cout << "dB" << endl;
}

void classEqTest(int sampleRate, int testVectorLength, int numberOfBands,
    const vector<double>& gains)
{
	printInputConfigurationData(sampleRate, testVectorLength, numberOfBands, gains);

	TestOrfanidisEq testEq;

	/* Input configuration. */
	FrequencyGrid fg;
	testEq.setFreqGrid(fg, numberOfBands);
	fg.printFrequencyGrid();

	/* Input data vector. */
	vector<eq_double_t> inVector(testVectorLength, 0);
	testEq.setUnitImpulse(inVector);

	Eq equalizer(fg, none);

	vector<eq_double_t> butterworthOut(testVectorLength, 0);
	equalizer.setEq(fg, butterworth);
	equalizer.setSampleRate(sampleRate);
	equalizer.changeGainsDb(gains);
	testEq.processEq(equalizer, inVector, butterworthOut);
	testEq.saveCSFile("butterworth.tstdat", butterworthOut);

	vector<eq_double_t> chebyshev1Out(testVectorLength, 0);
	equalizer.setEq(fg, chebyshev1);
	equalizer.setSampleRate(sampleRate);
	equalizer.changeGainsDb(gains);
	testEq.processEq(equalizer, inVector, chebyshev1Out);
	testEq.saveCSFile("chebyshev1.tstdat", chebyshev1Out);

	vector<eq_double_t> chebyshev2Out(testVectorLength, 0);
	equalizer.setEq(fg, chebyshev2);
	equalizer.setSampleRate(sampleRate);
	equalizer.changeGainsDb(gains);
	testEq.processEq(equalizer, inVector, chebyshev2Out);
	testEq.saveCSFile("chebyshev2.tstdat", chebyshev2Out);

	vector<eq_double_t> ellipticOut(testVectorLength, 0);
	equalizer.setEq(fg, elliptic);
	equalizer.setSampleRate(sampleRate);
	equalizer.changeGainsDb(gains);
	testEq.processEq(equalizer, inVector, ellipticOut);
	testEq.saveCSFile("elliptic.tstdat", ellipticOut);

}

void usage()
{
	cout << "usage keys:" << endl;
	cout << "-f - sampling frequency in Hz" << endl;
	cout << "-s - input vector size" << endl;
	cout << "-b - number of bands" << endl;
	cout << "-g - list of gains per band in db" << endl;
	cout << "NOTE: -g list length should be equal number of bands (-b option)" << endl;
	cout << "Example:" << endl;
	cout << "./eq -f 48000 -s 10000 -b 30 -g 0 0 0 0 3 0 5 0 0 16 0 0 0 0 0 0 0 0 0 0 0 -10 0 0 -5 0 0 -2 0 0" << endl;

	exit(1);
}

int main(int argc, char *argv[])
{
	int sampleRateHz;
	int vectorSize;
	int numberOfBands;
	vector<double> bandsGains;

	if (argc < 8)
		usage();

	/* Get sample rate in Hz. */
	if (strcmp(argv[1], "-f"))
		usage();

	sampleRateHz = atoi(argv[2]);
	if (sampleRateHz < 0 || sampleRateHz > 192000)
		usage();

	/* Get input vector size. */
	if (strcmp(argv[3], "-s"))
		usage();

	vectorSize = atoi(argv[4]);
	if (vectorSize < 0 || vectorSize > 192000)
		usage();

	/* Get number of bands. */
	if (strcmp(argv[5], "-b"))
		usage();

	numberOfBands = atoi(argv[6]);
	if (numberOfBands != 30) {
		cout << "ERROR: Only 30 bands eq now implemented" << endl;
		usage();
	}

	/* Get bands gains. */
	if (strcmp(argv[7], "-g") || argc < 7 + numberOfBands)
		usage();

	for (int i = 0, j = 8; i < numberOfBands; i++)
		bandsGains.push_back(atof(argv[j++]));

	classEqTest(sampleRateHz, vectorSize, numberOfBands, bandsGains);

	return 0;
}
