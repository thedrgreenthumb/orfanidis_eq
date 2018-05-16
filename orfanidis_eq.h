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

#pragma once

#include <cmath>
#include <vector>

namespace OrfanidisEq {

/*
 * Just version of implementation
 * for custom purposes.
 */
static const char* eq_version = "0.02";

/*
 * The float type usage here
 * could be cause the lack of precession.
 */
typedef double eq_double_t;

/*
 * Eq configuration constants.
 * The defaultEqBandPassFiltersOrder should be more then 2.
 */
static const eq_double_t defaultSampleFreqHz = 48000;
static const unsigned int defaultEqBandPassFiltersOrder = 4;

/* Default frequency values to get frequency grid. */
static const eq_double_t lowestGridCenterFreqHz = 31.25;
static const eq_double_t bandsGridCenterFreqHz = 1000;
static const eq_double_t lowestAudioFreqHz = 20;
static const eq_double_t highestAudioFreqHz = 20000;

/* Default gain values for every type of filter channel. */
static const eq_double_t eqGainRangeDb = 40;
static const eq_double_t eqGainStepDb = 1;
static const eq_double_t commonBaseGainDb = 3;
static const eq_double_t eqDefaultGainDb = 0;

/*
 * Allow to convert values between different representations.
 * Also, fast conversions between linear and logarithmic scales are included.
 */
class Conversions {
	int rangeDb;
	std::vector<eq_double_t> linGains;

	int linGainsIndex(eq_double_t x)
	{
		int inTx = (int)x;

		if ((x >= -rangeDb) && (x < rangeDb - 1))
			return rangeDb + inTx;

		return rangeDb;
	}

	Conversions();
public:
	Conversions(int range)
	{
		int step = -range;

		rangeDb = range;

		while (step <= range)
			linGains.push_back(db2Lin(step++));
	}

	eq_double_t fastDb2Lin(eq_double_t x)
	{
		int integerPart = (int)x;
		eq_double_t fracPart = x - integerPart;

		return linGains[linGainsIndex(integerPart)]*(1-fracPart) +
		    (linGains[linGainsIndex(integerPart + 1)])*fracPart;
	}

	eq_double_t fastLin2Db(eq_double_t x)
	{
		if ((x >= linGains[0]) && (x < linGains[linGains.size() - 1])) {
			for (unsigned int i = 0; i < linGains.size() - 2; i++) {
				if ((x >= linGains[i]) && (x < linGains[i + 1])) {
					int integerPart = i - rangeDb;
					eq_double_t fracPart = x - (int)(x);

					return integerPart + fracPart;
				}
			}
		}

		return 0;
	}

	static eq_double_t db2Lin(eq_double_t x)
	{
		return pow(10, x/20);
	}

	static eq_double_t lin2Db(eq_double_t x)
	{
		return 20*log10(x);
	}

	static eq_double_t rad2Hz(eq_double_t x, eq_double_t fs)
	{
		return 2*M_PI/x*fs;
	}

	static eq_double_t hz2Rad(eq_double_t x, eq_double_t fs)
	{
		return 2*M_PI*x/fs;
	}
};

/*
 * Filter frequency band representation.
 */
class Band {
	Band();
public:
	eq_double_t minFreq;
	eq_double_t centerFreq;
	eq_double_t maxFreq;

	Band(eq_double_t fl, eq_double_t fc, eq_double_t fh)
		: minFreq(fl), centerFreq(fc), maxFreq(fh) {}
};

/* Basic eq errors handling. */
typedef enum {
	no_error,
	invalid_input_data_error,
	processing_error
} eq_error_t;

/*
 * Frequency grid representation.
 * Allow to calculate all required frequences for bandpass filters
 * for equalizer construction.
 */
class FrequencyGrid {
	std::vector<Band> freqs;

public:
	FrequencyGrid()
	{}

	eq_error_t setBand(eq_double_t fl, eq_double_t fc, eq_double_t fh)
	{
		freqs.clear();
		return addBand(fl, fc, fh);
	}

	eq_error_t addBand(eq_double_t fl, eq_double_t fc, eq_double_t fh)
	{
		if (fl < fc && fc < fh) {
			freqs.push_back(Band(fl, fc, fh));
			return no_error;
		}

		return invalid_input_data_error;
	}

	eq_error_t addBand(eq_double_t fc, eq_double_t df)
	{
		if (fc >= df /2 ) {
			eq_double_t fmin = fc - df / 2;
			eq_double_t fmax = fc + df / 2;
			freqs.push_back(Band(fmin, fc, fmax));
			return no_error;
		}

		return invalid_input_data_error;
	}

	eq_error_t set5Bands(eq_double_t fc = bandsGridCenterFreqHz)
	{
		freqs.clear();

		if (lowestAudioFreqHz < fc && fc < highestAudioFreqHz) {
			/* Find lowest center frequency in the band. */
			eq_double_t lowestFc = fc;

			while (lowestFc > lowestGridCenterFreqHz)
				lowestFc /= 4.0;

			if (lowestFc < lowestGridCenterFreqHz)
				lowestFc *= 4.0;

			/* Calculate frequences. */
			eq_double_t f0 = lowestFc;
			for (unsigned int i = 0; i < 5 ; i++) {
				freqs.push_back(Band(f0 / 2, f0, f0 * 2));
				f0 *= 4;
			}

			return no_error;
		}

		return invalid_input_data_error;
	}

	eq_error_t set10Bands(eq_double_t fc = bandsGridCenterFreqHz)
	{
		freqs.clear();

		if (lowestAudioFreqHz < fc && fc < highestAudioFreqHz) {
			/* Find lowest center frequency in the band. */
			eq_double_t lowestFc = fc;

			while (lowestFc > lowestGridCenterFreqHz)
				lowestFc /= 2;
	
			if (lowestFc < lowestGridCenterFreqHz)
				lowestFc *= 2;

			/* Calculate frequences. */
			eq_double_t f0 = lowestFc;
			for (unsigned int i = 0; i < 10; i++) {
				freqs.push_back(
				    Band(f0 / pow(2, 0.5), f0, f0 * pow(2, 0.5)));

				f0 *= 2;
			}

			return no_error;
		}

		return invalid_input_data_error;
	}

	eq_error_t set20Bands(eq_double_t fc = bandsGridCenterFreqHz)
	{
		freqs.clear();

		if (lowestAudioFreqHz < fc && fc < highestAudioFreqHz) {
			/* Find lowest center frequency in the band. */
			eq_double_t lowestFc = fc;

			while (lowestFc > lowestAudioFreqHz)
				lowestFc /= pow(2, 0.5);

			if (lowestFc < lowestAudioFreqHz)
				lowestFc *= pow(2, 0.5);

			/* Calculate frequences. */
			eq_double_t f0 = lowestFc;
			for (unsigned int i = 0; i < 20; i++) {
				freqs.push_back(Band(f0 / pow(2, 0.25),
				    f0, f0 * pow(2, 0.25)));

				f0 *= pow(2, 0.5);
			}

			return no_error;
		}

		return invalid_input_data_error;
	}

	eq_error_t set30Bands(eq_double_t fc = bandsGridCenterFreqHz)
	{
		freqs.clear();

		if (lowestAudioFreqHz < fc && fc < highestAudioFreqHz) {
			/* Find lowest center frequency in the band. */
			eq_double_t lowestFc = fc;
			while (lowestFc > lowestAudioFreqHz)
				lowestFc /= pow(2.0, 1.0/3.0);

			if (lowestFc < lowestAudioFreqHz)
				lowestFc *= pow(2.0, 1.0/3.0);

			/* Calculate frequences. */
			eq_double_t f0 = lowestFc;
			for (unsigned int i = 0; i < 30; i++) {
				freqs.push_back(Band(f0 / pow(2.0, 1.0/6.0),
				    f0, f0 * pow(2.0, 1.0/6.0)));

				f0 *= pow(2,1.0/3.0);
			}

			return no_error;
		}

		return invalid_input_data_error;
	}

	unsigned int getNumberOfBands()
	{
		return freqs.size();
	}

	std::vector<Band> getFreqs()
	{
		return freqs;
	}

	unsigned int getFreq(unsigned int index)
	{
		if (index < freqs.size())
			return freqs[index].centerFreq;
		else
			return 0;
	}

	unsigned int getRoundedFreq(unsigned int index)
	{
		if (index < freqs.size()) {
			unsigned int freq = freqs[index].centerFreq;
			if (freq < 100) {
				return freq;
			} else if (freq >= 100 && freq < 1000) {
				unsigned int rest = freq % 10;
				if (rest < 5)
					return freq - rest;
				else
					return freq - rest + 10;
			} else if (freq >= 1000 && freq < 10000) {
				unsigned int rest = freq % 100;
				if (rest < 50)
					return freq - rest;
				else
					return freq - rest + 100;
			} else if (freq >= 10000) {
				unsigned int rest = freq%1000;
				if (rest < 500)
					return freq - rest;
				else
					return freq - rest + 1000;
			}
		}

		return 0;
	}
};

/*
 * Fourth order biquad section representation.
 */
class FOSection {
protected:
	eq_double_t b0, b1, b2, b3, b4;
	eq_double_t a0, a1, a2, a3, a4;

	eq_double_t numBuf[4];
	eq_double_t denumBuf[4];

	eq_double_t df1FOProcess(eq_double_t in)
	{
		eq_double_t out = 0;

		out+= b0*in;
		out+= (b1*numBuf[0] - denumBuf[0]*a1);
		out+= (b2*numBuf[1] - denumBuf[1]*a2);
		out+= (b3*numBuf[2] - denumBuf[2]*a3);
		out+= (b4*numBuf[3] - denumBuf[3]*a4);

		numBuf[3] = numBuf[2];
		numBuf[2] = numBuf[1];
		numBuf[1] = numBuf[0];
		*numBuf = in;

		denumBuf[3] = denumBuf[2];
		denumBuf[2] = denumBuf[1];
		denumBuf[1] = denumBuf[0];
		*denumBuf = out;

		return(out);
	}

public:
	FOSection()
	{
		b0 = 1; b1 = 0; b2 = 0; b3 = 0; b4 = 0;
		a0 = 1; a1 = 0; a2 = 0; a3 = 0; a4 = 0;

		for (unsigned int i = 0; i < 4; i++) {
			numBuf[i] = 0;
			denumBuf[i] = 0;
		}
	}

	virtual ~FOSection() {}

	virtual FOSection get()
	{
		return *this;
	}

	eq_double_t process(eq_double_t in)
	{
		return df1FOProcess(in);
	}
};

class ButterworthFOSection : public FOSection {
	ButterworthFOSection() {}
	ButterworthFOSection(ButterworthFOSection&) {}
public:
	ButterworthFOSection(eq_double_t beta,
	    eq_double_t s, eq_double_t g, eq_double_t g0,
	    eq_double_t D, eq_double_t c0)
	{
		b0 = (g*g*beta*beta + 2*g*g0*s*beta + g0*g0)/D;
		b1 = -4*c0*(g0*g0 + g*g0*s*beta)/D;
		b2 = 2*(g0*g0*(1 + 2*c0*c0) - g*g*beta*beta)/D;
		b3 = -4*c0*(g0*g0 - g*g0*s*beta)/D;
		b4 = (g*g*beta*beta - 2*g*g0*s*beta + g0*g0)/D;

		a0 = 1;
		a1 = -4*c0*(1 + s*beta)/D;
		a2 = 2*(1 + 2*c0*c0 - beta*beta)/D;
		a3 = -4*c0*(1 - s*beta)/D;
		a4 = (beta*beta - 2*s*beta + 1)/D;
	}

	FOSection get()
	{
		return *this;
	}
};

class ChebyshevType1FOSection : public FOSection {
	ChebyshevType1FOSection() {}
	ChebyshevType1FOSection(ChebyshevType1FOSection&) {}
public:
	ChebyshevType1FOSection(eq_double_t a,
	    eq_double_t c, eq_double_t tetta_b,
	    eq_double_t g0, eq_double_t s, eq_double_t b,
	    eq_double_t D, eq_double_t c0)
	{
		b0 = ((b*b + g0*g0*c*c)*tetta_b*tetta_b + 2*g0*b*s*tetta_b + g0*g0)/D;
		b1 = -4*c0*(g0*g0 + g0*b*s*tetta_b)/D;
		b2 = 2*(g0*g0*(1 + 2*c0*c0) - (b*b + g0*g0*c*c)*tetta_b*tetta_b)/D;
		b3 = -4*c0*(g0*g0 - g0*b*s*tetta_b)/D;
		b4 = ((b*b + g0*g0*c*c)*tetta_b*tetta_b - 2*g0*b*s*tetta_b + g0*g0)/D;

		a0 = 1;
		a1 = -4*c0*(1 + a*s*tetta_b)/D;
		a2 = 2*(1 + 2*c0*c0 - (a*a + c*c)*tetta_b*tetta_b)/D;
		a3 = -4*c0*(1 - a*s*tetta_b)/D;
		a4 = ((a*a + c*c)*tetta_b*tetta_b - 2*a*s*tetta_b + 1)/D;
	}

	FOSection get()
	{
		return *this;
	}
};

class ChebyshevType2FOSection : public FOSection {
	ChebyshevType2FOSection() {}
	ChebyshevType2FOSection(ChebyshevType2FOSection&) {}
public:
	ChebyshevType2FOSection(eq_double_t a,
	    eq_double_t c, eq_double_t tetta_b,
	    eq_double_t g, eq_double_t s, eq_double_t b,
	    eq_double_t D, eq_double_t c0)
	{
		b0 = (g*g*tetta_b*tetta_b + 2*g*b*s*tetta_b + b*b + g*g*c*c)/D;
		b1 = -4*c0*(b*b + g*g*c*c + g*b*s*tetta_b)/D;
		b2 = 2*((b*b + g*g*c*c)*(1 + 2*c0*c0) - g*g*tetta_b*tetta_b)/D;
		b3 = -4*c0*(b*b + g*g*c*c - g*b*s*tetta_b)/D;
		b4 = (g*g*tetta_b*tetta_b - 2*g*b*s*tetta_b + b*b + g*g*c*c)/D;

		a0 = 1;
		a1 = -4*c0*(a*a + c*c + a*s*tetta_b)/D;
		a2 = 2*((a*a + c*c)*(1 + 2*c0*c0) - tetta_b*tetta_b)/D;
		a3 = -4*c0*(a*a + c*c - a*s*tetta_b)/D;
		a4 = (tetta_b*tetta_b - 2*a*s*tetta_b + a*a + c*c)/D;
	}

	FOSection get()
	{
		return *this;
	}
};

class EllipticFOSection : public FOSection {
	EllipticFOSection() {}
	EllipticFOSection(EllipticFOSection&) {}
public:
	EllipticFOSection(eq_double_t a,
	    eq_double_t c, eq_double_t tetta_b,
	    eq_double_t g, eq_double_t s, eq_double_t b,
	    eq_double_t D, eq_double_t c0)
	{
		/*
		Initialize all coefficients of section.
		*/
	}

	FOSection get()
	{
		return *this;
	}
};

/*
 * Bandpass filter representation.
 */
class BPFilter {
public:
	BPFilter() {}
	virtual ~BPFilter() {}

	virtual eq_double_t process(eq_double_t in) = 0;
};

class ButterworthBPFilter : public BPFilter {
private:
	std::vector<FOSection> sections;

	ButterworthBPFilter() {}
public:
	ButterworthBPFilter(ButterworthBPFilter& f)
	{
		this->sections = f.sections;
	}
	
	ButterworthBPFilter(unsigned int N,
	    eq_double_t w0, eq_double_t wb,
	    eq_double_t G, eq_double_t Gb, eq_double_t G0)
	{
		/* Case if G == 0 : allpass. */
		if (G == 0 && G0 == 0) {
			sections.push_back(FOSection());
			return;
		}

		/* Get number of analog sections. */
		unsigned int r = N % 2;
		unsigned int L = (N - r) / 2;

		/* Convert gains to linear scale. */
		G = Conversions::db2Lin(G);
		Gb = Conversions::db2Lin(Gb);
		G0 = Conversions::db2Lin(G0);

		eq_double_t epsilon = pow(((eq_double_t)(G*G - Gb*Gb)) /
		    (Gb*Gb - G0*G0), 0.5);

		eq_double_t g = pow(((eq_double_t)G), 1.0 / ((eq_double_t)N));
		eq_double_t g0 = pow(((eq_double_t)G0), 1.0 / ((eq_double_t)N));
		eq_double_t beta = pow(((eq_double_t)epsilon),
		    -1.0 / ((eq_double_t)N)) * tan(wb / 2.0);

		eq_double_t c0 = cos(w0);
		if (w0 == 0)
			c0 = 1;

		if (w0 == M_PI /2 )
			c0 = 0;

		if (w0 == M_PI)
			c0 = -1;

		/* Calculate every section. */
		for (unsigned int i = 1; i <= L; i++) {
			eq_double_t ui = (2.0 * i - 1) / N;
			eq_double_t si = sin(M_PI * ui / 2.0);
			eq_double_t Di = beta*beta + 2*si*beta + 1;

			sections.push_back(
			    ButterworthFOSection(beta, si, g, g0, Di, c0));
		}
	}

	~ButterworthBPFilter() {}

	static eq_double_t computeBWGainDb(eq_double_t gain)
	{
		eq_double_t bwGain = 0;
		if (gain <= -6)
			bwGain = gain + commonBaseGainDb;
		else if (gain > -6 && gain < 6)
			bwGain = gain * 0.5;
		else if (gain >= 6)
			bwGain = gain - commonBaseGainDb;

		return bwGain;
	}

	virtual eq_double_t process(eq_double_t in)
	{
		eq_double_t p0 = in, p1 = 0;

		/* Process FO sections in serial connection. */
		for (unsigned int i = 0; i < sections.size(); i++) {
			p1 = sections[i].process(p0);
			p0 = p1;
		}

		return p1;
	}
};

class ChebyshevType1BPFilter : public BPFilter {
	std::vector<FOSection> sections;

	ChebyshevType1BPFilter() {}
public:
	ChebyshevType1BPFilter(unsigned int N,
	    eq_double_t w0, eq_double_t wb,
	    eq_double_t G, eq_double_t Gb, eq_double_t G0)
	{
		/* Case if G == 0 : allpass. */
		if(G == 0 && G0 == 0) {
			sections.push_back(FOSection());
			return;
		}

		/* Get number of analog sections. */
		unsigned int r = N % 2;
		unsigned int L = (N - r) / 2;

		/* Convert gains to linear scale. */
		G = Conversions::db2Lin(G);
		Gb = Conversions::db2Lin(Gb);
		G0 = Conversions::db2Lin(G0);

		eq_double_t epsilon = pow((eq_double_t)(G*G - Gb*Gb) /
		    (Gb*Gb - G0*G0), 0.5);

		eq_double_t g0 = pow((eq_double_t)(G0), 1.0 / N);
		eq_double_t alfa =
		    pow(1.0 / epsilon + pow(1 + pow(epsilon, -2.0), 0.5), 1.0 / N);

		eq_double_t beta =
		    pow(G / epsilon + Gb*pow(1 + pow(epsilon, -2.0), 0.5),1.0 / N);

		eq_double_t a = 0.5 * (alfa - 1.0 / alfa);
		eq_double_t b = 0.5*(beta - g0*g0*(1 / beta));
		eq_double_t tetta_b = tan(wb / 2);

		eq_double_t c0 = cos(w0);
		if (w0 == 0)
			c0 = 1;
		if (w0 == M_PI / 2)
			c0 = 0;
		if (w0 == M_PI)
			c0 = -1;

		/* Calculate every section. */
		for (unsigned int i = 1; i <= L; i++) {
			eq_double_t ui = (2.0*i - 1.0) / N;
			eq_double_t ci = cos(M_PI * ui / 2.0);
			eq_double_t si = sin(M_PI * ui / 2.0);
			eq_double_t Di = (a*a + ci*ci) * tetta_b * tetta_b +
			    2.0 * a * si * tetta_b + 1;

			sections.push_back(ChebyshevType1FOSection(a, ci,
			    tetta_b, g0, si, b, Di, c0));
		}
	}

	~ChebyshevType1BPFilter() {}

	static eq_double_t computeBWGainDb(eq_double_t gain)
	{
		eq_double_t bwGain = 0;
		if (gain <= -6)
			bwGain = gain + 1;
		else if (gain > -6 && gain < 6)
			bwGain = gain * 0.9;
		else if (gain >= 6)
			bwGain = gain - 1;

		return bwGain;
	}
	
	eq_double_t process(eq_double_t in)
	{
		eq_double_t p0 = in, p1 = 0;

		/* Process FO sections in serial connection. */
		for (unsigned int i = 0; i < sections.size(); i++) {
			p1 = sections[i].process(p0);
			p0 = p1;
		}

		return p1;
	}
};

class ChebyshevType2BPFilter : public BPFilter {
private:
	std::vector<FOSection> sections;

	ChebyshevType2BPFilter() {}
public:
	ChebyshevType2BPFilter(unsigned int N,
	    eq_double_t w0, eq_double_t wb,
	    eq_double_t G, eq_double_t Gb, eq_double_t G0)
	{
		/* Case if G == 0 : allpass. */
		if (G == 0 && G0 == 0) {
			sections.push_back(FOSection());
			return;
		}

		/* Get number of analog sections. */
		unsigned int r = N % 2;
		unsigned int L = (N - r) / 2;

		/* Convert gains to linear scale. */
		G = Conversions::db2Lin(G);
		Gb = Conversions::db2Lin(Gb);
		G0 = Conversions::db2Lin(G0);

		eq_double_t epsilon = pow((eq_double_t)((G*G - Gb*Gb) /
		    (Gb*Gb - G0*G0)), 0.5);

		eq_double_t g = pow((eq_double_t)(G), 1.0 / N);
		eq_double_t eu = pow(epsilon + sqrt(1 + epsilon*epsilon), 1.0 / N);
		eq_double_t ew = pow(G0*epsilon + Gb*sqrt(1 + epsilon*epsilon),
		    1.0 / N);

		eq_double_t a = (eu - 1.0 / eu) / 2.0;
		eq_double_t b = (ew - g*g / ew) / 2.0;
		eq_double_t tetta_b = tan(wb / 2);

		eq_double_t c0 = cos(w0);
		if (w0 == 0)
			c0 = 1;

		if (w0 == M_PI / 2)
			c0 = 0;

		if (w0 == M_PI)
			c0 = -1;

		/* Calculate every section. */
		for (unsigned int i = 1; i <= L; i++) {
			eq_double_t ui = (2.0 * i - 1.0) / N;
			eq_double_t ci = cos(M_PI * ui / 2.0);
			eq_double_t si = sin(M_PI * ui / 2.0);
			eq_double_t Di = tetta_b*tetta_b + 2*a*si*tetta_b +
			    a*a + ci*ci;

			sections.push_back(
			    ChebyshevType2FOSection(a, ci, tetta_b, g, si, b, Di, c0));
		}
	}

	~ChebyshevType2BPFilter(){}
	
	static eq_double_t computeBWGainDb(eq_double_t gain)
	{
		eq_double_t bwGain = 0;
		if (gain <= -6)
			bwGain = -commonBaseGainDb;
		else if (gain > -6 && gain < 6)
			bwGain = gain*0.3;
		else if (gain >= 6)
			bwGain = commonBaseGainDb;

		return bwGain;
	}

	eq_double_t process(eq_double_t in)
	{
		eq_double_t p0 = in, p1 = 0;

		/* Process FO sections in serial connection. */
		for (unsigned int i = 0; i < sections.size(); i++) {
			p1 = sections[i].process(p0);
			p0 = p1;
		}

		return p1;
	}
};

class EllipticType2BPFilter : public BPFilter {
private:
	std::vector<FOSection> sections;

	EllipticType2BPFilter() {}
public:
	EllipticType2BPFilter(unsigned int N,
	    eq_double_t w0, eq_double_t wb,
	    eq_double_t G, eq_double_t Gb, eq_double_t G0)
	{
		/* Calculate every section. */
		/*
			sections.push_back
		*/

	}

	~EllipticType2BPFilter(){}
	
	static eq_double_t computeBWGainDb(eq_double_t gain)
	{
		eq_double_t bwGain = 0;
		if (gain <= -6)
			bwGain = -commonBaseGainDb;
		else if (gain > -6 && gain < 6)
			bwGain = gain*0.3;
		else if (gain >= 6)
			bwGain = commonBaseGainDb;

		return bwGain;
	}

	eq_double_t process(eq_double_t in)
	{
		eq_double_t p0 = in, p1 = 0;

		/* Process FO sections in serial connection. */
		for (unsigned int i = 0; i < sections.size(); i++) {
			p1 = sections[i].process(p0);
			p0 = p1;
		}

		return p1;
	}
};

/*
 * The next filter types are supported.
 */
typedef enum {
	none,
	butterworth,
	chebyshev1,
	chebyshev2,
	elliptic
} filter_type;

/*
 * Representation of single precomputed equalizer channel
 * contain vector of filters for every band gain value.
 */
class EqChannel {
	eq_double_t f0;
	eq_double_t fb;
	eq_double_t samplingFrequency;
	eq_double_t gainRangeDb;
	eq_double_t gainStepDb;

	unsigned int currentFilterIndex;
	eq_double_t currentGainDb;

	std::vector<BPFilter*> filters;
	filter_type currentChannelType;

	EqChannel() {}

	unsigned int getFltIndex(eq_double_t gainDb)
	{
		unsigned int numberOfFilters = filters.size();
		eq_double_t scaleCoef = gainDb / gainRangeDb;

		return (numberOfFilters / 2) + (numberOfFilters / 2) * scaleCoef;
	}

	void cleanupFiltersArray()
	{
		for(unsigned int j = 0; j < filters.size(); j++)
			delete filters[j];
	}

public:
	EqChannel(filter_type ft,
	    eq_double_t fs, eq_double_t f0, eq_double_t fb,
	    eq_double_t gainRangeDb = eqGainRangeDb,
	    eq_double_t gainStepDb = eqGainStepDb)
	{
		samplingFrequency = fs;
		this->f0 = f0;
		this->fb = fb;
		this->gainRangeDb = gainRangeDb;
		this->gainStepDb = gainStepDb;
		currentGainDb = 0;
		currentFilterIndex = 0;
		currentChannelType = ft;

		setChannel(currentChannelType, samplingFrequency);
	}

	~EqChannel()
	{
		cleanupFiltersArray();
	}

	eq_error_t setChannel(filter_type ft, eq_double_t fs)
	{
		(void)fs;

		eq_double_t wb = Conversions::hz2Rad(fb, samplingFrequency);
		eq_double_t w0 = Conversions::hz2Rad(f0, samplingFrequency);

		for (eq_double_t gain = -gainRangeDb; gain <= gainRangeDb;
		    gain+= gainStepDb) {
			switch(ft) {
			case (butterworth): {
				eq_double_t bw_gain =
				    ButterworthBPFilter::computeBWGainDb(gain);

				ButterworthBPFilter* bf =
				    new ButterworthBPFilter(
				    defaultEqBandPassFiltersOrder, w0,
				    wb, gain, bw_gain, eqDefaultGainDb);

				filters.push_back(bf);
				break;
			}

			case (chebyshev1): {
				eq_double_t bwGain =
				    ChebyshevType1BPFilter::computeBWGainDb(gain);

				ChebyshevType1BPFilter* cf1 =
				    new ChebyshevType1BPFilter(
				    defaultEqBandPassFiltersOrder, w0,
				    wb, gain, bwGain, eqDefaultGainDb);

				filters.push_back(cf1);
				break;
			}

			case (chebyshev2): {
				eq_double_t bwGain =
				    ChebyshevType2BPFilter::computeBWGainDb(gain);

				ChebyshevType2BPFilter* cf2 =
				    new ChebyshevType2BPFilter(
				    defaultEqBandPassFiltersOrder, w0,
				    wb, gain, bwGain, eqDefaultGainDb);

				filters.push_back(cf2);
				break;
			}

			case (elliptic): {
				/*
					Allocate elliptic filter class object,
					and pass next parameters to constructor:
					(N, w0, wb, gain, bwGain, eqDefaultGainDb)

					filters.push_back
				*/
				break;
			}

			default: {
				currentChannelType = none;
				return invalid_input_data_error;
			}
			}
		}

		/* Get current filter index. */
		currentGainDb = 0;
		currentFilterIndex = getFltIndex(currentGainDb);

		return no_error;
	}

	eq_error_t setGainDb(eq_double_t db)
	{
		if (db > -gainRangeDb && db < gainRangeDb) {
			currentGainDb = db;
			currentFilterIndex = getFltIndex(db);

			return no_error;
		}

		return invalid_input_data_error;
	}

	eq_error_t SBSProcess(eq_double_t *in, eq_double_t *out)
	{
		*out = filters[currentFilterIndex]->process(*in);

		return no_error;
	}
};

static const char *getFilterName(filter_type type)
{
	switch(type) {
	case none:
		return "not initialized";
	case butterworth:
		return "butterworth";
	case chebyshev1:
		return "chebyshev1";
	case chebyshev2:
		return "chebyshev2";
	default:
		return "none";
	}
}

/*
 * Main equalizer class.
 */
class Eq {
	Conversions conv;
	eq_double_t samplingFrequency;
	FrequencyGrid freqGrid;
	std::vector<EqChannel*> channels;
	filter_type currentEqType;

	void cleanupChannelsArray()
	{
		for(unsigned int j = 0; j < channels.size(); j++)
			delete channels[j];
	}

public:
	Eq(FrequencyGrid &fg, filter_type eq_t) : conv(46)
	{
		samplingFrequency = defaultSampleFreqHz;
		freqGrid = fg;
		currentEqType = eq_t;
		setEq(freqGrid, currentEqType);
	}

	~Eq()
	{
		cleanupChannelsArray();
	}

	eq_error_t setEq(const FrequencyGrid& fg, filter_type ft)
	{
		cleanupChannelsArray();
		channels.clear();

		freqGrid = fg;
		currentEqType = ft;

		for (unsigned int i = 0; i < freqGrid.getNumberOfBands(); i++) {
			Band bFres = freqGrid.getFreqs()[i];

			EqChannel* eq_ch = new EqChannel(ft, samplingFrequency,
			    bFres.centerFreq, bFres.maxFreq - bFres.minFreq);

			channels.push_back(eq_ch);
			channels[i]->setGainDb(eqDefaultGainDb);
		}

		return no_error;
	}

	eq_error_t setEq(filter_type ft)
	{
		return setEq(freqGrid, ft);
	}

	eq_error_t setSampleRate(eq_double_t sr)
	{
		samplingFrequency = sr;

		return setEq(currentEqType);
	}

	eq_error_t changeGains(std::vector<eq_double_t> bandGains)
	{
		if (channels.size() == bandGains.size())
			for(unsigned int j = 0; j < channels.size(); j++)
				channels[j]->setGainDb(conv.fastLin2Db(bandGains[j]));
		else
			return invalid_input_data_error;

		return no_error;
	}

	eq_error_t changeGainsDb(std::vector<eq_double_t> bandGains)
	{
		if (channels.size() == bandGains.size())
			for(unsigned int j = 0; j < channels.size(); j++)
				channels[j]->setGainDb(bandGains[j]);
		else
			return invalid_input_data_error;

		return no_error;
	}

	eq_error_t changeBandGain(unsigned int bandNumber, eq_double_t bandGain)
	{
		if (bandNumber < channels.size())
			channels[bandNumber]->setGainDb(conv.fastLin2Db(bandGain));
		else
			return invalid_input_data_error;

		return no_error;
	}

	eq_error_t changeBandGainDb(unsigned int bandNumber, eq_double_t bandGain)
	{
		if (bandNumber < channels.size())
			channels[bandNumber]->setGainDb(bandGain);
		else
			return invalid_input_data_error;

		return no_error;
	}

	eq_error_t SBSProcessBand(unsigned int bandNumber, eq_double_t *in,
	    eq_double_t *out)
	{
		if (bandNumber < getNumberOfBands())
			channels[bandNumber]->SBSProcess(in, out);
		else
			return invalid_input_data_error; 

		return no_error;
	}

	eq_error_t SBSProcess(eq_double_t *in, eq_double_t *out)
	{
		eq_error_t err = no_error;
		eq_double_t inOut = *in;

		for (unsigned int i = 0; i < getNumberOfBands(); i++) {
			err = SBSProcessBand(i, &inOut, &inOut);
			if (err)
				return err;
		}

		*out = inOut;

		return no_error;
	}

	filter_type getEqType()
	{
		return currentEqType;
	}

	const char* getStringEqType()
	{
		return getFilterName(currentEqType);
	}

	unsigned int getNumberOfBands()
	{
		return freqGrid.getNumberOfBands();
	}
};

} //namespace orfanidis_eq
