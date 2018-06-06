#
# Copyright (c) 2018 Fedor Uporov
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import os
import numpy as np
from numpy import genfromtxt
from scipy import signal
import matplotlib.pyplot as plt


sample_rate_hz = 48000
test_data_vector_size = sample_rate_hz / 4;
number_of_bands = 30; # Do not change it

gains_db = "0, 0, 0, 0, 3, 0, 5, 0, 4, 16, \
            8, -16, 8, 0, 0, 0, 0, 0, 0, 0, \
            4, -10, -6, 0, -5, 0, 0, -2, 0, 0";

"""
gains_db = "0, 0, 1, 0, 2, 0, 3, 0, 4, 0, \
            5, 0, 6, 0, 8, 0, 16, 0, 0, 0, \
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0";
"""
"""
gains_db = "0, 0, 0, 16, 0, 0, 0, 0, 0, 0, \
            0, 0, 0, 0, 0, 0, -16, 0, 0, 0, \
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0";
"""

# Make and run
os.system("make")
os.system("./eq -f %d -s %d -b %d -g %s" % \
         (sample_rate_hz, test_data_vector_size, number_of_bands, gains_db))

# Plot frequency responses

butterworth_eq = genfromtxt('butterworth.tstdat', delimiter=',')
w, h = signal.freqz(butterworth_eq, 1, 65536, sample_rate_hz)
plt.title('Butterworth filter frequency response')
plt.semilogx((sample_rate_hz / 2)*(w / np.pi), 20*np.log10(abs(h)), 'b')
plt.ylabel('Amplitude [dB]', color='b')
plt.xlabel('Frequency [Hz]')
plt.grid()
plt.axis('tight')
plt.axis([20, 20000, -30, 30])
plt.show()

chebyshev1_eq = genfromtxt('chebyshev1.tstdat', delimiter=',')
w, h = signal.freqz(chebyshev1_eq, 1, 65536, sample_rate_hz)
plt.title('Chebyshev 1 filter frequency response')
plt.semilogx((sample_rate_hz / 2)*(w / np.pi), 20*np.log10(abs(h)), 'b')
plt.ylabel('Amplitude [dB]', color='b')
plt.xlabel('Frequency [Hz]')
plt.grid()
plt.axis('tight')
plt.axis([20, 20000, -30, 30])
plt.show()

chebyshev2_eq = genfromtxt('chebyshev2.tstdat', delimiter=',')
w, h = signal.freqz(chebyshev2_eq, 1, 65536, sample_rate_hz)
plt.title('Chebyshev 2 filter frequency response')
plt.semilogx((sample_rate_hz / 2)*(w / np.pi), 20*np.log10(abs(h)), 'b')
plt.ylabel('Amplitude [dB]', color='b')
plt.xlabel('Frequency [Hz]')
plt.grid()
plt.axis('tight')
plt.axis([20, 20000, -30, 30])
plt.show()

elliptic_eq = genfromtxt('elliptic.tstdat', delimiter=',')
w, h = signal.freqz(elliptic_eq, 1, 65536, sample_rate_hz)
plt.title('Elliptic filter frequency response')
plt.semilogx((sample_rate_hz / 2)*(w / np.pi), 20*np.log10(abs(h)), 'b')
plt.ylabel('Amplitude [dB]', color='b')
plt.xlabel('Frequency [Hz]')
plt.grid()
plt.axis('tight')
plt.axis([20, 20000, -30, 30])
plt.show()

# Cleanup
os.system("make clean")

