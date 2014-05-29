LPC_Speech_Anylisis_MathLab
===========================
On the project I perform analysis on a given speech signal (sx119.pcm speech in raw PCM format, with sampling frequency of 16 kHz, and samples represented as 16-bit integers shorts). The analysis was done using MATLAB. 
I will perform the first analysis according to the given specifications. 
* Window type: Hamming 
* Window length (frame length): 20 milliseconds 
* LPC analysis step: 5 milliseconds (i.e. 200 Hz analysis rate) 
* LPC method: autocorrelation 
* LPC prediction order: 16
For each frame of the signal I compute: 
* The vector a of the LPC coefficients 
* The prediction gain G 
* The prediction error sequence 
* The (truncated) impulse response corresponding to the LPC all-pole model 
* The LPC envelope (amplitude, in logarithmic scale, of the frequency response of the all-pole model) 

The second analysis will be Fundamental frequency (pitch) computation: 
* I write a procedure to compute, for each frame, the normalized autocorrelation r’(m)=r(m)/r(0), for 0<m<M (with M equivalent to 50 ms). 
* For the interval specified in above,I check graphically the position of the peak of r’(m) satisfying Mmin<m< Mmax (where Mmin and Mmax are equivalent to 400 Hz e 50 Hz, respectively) and state whether it may correspond to the pitch period. 

Finally we I will do analysis of the first two formants: I Consider again the interval specified in above: 
* I will Identify, by graphical analysis of the LPC envelope and of the Fourier Transform, the frequencies of the first and the second formants. adjust the order of the LPC model if needed. 
* Write a Matlab procedure which, starting from the LPC envelope, derives the frequencies of the first and the second formants. 
* Repeat the analysis of the steps on a new signal derived from the original one by filtering it with an FIR filter characterized by a transfer function H(z)=1-0.95z-1 



