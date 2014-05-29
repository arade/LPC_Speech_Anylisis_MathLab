LPC_Speech_Anylisis_MathLab
===========================
On the project I perform analysis on a given speech signal (sx119.pcm speech in raw PCM format, with sampling frequency of 16 kHz, and samples represented as 16-bit integers shorts). The analysis was done using MATLAB. 
I will perform the First  analysis according to the given specifications. 
•	Window type: Hamming 
•	Window length (frame length): 20 milliseconds 
•	LPC analysis step: 5 milliseconds (i.e. 200 Hz analysis rate) 
•	LPC method: autocorrelation 
•	LPC prediction order: 16
For each frame of the signal I compute: 
I.	the vector a of the LPC coefficients 
II.	The prediction gain G 
III.	The prediction error sequence 
IV.	The (truncated) impulse response corresponding to the LPC all-pole model 
V.	the LPC envelope (amplitude, in logarithmic scale, of the frequency response of the all-pole model) 

