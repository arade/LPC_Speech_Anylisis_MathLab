clear all;
clc;
fid=fopen('sx119.pcm','r');
[rw,lenght]=fread(fid,inf,'short');
fclose(fid);
p=16;                                   % Initialization of order of LPC 
ws=320;                                 %Initialization of window size of anylysis frame.Possible values 160, 320, 480, 640, 960, 1280
as=80;                                  % Initialization of anylisis stape
window=hamming(ws);                     % possible options: hamming(ws),hann(ws),blackman(ws),bartlett(ws),rectwin(ws)
rw(lenght+1:1:ceil(lenght/as)*as)=0;    %Zero padding in order to have the last anylysis frame
lenght=length(rw);                      % Here we calculate the lenght of our speach signal after padding
i=ceil(lenght/as);                     
max=ceil(lenght/as)+ceil(ws/as)-2;      %we find the number of our analysis frames
af=zeros(ws,max);                       % Initialization of matrix for each anlysis frames 
afh=zeros(ws,max);                      % Initialization of matrix for each anlysis frames after windowing
err=zeros(ws,max);                      % Initialization of matrix for each linear prediction error 
coff=zeros(p,max);                      % Initialization of matrix for LPC cofficents of each frame 
corr=zeros((2*p+1),max);                % Initialization of matrix for signal after corrolation 
trcorr=zeros(p,max);                    % Initialization of matrix for corrolated siglan after truncation
G=zeros(1,max);                         % Initialization of matrix for gain of each  anlysis frame 
h=zeros(1000,max);                      % Initialization of matrix for impuls responce of LPC of each  anlysis frame 
t=zeros(1000,max);                      % Initialization of matrix for time 
H=zeros(256,max);                       % Initialization of matrix for frequncey responce of LPC of each  anlysis frame 
w=zeros(256,max);                       % Initialization of matrix for omega 
afhfft=zeros(512,max);                  % Initialization of matrix for FFT of each anlysis frame
for k=1:max
     if (((k*as)/ws)<1)
        af((ws-(as*k)-(as/2)+1):ws,k)=rw(1:(as*k)+(as/2)) ;     
     end;   
    if (((k*as)/ws)>=1 && k<i)        
        af(1:ws,k)=rw(((as*k)+(as/2))-(ws-1):((as*k)+(as/2))); 
    end;
    if (k>=i)
        af(1:ws-(((k-i)*as)+(as/2)),k)=rw(((as*k)+(as/2))-(ws-1):lenght);
    end; 
        afh(1:ws,k)=af(1:ws,k).*window;
        afhfft(1:512,k)=fft(afh(1:ws,k),512);
        coff(1:(p+1),k)=lpc((afh(1:ws,k)),p);        
        corr(1:(2*p+1),k)=xcorr(afh(1:ws,k),p);
        trcorr(1:(p+1),k)=corr(p+1:2*p+1,k);
        G(1,k)=sqrt(trcorr(1:(p+1),k)'*coff(1:(p+1),k));
        err(1:ws,k)=filter(coff(1:(p+1),k),1,af(1:ws,k));
        [h(1:1000,k),t(1:1000,k)]=impz(1,coff(1:(p+1),k),1000);    
        [H(1:256,k),w(1:256,k)]=freqz(G(1,k),coff(1:(p+1),k), 256);          
end

%sptool;

%Windowed signal of 4 selected frames in the voiced and unvoiced part ofthe
%speech
k=[157 327 104 308];
figure;
subplot(2,2,1);
    plot(afh(1:ws,k(1)));
subplot(2,2,2);
    plot(afh(1:ws,k(2)));
subplot(2,2,3);
    plot(afh(1:ws,k(3)));
subplot(2,2,4);
    plot(afh(1:ws,k(4)));

%Frequency response of the LPC envelopes for 4 selected frames in the
%voiced and unvoiced part of the speech
figure;
k=[157 327 104 308];
subplot(2,2,1);
    plot(w/pi,20*log10(abs(H(1:256,k(1)))),'b');
    hold on;
    plot(w/pi,20*log10(abs(afhfft(1:256,k(1)))),'r');
    hold off;
subplot(2,2,2);
    plot(w/pi,20*log10(abs(H(1:256,k(2)))),'b');
    hold on;
    plot(w/pi,20*log10(abs(afhfft(1:256,k(2)))),'r');
    hold off;   
subplot(2,2,3);
    plot(w/pi,20*log10(abs(H(1:256,k(3)))),'b');
    hold on;
    plot(w/pi,20*log10(abs(afhfft(1:256,k(3)))),'r');
    hold off;
subplot(2,2,4);
    plot(w/pi,20*log10(abs(H(1:256,k(4)))),'b');
    hold on;
    plot(w/pi,20*log10(abs(afhfft(1:256,k(4)))),'r');
    hold off; 

%Impulse responses for the earlier selected frames speech signal 
k=[157 327 104 308];
figure;
subplot(2,2,1);
    plot(h(1:1000,k(1)));
subplot(2,2,2);
    plot(h(1:1000,k(2)));
subplot(2,2,3);
    plot(h(1:1000,k(3)));
subplot(2,2,4);
    plot(h(1:1000,k(4))); 
    k=[157 327 104 308];

%Prediction error for voiced and unvoiced part of the speech
k=[157 327 104 308];
figure;
subplot(2,2,1);
    plot(err(1:ws,k(1)));
subplot(2,2,2);
    plot(err(1:ws,k(2)));
subplot(2,2,3);
    plot(err(1:ws,k(3)));
subplot(2,2,4);
    plot(err(1:ws,k(4)));
   
%plot of the frequency response of LPC predictor filter and frequency
%response of truncated impulse responses
k=[327 104];
for i=1:length(k)
      IR=h(1:80,k(i));      
      IR2=h(1:320,k(i));
      IR3=h(1:320,k(i)).*hamming(320);
      HIR=fft(IR,512);
      HIR2=fft(IR2,512);
      HIR3=fft(IR3,512);     
      
      subplot(2,1,i);
      plot(w/pi,20*log10(HIR(1:256)),'g');
       hold on;
      plot(w/pi,20*log10(HIR2(1:256)),'k');
       hold on; 
      plot(w/pi,20*log10(HIR3(1:256)),'b');
       hold on
      plot(w/pi,20*log10(abs(H(1:256,k(i)))/G(1,k(i))),'r');   
       hold on;
end
 hold off;
 
% Code to check characteristic of any frame 
 k=[628 629 630 631];
for i=1:length(k)
    figure;
subplot(3,2,1);
    plot(af(1:ws,k(i)));
subplot(3,2,2);
    plot(afh(1:ws,k(i)));
subplot(3,2,3);
    plot(err(1:ws,k(i)));
subplot(3,2,4);
    plot(corr(1:2*p+1,k(i)));
subplot(3,2,5);
    plot(w/pi,20*log10(abs(H(1:256,k(i)))),'b');
    hold on;
    plot(w/pi,20*log10(abs(afhfft(1:256,k(i)))),'r');
    hold off;
subplot(3,2,6);
    plot(h(1:1000,k(i)));    
   
figure;  
      Hh=fft(h(1:80,k(i)),512);
      Hhh=fft(h(1:160,k(i)),512);
     plot(w/pi,20*log10(Hhh(1:256)),'b');
      hold on;
     plot(w/pi,20*log10(Hh(1:256)),'g');
      hold on;  
     plot(w/pi,20*log10(abs(H(1:256,k(i)))/G(1,k(i))),'r');   
      hold off;
end;

