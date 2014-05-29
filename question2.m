clear all;
clc;
fid=fopen('sx119.pcm','r');
[rw,lenght]=fread(fid,inf,'short');
fclose(fid);
p=16;                                   % Initialization of order of LPC 
ws=320;                                 %Initialization of window size of anylysis frame.Possible values 160, 320, 480, 640, 960, 1280
as=80;                                  % Initialization of anylisis stape
window(1:ws,1)=hamming(ws);
window(1:ws,2)=hann(ws); 
window(1:ws,3)=blackman(ws);
window(1:ws,4)=bartlett(ws); 
window(1:ws,5)=rectwin(ws);             % possible options: hamming(ws),hann(ws),blackman(ws),bartlett(ws),rectwin(ws)
rw(lenght+1:1:ceil(lenght/as)*as)=0;    %Zero padding in order to have the last anylysis frame
lenght=length(rw); 
i=ceil(lenght/as); 
max=ceil(lenght/as)+ceil(ws/as)-1;      %we find the number of our analysis frames
af=zeros(ws,max);                       % Initialization of matrix for each anlysis frames 
afh=zeros(ws,max);                      % Initialization of matrix for each anlysis frames after windowing
err=zeros(ws,max);                      % Initialization of matrix for each linear prediction error 
coff=zeros(p,max);                      % Initialization of matrix for LPC cofficents of each frame 
corr=zeros((2*p+1),max);                % Initialization of matrix for signal after corrolation 
trcorr=zeros(p,max);                    % Initialization of matrix for corrolated siglan after truncation
G=zeros(5,max);                         % Initialization of matrix for gain of each  anlysis frame 
h=zeros(1000,max);                      % Initialization of matrix for impuls responce of LPC of each  anlysis frame 
t=zeros(1000,max);                      % Initialization of matrix for time 
H=zeros(256,max);                       % Initialization of matrix for frequncey responce of LPC of each  anlysis frame 
w=zeros(256,max);                       % Initialization of matrix for omega 
afhfft=zeros(512,max);                  % Initialization of matrix for FFT of each anlysis frame
for j=1:5
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
        afh(1:ws,k)=af(1:ws,k).*window(1:ws,j);
        afhfft(1:512,k)=fft(afh(1:ws,k),512);
        coff(1:(p+1),k)=lpc((afh(1:ws,k)),p);        
        corr(1:(2*p+1),k)=xcorr(afh(1:ws,k),p);
        trcorr(1:(p+1),k)=corr(p+1:2*p+1,k);
        G(j,k)=sqrt(trcorr(1:(p+1),k)'*coff(1:(p+1),k));
        err(1:ws,k)=filter(coff(1:(p+1),k),1,af(1:ws,k));
        [h(1:1000,k),t(1:1000,k)]=impz(1,coff(1:(p+1),k),1000);    
        [H(1:256,k),w(1:256,k)]=freqz(G(j,k),coff(1:(p+1),k), 256);          
end
end
figure;
     plot(20*log10(abs(G(1,:))),'r');%hamming
     hold on;
     plot(20*log10(abs(G(2,:))),'g');%hann
     hold on;   
     plot(20*log10(abs(G(3,:))),'b');%blackman
     hold on;
     plot(20*log10(abs(G(4,:))),'c');%bartlett
     hold on;
     plot(20*log10(abs(G(5,:))),'y');%rectwin
     hold off;
   
     
    

