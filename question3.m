clear all;
clc;
fid=fopen('sx119.pcm','r');
[rw,lenght]=fread(fid,inf,'short');
fclose(fid);
for j=1:6
ws=[160 320 480 640 960 1280];              %Initialization of window size of anylysis frame.Possible values 160, 320, 480, 640, 960, 1280
as=80;                                      % Initialization of anylisis stape
window=hamming(ws(j));                      % possible options: hamming(ws),hann(ws),blackman(ws),bartlett(ws),rectwin(ws)
max=ceil(lenght/as)+ceil(ws(j)/as)-1;       %we find the number of our analysis frames
af=zeros(ws(j),max);                        % Initialization of matrix for each anlysis frames 
afh=zeros(ws(j),max);                       % Initialization of matrix for each anlysis frames after windowing
err=zeros(ws(j),max);                       % Initialization of matrix for each linear prediction error 
if (j==1)
p=16;                                       % Initialization of order of LPC 
rw(lenght+1:1:ceil(lenght/as)*as)=0;        %Zero padding in order to have the last anylysis frame
lenght=length(rw); 
i=ceil(lenght/as); 
coff=zeros(p,max);                      	% Initialization of matrix for LPC cofficents of each frame 
corr=zeros((2*p+1),max);                    % Initialization of matrix for signal after corrolation 
trcorr=zeros(p,max);                        % Initialization of matrix for corrolated siglan after truncation
G=zeros(6,max);                             % Initialization of matrix for gain of each  anlysis frame 
h=zeros(1000,max);                          % Initialization of matrix for impuls responce of LPC of each  anlysis frame 
t=zeros(1000,max);                          % Initialization of matrix for time 
H=zeros(256,max);                           % Initialization of matrix for frequncey responce of LPC of each  anlysis frame 
w=zeros(256,max);                           % Initialization of matrix for omega 
afhfft=zeros(512,max);                      % Initialization of matrix for FFT of each anlysis frame
end
for k=1:max
     if (((k*as)/ws(j))<1)
        af((ws(j)-(as*k)-(as/2)+1):ws(j),k)=rw(1:(as*k)+(as/2)) ;     
     end;   
    if (((k*as)/ws(j))>=1 && k<i)        
        af(1:ws(j),k)=rw(((as*k)+(as/2))-(ws(j)-1):((as*k)+(as/2)));        
    end;
    if (k>=i)
        af(1:ws(j)-(((k-i)*as)+(as/2)),k)=rw(((as*k)+(as/2))-(ws(j)-1):lenght);
    end; 
        afh(1:ws(j),k)=af(1:ws(j),k).*window;
        afhfft(1:512,k)=fft(afh(1:ws(j),k),512);
        coff(1:(p+1),k)=lpc((afh(1:ws(j),k)),p);        
        corr(1:(2*p+1),k)=xcorr(afh(1:ws(j),k),p);
        trcorr(1:(p+1),k)=corr(p+1:2*p+1,k);
        G(j,k)=sqrt(trcorr(1:(p+1),k)'*coff(1:(p+1),k));
        err(1:ws(j),k)=filter(coff(1:(p+1),k),1,af(1:ws(j),k));
        [h(1:1000,k),t(1:1000,k)]=impz(1,coff(1:(p+1),k),1000);    
        [H(1:256,k),w(1:256,k)]=freqz(G(j,k),coff(1:(p+1),k), 256);          
end
end
figure;
    plot(20*log10(abs(G(1,:))),'r');%160
     hold on
     plot(20*log10(abs(G(2,:))),'c');% 320
     hold on   
     plot(20*log10(abs(G(3,:))),'b');% 480
     hold on
     plot(20*log10(abs(G(4,:))),'m');% 640
     hold on
     plot(20*log10(abs(G(5,:))),'y');% 960
     hold on
     plot(20*log10(abs(G(6,:))),'g');% 1280
     hold off