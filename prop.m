% bcb / Nov. 2008
% Propagator matrix method from Aki and Richards
% Calculates the frequency dependent Reflection coefficient 
% for a layered half-space. This assumes that the angle of incidence is 0,
% so all waves are traveling vertically. This makes horizontal slowness
% zero, and vertical slowness (eta) = 1/Vs. 
% INPUT: model name (model formatted as Depth, Vs, and density)
% OUTPUT: R - result from Reflectivity Method in freq domain
%         RC - Reflection coefficient in time domain
%         RCS - RC convolved with ScS3

%w G rho B        freq, shear mod, density, Vs 
%v z P eta        wave number, depth, Prop Matrix, vertical slowness
%a b M            symbols from solution to reflection coefficient 
%c m n np model   counters, number of points, model name; i and j are used for complex numbers

%function [wavelet,w,R,RC,RCS] = prop(model)
%load values (z, B, rho) from model
clear all, close all
model = 'stix.txt'
data = load (model);
n=size(data,1)
N = 2^13
N2 = N/2

rho=zeros(1,n);B=zeros(1,n);z=zeros(1,n);G=zeros(1,n);eta=zeros(1,n);
%load model values into arrays
    for m = 1:1:n
        rho(m) = data(m,3);
        B(m) = data(m,2);
        z(m) = data(m,1);
        G(m) = B(m)^2*rho(m);
        eta(m) = 1.0/B(m);
    end

%calculate fourier frequencies
    w=zeros(1,N2);R=zeros(1,N2);
    dt=1.0;
    for c=1:1:N2
        w(c) = 2.0*pi*c/(N*dt);
    end

%loop over frequencies -- ends at nyquist freq (nyq*2pi)
    for c=1:1:N2
        P = [1,0;0,1];
        v=zeros(1,n);
        for m = 1:1:n
            v(m) = i*w(c)*eta(m); 
            
        end
        for m = 2:1:n-1
            P = P * [cosh(v(m)*(z(m)-z(m-1))), (1.0/(v(m)*G(m)))*sinh(v(m)*(z(m)-z(m-1)));...
                     v(m)*G(m)*sinh(v(m)*(z(m)-z(m-1))), cosh(v(m)*(z(m)-z(m-1)))];
        end
        a=i*w(c)*G(1)*eta(1);
        b=i*w(c)*G(m+1)*eta(m+1);
        M = a*( (P(1,1)+b*P(1,2)) / (P(2,1)+b*P(2,2)) );
        R(c) = (M-1.0) / (M+1.0);
    end % end freq loop
    
%create a conjugate symmetric array
%x(1) = DC, x(N2+1)=nyquest NEEDS TO BE REAL VALUED
%must be formated this way to get real valued results
    Rsym = zeros(1,N);
    Rsym(2:N2+1)=R(1:N2);
    Rsym(N2+2:N)=fliplr(conj(R(1:N2-1)));
    Rsym(N2+1)=real(Rsym(N2+1));
    RC = ifft(Rsym,'symmetric');
 
% can be used to test for Symmetry
%for n=1:1:N
%    if (Rsym(n) == conj(Rsym(mod(N-n+1,N)+1)))
%        tmp(n)=1;
%    else tmp(n)=0;
%    end
%end
   
%import ScS3 wavelet and pad with zeros
    wavelet = zeros(1,N);
    temp = load ('scs3.txt');
    for n=1:1:1500
        wavelet(n)=temp(n);
    end
    
%convolve wavelet with RC
    RCS = conv(RC,wavelet);

    phi=atan(imag(R)./real(R));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting stuff                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure, plot(w/(2*pi)*1000,phi)
%title('Phase')
%figure, plot(w/(2*pi)*1000,abs(R))
%title('Frequency Domain'),ylabel ('Amplitude'), xlabel ('Frequency (mHz)')
x1(1:N)=1:1:N;
%figure, plot(x1*dt,RC)
%title('Time Domain'),ylabel ('Amplitude'), xlabel ('Time (sec)')
%axis ([0 500 -1 1])
x2(1:(N*2)-1)=1:1:(N*2)-1;
%figure, plot(x2*dt,RCS)
%title('Reflection Coefficient'),ylabel ('Amplitude'), xlabel ('Time (sec)')
%axis ([500 1500 -1 1])
%figure, plot(x1*dt,wavelet)
%title('ScS3 Wavelet'),ylabel ('Amplitude'), xlabel ('Time (sec)')
%axis ([700 1300 -1 1])

f1=fopen('f.txt','w');
f2=fopen('t.txt','w');
f3=fopen('r.txt','w');
for c=1:1:N2
    fprintf(f1,'%f %f\n',w(c)/(2*pi)*1000,abs(R(c)));
end
for c=1:1:N
    fprintf(f2,'%f %f\n',x1(c)*dt,RC(c));
end
for c=1:1:(N*2)-1
    fprintf(f3,'%f %f\n',x2(c)*dt-935,RCS(c));
end
fclose(f1);fclose(f2);fclose(f3);