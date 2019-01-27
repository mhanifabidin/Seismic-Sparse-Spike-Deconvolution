% matlab script to demonstrate regularized inversion 
% to enhance timing accuracy (spiking) for 
% wedge response 
% -------------------------------------------------------------------------
%


clearvars
addpath('/Users/hanif/Documents/MATLAB/Sparse_decon','/Users/hanif/Documents/MATLAB/Sparse_decon/SegyMAT','/Users/hanif/Documents/MATLAB/Sparse_decon/SegyMAT/GUI');
seisfile = '610135_803007_h_NLmean_9_9_0_0015';
[dat,Trhead,SegyHead] = ReadSegy(strcat(seisfile,'.sgy'));
[nz nx]=size(dat);
% dat = dat(575:774,:);
dat = dat(575:674,520:799);
nz = 100;
nx = 280;
divider=1;
local_op = ceil(nx/divider);
smooth_op = ceil(nz/20);
% figure(3);
% colormap(gray);
% imagesc(dat);
% figure(4);
% plot(dat(1:nz,1));
% figure(5);
% plot(abs(ifft(dat(1:nz,1))));
if divider==1;
  wav = zeros(nz,1);
  for i = 1:nx;
    dat(1,i)=dat(1,i)*sin(0.0*0.5*pi); dat((nz-0),i)=dat((nz-0),i)*sin(0.0*0.5*pi);
    dat(2,i)=dat(2,i)*sin(0.1*0.5*pi); dat((nz-1),i)=dat((nz-1),i)*sin(0.1*0.5*pi);
    dat(3,i)=dat(3,i)*sin(0.2*0.5*pi); dat((nz-2),i)=dat((nz-2),i)*sin(0.2*0.5*pi);
    dat(4,i)=dat(4,i)*sin(0.3*0.5*pi); dat((nz-3),i)=dat((nz-3),i)*sin(0.3*0.5*pi);
    dat(5,i)=dat(5,i)*sin(0.4*0.5*pi); dat((nz-4),i)=dat((nz-4),i)*sin(0.4*0.5*pi);
    dat(6,i)=dat(6,i)*sin(0.5*0.5*pi); dat((nz-5),i)=dat((nz-5),i)*sin(0.5*0.5*pi);
    dat(7,i)=dat(7,i)*sin(0.6*0.5*pi); dat((nz-6),i)=dat((nz-6),i)*sin(0.6*0.5*pi);
    dat(8,i)=dat(8,i)*sin(0.7*0.5*pi); dat((nz-7),i)=dat((nz-7),i)*sin(0.7*0.5*pi);
    dat(9,i)=dat(9,i)*sin(0.8*0.5*pi); dat((nz-8),i)=dat((nz-8),i)*sin(0.8*0.5*pi);
    dat(10,i)=dat(10,i)*sin(0.9*0.5*pi); dat((nz-9),i)=dat((nz-9),i)*sin(0.9*0.5*pi);
    dat(11,i)=dat(11,i)*sin(1.0*0.5*pi); dat((nz-10),i)=dat((nz-10),i)*sin(1.0*0.5*pi);
    trace_fx = abs(ifft(dat(1:nz,i)));
    wav = wav + (trace_fx)/nx;
%     figure(17);
%     plot(wav);
  end
  wav = smooth(wav,smooth_op);
end
% figure(13);
% plot(wav);
% figure(12);
% plot(smooth(wav,11));
% nx = 64; nz = 128;               % nx: number of traces, nz number of time samples 
% [dat,wav] = make_wedge(nx,nz);   % create data and wavelet 
dat = dat./max(abs(dat(:)));     % scale data amplitude 

% w_t = fftshift(real(ifft(wav)))';% make the time domain wavelet 
% w_t = w_t./max(abs(w_t(:)));     % scale wavelet 
% w_s = wav./max(abs(wav(:)));     % scale frequency domain wavelet 
% colormap(gray);
% imagesc(dat);
dd = dat;                        % 
N = 2*nz;                        % N: number of time samples of output trace (N can be greater than nz for interpolation  
M = nz;                          % M: number time samples input trace 
dt1 = 1;  df1 = 1;               % time and frequency sampling input data for FT operator (for now take 1 and 1 
V1 = make_F(M,M,dt1,df1);        % make Fourier transform operator 
dt2 = M/N; df2 = 1;              % time and frequency sampling output data 
V2 = make_F(M,N,dt2,df2);        % make Fourier transform operator for interpolation 
D = make_D(N);                   % make the derivative operator for optional derivative wavelet deconvolution 
Psi = eye(N,N);                  % create the "data model": defines the basis in which the data is sparse (Identity matrix 
                                 % corresponds to spikes in the time domain
% Psi = D*Psi;                     % Derivative * Identity corresponds spiked derivative in the time domain 
for jx = 1:nx;                   % loop over traces 
  disp(['trace number ',num2str(jx)]); 
  if divider > 1
    dat(1,jx)=dat(1,jx)*sin(0.0*0.5*pi); dat((nz-0),jx)=dat((nz-0),jx)*sin(0.0*0.5*pi);
    dat(2,jx)=dat(2,jx)*sin(0.1*0.5*pi); dat((nz-1),jx)=dat((nz-1),jx)*sin(0.1*0.5*pi);
    dat(3,jx)=dat(3,jx)*sin(0.2*0.5*pi); dat((nz-2),jx)=dat((nz-2),jx)*sin(0.2*0.5*pi);
    dat(4,jx)=dat(4,jx)*sin(0.3*0.5*pi); dat((nz-3),jx)=dat((nz-3),jx)*sin(0.3*0.5*pi);
    dat(5,jx)=dat(5,jx)*sin(0.4*0.5*pi); dat((nz-4),jx)=dat((nz-4),jx)*sin(0.4*0.5*pi);
    dat(6,jx)=dat(6,jx)*sin(0.5*0.5*pi); dat((nz-5),jx)=dat((nz-5),jx)*sin(0.5*0.5*pi);
    dat(7,jx)=dat(7,jx)*sin(0.6*0.5*pi); dat((nz-6),jx)=dat((nz-6),jx)*sin(0.6*0.5*pi);
    dat(8,jx)=dat(8,jx)*sin(0.7*0.5*pi); dat((nz-7),jx)=dat((nz-7),jx)*sin(0.7*0.5*pi);
    dat(9,jx)=dat(9,jx)*sin(0.8*0.5*pi); dat((nz-8),jx)=dat((nz-8),jx)*sin(0.8*0.5*pi);
    dat(10,jx)=dat(10,jx)*sin(0.9*0.5*pi); dat((nz-9),jx)=dat((nz-9),jx)*sin(0.9*0.5*pi);
    dat(11,jx)=dat(11,jx)*sin(1.0*0.5*pi); dat((nz-10),jx)=dat((nz-10),jx)*sin(1.0*0.5*pi);
    wav = zeros(nz,1);
    for i = 1:local_op;
        if jx < ceil(local_op/2)
            trace_fx = abs(ifft(dat(1:nz,i)));
        elseif jx > (nx-floor(local_op/2))
            trace_fx = abs(ifft(dat(1:nz,(nx-local_op)+i)));
        else
            trace_fx = abs(ifft(dat(1:nz,(jx-ceil(local_op/2))+i)));
        end
        wav = wav + (trace_fx/local_op);
    end
    wav = smooth(wav,smooth_op);
  end
w_t = fftshift(real(ifft(wav)))';% make the time domain wavelet 
w_t = w_t./max(abs(w_t(:)));     % scale wavelet 
w_s = wav./max(abs(wav(:)));     % scale frequency domain wavelet 
% figure(17);
% plot(wav);
% figure(18);
% plot(w_t);
% figure(19);
V3 = diag(fftshift(w_s));        % make the wavelet convolution operator in the frequency domain (=weight function) 
V = V3*V2;                       % create the "measurement model" (wavelet convolution plus inverse FT) 
y = V1*dd(:,jx);                 % create the input data (Fourier transform of trace) 
% ---------------------------------------------------------------------
% set some of the parameter for the regularization 
k     = 1.;                      % order (k=1 corresponds to L_1, k=2 is L_2) 
la1   = .0;                     % weighting of the penalty function (amplitudes) for sparsity constraint 
la2   = .9;                      % weighting of the penalty function (derivative) for smoothness constraint 
tol   = 1e-7;                    % tolerance of the Conjugate Gradient scheme (matlab parameter) 
mitcg = 50;                      % max number of CG iterations 
mital = 70;                      % max number of total iterations 
e_tol = 1e-9;                    % tolerance of the overall scheme (min difference between subsequent iterations) 
iter  = 1;                       % iteration number 
e_up  = 1.;                      % change between two iterations 
Vp = V*Psi;                      % measurement plus data model 
f = zeros(N,1);                  % f output data 
gam = 0.999;                     % damping parameter 
% From here start iteration loop 
% calculate the matrices Eq. 6.6
while iter < mital  & e_up > e_tol
 %   disp(['iteration number: ' num2str(iter)]); 
 % create the Hessian matrix (Adjoint operator plus penalty functions) 
  H = 2*Vp'*Vp+ ... 
  k*la1.* make_LA1(f,k,N) + ... 
  k*la2*make_LA2(f,k,N,D); 
 % make the v vector; 
  v = (1.-gam)*(H*f) + gam*2*Vp'*y; 
 % solve conjugate gradient 
 [f_new,alarm] = pcg(H,v,tol,mitcg,[],[],f); 
 % calculate difference between iterations 
 e_up = (sum((abs(f_new)-abs(f)).^2)/sum(abs(f).^2+1e-7)).^0.5;
% disp([num2str(iter) ' ' num2str(e_up)])
 f = f_new; 
 iter = iter+1; 
end 
% output data 
xhat = Psi*f;
xn(1:N,jx) = real(xhat(1:N)'); 
end 
% save data 
%save('norne4d_diff_2006_2001_inline_1230_NLmeans_9_9_0_0015_local','dat','xn','wav');
save('la1_0_la2_5_D2','dat','xn','wav');
WriteSegyStructure(strcat(seisfile,'_spdc_testhanif.sgy'),SegyHead,Trhead,xn,'revision',0,'dsf',1);



