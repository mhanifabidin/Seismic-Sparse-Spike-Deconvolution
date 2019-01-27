clear,clc
clearvars;
fclose('all');

%change to your directory
addpath('C:\Matlab\Sparse_decon','C:\Matlab\Sparse_decon\SegyMAT','C:\Matlab\Sparse_decon\SegyMAT\GUI');
%parpool('local',4);


% Input Seismic Data

seisfile = '99_GNL91-112_PSTM_r6000_noAA_cdp3500-3800'; 
[dat] = ReadSegyFast(strcat(seisfile,'.sgy'));
[dat,Trhead,SegyHead] = ReadSegy(strcat(seisfile,'.sgy'));
dat = dat(:,101:132);
[nz nx]=size(dat);
oldnx = nx;
oldnz = nz;

%imagesc(dat);

%Smooth Data in the initial time to sabilize for wavelet extraction

divider=2;
local_op = ceil(nx/divider);
smooth_op = ceil(nz/20);
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
end


%create window selection for wavelet extraction (zero traces with 1 value
%in the time window

wav = zeros(nz,nx);
selectwin = zeros(nz,1);
selectwin(625:875,1)=1.;

% wavelet extraction

if divider==1;
  for i = 1:nx;
    trace_fx = abs(ifft((selectwin(1:nz,i)*(dat(1:nz,i)))));
    wav(:,1) = wav(:,1) + (trace_fx)/nx;
%     figure(17);
%     plot(wav);
  end
  wav(:,1) = smooth(wav(:,1),smooth_op);
  wav = repmat(wav(:,1),1,nx);
else
    for jx = 1:nx;
      for i2 = 1:local_op;
          if jx < ceil(local_op/2)
              trace_fx = abs(ifft((selectwin(1:nz,1).*dat(1:nz,i2))));
          elseif jx > (nx-floor(local_op/2))
              trace_fx = abs(ifft((selectwin(1:nz,1).*dat(1:nz,(nx-local_op)+i2))));
          else
              trace_fx = abs(ifft((selectwin(1:nz,1).*dat(1:nz,(jx-ceil(local_op/2))+i2))));
          end
          wav(:,jx) = wav(:,jx) + (trace_fx/local_op);
      end
      wav(:,jx) = smooth(wav(:,jx),smooth_op);
  end
end
