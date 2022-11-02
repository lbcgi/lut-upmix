function   process(hpl, hpr, framelen,nfft)
global msr
global analysisWnd
global synthesisWnd
 
% global u
% global v
% global x
% global y
global M_PI
global chmaskPtr
  n = length(hpl);
  hop = msr.hop;
%   ovplen = framelen - hop;
  w = hamming(framelen, 'periodic');
  e = sum(w .* w); 

  % find dimensions
  freqs = nfft/2 + 1; % two-sided + nyquist point
%   bloecke = floor((n-nfft)/hop) + 1;
%   X = zeros(freqs, bloecke);
  % do calculation
%   x_old = zeros(5,len1);
% ibdx = 1;
  old = zeros(msr.frameLen,5);
  ib = 1;
  for bdx = 0:hop:(n-framelen)
    % window time domain signal
     xr_win = hpr(bdx+1:bdx+framelen).*analysisWnd;
    ffts_r = fft(xr_win, nfft);
    xr_fft = ffts_r(1:freqs);

    xl_win = hpl(bdx+1:bdx+framelen).*analysisWnd;
    ffts_l = fft(xl_win, nfft);
    xl_fft = ffts_l(1:freqs);

    for i = 1:freqs
        lI = imag(xl_fft(i));
        lR = real(xl_fft(i));
        rI = imag(xr_fft(i));
        rR = real(xr_fft(i));
        magnitudeLeft = abs(xl_fft(i));
        magnitudeRight = abs(xr_fft(i));
        minMag = min(magnitudeLeft, magnitudeRight);
        ratMinLeft = minMag/(magnitudeLeft + eps);
         if (abs(ratMinLeft - msr.smoothedFunc(9,i)) < msr.smoothing)
             msr.smoothedFunc(9,i) = ratMinLeft; 
         elseif (ratMinLeft > msr.smoothedFunc(9,i))
             msr.smoothedFunc(9,i) = msr.smoothedFunc(9,i) + msr.smoothing;
         else
             msr.smoothedFunc(9,i) = msr.smoothedFunc(9,i) - msr.smoothing;
         end
         ratMinRight = minMag/(magnitudeRight + eps);
         if (abs(ratMinRight - msr.smoothedFunc(8,i)) < msr.smoothing)
             msr.smoothedFunc(8,i) = ratMinRight; 
         elseif (ratMinRight > msr.smoothedFunc(8,i))
             msr.smoothedFunc(8,i) = msr.smoothedFunc(8,i) + msr.smoothing;
         else
             msr.smoothedFunc(8,i) = msr.smoothedFunc(8,i) - msr.smoothing;
         end
         magnitudeSum = magnitudeLeft + magnitudeRight; 
         if(magnitudeSum<2.2204460492503131e-016)
             magDiff = 0;
         else
             magDiff = (magnitudeRight - magnitudeLeft)/ magnitudeSum;
         end
         magDiff = max(-1.0, min(1.0, magDiff));
         phaseL = atan2(lI,lR);
         phaseR = atan2(rI,rR);
         phaseDiff = abs(phaseR - phaseL);
         if (phaseDiff > M_PI) 
                phaseDiff = 2*M_PI - phaseDiff;
         end
         [x,y] = cartesianMap(msr.magPhaseDiff2Cartesian,magDiff,phaseDiff);
       
         leftSineTerm = sin(phaseL);
         leftCosineTerm = cos(phaseL);

         rightSineTerm = sin(phaseR);
         rightCosineTerm = cos(phaseR);

         centrePhase = atan2((lI + rI),(lR + rR));
         centreCosineTerm = cos(centrePhase); 
         centreSineTerm = sin(centrePhase); 
       
         [u,x] = map_to_grid(x, 21); 
         [v,y] = map_to_grid(y, 21);
         
          magSqrt = sqrt(magnitudeLeft^2 + magnitudeRight^2);
          
        opLFwd = ffts_l(i) - ffts_r(i) * msr.smoothedFunc(8,i); 
        opRFwd = ffts_r(i) - ffts_l(i)  * msr.smoothedFunc(9,i); 
        pmix = msr.mix;
        minusPmix = 1.0 - pmix;
         
        magnitude = lerpCompCplx(1,i,magSqrt,chmaskPtr,u,v,x,y);
        real_part = magnitude * leftCosineTerm;
        imag_part = magnitude * leftSineTerm;
        msr.freqDomainOut(1,i) = real_part + 1i*imag_part;

       magnitude = lerpCompCplx(2,i,magSqrt,chmaskPtr,u,v,x,y);
       real_part = magnitude * rightCosineTerm;
        imag_part = magnitude * rightSineTerm;
         msr.freqDomainOut(2,i) = real_part + 1i*imag_part;

            magnitude = lerpCompCplx(3,i,magSqrt,chmaskPtr,u,v,x,y);
           real_part = magnitude * centreCosineTerm;
		   imag_part = magnitude * centreSineTerm;
%             msr.freqDomainOut(3,bitRevFwd) = (real_part + imag_part) * 1.2;
% 			msr.freqDomainOut(3,bitRevSym) = (real_part - imag_part) * 1.2;
             msr.freqDomainOut(3,i) = real_part + 1i*imag_part;

            magnitude = lerpCompCplx(4,i,magSqrt,chmaskPtr,u,v,x,y);
            real_part = magnitude * leftCosineTerm;
			imag_part = magnitude * leftSineTerm;
%             msr.freqDomainOut(4,bitRevFwd) = (real_part + imag_part) * minusPmix + opLFwd * pmix;
% 			msr.freqDomainOut(4,bitRevSym) = (real_part - imag_part) * minusPmix + opLSym * pmix;
            msr.freqDomainOut(4,i) = (real_part + 1i*imag_part)*minusPmix + opLFwd * pmix;

            magnitude = lerpCompCplx(5,i,magSqrt,chmaskPtr,u,v,x,y);
            real_part = magnitude * rightCosineTerm;
			imag_part = magnitude * rightSineTerm;
%           msr.freqDomainOut(5,bitRevFwd) = (real_part + imag_part) * minusPmix + opRFwd * pmix;
% 			msr.freqDomainOut(5,bitRevSym) = (real_part - imag_part) * minusPmix + opRSym * pmix;
            msr.freqDomainOut(5,i) = (real_part + 1i*imag_part)*minusPmix + opRFwd * pmix;   
          
    end
            smoothedFunc_hist(ib) =  msr.smoothedFunc(1,10);
        ib = ib + 1;
     for j = 1:5
         fdm0 = msr.freqDomainOut(j,1:freqs);
         fdm1 = [fdm0 conj(fdm0(end-1:-1:2))]; % restore one-sided spectrum in reverse order
         timeTmp = real(ifft(fdm1)).';
         timeTmp = timeTmp(1:msr.frameLen).*synthesisWnd;
         % WOLA
         msr.timeDomainOut(bdx+1:bdx+msr.frameLen,j)...
             = [old(end - msr.ovpLen+1:end,j) + timeTmp(1:msr.ovpLen);timeTmp(msr.ovpLen+1:end)];
         old(:,j) = msr.timeDomainOut(bdx+1:bdx+msr.frameLen,j);
    end
%          msr.timeDomainOut(bdx+1:bdx+msr.frameLen,6) = (lp(:,1) + lp(:,2))*msr.crossBass;

 
  end
end