function   LLPAMSProcessNPR
global msr
global analysisWnd
global synthesisWnd
global M_PI
global chmaskPtr
global MAX_OUTPUT_BUFFERS
global mode
  % find dimensions
  freqs = msr.halfLen + 2; % two-sided + nyquist point
  TempLBuffer = zeros(1,msr.frameLen);
  TempRBuffer = zeros(1,msr.frameLen);
  for i = 0: msr.frameLen-1
      k = i + msr.mInputPos;
      if (k >= msr.frameLen)
			k = k - msr.frameLen;
      end
      w = analysisWnd(i+1);
      TempLBuffer(i+1) = msr.mInput(1,k+1) * w;
	  TempRBuffer(i+1) = msr.mInput(2,k+1) * w;
  end

       ffts_r = fft(TempRBuffer, msr.fftLen);
       xr_fft = ffts_r(1:freqs);
    
       ffts_l = fft(TempLBuffer, msr.fftLen);
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
             msr.freqDomainOut(3,i) = real_part + 1i*imag_part;

            magnitude = lerpCompCplx(4,i,magSqrt,chmaskPtr,u,v,x,y);
            real_part = magnitude * leftCosineTerm;
			imag_part = magnitude * leftSineTerm;
            msr.freqDomainOut(4,i) = (real_part + 1i*imag_part)*minusPmix + opLFwd * pmix;
            magnitude = lerpCompCplx(5,i,magSqrt,chmaskPtr,u,v,x,y);
            real_part = magnitude * rightCosineTerm;
			imag_part = magnitude * rightSineTerm;
%           msr.freqDomainOut(5,bitRevFwd) = (real_part + imag_part) * minusPmix + opRFwd * pmix;
% 			msr.freqDomainOut(5,bitRevSym) = (real_part - imag_part) * minusPmix + opRSym * pmix;
            msr.freqDomainOut(5,i) = (real_part + 1i*imag_part)*1 + opRFwd * 0;   
          
    end
       
     for j = 1:mode
         fdm0 = msr.freqDomainOut(j,1:freqs);
         fdm1 = [fdm0 conj(fdm0(end-1:-1:2))]; % restore one-sided spectrum in reverse order
         msr.timeDomainOut_right(j,1:msr.fftLen) = real(ifft(fdm1,msr.fftLen)).';
%          msr.timeDomainOut_right(j,1:msr.fftLen) = real(ifft(ffts_l,msr.fftLen)).';
     end
     	msr.mOutputBufferCount = msr.mOutputBufferCount + 1;

    	if (msr.mOutputBufferCount > MAX_OUTPUT_BUFFERS)
              'return';
%               return;
        end
%          msr.timeDomainOut(bdx+1:bdx+msr.frameLen,6) = (lp(:,1) + lp(:,2))*msr.crossBass;
    
%        k = i + msr.mInputPos + msr.smpShift;
% 		if (k >= msr.frameLen)
% 			k = k - msr.frameLen;
%         plot(msr.timeDomainOut_right(1,1:384));
        for kk = 0:msr.ovpLen-1  
            mm = kk + msr.mInputPos + msr.smpShift;
            if(mm >= msr.frameLen)
                mm = mm - msr.frameLen;
            end
          
%           msr.mOverlapStage2dash(1,kk+1) = msr.timeDomainOut_right(1,1+msr.smpShift + msr.ovpLen + kk);
            msr.mOutputBuffer(1,kk+1) = msr.mOverlapStage2dash(1,kk+1) + (msr.timeDomainOut_right(1,1+kk + msr.smpShift) * synthesisWnd(kk+1)) + msr.mInputSub(1,mm+1) * msr.minuscrossBass;%20%的低频混到左声道去
            msr.mOverlapStage2dash(1,kk+1) = msr.timeDomainOut_right(1,1+msr.smpShift + msr.ovpLen + kk) * synthesisWnd(1+kk + msr.ovpLen);
            msr.mOutputBuffer(2,kk+1) = msr.mOverlapStage2dash(2,kk+1) + (msr.timeDomainOut_right(2,1+kk + msr.smpShift) * synthesisWnd(kk+1)) + msr.mInputSub(1,mm+1) * msr.minuscrossBass;%20%的低频混到左声道去
            msr.mOverlapStage2dash(2,kk+1) = msr.timeDomainOut_right(2,1+msr.smpShift + msr.ovpLen + kk) * synthesisWnd(1+kk + msr.ovpLen);
            
           for ll = 3:5
                msr.mOutputBuffer(ll,kk+1) = msr.mOverlapStage2dash(ll,kk+1) + (msr.timeDomainOut_right(ll,1+kk + msr.smpShift) * synthesisWnd(kk+1));
                msr.mOverlapStage2dash(ll,kk+1) = msr.timeDomainOut_right(ll,1+msr.smpShift + msr.ovpLen + kk) * synthesisWnd(1+kk + msr.ovpLen);
           end
            msr.mOutputBuffer(6,kk+1) = (msr.mInputSub(1,mm+1) + msr.mInputSub(2,mm+1)) * msr.crossBass;
        end
     
        msr.mInputSamplesNeeded = msr.ovpLen;
  end
