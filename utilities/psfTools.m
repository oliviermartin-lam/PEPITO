classdef psfTools < handle
    % Tools facility for psf characterization
    
    methods (Static)
        
        function out = sr2wfe(SR,lambda)
            out = 1e9*sqrt(-log(SR))*lambda/2/pi;
        end
        
        function out = wfe2sr(wfe,lambda)
            out = exp(-(2*pi*wfe*1e-9/lambda)^2);
        end
        
        function [Flux,bg,ron,msk] = getFlux(psf)
            %Define the inner circle
            npx      = length(psf);
            x        = linspace(-1,1,npx);
            [X,Y]    = meshgrid(x);
            r        = hypot(X,Y);
            msk      = r>1;
            % Computing the residual background
            psfNoise = psf .*msk;
            bg       = median(psfNoise(msk));
            % Computing the read-out noise
            ron      = std(psfNoise(msk));
            %Computing the normalized flux            
            Flux     = sum(psf(:) - bg);
        end
        
        function [FWHMx,FWHMy,dFWHM,aRatio,theta,beta] = getFWHM(psf,pixelScale,rebin,method)
            
            % Gaussian and Moffat fitting are not really efficient on
            % anisoplanatic PSF. Prefer the coutour function in such a
            % case. The cutting method is not compliant to PSF not oriented
            % along x or y-axis.
            
            if nargin < 3
                rebin = 4;
            end
            if nargin < 4
                method = 'contour';
            end
            %Interpolation   
            im2     = pepitoTools.interpolateOtf(psf,rebin*size(psf,1));
            if strcmp(method,'cutting')
                % Brutal approach when the PSF is centered and aligned
                % x-axis FWHM
                imx     = im2(:,floor(end/2+1));
                idx     = imx >= (max(imx(:))/2.);
                w       = find(idx==1);
                FWHMx   = (max(w) - min(w))/rebin*pixelScale;%sqrt(4.*sum(idx)/pi)*pixelScale/rebin;
                % y-axis FWHM
                imy     = im2(floor(end/2+1),:);
                idx     = imy >= (max(imy(:))/2.);
                w       = find(idx==1);
                FWHMy   = (max(w) - min(w))/rebin*pixelScale;%sqrt(4.*sum(idx)/pi)*pixelScale/rebin;
                theta   = 0;
            elseif strcmp(method,'contour')
                % Contour approach~: something wrong about the ellipse
                % orientation
                C       = contourc(im2,max(im2(:))*[0.5 0.5]);
                if ~isempty(C)
                    % centering the ellispe
                    mx      = [max(C(1,2:end)),max(C(2,2:end))];
                    mn      = [min(C(1,2:end)),min(C(2,2:end))];
                    cent    = (mx+mn)/2;
                    wx      = C(1,2:end) - cent(1);
                    wy      = C(2,2:end) - cent(2);
                    % Get the module
                    r       = hypot(wx,wy)/rebin*pixelScale;
                    % Getting the FWHM
                    FWHMx   = 2*max(r);
                    FWHMy   = 2*min(r);
                    % Getting the ellipse orientation
                    xm      = wx(r == max(r));
                    ym      = wy(r == max(r));
                    theta   = mean(abs(cart2pol(xm,ym)*180/pi));%mean(180*atan(ym./xm)/pi);
                    % Angle are counted positively in the reverse clockwise
                else
                    FWHMx = 0;
                    FWHMy = 0;
                end
                % direction.
            elseif strcmp(method,'Gaussian')
                xdata   = pepitoTools.getFocalGrid(size(psf),pixelScale);
                x0      = [max(psf(:)),20,20,1,0,0,0];
                f       = @(x,xdata) pepitoTools.gaussianModel(x,xdata);
                beta    = lsqcurvefit(f,x0,xdata,psf);
                FWHMx   = beta(2)*2*sqrt(2*log(2));
                FWHMy   = beta(3)*2*sqrt(2*log(2));
                theta   = beta(4);
            elseif strcmp(method,'Moffat')
                xdata   = pepitoTools.getFocalGrid(size(psf),pixelScale);
                x0      = [max(psf(:)),20,20,1,0,0,0];
                f       = @(x,xdata) pepitoTools.moffatModel(x,xdata);
                beta    = lsqcurvefit(f,x0,xdata,psf);
                FWHMx   = 2*beta(2)*sqrt(2^(1./beta(4))-1);
                FWHMy   = 2*beta(3)*sqrt(2^(1./beta(4))-1);
                theta   = beta(5);
            end
            % Get Ellipticity
            dFWHM  = sqrt(2)*pixelScale/rebin/2;
            aRatio = max(FWHMx/FWHMy,FWHMy/FWHMx);
        end
        
        function [SR,dSR] = getStrehl(psf0,pupil,Samp)
            
           
            %1\ Get the Image OTF            
            psf     = pepitoTools.recenterPSF(psf0,4);
            otf     = pepitoTools.psf2otf(psf);
            otf     = otf/max(otf(:));
           
            %2\ Get the Diffraction-limit OTF
            notf    = size(psf0,2);
            if Samp >= 1
                otfDL   = pepitoTools.telescopeOtf(pupil,2*Samp);
                otfDL   = pepitoTools.interpolate(otfDL,notf,'spline');
            else
                otfDL   = pepitoTools.telescopeOtf(pupil,2);
                psfDL  = pepitoTools.otfShannon2psf(otfDL,2*Samp,notf);
                otfDL  = pepitoTools.psf2otf(psfDL);
            end
            otfDL   = otfDL/max(otfDL(:));
            
            
            %3\ Get the Strehl
            %-1+~mod(notf,2)/2/notf:2/notf:1-1/2/notf;
            u       = linspace(-1,1-~mod(notf,2)/notf,notf);
            Mdl       = trapz(u,trapz(u,otfDL));
            SR     = abs(trapz(u,trapz(u,otf)))/Mdl;
            
            %4\ Get the uncertainty from the precision on the maximum intensity value and the image sum         
            % note that it does not include the uncertainty du to dark subtraction
            [~,~,ron] = pepitoTools.getFlux(psf);                        
            Fim = sum(psf(:));
            Mim = max(psf(:));
	    
            dM = sqrt(ron^2 + Mim);
            dF = sqrt(notf^2*ron^2 + Fim);
            dSR = 5*SR*(dM/Mim + dF/Fim); % precision
        end
        
        function ee = getEncircledEnergy(psf)                                  
            [~,~,pr] = radial(psf);
            ee  =cumsum(pr)/sum(psf(:));                                        
        end
        
        function out = getFVU(xtrue,xest,nbox)
            if nargin > 2
                n   = length(xtrue);
                idx = floor(n/2+1-nbox/2):floor(n/2+nbox/2);
                xest = xest(idx,idx);
                xtrue= xtrue(idx,idx);
            end
            MSE = sum(sum((xest-xtrue).^2));
            VarX= sum(sum((xtrue - mean(xtrue(:))).^2));
            out = MSE/VarX;
        end
        
    end
end
