classdef pepito < handle
    properties
        % --------------------- Telescope configurations
        tel;       % telescope class
        xStars;
        yStars;
        % --------------------- Input frame
        SEframes;               % Observations
        LEframes;               % Initial guess-based reconstructed PSF or image
        SEpsf;                  % SE PSFs
        LEpsf;                  % LE PSFs
        LEpsf_noTT;             % TT-compensated LE PSF
        nref;
        psInMas;                % PSF pixel scale in mas/pixel
        Samp;             % PSF nyquist-sampling vector. Undersampling is accepted
        xPeaks;yPeaks;          % Center of extraction boxes
        dxPeaks;dyPeaks;        % 3-sigma uncertainties on cog in pixel
        fluxes;dfluxes          % Fluxes in ADU with 3-sigma rms uncertainty
        nBox;                   % PSF resolution in pixel
        % --------------------- OTFs/Model
        nSrc;
        nBaseline;
        otfDL;
        mskOtf;
        nOtf;                   % OTF sampling
        npt;                    % Phase sampling
        wvl;                    % Wavelength
        Dani_l;                 % Normalized anisoplanatic covariance matrices per layers
        % --------------------- Fitting process setup
        np_psf;
        ydata;
        normFactor;
        weightMap;
        weighting = true;
        ron;
        xdata;
        xinit;
        nLayer;                 % Number of atmospheric turbulent layers
        weightInit;             % Weight initial guess in m^(-5/3)
        heightInit;             % Heights initial guess in m
        
        % --------------------- Results
        r0;
        L0;
        vib;
        bg;
        Cn2h;
        dCn2h;
        dr0;
        dL0;
        dvib;
        dbg;
        im_fit;
        im_sky;
        xfinal;
        mqe;
    end
    
    
    
    methods
        
        function obj = pepito(SEframes,tel,psInMas,xStars,yStars,nBox,nLayer,varargin)
            inputs = inputParser;
            inputs.addRequired('SEframes',@isnumeric);
            inputs.addRequired('tel',@(x) isa(x,'telescope') );
            inputs.addRequired('psInMas',@isnumeric);
            inputs.addRequired('xStars',@isnumeric);
            inputs.addRequired('yStars',@isnumeric);
            inputs.addRequired('nBox',@isnumeric);
            inputs.addRequired('nLayer',@isnumeric);
            inputs.addParameter('npt',21,@isnumeric);
            inputs.addParameter('wvl',500e-9,@isnumeric);
            inputs.addParameter('nref',1,@isnumeric);
            inputs.addParameter('weightInit',ones(1,nLayer)/nLayer,@isnumeric);
            inputs.addParameter('heightInit',linspace(500,16e3,nLayer),@isnumeric);
            inputs.addParameter('ron',2,@isnumeric);
            
            inputs.parse(SEframes,tel,psInMas,xStars,yStars,nBox,nLayer,varargin{:});
            
            % Parsing inputs
            obj.SEframes   = SEframes;
            obj.tel           = tel;
            obj.psInMas  = psInMas;
            obj.xStars     = xStars;
            obj.yStars     = yStars;
            obj.nBox       = nBox;
            obj.nLayer     = nLayer;
            obj.ron          = inputs.Results.ron;
            obj.npt          = inputs.Results.npt;
            obj.wvl          = inputs.Results.wvl;
            obj.nref          = inputs.Results.nref;
            obj.weightInit = inputs.Results.weightInit;
            obj.heightInit = inputs.Results.heightInit;
            
            obj.nSrc       = length(xStars);
            obj.Samp     = 1e3*180*3600/pi*obj.wvl/obj.tel.D/2/obj.psInMas;
            obj.nOtf       = 2*round(nBox/obj.Samp);
        end
        
        function obj = runPepito(obj,initParam)
                      
             
            % ------------------ STEP 1: SEEING ESTIMATION -------------------- %
            
            % 2.1 STACKING
            obj.LEframes = pepito.stackingFrames(obj.SEframes,3);
            
            % 2.2. CROPING
            obj.SEpsf = pepito.extractPSF(obj.SEframes,obj.xStars,obj.yStars,obj.nBox);
            
            % 2.3. STACKING
            obj.LEpsf = pepito.stackingFrames(obj.SEpsf,3);
            
            % 2.4. SEEING/L0 ESTIMATION
            obj = obj.measureSeeing(initParam);
            
            
            % ----------------- STEP 2: CN2 RETRIEVAL ------------------ %
            
            %3.1. TT SUBTRACTION
            obj.LEpsf_noTT = obj.compensateTipTilt(obj.SEpsf,obj.nref);
            
            %3.2. MODEL INSTANTIATION
            obj = obj.instantiateModel();
            
            %3.3. ITERATIVE NL LS MINIMIZATION
            obj = obj.measureCn2h();
            
        end
        
        
        %% SEEING ESTIMATION
        
        function obj = measureSeeing(obj,LEim,varargin)
            inputs = inputParser;
            inputs.addRequired('LEim',@isnumeric);
            inputs.addParameter('r0Init',0.15,@isnumeric);
            inputs.addParameter('L0Init',25,@isnumeric);
            inputs.addParameter('vibInit',[10,10,-pi/4],@isnumeric);
            inputs.addParameter('fInit',ones(1,size(LEim,3)),@isnumeric);
            inputs.addParameter('xInit',zeros(1,size(LEim,3)),@isnumeric);
            inputs.addParameter('yInit',zeros(1,size(LEim,3)),@isnumeric);
            inputs.addParameter('weighting',false,@islogical);
            inputs.parse(LEim,varargin{:});
            
            obj.im_sky = LEim;
            r0Init = inputs.Results.r0Init;
            L0Init = inputs.Results.L0Init;
            vibInit = inputs.Results.vibInit;
            fInit = inputs.Results.fInit;
            xInit = inputs.Results.xInit;
            yInit = inputs.Results.yInit;
            obj.weighting = inputs.Results.weighting;
            
            
            % -------------- 1.DEFINE SPATIAL FREQUENCY     ----------------- %
            % Frequency
            [obj.nBox(1),obj.nBox(2),nS]  = size(LEim);
            [nIm_x,nIm_y]    = size(obj.tel.pupil);
            [Fx,Fy]                     = freqspace([nIm_x,nIm_y],'meshgrid'); %from -1 to 1 m-1
            obj.xdata                 = cell(1,2);
            obj.xdata{1}           = Fx*nIm_x/2/(obj.tel.D*obj.Samp);
            obj.xdata{2}          = Fy*nIm_y/2/(obj.tel.D*obj.Samp);
            % Angular positions
            xx                         = (-nIm_x/2:1:nIm_x/2-1);
            yy                          = (-nIm_y/2:1:nIm_y/2-1);
            [obj.xdata{3},obj.xdata{4}]  = meshgrid(xx,yy);
            
            % -------------- 2. DEFINE TELESCOPE OTF     ----------------- %
            P = obj.tel.pupil;
            obj.otfDL = tools.interpolateOtf(tools.pupil2otf(P,0*P,1),size(P));
            obj.otfDL = obj.otfDL/max(obj.otfDL(:));
            
            % ------------------- 3. PSF FITTING ----------------- %
            
            %1. Data
            obj.ydata =  obj.im_sky;
            obj.normFactor = zeros(1,nS);
            xBg = obj.normFactor;
            
            for iFrame = 1:nS
                tmp =  obj.im_sky(:,:,iFrame);
                xBg(iFrame) = std(tmp(:));
                % Normalization
                obj.normFactor(iFrame) = sum(tmp(:));
                obj.ydata(:,:,iFrame) =  tmp/obj.normFactor(iFrame);
                % Weighting matrix
                wtmp         = 1./sqrt((max(tmp,0)+obj.ron^2));
                obj.weightMap(:,:,iFrame) = wtmp;
            end
            
            if obj.weighting
                obj.weightMap = obj.weightMap.*(obj.ydata>0);
            else
                obj.weightMap = ones(size(obj.ydata)).*(obj.ydata>0);
            end
            obj.ydata = obj.ydata.*obj.weightMap;
            
            %2. Model
            PSFMOD = @(x,xdata) pepito.psfKolmo(x,xdata,obj.otfDL,obj.tel.D,obj.Samp,obj.nBox);
            FUN = @(x,xdata) pepito.multiplePSF(PSFMOD(x(1:obj.np_psf),xdata),x(obj.np_psf+1:end),obj.nBox).*obj.weightMap;
            
            %3. Fitting options and initial guess
            lX    = [xInit,yInit] - 2;
            uX    = [xInit,yInit] + 2;
            lF    = zeros(1,nS);
            uF    = 10*ones(1,nS);
            if ~isempty(vibInit)
                lV     = [0,0,-pi];
                uV   = [10,10,pi];
                nVib = 3;
            else
                lV = [];
                uV = [];
                nVib = 0;
            end
            lb    = [0.05,1e-2,lV,lX,lF,-5*xBg];
            ub    = [50,1e3,uV,uX,uF,5*xBg];
            
            obj.xinit = [r0Init,L0Init,vibInit,xInit,yInit,fInit,zeros(1,nS)];
            obj.np_psf = 2+nVib;
            
            opt = optimoptions(@lsqcurvefit,'MaxIter',1e3,...
                'TolX',1e-13,'TolFun',1e-13,'MaxFunEvals',3e3,...
                'InitDamping',1,'Display','iter');
            
            %3. Fitting procedure
            [beta,~,R,~,~,~,J] = lsqcurvefit(FUN,obj.xinit,obj.xdata,obj.ydata,lb,ub,opt);
            
            %4. Unpacking results
            obj.r0     = beta(1);
            obj.L0     = beta(2);
            if nVib>1
                obj.vib     = beta(3:2+nVib);
            else 
                obj.vib = [];
            end
            obj.xPeaks = beta(3+nVib:2+nVib+nS);
            obj.yPeaks = beta(nS+3+nVib:2+nVib+2*nS);
            obj.fluxes = beta(2*nS+3+nVib:2+nVib+3*nS).*obj.normFactor;
            obj.bg = beta(3+nVib+3*nS:end);
            
            %5. Measurements uncertainties
            dbeta       = diff(nlparci(beta,R,'jacobian',J),1,2);
            for i = 1:numel(dbeta)
                dbeta(i) = diff(nlparci(beta(i),R,'jacobian',J(:,i)),1,2);
            end
            obj.dr0     = dbeta(1);
            obj.dL0     = dbeta(2);
            if nVib>1
                obj.dvib     = dbeta(3:2+nVib);
            else
                obj.dvib = [];
            end
            obj.dxPeaks = dbeta(3+nVib:2+nVib+nS);
            obj.dyPeaks = dbeta(nS+3+nVib:2+nVib+2*nS);
            obj.dfluxes = dbeta(2*nS+3+nVib:2+nVib+3*nS).*obj.normFactor';
            obj.dbg = dbeta(3+nVib+3*nS:end);
            
            %6. Final image
            obj.xfinal = beta;
            obj.xfinal = [obj.r0,obj.L0,obj.vib,obj.xPeaks,obj.yPeaks,obj.fluxes,obj.bg];
            w = obj.weightMap;
            obj.weightMap = 1;
            obj.im_fit = FUN(obj.xfinal,obj.xdata);
            obj.weightMap = w;
            
            %7. Mean Residual Error
            obj.mqe = zeros(1,obj.nSrc);
            for iFrame = 1:obj.nSrc
                ims = obj.im_sky(:,:,iFrame);
                imr = ims - obj.im_fit(:,:,iFrame);
                obj.mqe(iFrame) = sqrt(sum(imr(:).^2))/sum(ims(:));
            end
            
        end
        
        %% CN2H ESTIMATION
        
        function obj = instantiateModel(obj)
            
            % Separation
            sepInArcsec   = [obj.xStars - obj.xStars(obj.nref); obj.yStars - obj.yStars(obj.nref)];
            srcBase       = sepInArcsec(:,1:obj.nSrc~=obj.nref);
            obj.nBaseline = size(srcBase,2);
            
            % Instantiate the Diffraction-limited OTF
            obj.otfDL  = pepito.zonalCovarianceToOtf(zeros(obj.npt^2),obj.nOtf,obj.tel.D,obj.tel.D/(obj.npt-1),true(obj.npt));
            obj.mskOtf = obj.otfDL>1e-5;
            
            % Instantiate the normalized SF function
            obj.Dani_l = zeros(obj.npt^2,obj.npt^2,obj.nLayer,obj.nBaseline);
            for iSrc = 1:obj.nBaseline
                obj.Dani_l(:,:,:,iSrc) = pepito.anisoplanaticCovarianceModel(obj.tel.D,1,obj.L0Init,...
                    srcBase(:,iSrc),ones(1,obj.nLayer),obj.heightInit,obj.wvl,obj.npt);
            end
            
        end
        
        function [PSF,J] = forwardPSFmodel(obj,Cn2,otf0)
            
            % 1/ Forward model
            %otf0 is obtained from the image at the obj.psInMas resolution
            PSF = zeros(obj.nBox,obj.nBox,obj.nBaseline);
            
            if nargout > 1
                J = zeros(obj.nBox,obj.nBox,obj.nBaseline,obj.nLayer);
            end
            
            % Get the Aniso SF
            Dani = zeros(obj.npt^2,obj.npt^2,obj.nBaseline);
            for l=1:obj.nLayer
                Dani   = Dani + squeeze(obj.Dani_l(:,:,l,:)*Cn2(l));
            end
            
            % Get the aniso filter at a Shanon sampling
            for iSrc = 1:obj.nBaseline
                Kani = pepito.zonalCovarianceToOtf(-0.5*Dani(:,:,iSrc),obj.nOtf,obj.tel.D,obj.tel.D/(obj.npt-1),true(obj.npt));
                Kani(obj.mskOtf) = Kani(obj.mskOtf)./obj.otfDL(obj.mskOtf);
                Kani = pepito.interpolateOtf(Kani,obj.nBox);
                % Get the PSF
                tmp           = pepito.otf2psf(Kani.*otf0);
                F0            = sum(tmp(:));
                PSF(:,:,iSrc) = tmp/F0.*obj.weightsMap;
                
                % 2/ Gradients
                if nargout > 1
                    for l=1:obj.nLayer
                        dKani_l= pepito.zonalCovarianceToOtf(-0.5*obj.Dani_l(:,:,l,iSrc),obj.nOtf,obj.tel.D,obj.tel.D/(obj.npt-1),true(obj.npt));
                        dKani_l(obj.mskOtf) = dKani_l(obj.mskOtf)./obj.otfDL(obj.mskOtf);
                        dKani_l     = -0.5*Kani.*pepito.interpolateOtf(dKani_l,obj.nBox);
                        tmp         = pepito.otf2psf(dKani_l.*otf0);
                        J(:,:,iSrc,l) = tmp/F0;
                    end
                end
            end
            
            PSF = reshape(PSF,obj.nBox,obj.nBox*obj.nBaseline);
            if nargout > 1
                J = reshape(J,obj.nBox*obj.nBox*obj.nBaseline,obj.nLayer);
            end
            
        end
        
        function obj = measureCn2h(obj)
            
            % Observations
            idOff = 1:obj.nSrc~=obj.nref;
            obs_w = reshape(bsxfun(@mtimes,obj.LEpsf_noTT(:,:,idOff),obj.weightMap),obj.nBox,[]);
            
            % Reference OTF
            otf0  = pepito.psf2otf(obj.LEpsf_noTT(:,:,obj.nref));
            
            % PSF Model
            FUN   = @(x,xdata) obj.forwardPSFmodel(x,xdata);
            
            % Fitting options
            lb  = obj.weightInit*0;
            ub  = obj.weightInit*10;
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',3e2,...
                'InitDamping',1,'Display','iter','Jacobian','off');
            
            % Run the fitting process
            [beta,~,R,~,~,~,J] = lsqcurvefit(FUN,obj.weightInit,otf0,obs_w,lb,ub,opt);
            
            % Get results
            obj.Cn2h  = beta;
            obj.dCn2h = diff(nlparci(real(beta),R,'jacobian',J),1,2);
            obj.r0    = sum(beta)^(-3/5);
            obj.dr0   = obj.r0*3/5*sqrt(sum(obj.dCn2h.^2))/sum(obj.Cn2h);
        end
        
    end
    
    methods (Static)
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     IMAGE MANIPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function LEframes = stackingFrames(SEframes,ndim)
            LEframes = squeeze(mean(SEframes,ndim));
        end
        
        function im_ = extractPSF(frame,xP,yP,nBox)
            
            nD   = ndims(frame);
            if nD>2
                nExp = size(frame,3);
            else
                nExp = 1;
            end
            
            nP   = length(xP);
            im_  = zeros(nBox,nBox,nExp,nP);
            
            for k=1:nExp
                for p=1:nP
                    %indexes
                    idx = round(xP(p) - nBox/2 +1: xP(p) + nBox/2);
                    idy = round(yP(p) - nBox/2 +1: yP(p) + nBox/2);
                    %cropping
                    tmp = frame(idx,idy,k);
                    im_(:,:,k,p) =  tmp - median(tmp(:));
                end
            end
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     TT COMPENSATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = compensateTipTilt(frame)
            
            % Measure tip-tilt
            [xOn,yOn] = pepito.getBarycenter(frame);
            % Compensate tip-tilt
            SEframes_c = pepito.recenterFrames(frame,xOn,yOn);
            % Long-exposure
            out        = pepito.stackingFrames(SEframes_c,3);
        end
        
        function [xOn,yOn] = getBarycenter(frame,varargin)
            inputs = inputParser;
            inputs.addRequired('frame',@isnumeric);
            inputs.addParameter('thresh',1,@isnumeric);
            inputs.parse(frame,varargin{:});
            
            thresh = inputs.Results.thresh;
            
            % Instantiating outputs
            [~,~,nExp,nSrc] = size(frame);
            xOn   = zeros(nExp,nSrc);
            yOn   = zeros(nExp,nSrc);
            
            % loop on exposures and sources
            for k=1:nExp
                for j=1:nSrc
                    [xOn(k,j),yOn(k,j)] = tools.cog(squeeze(frame(:,:,k,j)),thresh);
                end
            end
        end
        
        function SEframes_c = recenterFrames(SEframes,xOn,yOn)
            
            % Get dimension and number of baselines
            [nX,nY,nExp,nSrc]  = size(SEframes);
            % Define the ouput matrix
            SEframes_c = zeros(nX,nY,nExp,nSrc,nSrc);
            
            for k=1:nExp                                          %loop on exposures
                for i=1:nSrc                                        % loop on sources
                    tmp = squeeze(SEframes(:,:,k,i));
                    for j=1:nSrc                                    % loop on references
                        dX = [nX/2-xOn(k,j),nY/2-yOn(k,j)];% Get the amount of angular shifts
                        SEframes_c(:,:,k,i,j) = tools.translateImage(tmp,dX);
                        SEframes_c(:,:,k,i,j) = SEframes_c(:,:,k,i,j)/sum(sum(squeeze(SEframes_c(:,:,k,i,j))));
                    end
                end
            end
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     SEEING-LIMITED PSF MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = psfKolmo(x,xdata,otfDL,D,Samp,npsf)
           
            %1\ Grab inputs
            fx = xdata{1}; fy = xdata{2}; f = hypot(fx,fy);           
            r0 = x(1); L0 = x(2);
                       
            %2\ Atmospheric PSD
            k   = (24.*gamma(6./5)./5).^(5./6).*(gamma(11./6).^2./(2.*pi.^(11./3)));
            psd = k*r0.^(-5./3)*(f.^2 + 1./L0.^2).^(-11./6);
            % Remove piston
            psd = psd.*(1 - 4*tools.sombrero(1,pi*D*f).^2);
                              
            % Get the atmospheric OTF
            otfAtm = fftshift(tools.psd2otf(ifftshift(psd),1/(D*Samp)^2));
            
            %3\ Multiply with the telescope OTF
            otf = otfAtm.*otfDL;
            
            %4\ Convolve with a Gaussian
            if length(x)>2
                X = xdata{3}; Y = xdata{4};
                G = pepito.gaussianModel(x(3:5),[X,Y]);
                otf = otf.*G;
            end
            
            %5\ Get the PSF
            out = tools.otfShannon2psf(otf,Samp,npsf);
        end
        
        function out = gaussianModel(x,xdata)
            % ------- Grabbing parameters ---------%
            ax = x(1); %x spreading
            ay = x(2); %y-spreading
            th = x(3);  %rotation
            
            % ------- Including shifts ---------
            n     = size(xdata,2);
            X     = xdata(:,1:n/2);
            Y     = xdata(:,n/2+1:end);
            %Rotation
            Xr    = X.*cos(th) + Y.*sin(th);
            Yr    = Y.*cos(th) - X.*sin(th);
            
            % Gaussian expression
            out = exp(-0.5*((Xr./ax).^2 + (Yr./ay).^2) );
            
        end
        
        %function out = psfAniso(x,pupil,Samp);
            
        %end
        
        function out = multiplePSF(psfModel,pStars,nBox)
            
            nBox_x = nBox(1);
            nBox_y = nBox(2);
            
            % Unpack input
            nS      = length(pStars)/4;
            xSrc = pStars(1:nS);
            ySrc = pStars(1+nS:2*nS);
            fluxS= pStars(1+2*nS:3*nS);
            dX   = [xSrc;ySrc];
            
            
            out = zeros(nBox_x,nBox_y,nS);
            
            for iS = 1:nS
                if size(psfModel,3) == 1
                    psfi = psfModel;
                else
                    psfi = psfModel(:,:,iS);
                end
                % Translate PSFs at sources position
                tmp  = fluxS(iS)*tools.translateImage(psfi,dX(:,iS));
                % Concatenate frames
                out(:,:,iS) =  tmp + pStars(end-nS+iS);
            end
        end
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     ANISOPLANATISM MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function Dani = anisoplanaticCovarianceModel(D,r0,L0,sepInArcsec,weights,heights,wvl,npt,varargin)
            inputs = inputParser;
            inputs.addRequired('D',@isnumeric);
            inputs.addRequired('r0',@isnumeric);
            inputs.addRequired('L0',@isnumeric);
            inputs.addRequired('sepInArcsec',@isnumeric);
            inputs.addRequired('weights',@isnumeric);
            inputs.addRequired('heights',@isnumeric);
            inputs.addRequired('wvl',@isnumeric);
            inputs.addRequired('npt',@isnumeric);
            inputs.addParameter('flagToeplitz',false,@islogical);
            inputs.parse(D,r0,L0,sepInArcsec,weights,heights,wvl,npt,varargin{:});
            
            flagToeplitz = inputs.Results.flagToeplitz;
            
            % Get inputs
            nLayer  = length(weights);
            f0      = 2*pi/L0;
            fracR0  = weights.*0.06*wvl^2*r0^(-5/3);
            % Phase sample locations in the pupil
            x       = -D/2:D/(npt-1):D/2;
            [x1,y1] = meshgrid(x);
            X1      = (x1(:)*ones(1,npt^2))';
            Y1      = repmat(y1,[npt,npt]);
            % Samples separation in the pupil
            rhoX    = bsxfun(@minus,X1,x1(:));
            rhoY    = bsxfun(@minus,Y1,y1(:)');
            % Tip-tilt filter
            TT      = [x1(:),y1(:)];
            FTT     = 1;%TT*pinv(TT);
            % Instantiation
            Ialpha  = @(x,y) mcDonald(f0*hypot(x,y));
            I1      = Ialpha(rhoX,rhoY);
            cte     = 0.12184*(2*pi/wvl)^2;
            Dani    = zeros(npt^2,npt^2,nLayer);
            % Anisoplanatism Structure Function
            thx     = sepInArcsec(1)/206264.8;
            thy     = sepInArcsec(2)/206264.8;
            
            
            for l = 1:nLayer
                zl   = heights(l);
                I2    = Ialpha(rhoX+zl*thx,rhoY+zl*thy);
                I3    = Ialpha(zl*thx,zl*thy);
                I4    = Ialpha(rhoX-zl*thx,rhoY-zl*thy);
                Dani(:,:,l)  = Dani(:,:,l) + cte*L0^(5/3.)*fracR0(l).*FTT*(1.2 - 2*I1 + I2 - 2*I3  + I4)*FTT';
            end
            
            function out = mcDonald(x)
                out = x.^(5/6.).*besselk(5./6,x)./(2^(5/6.)*gamma(11/6.)) ;
                out(find(x == 0)) = 3/5.;
            end
            
            if flagToeplitz
                tmp = zeros(2*npt-1,2*npt-1,nLayer);
                for l=1:nLayer
                    tmp(:,:,l) = pepito.covMatrix2Map(Dani(:,:,l),npt,npt);
                end
                Dani = tmp;
            end
        end
        
        function [map,num] = covMatrix2Map(mat,n1,n2,varargin)
            %% convert covariance Matrix to covariance Map
            
            inputs = inputParser;
            inputs.addRequired('mat',@ismatrix);
            inputs.addRequired('n1',@isnumeric);
            inputs.addRequired('n2',@isnumeric);
            inputs.addParameter('mask1', [], @islogical );
            inputs.addParameter('mask2', [], @islogical );
            inputs.parse(mat,n1,n2,varargin{:});
            
            mask1 = inputs.Results.mask1;
            if isempty(mask1)
                mask1 = true(n1);
            end
            mask2 = inputs.Results.mask2;
            if isempty(mask2)
                mask2 = true(n2);
            end
            
            n = n1+n2-1;
            T = cell(n,1);
            for k=1:n
                T{k} = pepito.toeplitz(n2:-1:1,n2:n)+n*(k-1);
            end
            T = cell2mat(pepito.toeplitz(T(n2:-1:1),T(n2:n)));
            T(~mask2(:),:) = [];
            T(:,~mask1(:)) = [];
            map = zeros(n);
            num = zeros(n);
            for k=1:length(T(:))
                map(T(k))=map(T(k))+mat(k);
                num(T(k))=num(T(k))+1;
            end
            num(num==0)=1;
            map=map./num;
            
        end
        
        function t = toeplitz(c,r)
            % this version works on numeric vector as well as on cells vector
            r = r(:);                               % force column structure
            p = length(r);
            m = length(c);
            x = [r(p:-1:2) ; c(:)];                 % build vector of user data
            cidx = uint16(0:m-1)';
            ridx = uint16(p:-1:1);
            subscripts = cidx(:,ones(p,1)) + ridx(ones(m,1),:);  % Toeplitz subscripts
            t = x(subscripts);                                   % actual data
        end
                        
    end                    
end