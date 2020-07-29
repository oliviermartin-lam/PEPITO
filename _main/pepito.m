classdef pepito < handle
    properties
        % --------------------- Experiment configurations
        tel;                    % Telescope structure
        atm;                    % Atmosphere structure
        xStars;
        yStars;
        nStars;
        psInMas;                % PSF pixel scale in mas/pixel
        samp;                   % PSF nyquist-sampling vector. Undersampling is accepted
        k_;
        samp_;        
        nBox;                   % PSF resolution in pixel
        wvl;                    % Wavelength
        ron;
        % --------------------- Input frame
        cube_SE;               % Cube of short-exposure detector frames
        im_LE;                 % Stacked detector frames
        psf_SE;                % Short-exposure sub-images of each star
        psf_LE;                % Long-exposure sub-images of each star
        psf_LE_noTT;           % TT-compensated Long-exposure sub-images of each star
        FWHM;
        tipTilt;
        % --------------------- OTFs/Model
        f_;
        X_;
        mod;
        psf_mod
        nBaseline;
        otfDL;
        mskOtf;
        npt;                    % Phase sampling
        Dani_l;                 % Normalized anisoplanatic covariance matrices per layers
        nref;
        % --------------------- Fitting process setup
        np_psf;
        im_sky;
        ydata;
        kappa;
        weightMap;        
        xdata;
        xinit;
        % --------------------- Results
        r0;dr0;        
        L0;dL0;                
        Cn2h;dCn2h;        
        theta;dtheta;                
        vib;dvib;                                             
        xPeaks;yPeaks;          % Center of extraction boxes
        dxPeaks;dyPeaks;        % 3-sigma uncertainties on cog in pixel
        fluxes;dfluxes          % Fluxes in ADU with 3-sigma rms uncertainty
        bg;dbg;
        im_fit;        
        xfinal;
        xprec;
        mqe;
        % flags
        fitCn2;
        fitL0;
        fitBg;
        fitVib;
        flagSeeing;
        flagToeplitz;
    end
    
    properties (Access=private)
        pistonFilter_;
        psdKolmo_;
        xStars_;
        yStars_;      
        beta_;
        R_;
        J_;
        nFovFit_;
        nOtf_;
        otf_baseline_;
        validPhase_;
        oDL_;
        mskOtf_;
        dk_;
    end
    
    methods
        
        function obj = pepito(parFile,varargin)
            inputs = inputParser;
            inputs.addRequired('parFile', @ischar );           
            inputs.addParameter('fitCn2', true,@islogical );           
            inputs.addParameter('flagToeplitz', true,@islogical );           
            inputs.parse(parFile,varargin{:});
            fitCn2 = inputs.Results.fitCn2;
            obj.flagToeplitz = inputs.Results.flagToeplitz;
            
            %% 1\ Read the parameter file       
            eval(parFile)
            
            % Telescope configuration
            tel = [];tel.D = D;tel.cobs = cobs;tel.resolution = nPup;
            
            Dpix = nPup/2;
            x  = 1:nPup;
            [X,Y] = meshgrid(x);
            X = X - (nPup/2+1);
            Y = Y - (nPup/2+1);
            R = hypot(X,Y);
            tel.pupil = (R < Dpix) .* (R >= Dpix*cobs);

            % A priori atmosphere
            atm = [];
            atm.r0 = r0;
            atm.L0 = L0;
            atm.weight = weight;
            atm.height = height;
            atm.nLayer = numel(height);
            atm.windSpeed = windSpeed;
            atm.windDir = windDir;
            
            % Update the structure
            obj.tel     = tel;
            obj.atm     = atm;
            obj.psInMas = psInMas;
            obj.xStars_ = xStars_0;
            obj.yStars_ = yStars_0;
            obj.nBox    = nPSF;
            obj.ron     = ron;
            obj.wvl     = wvl;                        
            obj.nStars  = length(xStars_0);
            
            
            %% 2\ Define the sampling
            obj.samp = 1e3*180*3600/pi*obj.wvl/obj.tel.D/obj.psInMas;            
            if obj.samp >=2
                obj.k_ = 1;
            else
                obj.k_ = fix(ceil(2.0/obj.samp)); % works for oversampling
            end
            obj.samp_ = obj.k_ * obj.samp;          
            obj.nFovFit_ = obj.nBox + 10;
            obj.nOtf_  = obj.nFovFit_*obj.k_/2;
            
            %% 3\ Define the frequencies vector
            
            % Spatial Frequencies
            tmp = pepito.freq_array(obj.nOtf_,obj.samp_,obj.tel.D);
            fx = tmp{1};fy = tmp{2};    
            obj.f_ = hypot(fx,fy);
             
            % Atmospheric PSD and piston filter
            obj.pistonFilter_ = (1 - 4*pepitoTools.sombrero(1,pi*obj.tel.D*obj.f_).^2);
            cte = (24.*gamma(6./5)./5).^(5./6).*(gamma(11./6).^2./(2.*pi.^(11./3)));
            obj.psdKolmo_ = cte*(obj.f_.^2 + 1./obj.atm.L0.^2).^(-11./6);
            obj.psdKolmo_ = obj.psdKolmo_ .* obj.pistonFilter_;
            
            % Angular frequencies             
            x = ((1:obj.nOtf_) - obj.nOtf_/2+1)/obj.nOtf_ * obj.samp_ * (2*obj.tel.D/obj.wvl/constants.radian2arcsec);
            [X,Y] = meshgrid(x); 
            obj.X_ = {X,Y};  
           
                      
            % Diffraction-limited telescope
            obj.otfDL = pepitoTools.interpolateOtf(pepitoTools.pupil2otf(...
                obj.tel.pupil,0*obj.tel.pupil,obj.samp_/2),obj.nOtf_);
            
            %% 4 Define the anisoplanatism model
            obj.fitCn2 = fitCn2;
            
            if obj.fitCn2
                obj = obj.instantiateAnisoplanatism(inputs.Results.flagToeplitz);
            end
        end
        
        function obj = runPepito(obj,cube_SE,varargin)
            inputs = inputParser;
            inputs.addRequired('cube_SE',@isnumeric);
            inputs.addParameter('method','PSF',@ischar);
            inputs.addParameter('ron',0,@isnumeric);
            inputs.addParameter('fitL0',false,@islogical);
            inputs.parse(cube_SE,varargin{:});
             
            
            % --------------- STEP 1 : STARS DETECTION ------------------ %
                        
            % 1.1 STACKING
            obj.cube_SE = cube_SE; % Cubes of short exposure detector frames
            obj.im_LE = squeeze(sum(obj.cube_SE,3)); % Stacked image
            
            % 1.2. DETERMINING THE STARS POSITION
            [obj.psf_LE, obj.xStars, obj.yStars] = pepito.findPSF(obj.im_LE,obj.xStars_,obj.yStars_,obj.nBox);
                                                
            % -------------- STEP 2 : SEEING ESTIMATION ----------------- %                      
            obj.flagSeeing = false;
            obj = obj.measureSeeing(obj.psf_LE,'ron',inputs.Results.ron,...
                'fitL0',inputs.Results.fitL0,'method',inputs.Results.method);
            
            % Note : the seeing and Cn2(h) estimated are kept separated in
            % order to assess the outer scale first and instantiate the
            % anisoplanatism model secondly. This saves a lot of computing
            % time, although estimating an outer scale profile will be
            % investigated in the next stage of this work
            
            % ----------------- STEP 3 : CN2 RETRIEVAL ------------------ %
            if obj.fitCn2
                %3.1 GETTING THE SHORT-EXPOSURE PSF
                obj.psf_SE = pepito.extractPSF(obj.cube_SE,obj.xStars,obj.yStars,obj.nBox);
                
                %3.2. TT SUBTRACTION
                tmp = pepito.compensateTipTilt(obj.psf_SE);                                               
                
                %3.3. GET THE BASELINE PSFs/OTFs 
                [obj.psf_LE_noTT,obj.otf_baseline_] = pepito.baselinePSF(tmp,obj.nBaseline,obj.nOtf_);
                
                %3.4. ITERATIVE NL LS MINIMIZATION
                obj = obj.measureCn2h(obj.psf_LE_noTT);
            end
        end
               
        
        %% SEEING/CN2H ESTIMATION
        
        function obj = measureSeeing(obj,LEim,varargin)
            inputs = inputParser;
            inputs.addRequired('LEim',@isnumeric);           
            inputs.addParameter('method','PSF',@ischar); 
            inputs.addParameter('fitL0',false,@islogical);
            inputs.addParameter('fitVib',true,@islogical);
            inputs.addParameter('fitBg',true,@islogical);
            inputs.addParameter('ron',0,@isnumeric);
            inputs.parse(LEim,varargin{:});
            
            obj.im_sky  = LEim;
            obj.fitL0   = inputs.Results.fitL0;
            obj.fitVib  = inputs.Results.fitVib;
            obj.fitBg  = inputs.Results.fitBg;
            obj.ron     = inputs.Results.ron;
                                                        
                    
            % fitted parameters
            
            if obj.fitL0
                lL0 = obj.tel.D;uL0 = 1e3;
                L0Init = obj.atm.L0;
            else
                lL0 = [];uL0 = [];L0Init = [];
            end
            if obj.fitVib
                lV   = [0,0,0];uV   = [1e2,1e2,1e2];vibInit = [.1,.1,0];
            else
                lV = [];uV = [];vibInit = [];
            end
            
            if obj.fitBg
                lB = -10;uB = 20;bgInit= 0;
            else
                lB = [];uB = [];bgInit = [];
            end
                
            
            if strcmp(inputs.Results.method,'FWHM')
            
                %1\ ESTIMATE PSF FWHM IN MAS
                obj.FWHM = zeros(2,obj.nStars);
                for iS=1:obj.nStars
                    tmp = obj.psf_LE(:,:,iS);
                    [~,bb] = pepitoTools.getFlux(tmp);
                    [obj.FWHM(1,iS),obj.FWHM(2,iS)] = pepitoTools.getFWHM(tmp-bb,obj.psInMas/1e3,4,'contour');
                end
                
               
                %3\ FITTING OPTIONS               
                
                obj.xinit = [0.15,vibInit(1:2)];       
                lb = [0.0,lL0,lV(1:2)];
                ub = [1e2,uL0,uV(1:2)];
                obj.xinit = [0.15,L0Init,vibInit(1:2)];
                opt = optimoptions(@lsqcurvefit,'MaxIter',1e3,'TolX',1e-20,...
                    'TolFun',1e-20,'MaxFunEvals',5e3,'Display','iter');
                
                %3\ FWHM MODEL
                if obj.fitL0 
                    if obj.fitVib
                    obj.mod = @(x,xdata) [hypot(0.976*3600*180/pi*obj.wvl/x(1) * sqrt(1 - 2.183*(x(1)/x(2))^(0.356)),x(3))*ones(1,obj.nStars);...
                        hypot(0.976*3600*180/pi*obj.wvl/x(1) * sqrt(1 - 2.183*(x(1)/x(2))^(0.356)),x(4))*ones(1,obj.nStars)];
                    else
                        obj.mod = @(x,xdata) 0.976*3600*180/pi*obj.wvl/x(1) * sqrt(1 - 2.183*(x(1)/x(2))^(0.356)) *ones(2,obj.nStars);
                    end
                else
                    if obj.fitVib
                        obj.mod = @(x,xdata) [hypot(0.976*3600*180/pi*obj.wvl/x(1),x(2))*ones(1,obj.nStars);...
                            hypot(0.976*3600*180/pi*obj.wvl/x(1),x(3))*ones(1,obj.nStars)];
                    else
                        obj.mod = @(x,xdata) 0.976*3600*180/pi*obj.wvl/x(1) * ones(2,obj.nStars);
                    end
                end
                
                %4\ FIT
                [beta,~,R,~,~,~,J] = lsqcurvefit(obj.mod,obj.xinit,[],obj.FWHM,lb,ub,opt);
                
                %5\ GETTING RESULTS
                dbeta = diff(nlparci(beta,R,'jacobian',J),1,2);
                for i = 1:numel(dbeta)
                    dbeta(i) = diff(nlparci(beta(i),R,'jacobian',J(:,i)),1,2);
                end
                
                obj.r0 = beta(1);
                obj.dr0 = dbeta(1);
                if obj.fitL0
                    obj.L0 = beta(2);
                    obj.dL0 = dbeta(2);
                    nL0 = 1;
                else
                    obj.L0 = [];
                    obj.dL0 = [];
                    nL0 = 0;
                end
                if obj.fitVib
                    obj.vib = beta(2 + nL0:end);
                    obj.dvib = dbeta(2+nL0:end);
                end
            
            elseif strcmp(inputs.Results.method,'PROFILE')
                
                %1\ INPUT DATA
                obj.ydata = zeros(obj.nStars,obj.nBox/2);
                for k=1:obj.nStars
                    tmp = obj.im_sky(:,:,k);
                    [m,n] = find(tmp == max(tmp(:)));
                    obj.ydata(k,:) = radial(tmp,m,n);
                end
                obj.ydata = median(obj.ydata,1);
                obj.ydata = obj.ydata/max(obj.ydata(:));
                
               
                %2\ FITTING OPTIONS                              
                cte = (24.*gamma(6./5)./5).^(5./6).*(gamma(11./6).^2./(2.*pi.^(11./3)));
                obj.psdKolmo_ = cte*(obj.f_.^2 + 1./obj.atm.L0.^2).^(-11./6);
                obj.psdKolmo_ = obj.psdKolmo_ .* obj.pistonFilter_;
                                                
                lb = [0.9^(-5/3),lL0,lV,0,lB];
                ub = [0.03^(-5/3),uL0,uV,1e5,uB];
                obj.xinit = [0.15^(-5/3),L0Init,vibInit,1,bgInit];
                opt = optimoptions(@lsqcurvefit,'MaxIter',1e3,'TolX',1e-18,...
                    'TolFun',1e-18,'MaxFunEvals',5e3,'Display','iter');
                  
                %3\ MODEL
                obj.xdata = obj.X_;
                if obj.fitBg
                    obj.mod = @(x,xdata) x(end-1)*pepitoTools.normalize(radial(pepitoTools.crop(obj.psfModel_seeing(x(1:end-2),xdata),obj.nBox)) + x(end),'max');
                else
                    obj.mod = @(x,xdata) x(end)*pepitoTools.normalize(radial(pepitoTools.crop(obj.psfModel_seeing(x,xdata),obj.nBox)),'max');
                end
                
                
                %4\ FIT : non unicity of the solution
                [beta,~,R,~,~,~,J] = lsqcurvefit(obj.mod,obj.xinit,obj.xdata,obj.ydata,lb,ub,opt);
                
                %5\ GETTING RESULTS
                dbeta = diff(nlparci(beta,R,'jacobian',J),1,2);
                for i = 1:numel(dbeta)
                    dbeta(i) = diff(nlparci(beta(i),R,'jacobian',J(:,i)),1,2);
                end
                obj.xfinal = beta;
                obj.r0 = beta(1)^(-3/5);
                obj.dr0 = 3/5*dbeta(1)/beta(1).^(8/5);
                
                if obj.fitL0
                    obj.L0 = beta(2);
                    obj.dL0 = dbeta(2);
                    nL0 = 1;
                else
                    obj.L0 = [];
                    obj.dL0 = [];
                    nL0 = 0;
                end
                
                if obj.fitVib
                    obj.vib = beta(2+nL0:4+nL0);
                    obj.dvib = dbeta(2+nL0:4+nL0);
                end
                
                if obj.fitBg
                    obj.bg = beta(end);
                    obj.dbg = dbeta(end);
                end
                
                
            elseif strcmp(inputs.Results.method,'PSF')
                
                % --------------------- 1. PSF FITTING ---------------------- %
                
                %1. Data
                % do not consider pixels strictly equal to zero
                obj.weightMap = obj.im_sky ~= 0;
                %Account for photon and read-out noise
                if obj.ron ~= 0
                    obj.weightMap = obj.weightMap .* 1./sqrt(max(obj.im_sky,0) + obj.ron^2);
                end
                obj.ydata = obj.im_sky .* obj.weightMap;
                %normalization
                obj.kappa = abs(sum(obj.im_sky(:)));
                obj.ydata = obj.ydata/obj.kappa;
                
                %2. Fitting options and initial guess
                cte = (24.*gamma(6./5)./5).^(5./6).*(gamma(11./6).^2./(2.*pi.^(11./3)));
                obj.psdKolmo_ = cte*(obj.f_.^2 + 1./obj.atm.L0.^2).^(-11./6);
                obj.psdKolmo_ = obj.psdKolmo_ .* obj.pistonFilter_;
                
                dMax = 10;
                lX = zeros(1,2*obj.nStars) - dMax;
                uX = zeros(1,2*obj.nStars) + dMax;
                lF = zeros(1,obj.nStars);
                uF = 1e2*ones(1,obj.nStars);
                           
                lb = [1^(-5/3),lL0,lV,lX,lF,lB];
                ub = [0.03^(-5/3),uL0,uV,uX,uF,uB];
                
                obj.xinit = [0.1^(-5/3),L0Init,vibInit,zeros(1,2*obj.nStars),ones(1,obj.nStars)/obj.nStars,bgInit];
                obj.np_psf = 1 + numel(L0Init) + numel(vibInit);
                obj.xdata = obj.X_;
                
                opt = optimoptions(@lsqcurvefit,'MaxIter',1e3,'TolX',1e-14,...
                    'TolFun',1e-18,'MaxFunEvals',1e4,'Display','iter');
                
                %3. Model
                obj.psf_mod = @(x,xdata) obj.psfModel_seeing(x,xdata);
                obj.mod = @(x,xdata) pepito.multiplePSF(obj.psf_mod(x(1:obj.np_psf),xdata),x(obj.np_psf+1:end),obj.nBox).*obj.weightMap;
                
                %4. Fitting procedure
                [beta,~,R,~,~,~,J] = lsqcurvefit(obj.mod,obj.xinit,obj.xdata,obj.ydata,lb,ub,opt);
                
                %5. Unpacking results
                obj.flagSeeing = true;
                obj.updateResults(beta,R,J);
            end
            obj.flagSeeing = true;
        end
                       
        function obj = measureCn2h(obj,LEim,varargin)
            inputs = inputParser;
            inputs.addRequired('LEim',@isnumeric);
            inputs.addParameter('ron',0,@isnumeric);
            inputs.parse(LEim,varargin{:});
            obj.ron = inputs.Results.ron;
            
            %1\ Data                                                      
            obj.im_sky = LEim;
            % do not consider pixels strictly equal to zero
            obj.weightMap = obj.im_sky ~= 0;
            % Account for photon and read-out noise
            if obj.ron ~= 0
                obj.weightMap = obj.weightMap .* 1./sqrt(max(obj.im_sky,0) + obj.ron^2);
            end
            obj.ydata = obj.im_sky .* obj.weightMap;
            % Normalization
            obj.kappa = abs(sum(obj.im_sky(:)));
            obj.ydata = obj.ydata/obj.kappa;
            
               
            %2\ Fitting options
            w = obj.atm.weight(obj.atm.height ~=0);
            nL = nnz(obj.atm.height);
            
            lX = zeros(1,2*obj.nBaseline) - 10;
            uX = zeros(1,2*obj.nBaseline) + 10;
            lF = zeros(1,obj.nBaseline);
            uF = 1e5*ones(1,obj.nBaseline);
            lb = [zeros(1,nL),lX,lF,-1e10*ones(1,obj.nBaseline)];
            ub = [1e3*ones(1,nL),uX,uF,1e10*ones(1,obj.nBaseline)];
            
            obj.xinit = [0.15^(-5/3) * w/sum(w(:)),zeros(1,2*obj.nBaseline),ones(1,obj.nBaseline),zeros(1,obj.nBaseline)];
            obj.np_psf = nL;
                                 
            % Reference OTF
            obj.xdata = obj.otf_baseline_;
            
            % PSF Model
            obj.psf_mod = @(x,xdata) obj.psfModel_aniso(x,xdata);
            obj.mod = @(x,xdata) pepito.multiplePSF(obj.psf_mod(x(1:obj.np_psf),xdata),x(obj.np_psf+1:end),obj.nBox).*obj.weightMap;

            % Fitting options                        
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',5e3,...
                'Display','iter','Jacobian','off');
            
            % Run the fitting process
            [obj.beta_,~,obj.R_,~,~,~,obj.J_] = lsqcurvefit(obj.mod,obj.xinit,obj.xdata,obj.ydata,lb,ub,opt);
            
            % Get results
            obj.flagSeeing = false;
            obj.updateResults(obj.beta_,obj.R_,obj.J_);
        end
                 
        %% PSFS MODEL
        
        function [out,otf] = psfModel_seeing(obj,x,xdata)
            
            %1\ Grab inputs            
            obj.atm.r0 = x(1)^(-3/.5); 
            if obj.fitL0
                obj.atm.L0 = x(2);
                nL0 = 1;
                cte = (24.*gamma(6./5)./5).^(5./6).*(gamma(11./6).^2./(2.*pi.^(11./3)));
                obj.psdKolmo_ = cte*(obj.f_.^2 + 1./obj.atm.L0.^2).^(-11./6);  
                obj.psdKolmo_ = obj.psdKolmo_ .* obj.pistonFilter_;
            else
                nL0 = 0;
            end
            
            %2\ Get the autocovariance map
            Bphi = fft2(fftshift(x(1) * obj.psdKolmo_)) / (obj.tel.D * obj.samp_)^2;            
                                 
            %3\ Multiply with the telescope OTF
            otf = obj.otfDL .* exp(-0.5*fftshift(real(2*(max(Bphi(:)) - Bphi))));
            
            %4\ Convolve with a Gaussian jitter kernel
            if numel(x) > 1 + nL0
                % Get data
                var_x = x(2+nL0); var_y = x(3+nL0); var_xy = x(4+nL0);                       
                    X = xdata{1}; Y = xdata{2};                                       
                    % Gaussian expression
                    G = exp(-0.5*(var_x * X.^2 + var_y * Y.^2 + var_xy * X.*Y ) );
                    G = G/sum(G(:));
                    otf = otf .* G;               
            end
            otf = otf/max(otf(:));
            
            %5\ Get the PSF
            out = real(fftshift(ifft2(fftshift( otf ))));
            out = pepitoTools.interpolateOtf(out,obj.nFovFit_);
            out = out/sum(out(:));
        end
                          
        function obj = instantiateAnisoplanatism(obj,flagToeplitz)
                                     
            % Baselines
            pix2rad = obj.psInMas/1e3 * pi/180/3600;
            obj.nBaseline = obj.nStars * (obj.nStars - 1)/2;            
            
            % AA filter
            obj.dk_ = obj.nOtf_/2+1;
            x = (-1+1/obj.dk_:2/obj.dk_:1-1/obj.dk_);
            [X,Y] = meshgrid(x,x);
            AA    = [X(:),Y(:)];
            FAA   = AA*pinv(AA);
                     
            % Structure function
            atm_oomao = atmosphere(obj.wvl,0.15,'L0',obj.atm.L0,'altitude',obj.atm.height,...
                'fractionnalR0',obj.atm.weight,'windSpeed',obj.atm.windSpeed,'windDir',obj.atm.windDir);
            [aS,zS] = cart2pol((obj.xStars_-mean(obj.xStars_))*pix2rad,(obj.yStars_ - mean(obj.yStars_))*pix2rad);
            src = source('zenith',zS,'azimuth',aS);
            
            % Loops on baselines and layers
            nLayer = nnz(obj.atm.height);
            if flagToeplitz
                obj.Dani_l = zeros(obj.nOtf_,obj.nOtf_,nLayer,obj.nBaseline);
            else
                obj.Dani_l = zeros(obj.dk_^2,obj.dk_^2,nLayer,obj.nBaseline);
                obj.validPhase_ = true(obj.dk_);
                obj.oDL_ = tools.zonalCovarianceToOtf(zeros(obj.dk_^2),obj.nOtf_,obj.tel.D,obj.tel.D/(obj.dk_-1),obj.validPhase_);
                obj.oDL_ = obj.oDL_/max(obj.oDL_(:));
                obj.mskOtf_ = obj.oDL_ > 1e-5;
            end
            
            ll = 0;
            for l = find(obj.atm.height~=0)
                iB = 0;
                ll = ll+1;
                f0 = (atm_oomao.r0^(-5/3)*atm_oomao.layer(l).fractionnalR0);
                for iS = 1:obj.nStars
                    for jS = iS+1:obj.nStars
                        iB = iB + 1;
                        [Css,Cgs] = phaseStats.spatioAngularCovarianceMatrix(obj.dk_,obj.tel.D,atm_oomao.slab(l),src(jS),'srcCC',src(iS));
                        covMat    = FAA*(2*Css - Cgs{1} - Cgs{1}')*FAA'/f0;
                        if flagToeplitz
                            covMap    = pepitoTools.interpolateOtf(pepitoTools.covMatrix2Map(covMat,obj.dk_,obj.dk_),obj.nOtf_);
                            obj.Dani_l(:,:,ll,iB) = 2 * (max(covMap(:)) - covMap);
                        else
                            obj.Dani_l(:,:,ll,iB) = covMat;
                        end
                    end
                    
                end
            end                      
        end
                
        function out = psfModel_aniso(obj,Cn2,otf0)
                                                           
            out = zeros(obj.nOtf_,obj.nOtf_,obj.nBaseline);
            for iB = 1:obj.nBaseline
                %total phase structure function
                Dani = squeeze(sum(bsxfun(@times, obj.Dani_l(:,:,:,iB) , reshape(Cn2,1,1,[])), 3));
                
                % Aniso filter
                if obj.flagToeplitz
                    Kaniso = exp(-0.5 * Dani );
                else
                    Kaniso = obj.mskOtf_ .* tools.zonalCovarianceToOtf(Dani,obj.nOtf_,obj.tel.D,obj.tel.D/(obj.dk_-1),obj.validPhase_);
                    Kaniso(obj.mskOtf_) = Kaniso(obj.mskOtf_)./obj.oDL_(obj.mskOtf_);
                end
                % Get the PSF
                tmp = real(fftshift(ifft2(fftshift( otf0(:,:,iB) .* Kaniso ))));
                out(:,:,iB) = tmp/sum(tmp(:));
            end                       
        end
        
        %% RESULTS
        
        function displayResults(obj,varargin)
            inputs = inputParser;
            inputs.addParameter('fontsize',24,@isnumeric);
            inputs.parse(varargin{:});
            fontsize = inputs.Results.fontsize;
            
            close all;
            if obj.flagSeeing
                nn = obj.nStars;
            else
                nn = obj.nBaseline;
            end
            if ~isempty(obj.bg)
                bgMaps = reshape(bsxfun(@times,obj.bg,ones(obj.nBox^2,nn)),obj.nBox,obj.nBox,nn);
                ims = reshape(obj.im_sky - bgMaps,obj.nBox,[]);
                imf = reshape(obj.im_fit - bgMaps,obj.nBox,[]);
            else
                bgMaps = 0;
                ims = reshape(obj.im_sky,obj.nBox,[]);
                imf = reshape(obj.im_fit,obj.nBox,[]);
            end
            
            imr = imf - ims;
            A = abs([ims;imf;imr]);
            
            % 2D PSF
            figure;
            imagesc(asinh(A))
            pbaspect([nn,3,1]);
            cb = colorbar();
            cb.FontSize = fontsize;
            cb.TickLabelInterpreter = 'latex';
            set(gca,'XTick',[],'YTick',[]);

            
            % PSF profile
            prof_sky = zeros(nn,obj.nBox/2);
            prof_mod = zeros(nn,obj.nBox/2);
            for k=1:nn
                prof_sky(k,:) = radial(obj.im_sky(:,:,k) - bgMaps(k));
                prof_mod(k,:) = radial(obj.im_fit(:,:,k) - bgMaps(k));                
            end
                                
            prof_sky = median(prof_sky,1);
            prof_mod = median(prof_mod,1);
            
            figure;
            x = linspace(0,obj.nBox/2,obj.nBox/2)*obj.psInMas/1e3;
            semilogy(x,prof_sky,'ks','MarkerFaceColor','k','MarkerSize',7);
            hold on;
            semilogy(x,prof_mod,'b--');
            %semilogy(x,1e2*abs(prof_mod-prof_sky),'b--');
            xlabel('Angular separation (arcsec)','interpreter','latex','fontsize',fontsize);
            ylabel('PSF azimuthal profile (ADU)','interpreter','latex','fontsize',fontsize);
            legend({'Average PSF profile','Modelled profile'},'interpreter','latex','fontsize',fontsize);
            pbaspect([1,1,1]);
            set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
            xlim([0,max(x)])
        end
      
        function updateResults(obj,beta,R,J)
                                    
            % Measurements uncertainties
            dbeta = diff(nlparci(beta,R,'jacobian',J),1,2);
            for i = 1:numel(dbeta)
                dbeta(i) = diff(nlparci(beta(i),R,'jacobian',J(:,i)),1,2);
            end
            
                
            if obj.flagSeeing
                % SEEING ESTIMATION CASE                                
                
                % r0
                obj.r0 = beta(1)^(-3/5);
                obj.dr0 = 3/5*dbeta(1)/beta(1).^(8/5);
                obj.atm.r0 = obj.r0;
                
                % outer scale
                if obj.fitL0
                    obj.L0  = beta(2);
                    obj.dL0 = dbeta(2);
                    nL0 = 1;
                else
                    obj.L0 = [];
                    obj.dL0 = [];
                    nL0 = 0;
                end
                
                %jitter
                if obj.fitVib
                    nVib = 3;
                    obj.vib     = beta(2 + nL0 : 1 + obj.fitL0 + nVib);
                    obj.dvib = dbeta(2 + obj.fitL0 : 1 + obj.fitL0 + nVib)';
                else
                    nVib = 0;
                    obj.vib = [];
                    obj.dvib = [];
                end
                
                               
                % Photometry/astrometry
                obj.xPeaks = beta(2 + nL0 + nVib                : 1 + nL0 + nVib + obj.nStars);
                obj.yPeaks = beta(2 + nL0 + nVib + obj.nStars   : 1 + nL0 + nVib + 2*obj.nStars);
                obj.fluxes = beta(2 + nL0 + nVib + 2*obj.nStars : 1 + nL0 + nVib + 3*obj.nStars).*obj.kappa;
                obj.dxPeaks = dbeta(2 + nL0 + nVib                : 1 + nL0 + nVib + obj.nStars)';
                obj.dyPeaks = dbeta(2 + nL0 + nVib + obj.nStars   : 1 + nL0 + nVib + 2*obj.nStars)';
                obj.dfluxes = dbeta(2 + nL0 + nVib + 2*obj.nStars : 1 + nL0 + nVib + 3*obj.nStars)'.*obj.kappa;
                
                 % background
                if obj.fitBg
                    obj.bg  = beta(2 + nL0 + nVib + 3*obj.nStars : end)*obj.kappa;
                    obj.dbg     = dbeta(2 + nL0 + nVib + 3*obj.nStars : end)'*obj.kappa;
                else
                    obj.bg = [];
                    obj.dbg = [];
                end
                                          
                
                
                obj.xfinal = [obj.r0^(-5/3),obj.L0,obj.vib,obj.xPeaks,obj.yPeaks,obj.fluxes,obj.bg];
                obj.xprec = [obj.dr0,obj.dL0,obj.dvib,obj.dxPeaks,obj.dyPeaks,obj.dfluxes,obj.dbg];
                
                nn = obj.nStars;
                obj.atm.r0 = beta(1);
                if numel(obj.bg) == 1
                    obj.bg = obj.bg*ones(1,obj.nStars);
                end
                
            else
                
                % CN2H ESTIMATION CASE
                nL        = nnz(obj.atm.height);
                obj.Cn2h  = beta(1:nL);
                obj.dCn2h = dbeta(1:nL)';                
                
                %update the atmosphere
                tmp       = [obj.r0^(-5/3) - sum(obj.Cn2h),obj.Cn2h];
                obj.atm.weight = tmp/sum(tmp);
                obj.atm.r0 = obj.r0;
                
                % isoplanatic angle calculation
                cst = 3600*180/pi * ((24*gamma(6/5)/5)^(-5/6))^(3/5)* obj.atm.r0;
                hb  = sum( obj.atm.weight .* obj.atm.height.^(5/3) )^(3/5);
                obj.theta = cst/hb;
                
                dW = obj.dCn2h .* [obj.Cn2h*sum(obj.Cn2h)^(-2) + sum(obj.Cn2h)^(-1)];
                % I do not understand this partial derivative here, it's
                % not homogeneous
                dhb = sum( 3/5 * sum( obj.atm.weight .* obj.atm.height.^(5/3) ).^(-2/5) .* obj.atm.height.^(5/3) .* [0,dW] );                
                obj.dtheta = obj.theta*(obj.dr0/obj.r0 + dhb/hb );
                               
                % results
                obj.xPeaks = beta(nL + 1 : nL + obj.nBaseline);
                obj.yPeaks = beta(nL + 1 + obj.nBaseline : nL + 2*obj.nBaseline);
                obj.fluxes = beta(nL + 1 + 2*obj.nBaseline : nL+ 3*obj.nBaseline).*obj.kappa;
                obj.bg     = beta(nL + 1 + 3*obj.nBaseline : end)*obj.kappa;
               
                if numel(obj.bg) == 1
                    obj.bg = obj.bg*ones(1,obj.nBaseline);
                end
                
                %uncertainties
                obj.dxPeaks = dbeta(nL + 1 : nL  + obj.nBaseline)';
                obj.dyPeaks = dbeta(nL + 1 + obj.nBaseline : nL + 2*obj.nBaseline)';
                obj.dfluxes = dbeta(nL + 1 + 2*obj.nBaseline : nL + 3*obj.nBaseline)'.*obj.kappa;
                obj.dbg     = dbeta(nL + 1 + 3*obj.nBaseline : end)'*obj.kappa;
                
                obj.xfinal = [obj.Cn2h,obj.xPeaks,obj.yPeaks,obj.fluxes,obj.bg];
                obj.xprec   = [obj.dCn2h,obj.dxPeaks,obj.dyPeaks,obj.dfluxes,obj.dbg];
                
                nn = obj.nBaseline;
            end
            
            
            %6. Final image
            obj.mod = @(x,xdata) pepito.multiplePSF(obj.psf_mod(x(1:obj.np_psf),xdata),x(obj.np_psf+1:end),obj.nBox);
            obj.im_fit = obj.mod(obj.xfinal,obj.xdata);
            
            %7. Mean Residual Error
            obj.mqe = zeros(1,nn);
            for iFrame = 1:nn
                ims = obj.im_sky(:,:,iFrame);
                imr = ims - obj.im_fit(:,:,iFrame);
                obj.mqe(iFrame) = 1e2*sqrt(sum(imr(:).^2))/sum(ims(:));
            end
        end
    end
    
    methods (Static)
                       
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
        
        function [psf_LE,xStars,yStars] = findPSF(im_LE,xStars_,yStars_,nBox)
            
            nStars = numel(xStars_);            
            psf_LE = zeros(nBox,nBox,nStars); %short exposure observed PSFs
            xStars = zeros(1,nStars);
            yStars = zeros(1,nStars);
            
            for j=1:nStars
                % First crop to detect the peaks
                idxj = xStars_(j) - nBox + 1 : xStars_(j) + nBox;
                idyj = yStars_(j) - nBox + 1 : yStars_(j) + nBox;
                % Peak detection
                im_j    = im_LE(idxj,idyj);
                tmp     = im_j;
                tmp(:)  = medfilt1(im_j(:));
                [xm,ym] = find(tmp == max(tmp(:)));
                p = zeros(1,numel(xm));
                for n=1:numel(xm)
                    p(n) = im_j(xm(n),ym(n));
                end
                nm = find(p == max(p));
                nm = nm(1);
                % crop
                idxj = xm(nm) - nBox/2 + 1 : xm(nm) + nBox/2;
                idyj = ym(nm) - nBox/2 + 1 : ym(nm) + nBox/2;
                psf_LE(:,:,j) = im_j(idxj,idyj);
                % Get stars position
                yStars(j) = ym(nm) + yStars_(j) - nBox;
                xStars(j) = xm(nm) + xStars_(j) - nBox;                
            end
        end
              
        function [out,tipTilt] = compensateTipTilt(frame)
            
            % Measure tip-tilt
            [xOn,yOn] = pepito.getBarycenter(frame);            
            
            [nStars,nFrames] = size(xOn);
            tipTilt = zeros(2,nStars,nFrames);
            tipTilt(1,:,:) = xOn - mean(xOn,2);
            tipTilt(2,:,:) = yOn - mean(yOn,2);   
            % Compensate tip-tilt
            SEframes_c = pepito.recenterFrames(frame,xOn,yOn);
            % Long-exposure
            out        = squeeze(sum(SEframes_c,3));
            
        end
             
        function [xOn,yOn] = getBarycenter(frame,varargin)
            inputs = inputParser;
            inputs.addRequired('frame',@isnumeric);
            inputs.addParameter('thresh',0,@isnumeric);
            inputs.parse(frame,varargin{:});
            
            thresh = inputs.Results.thresh;
            
            % Instantiating outputs
            [~,~,nExp,nSrc] = size(frame);
            xOn   = zeros(nSrc,nExp);
            yOn   = zeros(nSrc,nExp);
            
            % loop on exposures and sources
            for k=1:nExp
                for j=1:nSrc
                    tmp = squeeze(frame(:,:,k,j));
                    [~,~,ron] = pepitoTools.getFlux(tmp);
                    [xOn(j,k),yOn(j,k)] = tools.cog(tmp,ron*3);
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
                    tmp(tmp~=tmp) = 0;
                    for j=1:nSrc                                    % loop on references
                        dX = -[xOn(j,k) - nY/2 , yOn(j,k) - nX/2];% Get the amount of angular shifts
                        tmp2 = pepitoTools.translateImage( tmp , dX );      
                        %tmp2 = tools.shift(tmp,dX(1),dX(2));
                        SEframes_c(:,:,k,i,j) = tmp/max(abs(tmp2(:)));
                    end
                end
            end
            
        end             
        
        function [psf,otf] = baselinePSF(psfLE,nBaseline,nOtf)
            
            nBox= size(psfLE,1);
            nStars = size(psfLE,3);
            
            psf = zeros(nBox,nBox,nBaseline);
            otf = zeros(nOtf,nOtf,nBaseline);
            
            x = 2*((1:nBox) - nBox/2+1)/nBox;
            [X,Y] = meshgrid(x); 
            X_lr = {X,Y};  
            
            
            iB = 0;
            p = 0;
            for iS = 1:nStars
                tmp = psfLE(:,:,iS,iS);
                [~,bg] = pepitoTools.getFlux(tmp);
                tmp = pepitoTools.psf2otf(tmp - bg);
                tmp = tmp/max(tmp(:));
                % fit moffat
                otf_moff = pepitoTools.fitImage(tmp,X_lr,'moffat');
                otf0 = pepitoTools.interpolateOtf(otf_moff,nOtf);
                %otf0 = pepitoTools.moffat(beta,X_hr);
                otf0 = otf0/max(otf0(:));
                for jS = iS+1:nStars
                    iB = iB+1;
                    %image
                    psf(:,:,iB) = 0.5*(psfLE(:,:,jS,iS) + psfLE(:,:,iS,jS));
                    %otf
                    otf(:,:,iB) = otf0;
                end
            end
        end
        
        function fxfy= freq_array(nX,samp,D)
            pix2freq = 1.0/(D*samp);
            x = 1:nX;
            [fx,fy] = meshgrid(x);
            %central frequency set to 0
            fx = fx - (nX/2+1);
            fy = fy - (nX/2+1);
            fxfy = {fx*pix2freq,fy*pix2freq};
        end
        
        function out = multiplePSF(psfModel,pStars,nBox)
            
            if isscalar(nBox)
                nBox_x = nBox;
                nBox_y = nBox;
            else
                nBox_x = nBox(1);
                nBox_y = nBox(2);
            end
            
            % Unpack input
            nS   = length(pStars);
            if mod(nS,4) == 0
                nS = nS/4;                
                bg   = pStars(3*nS+1:end);                
            elseif mod(nS,3) == 0
                nS = 1;
                bg = 0;
            else
                nS = (nS-1)/3;                
                bg   = pStars(end)*ones(1,nS);        
            end
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
                tmp  = fluxS(iS)*pepitoTools.crop(pepitoTools.translateImage(psfi,dX(:,iS)),[nBox_x,nBox_y]);
                %tmp = fluxS(iS) * pepitoTools.crop(tools.shift(psfi,dX(1,iS),dX(2,iS)),[nBox_x,nBox_y]);
                % Concatenate frames
                out(:,:,iS) =  tmp + bg(iS);
            end
        end
                             
    end                    
end




        