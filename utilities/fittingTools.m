classdef fittingTools < handle
    
    methods (Static)
        
        % PSF-fitting
        function [out,imFit,beta] = findStellarParameters(im,psfModel,initParam)
            %1. Model
            nIm          = size(im);                     
            normFactor   = sum(im(:));
            y = im/normFactor;
            
            %weights
            %[F,bg,ron] = pepitoTools.getFlux(y);
            %im_tmp = im_tmp-bg;
            %w = 1./sqrt(max(im_tmp,ron) + ron^2);
            %y = im_tmp.*w.*(im_tmp>0);
           
            stellarModel = @(x,xdata) pepitoTools.stellarFieldModel(x,xdata,nIm);
            nParam = length(initParam);
            %2. Fitting options and initial guess
            if mod(nParam,3)
                % maybe a background estimation
                if mod(nParam-1,3)
                    fprintf('You must provide a valid initial guess vector, which should size 3*nS or 3*nS +1 with nS the number of stars\n');
                    return
                    out = [];
                    imFit = [];
                    beta = [];
                else
                    nS    = (length(initParam)-1)/3;
                    lX    = initParam(1:2*nS) - 5;
                    uX    = initParam(1:2*nS) + 5;
                    lF    = zeros(1,nS);
                    uF    = 10*ones(1,nS);
                    lb    = [lX,lF,-5*std(y(:))];
                    ub    = [uX,uF,5*std(y(:))];
                end
            else
                % no background estimation
                nS    = (length(initParam))/3;
                lX    = initParam(1:2*nS) - 5;
                uX    = initParam(1:2*nS) + 5;
                lF    = zeros(1,nS);
                uF    = 10*ones(1,nS);
                lb    = [lX,lF];
                ub    = [uX,uF];
            end
            
            opt = optimoptions(@lsqcurvefit,'MaxIter',3e2,...
                'TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',3e2,...
                'InitDamping',1,'Display','iter');
            
            %3. Fitting procedure
            [beta,~,R,~,~,~,J] = lsqcurvefit(stellarModel,initParam,psfModel,...
                y,lb,ub,opt);
            
            %4. Unpacking results
            xS = beta(1:nS);
            yS = beta(nS+1:2*nS);
            fS = beta(2*nS+1:3*nS)*normFactor;
            
            %5. Measurements uncertainties
            dbeta = diff(nlparci(beta,R,'jacobian',J),1,2);
            dX    = dbeta(1:nS);
            dY    = dbeta(nS+1:2*nS);
            dF    = dbeta(2*nS+1:3*nS)*normFactor;
            
            % Concatenating results
            
            out = zeros(nParam,nS);
            for iS=1:nS
                out(1,iS) = xS(iS);
                out(2,iS) = yS(iS);
                out(3,iS) = fS(iS);
                out(4,iS) = dX(iS);
                out(5,iS) = dY(iS);
                out(6,iS) = dF(iS);
                if mod(nParam,3)
                    out(7,iS) = beta(end)*normFactor;
                end
            end
            
            % Fitted model
            if ~mod(nParam,3)
                res = [xS,yS,fS];
            else
                res = [xS,yS,fS,beta(end)*normFactor];
            end
            
            imFit =  pepitoTools.stellarFieldModel(res,psfModel,nIm);
        end
        
        function out = stellarFieldModel(pStars,psfModel,nIm)
            
            nParam = length(pStars);
            if ~mod(nParam,3)                
                nS = nParam/3;
                xS = pStars(1:nS);
                yS = pStars(1+nS:2*nS);
                fS = pStars(1+2*nS:3*nS);
                bg = 0;
            else
                nS = (nParam-1)/3;
                xS = pStars(1:nS);
                yS = pStars(1+nS:2*nS);
                fS = pStars(1+2*nS:3*nS);
                bg = pStars(end);
            end
            
            out = 0*psfModel;
            for iS = 1:nS
                % Translate the PSF
                psf_i = pepitoTools.translateImage(psfModel,[xS(iS),yS(iS)]);
                % Update the image and Flux scaling:
                out = out + psf_i*fS(iS);
            end
            out = pepitoTools.crop(out,nIm) + bg;
        end
        
        function out = multipleImage(x,im)
            %Grab locations
            n       = length(x)/3;
            [nx,ny] = size(im);
            xloc    = x(1:n);
            yloc    = x(n+1:2*n);
            flux    = x(2*n+1:end);
            
            if any(abs(xloc) > nx) || any(abs(yloc) > ny)
                id = abs(xloc) > nx | abs(yloc) > ny;
                nn = sum(id);
                if nn == 1
                    fprintf('Warning: 1 PSF are out the field\n');
                else
                    fprintf('Warning: %d PSFs are out the field\n',nn);
                end
                flux(id) = 0;
            end
            
            % For xloc=yloc=0 and flux=1, the procedure does modify the PSF
            % while it mustn't. To be check.
            %nx      = 2*nx;
            %ny      = 2*ny;
            out     = zeros(nx,ny);
            otf     = pepitoTools.psf2otf(im);
            otf     = otf/max(otf(:));
            [u,v]   = freqspace([nx ny],'meshgrid');
            for i=1:n
                yi    = yloc(i);
                xi    = xloc(i);
                fftPhasor = exp(-1i*pi*(u*xi + v*yi));
                map_i = pepitoTools.otf2psf(otf.*fftPhasor);
                % Adding images
                out = out + flux(i)*map_i/sum(map_i(:));
            end
        end
        
        % Gaussian/Moffat-fitting
        function out = gaussian(x,xdata)
            % ------- Grabbing parameters ---------%
            I0 = x(1);          %Amplitude
            ax = x(2);          %x spreading
            ay = x(3);          %y-spreading
            th = x(4)*pi/180.;  %rotation
            x0 = x(5);          %x-shift
            y0 = x(6);          %y-shift
            
            % ------- Including shifts ---------
            X     = xdata{1};
            Y     = xdata{2};
            %Shifts
            X     = X - x0;
            Y     = Y - y0;
            %Rotation
            Xr    = X.*cos(th) + Y.*sin(th);
            Yr    = Y.*cos(th) - X.*sin(th);
            
            % Gaussian expression
            out = I0.*exp(-0.5*((Xr./ax).^2 + (Yr./ay).^2) );
        end
        
        function out = moffat(x,xdata,norm)
            
            % ------- Grabbing parameters ---------%
            I0 = x(1);          %Amplitude
            ax = x(2);          %x spreading
            ay = x(3);          %y-spreading
            be = x(4);          %center slope
            th = x(5);          %rotation
            x0 = x(6);          %x-shift
            y0 = x(7);          %y-shift
            
            % ------- Including shifts ---------
            X     = xdata{1};
            Y     = xdata{2};
            %Shifts
            X     = X - x0;
            Y     = Y - y0;
            %Rotation
            Xr    = X.*cos(th) + Y.*sin(th);
            Yr    = Y.*cos(th) - X.*sin(th);
            
            % Moffat expression
            out = I0.*(1. + (Xr./ax).^2 + (Yr./ay).^2).^(-be);
            %Normalization
            if exist('norm','var')
                out = out*(be-1)/pi/ax/ay;
            end
        end
        
          function out = multiAnalyticIsoplanatic(x,xdata,type)
            
            if strcmp(type,'gaussian')
                f = @(x,xdata) pepitoTools.gaussian(x,xdata);
                nS = (length(x)-3)/3;
                xfix = x(1+nS:3+nS);
                fS   =  x(1:nS);
                xS = x(4+nS:3+2*nS);
                yS = x(4+2*nS:end);
            else
                f = @(x,xdata) pepitoTools.moffat(x,xdata);
                nS = (length(x)-4)/3;
                xfix = x(1+nS:4+nS);
                fS   =  x(1:nS);
                xS = x(5+nS:4+2*nS);
                yS = x(5+2*nS:end);
            end
            
            
            out = zeros(size(xdata,1));
            for k=1:nS
                out = out + f([fS(k) xfix,xS(k) yS(k)],xdata);
            end
          end
        
        function out = multiAnalyticAniso(x,xdata,type)
            
            if strcmp(type,'gaussian')
                f = @(x,xdata) pepitoTools.gaussian(x,xdata);
                nG  = length(x)/6;
            else
                f = @(x,xdata) pepitoTools.moffat(x,xdata);
                nG  = length(x)/7;
            end
            
            
            out = zeros(size(xdata,1));
            for i=1:nG
                out = out + f(x(1+(i-1)*6:i*6),xdata);
            end
        end
        
        % Data fitting
        function p = polyfitweighted(x,y,n,w)
            % polyfitweighted.m
            % -----------------
            %
            % Find a least-squares fit of 1D data y(x) with an nth order
            % polynomial, weighted by w(x).
            %
            % By S.S. Rogers (2006), based on polyfit.m by The MathWorks, Inc. - see doc
            % polyfit for more details.
            %
            % Usage
            % -----
            %
            % P = polyfitweighted(X,Y,N,W) finds the coefficients of a polynomial
            % P(X) of degree N that fits the data Y best in a least-squares sense. P
            % is a row vector of length N+1 containing the polynomial coefficients in
            % descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1). W is
            % a vector of weights.
            %
            % Vectors X,Y,W must be the same length.
            %
            % Class support for inputs X,Y,W:
            %    float: double, single
            %
            
            % The regression problem is formulated in matrix format as:
            %
            %    yw = V*p    or
            %
            %          3    2
            %    yw = [x w  x w  xw  w] [p3
            %                            p2
            %                            p1
            %                            p0]
            %
            % where the vector p contains the coefficients to be found.  For a
            % 7th order polynomial, matrix V would be:
            %
            % V = [w.*x.^7 w.*x.^6 w.*x.^5 w.*x.^4 w.*x.^3 w.*x.^2 w.*x w];
            
            if ~isequal(size(x),size(y),size(w))
                error('X and Y vectors must be the same size.')
            end
            
            x = x(:);
            y = y(:);
            w = w(:);
            
            
            % Construct weighted Vandermonde matrix.
            %V(:,n+1) = ones(length(x),1,class(x));
            V(:,n+1) = w;
            for j = n:-1:1
                V(:,j) = x.*V(:,j+1);
            end
            
            % Solve least squares problem.
            [Q,R] = qr(V,0);
            ws = warning('off','all');
            p = R\(Q'*(w.*y));    % Same as p = V\(w.*y);
            
            warning(ws);
            if size(R,2) > size(R,1)
                warning('polyfitweighted:PolyNotUnique', ...
                    'Polynomial is not unique; degree >= number of data points.')
            elseif condest(R) > 1.0e10
                if nargout > 2
                    warning('polyfitweighted:RepeatedPoints', ...
                        'Polynomial is badly conditioned. Remove repeated data points.')
                else
                    warning('polyfitweighted:RepeatedPointsOrRescale', ...
                        ['Polynomial is badly conditioned. Remove repeated data points\n' ...
                        '         or try centering and scaling as described in HELP POLYFIT.'])
                end
            end
            p = p.';          % Polynomial coefficients are row vectors by convention.
        end
        
    end
end