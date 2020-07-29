classdef covarianceTools < handle
    
    methods (Static)
        
        function varargout = spatioAngularCovarianceMatrix(sampling,range,atm,srcAC,varargin)
            %% SPATIOANGULARCOVARIANCEMATRIX Phase spatio-angular covariance meta matrix
            %
            % S = spatioAngularCovarianceMatrix(sampling,range,atm,src1)
            % computes the spatio-angular auto-correlation meta-matrix of the
            % wavefront between all the sources srcAC. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % C = spatioAngularCovarianceMatrix(sampling,range,atm,src1,src2)
            % computes the spatio-angular cross-correlation meta-matrix of the
            % wavefront between all src2 and src1. The phase is
            % sampling with sampling points in the given range and
            % propagates through the atmosphere defined by the object atm
            %
            % [S,C] = spatioAngularCovarianceMatrix(...) computes both
            % auto- and cross-correlation meta-matrix
            %
            % ... = spatioAngularCovarianceMatrix(...,'mask',mask) restrict
            % the sampling within the mask
            
            inputs = inputParser;
            inputs.addRequired('sampling',@isnumeric);
            inputs.addRequired('range',@isnumeric);
            inputs.addRequired('atm',@isstruct);
            inputs.addRequired('srcAC',@isstruct);
            inputs.addOptional('srcCC',[],@isstruct);
            inputs.addOptional('tipTilt',false,@islogical);
            inputs.addParameter('mask',true(sampling),@islogical);
            inputs.addParameter('lag',0,@isnumeric);
            inputs.addParameter('xyOutput',[],@isnumeric);
            inputs.parse(sampling,range,atm,srcAC,varargin{:});
            
            m_srcCC = inputs.Results.srcCC;
            tipTilt = inputs.Results.tipTilt;
            m_mask  = inputs.Results.mask;
            m_tau   = inputs.Results.lag;
            xyOutput= inputs.Results.xyOutput;
            
            [m_x,m_y] = meshgrid( linspace(-1,1,sampling)*range/2 );
            m_nGs   = length(srcAC);
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            m_cst      = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                L0r0ratio;
            m_cstL0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            m_cstr0 = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                atm.r0.^(-5./3);
            m_L0      = atm.L0;
            
            m_nLayer   = atm.nLayer;
            m_altitude = [atm.heights];
            m_fr0      = [atm.weights];
            [m_windVx,m_windVy] = pol2cart([atm.windDirection],[atm.windSpeed]);
            m_srcACdirectionVector = cat(2,[srcAC.x;srcAC.y]);
            m_srcACheight          = [srcAC.height];
            
            if nargout==2
                
                tic
                varargout{1} = autoCorrelation(m_x,m_y,m_mask,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cstL0,m_cst);
                toc
                
                m_xAC = m_x(m_mask);
                m_yAC = m_y(m_mask);
                if isempty(xyOutput)
                    m_xCC = m_xAC;
                    m_yCC = m_yAC;
                else
                    m_xCC = xyOutput(:,1);
                    m_yCC = xyOutput(:,2);
                end
                
                tic
                varargout{2} = crossCorrelation(m_srcCC,m_xAC,m_yAC,m_xCC,m_yCC,...
                    m_nGs,m_srcACdirectionVector,m_srcACheight,...
                    m_nLayer,m_altitude,m_fr0,...
                    m_L0,m_cstL0,m_cst,m_windVx,m_windVy,m_tau);
                toc
                
            else
                
                if isempty(m_srcCC)
                    
                    varargout{1} = autoCorrelation(m_x,m_y,m_mask,...
                        m_nGs,m_srcACdirectionVector,m_srcACheight,...
                        m_nLayer,m_altitude,m_fr0,...
                        m_L0,m_cstL0,m_cst);
                    
                else
                    
                    if tipTilt
                        varargout{1} = crossCorrelationTT(m_srcCC,m_x,m_y,m_mask,...
                            m_nGs,m_srcACdirectionVector,m_srcACheight,...
                            m_nLayer,m_altitude,m_fr0,...
                            m_L0,m_cstr0,range);
                    else
                        m_xAC = m_x(m_mask);
                        m_yAC = m_y(m_mask);
                        if isempty(xyOutput)
                            m_xCC = m_xAC;
                            m_yCC = m_yAC;
                        else
                            m_xCC = xyOutput(:,1);
                            m_yCC = xyOutput(:,2);
                        end
                        varargout{1} = crossCorrelation(m_srcCC,m_xAC,m_yAC,m_xCC,m_yCC,...
                            m_nGs,m_srcACdirectionVector,m_srcACheight,...
                            m_nLayer,m_altitude,m_fr0,...
                            m_L0,m_cstL0,m_cst,m_windVx,m_windVy,m_tau);
                    end
                    
                end
                
            end
            
            function C = crossCorrelationTT(srcTT,x,y,mask,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstr0,D)
                
                fprintf(' -->> Cross-correlation meta-matrix (phase/tip-tilt)!\n')
                
                nSs = length(srcTT);
                srcTTdirectionVector = cat(2,srcTT.directionVector);
                C = cellfun( @(x) zeros(sum(mask(:))) , cell(nSs,nGs) , 'UniformOutput' , false );
                x = x(mask);
                y = y(mask);
                f02 = 1./L0.^2;
                
                for k=1:nSs*nGs
                    
                    [kSs,iGs] = ind2sub([nSs,nGs],k);
                    buf = 0;
                    
                    for kLayer=1:nLayer
                        
                        beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcACheight(iGs);
                        iZ = complex( x*scale + beta(1) , y*scale + beta(2) );
                        
                        betaTt = srcTTdirectionVector(:,kSs)*altitude(kLayer);
                        zTt = iZ.' - complex(betaTt(1),betaTt(2));
                        rTt = abs(zTt);
                        oTt = angle(zTt);
                        
                        Inm = @(r) quadgk( @(f) (f.^2 + f02).^(-11./6).*...
                            besselj(2,pi.*f.*D).*...
                            besselj(1,2.*pi.*f.*r),0,Inf);
                        out = fr0(kLayer)*cstr0.*arrayfun(Inm,rTt).*8/D/scale;
                        
                        buf = buf + ...
                            [ out.*cos(oTt) ; out.*sin(oTt) ];
                        
                    end
                    
                    C{k} = buf;
                    
                end
                
                buf = C;
                C = cell(nSs,1);
                for k=1:nSs
                    C{k} = cell2mat(buf(k,:));
                end
                
            end
            
            
            function C = crossCorrelation(srcCC,xAC,yAC,xCC,yCC,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstL0,cst,windVx,windVy,tau)
                
                fprintf(' -->> Cross-correlation meta-matrix!\n')
                
                nSs = length(srcCC);
                srcCCdirectionVector = cat(2,[srcCC.x;srcCC.y]);
                srcCCheight          = [srcCC.height];
                C = cellfun( @(x) zeros(length(xCC),length(xAC)) , cell(nSs,nGs) , 'UniformOutput' , false );
                cstL0CC = cstL0*ones(length(xCC),length(xAC));
                
                for k=1:nSs*nGs
                    
                    [kSs,iGs] = ind2sub([nSs,nGs],k);
                    buf = 0;
                    
                    for kLayer=1:nLayer
                        
                        
                        beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcACheight(iGs);
                        iZ = complex( xAC*scale + beta(1) , yAC*scale + beta(2) );
                        
                        betaSs = srcCCdirectionVector(:,kSs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcCCheight(kSs);
                        zSs = complex( ...
                            xCC*scale + betaSs(1) + windVx(kLayer)*tau, ...
                            yCC*scale + betaSs(2) + windVy(kLayer)*tau );
                        
                        rho   = abs(bsxfun(@minus,zSs,iZ.'));
                        out   = cstL0CC;
                        index = rho~=0;
                        u          = 2.*pi.*rho(index)./L0;
                        out(index) = cst.*u.^(5./6).*besselk(5./6,u);
                        
                        buf = buf + fr0(kLayer)*out;
                        
                    end
                    
                    C{k} = buf;
                    
                end
                
                buf = C;
                C = cell(nSs,1);
                for k=1:nSs
                    C{k} = cell2mat(buf(k,:));
                end
                
            end
            
            function S = autoCorrelation(x,y,mask,...
                    nGs,srcACdirectionVector,srcACheight,...
                    nLayer,altitude,fr0,...
                    L0,cstL0,cst)
                
                fprintf(' -->> Auto-correlation meta-matrix!\n')
                
                kGs = reshape( triu( reshape(1:nGs^2,nGs,nGs) , 1) , 1, []);
                kGs(1) = 1;
                kGs(kGs==0) = [];
                S = cellfun( @(x) zeros(sum(mask(:))) , cell(1,length(kGs)) , 'UniformOutput' , false );
                cstL0AC = cstL0*ones(sampling);
                
                for k=1:length(kGs)
                    
                    [iGs,jGs] = ind2sub( [nGs,nGs] , kGs(k) );
                    buf = 0;
                    
                    for kLayer=1:nLayer
                        
                        %                         atmSlab = slab(atm,kLayer);
                        
                        beta = srcACdirectionVector(:,iGs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcACheight(iGs);
                        iZ = complex( x*scale + beta(1) , y*scale + beta(2) );
                        
                        beta = srcACdirectionVector(:,jGs)*altitude(kLayer);
                        scale = 1 - altitude(kLayer)/srcACheight(jGs);
                        jZ  = complex( x*scale + beta(1) , y*scale + beta(2) );
                        
                        z1 = iZ;
                        z2 = jZ;
                        [nz,mz] = size( z1 );
                        % First Row
                        %                         r  =  phaseStats.covariance( abs(z2-z1(1)) , atm );
                        rho = abs(z2-z1(1));
                        r   = cstL0AC;
                        index         = rho~=0;
                        u             = 2.*pi.*rho(index)./L0;
                        r(index) = cst.*u.^(5./6).*besselk(5./6,u);
                        
                        r   = mat2cell( fr0(kLayer)*r  , nz , ones(mz,1));
                        % First column in first blocks fow
                        %                         c  =  phaseStats.covariance( abs( bsxfun( @minus , z2(1:nz:nz^2), z1(1:nz).' ) ) , atm );
                        rho = abs( bsxfun( @minus , z2(1:nz:nz^2), z1(1:nz).' ) );
                        c   = cstL0AC;
                        index         = rho~=0;
                        u             = 2.*pi.*rho(index)./L0;
                        c(index) = cst.*u.^(5./6).*besselk(5./6,u);
                        
                        c   = mat2cell( fr0(kLayer)*c , nz , ones(mz,1) );
                        % First block rows
                        rr   = cellfun( @(x,y) pepitoTools.toeplitz(x,y) , c , r , 'UniformOutput' , false);
                        
                        % First Column
                        %                         c  =  phaseStats.covariance( abs(z1-z2(1)) , atm );
                        rho = abs(z1-z2(1));
                        c   = cstL0AC;
                        index         = rho~=0;
                        u             = 2.*pi.*rho(index)./L0;
                        c(index) = cst.*u.^(5./6).*besselk(5./6,u);
                        
                        c   = mat2cell( fr0(kLayer)*c  , nz , ones(mz,1));
                        % First row in first blocks column
                        %                         r  =  phaseStats.covariance( abs( bsxfun( @minus , z1(1:nz:nz^2), z2(1:nz).' ) ) , atm );
                        rho = abs( bsxfun( @minus , z1(1:nz:nz^2), z2(1:nz).' ) );
                        r   = cstL0AC;
                        index         = rho~=0;
                        u             = 2.*pi.*rho(index)./L0;
                        r(index) = cst.*u.^(5./6).*besselk(5./6,u);
                        
                        r   = mat2cell( fr0(kLayer)*r , nz , ones(mz,1) );
                        % First blocks column
                        cc   = cellfun( @(x,y) pepitoTools.toeplitz(x,y) , c , r , 'UniformOutput' , false);
                        
                        out = cell2mat( pepitoTools.toeplitz(cc,rr) );
                        out(~mask,:) = [];
                        out(:,~mask) = [];
                        
                        buf = buf + out;
                        
                    end
                    
                    S{k} = buf;
                    
                end
                
                buf = S;
                S = cellfun( @(x) zeros(sum(mask(:))) , cell(nGs) , 'UniformOutput' , false );
                S(kGs) = buf;
                S(1:nGs+1:nGs^2) = S(1,1);
                S = cell2mat(S);
                S = triu(S,1)+triu(S)';
                
            end
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
                T{k} = pepitoTools.toeplitz(n2:-1:1,n2:n)+n*(k-1);
            end
            T = cell2mat(pepitoTools.toeplitz(T(n2:-1:1),T(n2:n)));
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
        
        function mat = covMap2Matrix(map,n1,n2,varargin)
            %% convert covariance Map to covariance Matrix
            
            inputs = inputParser;
            inputs.addRequired('mat',@ismatrix);
            inputs.addRequired('n1',@isnumeric);
            inputs.addRequired('n2',@isnumeric);
            inputs.addParameter('mask1', [], @islogical );
            inputs.addParameter('mask2', [], @islogical );
            inputs.parse(map,n1,n2,varargin{:});
            
            mask1 = inputs.Results.mask1;
            if isempty(mask1)
                mask1 = true(n1);
            end
            mask2 = inputs.Results.mask2;
            if isempty(mask2)
                mask2 = true(n2);
            end
            
            n = n1+n2-1;
            T = mat2cell(map,n,ones(1,n));
            c = cellfun(@(x) x(n2:-1:1),T,'UniformOutput',false);
            r = cellfun(@(x) x(n2:n),T,'UniformOutput',false);
            T=cellfun(@(x,y) pepitoTools.toeplitz(x,y),c,r,'UniformOutput',false);
            T=pepitoTools.toeplitz(T(n2:-1:1),T(n2:n));
            mat = cell2mat(T);
            
            mat(~mask2(:),:) = [];
            mat(:,~mask1(:)) = [];
            
        end
        
        function out = mcDonald(x)
            out = x.^(5/6.).*besselk(5./6,x)./(2^(5/6.)*gamma(11/6.)) ;
            out(find(x == 0)) = 3/5.;
        end
        
         function out = sombrero(n,x)
            %% SOMBRERO Order n sombrero function
            %
            % out = sombrero(n,x) computes besselj(n,x)/x
            
            if n==0
                out = besselj(0,x)./x;
            else
                if n>1
                    out = zeros(size(x));
                else
                    out = 0.5*ones(size(x));
                end
                u = x~=0;
                x = x(u);
                out(u) = besselj(n,x)./x;
            end
        end
    end
end
