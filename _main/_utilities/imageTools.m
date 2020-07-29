classdef imageTools < handle
    
    % Tools facility for image manipulation
    
    methods (Static)
          
        function out = normalize(input,method)
            if strcmp(method,'max')
                out = input/max(input(:));
            else
                out = input/sum(input(:));
            end
        end
        function makeAxisSquared(n)
            for i=1:numel(n)
                axesHandles = findobj(get(figure(n(i)),'Children'), 'flat','Type','axes');
                % Set the axis property to square
                axis(axesHandles,'square');
            end
        end
                        
        function out = crop(input, ncrop)
            
            if isvector(input)
                out = input(1:ncrop);
            else
                nFrames = size(input,3);
                
                if isscalar(ncrop)
                    ncrop = [ncrop ncrop];
                end
                
                dim = size(input);
                out = zeros(ncrop(1), ncrop(2), nFrames);
                for iFrame = 1:nFrames
                    if all(ncrop<dim)
                        deb = round((dim - ncrop) / 2 + 1);
                        fin = round((dim + ncrop) / 2);
                        out(:,:,iFrame) = input(deb(1):fin(1), deb(2):fin(2), iFrame);
                    else
                        deb = round((ncrop-dim) / 2 + 1);
                        fin = round((ncrop+dim) / 2);
                        out(deb(1):fin(1), deb(1):fin(1), iFrame) = input(:,:,iFrame);
                    end
                end              
            end
        end
        
        function frame = createDeadPixFrame(badPixelMap)
            %dpframe = createDeadPixFrame(badPixelMap)
            %badPixelMap is the map of dead pixels
            %frame  is the image to be corrected
            
            %The dead pixel is replaced by a weighted average of the neighbours,
            %1 2 1
            %2 X 2
            %1 2 1
            %when they are "available". "Available" means that the sum of the
            %weights of neighbouring pixels must exceeds 4.
            
            %If no neighbouring pixel is available, the dead pixel is not
            %corrected, but a new "dead pixel map" is created, and the function is
            %called once again (recursive calls).
            
            % Get the number of dead pixels
            [sx,sy]      = size(badPixelMap);
            npixnoncorr  = 0;
            [nnx,nny]    = find(badPixelMap);
            nn1D         = find(badPixelMap(:));
            nDeadPix     = length(nn1D);
            %Instantiation
            tmp          = badPixelMap*0;
            frame        = zeros(nDeadPix,10,2); %2nd row: #pixel (one pixel + 8 neighbors)
            frame(:,:,1) = 1;                    %3rd row: adresses
            
            %loop on Pixel
            for i=1:nDeadPix
                nb = 2;
                frame(i,1,1) = nn1D(i);  % 1st row = bad pixel
                frame(i,2,1) = 0;        % number of used neighbour pixel for correction
                x            = nnx(i);
                y            = nny(i);
                wcum         = 0;
                
                % Edges neighbours
                if x>0 && x<=sx && y+1>0 && y+1<=sy
                    if ~badPixelMap(x,y+1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i) + sx;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                
                if x>0 && x<=sx && y-1>0 && y-1<=sy
                    if ~badPixelMap(x,y-1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-sx;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                if x+1>0 && x+1<=sx && y>0 && y<=sy
                    if ~badPixelMap(x+1,y)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)+1;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                if x-1>0 && x-1<=sx && y>0 && y<=sy
                    if ~badPixelMap(x-1,y)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-1;
                        frame(i,nb,2) = 2;
                        wcum          = wcum + 2;
                    end
                end
                
                % Diagonal neighbours
                if x+1>0 && x+1<=sx && y+1>0 && y+1<=sy
                    if ~badPixelMap(x+1,y+1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)+1+sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                if x-1>0 && x-1<=sx && y+1>0 && y+1<=sy
                    if ~badPixelMap(x-1,y+1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-1+sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                if x+1>0 && x+1<=sx && y-1>0 && y-1<=sy
                    if~badPixelMap(x+1,y-1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)+1-sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                if x-1>0 && x-1<=sx && y-1>0 && y-1<=sy
                    if~badPixelMap(x-1,y-1)
                        nb            = nb+1;
                        frame(i,nb,1) = nn1D(i)-1-sx;
                        frame(i,nb,2) = 1;
                        wcum          = wcum + 1;
                    end
                end
                
                % Take decision regarding the number of avalaible neighbours
                if wcum<4   %not enough neigbours
                    npixnoncorr      = npixnoncorr + 1;
                    tmp(x,y)         = tmp(x,y) + 1;
                    frame(i,3:end,1) = 1;    % pixel adresses set to 1
                    frame(i,:,2)     = 0;    % weights set to 0
                else
                    frame(i,2,1)     = nb;    %number of correcting pixels
                end
            end
            
            if npixnoncorr
                frame_suppl = pepitoTools.createDeadPixFrame(tmp);
                nSup        = size(frame_suppl,1);
                %Frame concatenation
                new_frame                     = zeros(nDeadPix+nSup,10,2);
                new_frame(1:nDeadPix,:,:)     = frame;
                new_frame(1+nDeadPix:end,:,:) = frame_suppl;
                frame                         = new_frame;
            end
            
        end
        
        function imCor = corrDeadPixFrame( frame, im )
            % imCor = corrDeadPixFrame(frame,im)
            %Correcting the bad pixel on a image from the bad pixel map
            
            npixd = size(frame,1);
            imCor = im;
            for i=1:npixd
                w =  sum(frame(i,3:end,2));
                if w~=0
                    imCor(frame(i,1,1)) = sum( im(frame(i,3:end,1)) .* frame(i,3:end,2)) / double(w);
                end
            end
        end
        
        function imagerot = rotateIm(im_origin,angleInDegree)
            switch mod(angleInDegree, 360)
                % Special cases
                case 0
                    imagerot = im_origin;
                case 90
                    imagerot = rot90(im_origin);
                case 180
                    imagerot = im_origin(end:-1:1, end:-1:1);
                case 270
                    imagerot = rot90(im_origin(end:-1:1, end:-1:1));
                otherwise
                    % Zero-pad the image
                    [Rows, Cols] = size(im_origin);
                    Diagonal = hypot(Rows,Cols);
                    RowPad = ceil(Diagonal - Rows) + 2;
                    ColPad = ceil(Diagonal - Cols) + 2;
                    im_pad = zeros(Rows+RowPad, Cols+ColPad);
                    im_pad(ceil(RowPad/2):(ceil(RowPad/2)+Rows-1),ceil(ColPad/2):(ceil(ColPad/2)+Cols-1)) = im_origin;
                    
                    % Convert to radians and create transformation matrix
                    imagerot=zeros(size(im_pad)); % midx and midy same for both
                    angleInRad = angleInDegree*pi/180;
                    midx=ceil((size(im_pad,1)+1)/2);
                    midy=ceil((size(im_pad,2)+1)/2);
                    
                    for i=1:size(im_pad,1)
                        for j=1:size(im_pad,2)
                            
                            x= (i-midx)*cos(angleInRad)+(j-midy)*sin(angleInRad);
                            y=-(i-midx)*sin(angleInRad)+(j-midy)*cos(angleInRad);
                            x=round(x)+midx;
                            y=round(y)+midy;
                            
                            if (x>=1 && y>=1 && x<=size(im_pad,2) && y<=size(im_pad,1))
                                imagerot(i,j)=im_pad(x,y); % k degrees rotated image
                            end
                        end
                    end
                    imagerot = pepitoTools.crop(imagerot,size(im_origin));
            end
        end
        
        function out = enlargePupil(P,n)
            % TO BE MERGED WITH enlargeOtf
            [nx,ny,nf] = size(P);
            %if n*nx ~= round(n*nx)
            %    fprintf('Warning: non integer support, ratio adjusted to %f\n',round(n*nx)/nx);
            %end
            xi      = floor((nx*n-nx)/2+1);
            xf      = floor((nx*n+nx)/2);
            yi      = floor((ny*n-ny)/2+1);
            yf      = floor((ny*n+ny)/2);
            if length(size(P)) == 2
                out  = zeros(round(n*nx),round(ny*n));
                out(xi:xf,yi:yf) = P;
            elseif length(size(P)) == 3
                out = zeros(round(n*nx),round(ny*n),nf);
                out(xi:xf,yi:yf,:) = P;
            end
        end
        
        function out = enlargeOtf(otf,n)
            if true
                out = padarray(otf,floor((n-1)*size(otf)/2),'both');
            else
                % Otf sizes
                nx  = size(otf,1);
                nx2 = round(n*nx);
                out = zeros(nx2);
                
                % Zero-padding
                if mod(nx2,2) == 0
                    idx = floor(0.5*(nx2-nx) + 1): floor(0.5*(nx2 + nx));
                else
                    idx = floor(0.5*(nx2-nx)+1: 0.5*(nx2+nx));
                end
                idx = floor(nx2/2 + 1 - nx/2): floor(nx2/2 + nx/2);
                out(idx,idx) = otf;
            end
        end
        
        function out = interpolateOtf(otf,nRes,method)
            if nargin < 3
                method = 'spline';
            end
            
            if numel(nRes) == 1
                nRes = [nRes nRes];
            end
            
            % Define input meshgrid
            notf_x = size(otf,1);
            notf_y = size(otf,2);
            if mod(notf_x,2)==0
                u1D = (-notf_x/2:1:notf_x/2-1)*2/notf_x;
            else
                u1D = (-floor(notf_x/2):1:floor(notf_x/2))*2/notf_x;
            end
            if mod(notf_y,2)==0
                v1D = (-notf_y/2:1:notf_y/2-1)*2/notf_y;
            else
                v1D = (-floor(notf_y/2):1:floor(notf_y/2))*2/notf_y;
            end
            [Uxi,Uyi]= meshgrid(u1D,v1D);
            
            % Define output meshgrid
            nRes_x = nRes(1);
            nRes_y = nRes(2);
            
            if mod(nRes_x,2)==0
                u1D2 = (-nRes_x/2:1:nRes_x/2-1)*2/nRes_x;
            else
                u1D2 = (-floor(nRes_x/2):1:floor(nRes_x/2))*2/nRes_x;
            end
            
            if mod(nRes_y,2)==0
                v1D2 = (-nRes_y/2:1:nRes_y/2-1)*2/nRes_y;
            else
                v1D2 = (-floor(nRes_y/2):1:floor(nRes_y/2))*2/nRes_y;
            end
            [Uxo,Uyo]= meshgrid(u1D2,v1D2);
            
            % Interpolation
            out = interp2(Uxi,Uyi,otf,Uxo,Uyo,method);
        end
        
        function out = interpolate(P,nRes,method)
            
            if nargin<3
                method='linear';
            end
            
            if numel(nRes) == 1
                nRes = [nRes,nRes];
            end
            
            nP_x       = size(P,1);
            nP_y       = size(P,2);
            xi       = linspace(-1,1,nP_x);
            yi       = linspace(-1,1,nP_y);
            xo       = linspace(-1,1,nRes(1));
            yo       = linspace(-1,1,nRes(2));
            [Xi,Yi]  = meshgrid(xi,yi);
            [Xo,Yo]  = meshgrid(xo,yo);
            out      = interp2(Xi,Yi,P,Xo,Yo,method);
        end
        
        function im_ = translateImage(image,dX,val)
            if nargin < 3
                val = 0;
            end
            
            dx = dX(1);
            dy = dX(2);
            
            % Zero-pad the suppot to avoid fft buffer circulation effect
            [nx,ny] = size(image);
            % nx -> rows
            % ny -> columns
            im_ = padarray(image,floor([nx,ny]/2),'both');
            [ny_,nx_] = size(im_);
            
            % Derives the fourier phasor            
            [u,v] = freqspace([ny_,nx_],'meshgrid');
            phasor = exp(-1i*pi*(u*dy+v*dx));
          
            % Translates the image
            im_ = real( ifft2( fft2( im_ ) .* fftshift(phasor) ) );
            %otf = pepitoTools.psf2otf(im_);
            %otf = otf/max(otf(:));
            %im_ = pepitoTools.otf2psf( otf .* phasor );
            
            %0.0104    0.0173    0.0075    0.0094
            % Empty zone filling
            if any(size(im_)>size(image))
                ixi = round(1:nx/2 + dx);
                ixf = round(3*nx/2+dx+1):size(im_,1);
                iyi = round(1:ny/2 + dy);
                iyf = round(3*ny/2+dy+1):size(im_,2);
                im_(ixi,:) = val;
                im_(ixf,:) = val;
                im_(:,iyi) = val;
                im_(:,iyf) = val;
                % Cropping
                im_ = pepitoTools.crop(im_,[nx,ny]);
            end
            
            im_ = im_/sum(im_(:))*sum(image(:));

            
%             % Zero-pad the support 
%             [ny,nx] = size(image);
%             % nx -> rows
%             % ny -> columns
%             im_ = padarray(image,floor([ny,nx]/2),'both');
%             [ny_,nx_] = size(im_);
%             
%             % Derives the fourier phasor
%             dx = dX(1);
%             dy = dX(2);            
%             py = exp(-2i*pi*dy*(0:ny_-1)/ny_);
%             px = exp(-2i*pi*dx*(0:nx_-1)/nx_);
%             phasor = py' * px;
%             
%             % Translates the image
%             %otf = fft2( ( im_ ) );
%             im_ = real( fft2( fft2( im_ ) .* phasor ) );
%                       
%   
%             % Empty zone filling
%             if any(size(im_)>size(image))
%                 ixi = round(1:nx/2 + dx);
%                 ixf = round(3*nx/2+dx+1):size(im_,2);
%                 iyi = round(1:ny/2 + dy);
%                 iyf = round(3*ny/2+dy+1):size(im_,1);
%                 im_(iyi,:) = val;
%                 im_(iyf,:) = val;
%                 im_(:,ixi) = val;
%                 im_(:,ixf) = val;
%                 % Cropping
%                 im_ = pepitoTools.crop(im_,[ny,nx]);
%             end
%             
%             % Normalization
%             im_ = im_/sum(im_(:))*sum(image(:));
        end
        
        function [imCor,otf_lr] = recenterPSF(psf,overSampling)
            flux          = sum(psf(:));
            [npsfx,npsfy] = size(psf);
            npsfx2        = npsfx*overSampling;
            npsfy2        = npsfy*overSampling;
            
            % Get the high-resolution PSF
            if overSampling > 1
                psf_hr = pepitoTools.interpolateOtf(psf,npsfx2);
            else
                psf_hr = psf;
            end
            % Get the max value
            mx        = max(psf_hr(:));
            [idx,idy] = find(psf_hr == mx);
            idx = idx(1);
            idy = idy(1);
            dx        = floor(npsfx2/2-idx)+1;
            dy        = floor(npsfy2/2-idy)+1;
            if (dx~=0) | (dy~=0)
                % Get the OTF
                otf_hr = pepitoTools.psf2otf(psf_hr);
                % Apply the Phasor
                [u,v]     = freqspace(length(otf_hr),'meshgrid');
                fftPhasor = exp(-1i.*pi.*(u*dy+v*dx));
                otf_hr    = otf_hr.*fftPhasor;
                % Get the PSF low-resolution
                imCor  = pepitoTools.otf2psf(otf_hr);
                imCor  = pepitoTools.interpolateOtf(imCor,npsfx);
                imCor  = flux*imCor/sum(imCor(:));
                otf_lr = pepitoTools.psf2otf(imCor);
            else
                imCor = psf;
                otf_lr = pepitoTools.psf2otf(imCor);
            end
        end
        
        function [x0, y0] = cog(im, threshold)
            
            %% Check for Threshold
            if nargin == 1
            else
                ind = im < threshold;
                im(ind) = 0.0;
            end
            % Calculate Centre of Gravity
            [rc, cc] = ndgrid(1:size(im, 1), 1:size(im, 2));
            Mt = sum(im(:));
            x0 = sum(im(:) .* rc(:))/Mt;
            y0 = sum(im(:) .* cc(:))/Mt;
        end
    end
end