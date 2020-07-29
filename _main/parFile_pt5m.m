% SOURCES
xStars_0 = [122 158 449 524];% 52 449];
yStars_0 = [703 708 464 325];% 272 46];

% ATM
if fitCn2
    r0      = 0.15;
    L0      = 16;
    height  = [0 16e3];%0.5 1 2 4 8 16]*1e3;
    weight  = [0.7 0.3];%0.1 0.1 0.05 0.15 0.1 0.1];
    windSpeed = [5 15];%5.5 6 10 8 25 15];
    windDir  = [0, pi/6];%, pi/2, pi/3, -pi/2, pi/6, -pi];
else
    r0      = 0.1;
    L0      = 16;
    height  = [0 10]*1e3;
    weight  = [0.7 0.3];
    windSpeed = [5 20];
    windDir  = [0,pi/4];
end

% TELESCOPE
D       = 0.5;
cobs    = 0.47;

% CAMERA
wvl     = 500e-9;
psInMas = 495.6;
%495.6;
Samp    = constants.radian2mas*wvl/D/psInMas;
ron     = 7;
gain    = 1.3;

% RESOLUTION
nPSF   = 26; %PSF
nPup   = fix(ceil(2.0/Samp)) * nPSF; % pupil


%0.0104    0.0173    0.0075    0.0094