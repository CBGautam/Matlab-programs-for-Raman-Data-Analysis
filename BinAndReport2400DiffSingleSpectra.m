clear all;

%% Datafiles

%load('C:\Users\j_r847\Desktop\3percent\SynthRefference3Vars','H','M');
%load('C:\Users\j_r847\Desktop\RajuSamples\Dia_OnSi120218Vars','H','M');
%load('U:\HoltzGroup\Raju\StressMeasurementsJRA\RajuDiamondShotsJRA\3percent\D170403Cross300s_fullCCD_6_1Vars','H','M');
%load('U:\HoltzGroup\Raju\StressMeasurementsJRA\RajuDiamondShotsJRA\4_5Percent\D1704010Cross300s_fullCCD_1_1Vars','H','M');
%load('U:\HoltzGroup\Raju\StressMeasurementsJRA\RajuDiamondShotsJRA\4_5Percent\D1704010Cross120s_fullCCD_2_1Vars','H','M');
%load('U:\HoltzGroup\Raju\StressMeasurementsJRA\RajuDiamondShotsJRA\4_5Percent\D1704010Plane300s_fullCCD_9_1Vars','H','M');
%load('U:\HoltzGroup\Raju\StressMeasuremntsRetest\4_5_8hr\TextFiles\4_5_8hr_120sFullChipCross_0_1Vars','H','M');
%load('U:\HoltzGroup\Raju\StressMeasuremntsRetest\4_5_8hr\TextFiles\4_5_8hr_120sFullChipPlane_0_1Vars','H','M');
%load('U:\HoltzGroup\Raju\StressMeasuremntsRetest\4_5_8hr\TextFiles\4_5_8hr_240sFullChipCross_0_1Vars','H','M');
load('U:\HoltzGroup\Raju\StressMeasuremntsRetest\4_5_8hr\TextFiles\4_5_8hr_40sSpectra_86_1Vars','H','M');


%% Physical Constants

totpixels=400; %total amount of strips
CCDDistance=60; %image size in Microns;
pixsize= CCDDistance/totpixels;

%% 

figure

waterfall(H);


Curve2=[];



%for strip=10:399
strip = 320;
aproxX=1336;

step=0.6155;     %CCD datapoint separation
t=1709.952;     %Max Wavenumber
for n= 1:901
  
XX(n)=t;  %X axis values
t=t-step;

end





%% Binning  Setup

numbin=1;  %Number of Spectra per Bin
DatSetTot=1;

startpixel=1;
endpixel= 2;

range=endpixel-startpixel;
jumpstep=(range-mod(range,DatSetTot))/DatSetTot;

dataYSEP=70




%% Main Data Analysis

figure( 'Name', 'Representative Spectra' );hold on

A = H(1:901);     
%for n=startpixel:endpixel


[fitresult,gof]=FitLorentz(XX, A, aproxX);         %lorentz fit
 

[fitresultII, gofy]= createSpline(XX, A);         %Raw Data SplineFit

   

% %Proof Figure ERASE LATER
% plot(XX,B+dataYSEP*b)
% A=B;


backAss=933;
windowStart=1160;
windowEnd=1709;

 
  for d=windowStart:windowEnd                         %Working Arrays loading
    e=1710-d;                                       
AA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                            %Subtract Lorentzian from Raw store in CC
    end

   
background=fitresult.back;

    for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       
XXX(e)=d;
BB(e)=fitresultII(d);              %Remove Background from Raw Spline and Store result in BB 
CC(e)=(BB(e)-AA(e));                              %Subtract Lorentzian from Raw store in CC
    end

   for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       
XXX(e)=d;
BB(e)=fitresultII(d);              %Remove Background from Raw Spline and Store result in BB 
CC(e)=(BB(e)-AA(e))+fitresult.back;                              %Subtract Lorentzian from Raw store in CC
    end
    
    
[fitresultNDC, gofII] = createFourFit(CC);        %Fit Fourier Series to CC

figure
hold on
scatter(XXX,CC,4)

    for z=1:550
DD(z)=fitresultNDC(z);                            %Store NDC fit in DD


    end
plot(XXX,DD);
hold off;
    
BBB=(BB-DD);
    
% length asignation (Z)
%Z(b)=b*jumpstep*pixsize;

[fitresult,gof]=FitLorentz(XXX, BBB, aproxX); 

for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

AAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end


for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

CCC(e)=(BB(e)-AAA(e));                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end





figure;
hold on
scatter(XXX,CCC,4);

[fitresultNDC, gofII] = createFourFit(CCC); 

 for z=1:550
DDD(z)=fitresultNDC(z);                            %Store NDC fit in DD


 end
    
 plot(XXX,DDD);
 hold off;

BBBB=(BB-DDD);
[fitresult,gof]=FitLorentz(XXX, BBBB, aproxX); 

for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

AAAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end

for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

CCCC(e)=(BB(e)-AAAA(e));                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end

figure;
hold on;
scatter(XXX,CCCC,4);

[fitresultNDC, gofII] = createFourFit(CCCC); 

for z=1:550
DDDD(z)=fitresultNDC(z);                            %Store NDC fit in DD


end

plot(XXX,DDDD); 
hold off;

for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

AAshow(e)=AAAA(e)+min(BB);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end

EEEE=(AAAA+DDDD)-BB;

figure;
plot(EEEE);

[fitresultRipple, gofIV] = createSinFit(XXX,EEEE); 


for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

Ripple(e)=fitresultRipple(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end




for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

Ripple(e)=fitresultRipple(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end

aproxhump=1340;

[fitresultHump,gof]=FitLorentz(XXX, EEEE, aproxhump);

for d=1160:1709                               %Working Arrays loading
    e=1710-d;                                       

Hump(e)=fitresultHump(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end

figure
hold on
plot(XXX,AAAA+DDDD+Hump+Ripple);
plot(XXX,BB);
hold off;




%% ERROR asignation
error=(gof.adjrsquare);
error2=(gofII.adjrsquare);
%error3=.3*sqrt((error(b)/error2(b))/2);


%% Plot All Curves
figure
hold on;
plot(XXX,AAshow);
plot(XXX,BB);
%plot(CC);
plot(XXX,DDDD);
hold off;
figure
hold on
scatter(XXX,BB,4);
plot(XXX,AAAA+DDDD,'LineWidth',1.5);
hold off;




%Numerical Integration of fits
Phase=trapz(AA)/(trapz(AAAA)+trapz(DDDD))
% PeakPosition(b)=fitresult.PP
% 
% Width(b)=fitresult.w
% Area(b)= trapz(AA);
Width=fitresult.w

%Plot Result

% figure( 'Name', 'Diamond Phase' );
% errorbar(Z,Phase,error3);
% figure( 'Name', 'Diamond Phase as function o Z in microns' );
% plot(Z,Phase);
% figure( 'Name', 'Peak Position v Position' );
% plot(Z,PeakPosition);
% figure( 'Name', 'Diamond Line Width Vs Position' );
% plot(Z,Width);


%% Data Export

filename='U:\HoltzGroup\Raju\StressMeasuremntsRetest\4_5_8hr\TextFiles\4_5_8hr_40sSpectra_86_1Analysis'
save(filename,'AAAA','DDDD','BB','Width','XXX');






%% FUNCTIONS
function [fitresult, gof] = FitLorentz(XX, A, aproxX)
%CREATEFIT(XX,A)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : XX
%      Y Output: A
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 04-Feb-2019 11:25:46


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( XX, A );

% Set up fittype and options.
% ft = fittype( '(a/(p*((x-PP)^2 + 1))+b)+background', 'independent', 'x','dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [aproxX 0.99 0.78 360 4];

% ft = fittype( '(a/(p*((x-PP)^2 + 1)))+back', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [aproxX 0.978732191344639 10 4];
ft = fittype( 'back+(2*a/3.1416)*(w/(4*(x-PP)^2+w^2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'OFF';
opts.Lower = [-Inf 0 0 2];
opts.StartPoint = [aproxX 20 350 8];
opts.Upper = [Inf Inf Inf Inf];
%            [pp   A  back  width]




% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.

figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'A vs. XX', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel XX
ylabel A
grid on


end



function [fitresult, gof] = createSpline(XX, A)
%CREATEFIT(XX,A)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : XX
%      Y Output: A
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 12-Mar-2019 12:24:46


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( XX, A );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
opts.SmoothingParam = 0.99999999;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'A vs. XX', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel XX
% ylabel A
% grid on
end

function [fitresult, gof] = createFourFit(CC)
%CREATEFIT(CC)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      Y Output: CC
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 12-Mar-2019 13:06:15


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( [], CC );

% Set up fittype and options.
ft = fittype( 'fourier7' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.001];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'CC', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% ylabel CC
% grid on


end

function [fitresult, gof] = createSinFit(XXX, EEEE)
%CREATEFIT(XXX,EEEE)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : XXX
%      Y Output: EEEE
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 01-Apr-2019 10:14:06


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( XXX, EEEE );

% Set up fittype and options.
ft = fittype( 'sin2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf];
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
opts.StartPoint = [11.4710762812716 12.7312804501037 2.25720593268414 7.36512926531848 10.9125261000889 1.13539654538716];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'EEEE vs. XXX', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel XXX
ylabel EEEE
grid on

end
