clear all;
Peaks=[];
%AA=[];
K=[];
FF=[];
Curve1=[];

%%   Load Files 


path='C:\Users\gauta\Desktop\Friday_Small_1180302B\5by5PerfectlyCoalescedFriday_1180302B431' %file path minus extension
LoadExtension='Vars';
SaveExtension='Analysis';

filename1=strcat(path,LoadExtension)
load(filename1,'H','M','MaxWave','CCDStep');

filename2=strcat(path,SaveExtension)

%% Expected Peak Locations

AproxX=[567];
correction=0;
%% window setup
winstart=500
winend=600

%%  CCD Setup %%


 
figure
waterfall(H);

step=CCDStep;       %CCD datapoint separation
t=MaxWave;          %Max Wavenumber

   
for n= 1:1340
  
XX(n)=t;  %X axis values
t=t-(step+correction);

end
%%  Main Analysis
max=size(AproxX);   % stablishes the total number of peaks to find and fit 

figure;
hold on;
for n=1:399         %First Loop runs through all  the pixel strips
    A=H(:,n);
  Z(n)=n*.1;
for peakNum=1:max(2)        %This Loop runs through all the peaks per strip and fits them 
 
 
    [fitresult, gof] = LorentzFit(XX, A, (AproxX(peakNum))); 
    [fitresultII, gofy]= createSpline(XX, A);         %Raw Data SplineFit
    

 for d=winstart:winend                               %Working Arrays loading
    e=winend+1-d;                                     
    AA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                            %Subtract Lorentzian from Raw store in CC
    end

   
background=fitresult.back;

for d=winstart:winend                               %Working Arrays loading
    e=winend+1-d;                                      
    XXX(e)=d;
    BB(e)=fitresultII(d);              %Remove Background from Raw Spline and Store result in BB 
    CC(e)=(BB(e)-AA(e));                              %Subtract Lorentzian from Raw store in CC
end

   [fitresultNDC, gofII] = createFourFit(CC);   
    
for z=1:(winend-winstart+1)
   DD(z)=fitresultNDC(z);                            %Store NDC fit in DD


end
   
% length asignation (Z)
%Z(b)=(b*jumpstep*pixsize)-1.5;

%% Details

BBB=(BB-DD);
    
% length asignation (Z)
%Z(b)=b*jumpstep*pixsize;

[fitresult,gof]=LorentzFit(XXX, BBB, AproxX(peakNum));

 
for d=winstart:winend                               %Working Arrays loading
    e=winend+1-d;                                     

    AAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end


for d=winstart:winend                               %Working Arrays loading
    e=winend+1-d;                                      

    CCC(e)=(BB(e)-AAA(e));                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end





% figure;
% hold on
% scatter(XXX,CCC,4);     

%% 
[fitresultNDC, gofII] = createFourFit(CCC); 

 for z=1:(winend-winstart+1)
     DDD(z)=fitresultNDC(z);                            %Store NDC fit in DD


 end
    
%  plot(XXX,DDD);
%  hold off;

BBBB=(BB-DDD);
[fitresult,gof]=LorentzFit(XXX, BBBB, AproxX(peakNum)); 

for d=winstart:winend                               %Working Arrays loading
    e=winend+1-d;                                       

    AAAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
end

%for d=1200:1799                               %Working Arrays loading
 %   e=1800-d;                                        

%CCCC(e)=(BB(e)-AAAA(e));                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
%end

% figure;
% hold on;
% scatter(XXX,CCCC,4);

%[fitresultNDC, gofII] = createFourFit(CCCC); 

%for z=1:600
%DDDD(z)=fitresultNDC(z);                            %Store NDC fit in DD

%end
%% Data Load onto export variables

    errors=confint(fitresult);
    %% probable modification for data selection
    % If Intensities(peaknum,n)< Threshold
    % PeakPositions(peakNum,n)=terminationValue;   
    % else
    %PeakPositions(peakNum,n)=fitresult.PP;
    %end
    
    PeakPositions(peakNum,n)=fitresult.PP;
    PPerr(peakNum,n)=(errors(2,1)-errors(1,1))/2;

    Widths(peakNum,n)=fitresult.w;
    Werr(peakNum,n)=(errors(2,4)-errors(1,4))/2;
    
    Intensities(peakNum,n)=fitresult.a;
    Ierr(peakNum,n)=(errors(2,2)-errors(1,2))/2;
    
       for i = 500:600

        L(i)=(fitresult(i)-fitresult.back);
        
        end
    K{peakNum}=L;

    end
    
    KK=K{1};
    if max(2)>1
      
        for n=2:max(2)
           KK=KK+K{n};
        end
    end
    
FF=cat(1,FF,KK);

Curve1=cat(1,Curve1,AAAA);

plot(KK);

end
%figure;
%waterfall(FF);
%figure;
%waterfall(Curve1);

for peakNum=1:max(2)
    
   figure('name',"Peak position as function of position for Peak at "+AproxX(peakNum)+" wavenumber")
   errorbar(Z,PeakPositions(peakNum,:),PPerr(peakNum,:));
   
   figure('name',"Linewidth as function of position for Peak at "+AproxX(peakNum)+" wavenumber")
   errorbar(Z,Widths(peakNum,:),Werr(peakNum,:));
   
   figure('name',"Intensity as function of position for Peak at "+AproxX(peakNum)+" wavenumber")
   errorbar(Z,Intensities(peakNum,:),Ierr(peakNum,:));
end

filename2
save(filename2,'AAAA','XX','PeakPositions','PPerr','Widths','Werr','Intensities','Ierr')

%% FUNCTIONS
function [fitresult, gof] = LorentzFit(XX, A, aproxX)
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
% opts.StartPoint = [aproxX 0.978732191344639 370 4];
ft = fittype( 'back+(2*a/3.1416)*(w/(4*(x-PP)^2+w^2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
pplow=aproxX-5;
pphigh=aproxX+5;
opts.Lower = [pplow 5 100 1];
opts.StartPoint = [aproxX 10 206 8];
opts.Upper = [pphigh Inf Inf 25];
%            [pp   A  back  width]




% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

%Plot fit with data.

% figure( 'Name', 'lorentz fit' );
% h = plot( fitresult, xData, yData );
% legend( h, 'A vs. XX', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel XX
% ylabel A
% grid on


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
opts.SmoothingParam = 0.9999999;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

%Plot fit with data.
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
ft = fittype( 'fourier6' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0 0 0.012];

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