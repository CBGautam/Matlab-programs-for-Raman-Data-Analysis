clear all;
Peaks=[];
AAAA=[];
K=[];
%% Datafiles

path='C:\Users\c_g640\Desktop\October\oct4_newtrial_120sec823'%file path minus extension
LoadExtension='Vars';
SaveExtension='Analysis';
filename1=strcat(path,LoadExtension)
load(filename1,'H','M','MaxWave','CCDStep');
filename2=strcat(path,SaveExtension)

%% Physical Constants

% totpixels=400; %total amount of strips
% CCDDistance=35; %image size in Microns;
% pixsize= CCDDistance/totpixels;
aproxX=567;  %approximate peak location

%%  CCD Setup

figure

waterfall(H);

 Curve1=[];
% Curve2=[];
 Curve3=[];
% Error3=[];
% RawData=[];

%for strip=10:399
%strip = 320;

step=CCDStep; %CCD datapoint separation
t=MaxWave; %Max Wavenumber
%step=0.61;     %CCD datapoint separation
%t=1831.793;     %Max Wavenumber

for n= 1:1340
  
XX(n)=t;  %X axis values
t=t-step;

end
%% Window Setting
winstart=500
winend=600
%% Binning  Setup
startpixel=100; %starting Strip
endpixel= 280;
numbin=10;  %Number of Spectra per Bin
DatSetTot=(endpixel-startpixel)/numbin; % Total Number of Bins

range=endpixel-startpixel;
jumpstep=(range-mod(range,DatSetTot))/DatSetTot;

%dataYSEP=60
%HH=H(:,startpixel:endpixel);
%theResult=sepblockfun(HH, [1,numbin],'sum'); %max, min, sum, prod, mean

%% Main Data Analysis

figure( 'Name', 'Representative Spectra' );hold on
for b=1:DatSetTot
    for num=1:numbin    
B = H(1:1340,startpixel);     
%for n=startpixel:endpixel
n=startpixel+num+(b*jumpstep)

   A = H(1:1340,n);                               %Spectra Harvesting from file array
B=A+B;

[fitresult,gof]=createFit(XX, A, aproxX);         %lorentz fit
%[fitresultII, gofy]= createSpline(XX, A);         %Raw Data SplineFit

   % end

%Proof Figure ERASE LATER
% plot(XX,B+dataYSEP*b)
% A=B;


 %for d=1200:1799                               %Working Arrays loading
  %  e=1800-d;                                     
%AA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                            %Subtract Lorentzian from Raw store in CC
   % end

   
background=fitresult.back;

  %  for d=1200:1799                               %Working Arrays loading
   % e=1800-d;                                       
%XXX(e)=d;
%BB(e)=fitresultII(d);              %Remove Background from Raw Spline and Store result in BB 
%CC(e)=(BB(e)-AA(e));                              %Subtract Lorentzian from Raw store in CC
  %  end

 %[fitresultNDC, gofII] = createFourFit(CC);   
    
    %for z=1:600
%DD(z)=fitresultNDC(z);                            %Store NDC fit in DD


   % end
   
% length asignation (Z)
%Z(b)=(b*jumpstep*pixsize)-1.5;



%% Details

%BBB=(BB-DD);
    
% length asignation (Z)
%Z(b)=b*jumpstep*pixsize;

%[fitresult,gof]=createFit(XXX, BBB, aproxX);

 
%for d=1200:1799                               %Working Arrays loading
 %   e=1800-d;                                       

%AAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
%end


%for d=1200:1799                               %Working Arrays loading
 %   e=1800-d;                                       

%CCC(e)=(BB(e)-AAA(e));                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
%end





% figure;
% hold on
% scatter(XXX,CCC,4);     

%% 
%[fitresultNDC, gofII] = createFourFit(CCC); 

 %for z=1:600
%DDD(z)=fitresultNDC(z);                            %Store NDC fit in DD


 %end
    
%  plot(XXX,DDD);
%  hold off;

%BBBB=(BB-DDD);
%[fitresult,gof]=createFit(XXX, BBBB, aproxX); 

%for d=1200:1799                               %Working Arrays loading
  %  e=1800-d;                                        

%AAAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                         %Subtract Lorentzian from Raw store in CC
%end

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

%plot(XXX,DDDD); 
%hold off;

%for d=1200:1799                               %Working Arrays loading
 %   e=1800-d;                                         

%AAshow(e)=AAAA(e)+min(BB); %Remove Background from Lorentzian and store in AA
%Dianobak(e)=AAAA(e)+(0-min(AAAA));
%NDCnobak(e)=DDDD(e)+(0-min(DDDD));  %Subtract Lorentzian from Raw store in CC

%PercentErr(e)=(100*abs((BB(e)-(AAAA(e)+DDDD(e)))/((AAAA(e)+DDDD(e)))));
%end

%EEEE=(AAAA+DDDD)-BB;

% figure;
% plot(EEEE);



%% ERROR asignation
%error(b)=(gof.adjrsquare);
%error2(b)=(gofII.adjrsquare);
%error3(b)=.3*sqrt((error(b)/error2(b))/2);
%TotpercentE(b)=mean(PercentErr);
%sigma=confint(fitresult);
%PPuncertainty(b)=fitresult.PP-sigma(1,1);
%Wuncertainty(b)=fitresult.w-sigma(1,4);
errors=confint(fitresult);
    
    PeakPositions=fitresult.PP;
    PPerr=(errors(2,1)-errors(1,1))/3.92;

    Widths=fitresult.w;
    Werr=(errors(2,4)-errors(1,4))/3.92;
    
    Intensities=fitresult.a;
    Ierr=(errors(2,2)-errors(1,2))/3.92;
     for i = 495:1460

        L(i)=(fitresult(i)-fitresult.back);
        
        end
    K=L;

    end
    
    KK=K;
%     if max(2)>1
%       
%         for n=2:max(2)
%             KK=KK+K{n};
%         end
%     end
    
AAAA=cat(1,AAAA,KK);


plot(KK);

end
figure;
waterfall(AAAA);


%% Plot All Curves
% figure
% hold on;
% 
% plot(XXX,AA+b*dataYSEP);
% plot(XXX,BB+b*dataYSEP);
% plot(CC);
% plot(XXX,DD+b*dataYSEP);
% waterfall(Curve2);
% hold off;

% figure
% hold on;
% plot(XXX,AAshow);
% plot(XXX,BB);
% %plot(CC);
% plot(XXX,DDDD);
% hold off;
% figure( 'Name', 'Plot No'+"b" );
% hold on
% scatter(XXX,BB,4);
% plot(XXX,AAAA+DDDD,'LineWidth',1.5);
% hold off;

% Curve1=cat(1,Curve1,DDDD);
% Curve2=cat(1,Curve2,AAAA+DDDD);
Curve3=cat(1,Curve3,AAAA);
% Error3=cat(1,Error3,PercentErr);
% RawData=cat(1,RawData,BB);

% waterfall(Curve2);

%Numerical Integration of fits
% Phase(b)=trapz(AAAA)/(trapz(AAAA)+trapz(DDDD));
% PeakPosition(b)=fitresult.PP

% Width(b)=fitresult.w
% AreaDiamond(b)= trapz(Dianobak);
% AreaNDC(b)=trapz(NDCnobak);
%end

%Plot Result
figure
waterfall(Curve1);
figure
waterfall(Curve3);
figure( 'Name', 'Raw' );
waterfall(RawData);
figure( 'Name', 'Total Fit Error' );
waterfall(Error3);

% figure( 'Name', 'Diamond Phase' );
% errorbar(Z,Phase,error3);
% figure( 'Name', 'Diamond Phase as function o Z in microns' );
% plot(Z,Phase);
figure( 'Name', 'Peak Position v Position' );
plot(Z,PeakPosition);
figure( 'Name', 'Diamond Line Width Vs Position' );
plot(Z,Width);

% figure( 'Name', 'Mean Percent Error' );
% plot(TotpercentE);




%% Data Save
filename2
%filename='C:\Users\gauta\Desktop\CB\Caliberation\Analysis'
save(filename2,'AAAA','XX','PeakPositions','PPerr','Widths','Werr','Intensities','Ierr');




%% FUNCTIONS
function [fitresult, gof] = createFit(XX, A, aproxX)
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
opts.Lower = [-Inf 0.1 0 1];
opts.StartPoint = [aproxX 4 200 4];
opts.Upper = [Inf Inf Inf 40];
%            [pp   A  back  width]




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
opts.SmoothingParam = 0.999999999;

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