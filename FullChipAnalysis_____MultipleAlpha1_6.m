clear all;


%%  Array Declaration

Curve10=[];
Curve20=[];
Curve30=[];
Error30=[];
RawData0=[];

Curve190=[];
Curve290=[];
Curve390=[];
Error390=[];
RawData90=[];


PP0=[];
PP0Err=[];
I0=[];
I0Err=[];

PP90=[];
PP90Err=[];
I90=[];
I90Err=[];




%% Datafiles
path='C:\Users\j_r847\Desktop\New folder (4)\fullchip_across_section_fullchip_112' %file path minus extension
LoadExtension='Vars';
SaveExtension='Analysis';

filename1=strcat(path,LoadExtension)
load(filename1,'H','M','MaxWave','CCDStep');

filename2=strcat(path,SaveExtension)


correction =70;
stepSize=CCDStep;


%% Multi peak setup  %%

PeakLocations=[1330, 525]



%% Physical Constants

totpixels=400; %total amount of strips
CCDDistance=35; %image size in Microns;
pixsize= CCDDistance/totpixels;

aproxX=PeakLocations(1);  %approximate peak location


%%  CCD Setup


 
figure
waterfall(H);

step=CCDStep;       %CCD datapoint separation
t=MaxWave;          %Max Wavenumber

   
for n= 1:1340
  
XX(n)=t+correction;  %X axis values
t=t-step;

end

astep=0;

%% Binning  Setup



startpixel=150; %starting Strip
endpixel= 350;


numbin=1;  %Number of Spectra per Bin
DatSetTot=200;% Total Number of Bins

range=endpixel-startpixel;
jumpstep=(range-mod(range,DatSetTot))/DatSetTot;

dataYSEP=60;


%% Data Window Setup

WinStart=500%cast((MaxWave-(CCDStep*1340)),int8);
WinEnd=1400%cast(MaxWave,int8);
showpix = 351;


for showpix=startpixel:endpixel

            %% EVEN Analysis
stringg="Even"
showpix 

%             figure( 'Name', '0 Degree Analyzer Laurentzian Fits' );hold on
for b=1:DatSetTot;
B = H(1:1340,startpixel+(b*jumpstep));   
   for num=1:numbin    
   

n=startpixel+num+(b*jumpstep)

   A = H(1:1340,n);                               %Spectra Harvesting from file array
B=A+B;
end

[fitresult,gof]=createFit(XX, B, aproxX);         %lorentz fit
[fitresultII, gofy]= createSpline(XX, B);  

         for d=WinStart:(WinEnd-1)                               %Working Arrays loading
            e=WinEnd-d;                                     
        AA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                                    %Subtract Lorentzian from Raw store in CC
            end


        background=fitresult.back;

            for d=WinStart:(WinEnd-1)                               %Working Arrays loading
            e=WinEnd-d;                                       
        XXX(e)=d;
        BB(e)=fitresultII(d);              %Remove Background from Raw Spline and Store result in BB 
        CC(e)=(BB(e)-AA(e));                              %Subtract Lorentzian from Raw store in CC
            end

         [fitresultNDC, gofII] = createFourFit(CC);   

            for z=1:(WinEnd-WinStart)
                DD(z)=fitresultNDC(z);                            %Store NDC fit in DD


            end

        % length asignation (Z)
        Z(b)=(b*jumpstep*pixsize)-1.5;



        %% Details

        BBB=(BB-DD);

        % length asignation (Z)
        %Z(b)=b*jumpstep*pixsize;

        [fitresult,gof]=createFit(XXX, BBB, aproxX);


        for d=WinStart:(WinEnd-1)                               %Working Arrays loading
            e=WinEnd-d;                                       

        AAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                                 %Subtract Lorentzian from Raw store in CC
        end


        for d=WinStart:(WinEnd-1)                               %Working Arrays loading
            e=WinEnd-d;                                  

        CCC(e)=(BB(e)-AAA(e));                %Remove Background from Lorentzian and store in AA
                                 %Subtract Lorentzian from Raw store in CC
        end





        %% 
        [fitresultNDC, gofII] = createFourFit(CCC); 

         for z=1:(WinEnd-WinStart)
        DDD(z)=fitresultNDC(z);                            %Store NDC fit in DD


         end

        %  plot(XXX,DDD);
        %  hold off;

        BBBB=(BB-DDD);
        [fitresult,gof]=createFit(XXX, BBBB, aproxX); 

        for d=WinStart:(WinEnd-1)                               %Working Arrays loading
            e=WinEnd-d;                                                    

        AAAA(e)=fitresult(d);                %Remove Background from Lorentzian and store in AA
                                 %Subtract Lorentzian from Raw store in CC
        end

        for d=WinStart:(WinEnd-1)                               %Working Arrays loading
            e=WinEnd-d;                                                    

        CCCC(e)=(BB(e)-AAAA(e));                %Remove Background from Lorentzian and store in AA
                                 %Subtract Lorentzian from Raw store in CC
        end

 

        [fitresultNDC, gofII] = createFourFit(CCCC); 

        for z=1:(WinEnd-WinStart)
        DDDD(z)=fitresultNDC(z);                            %Store NDC fit in DD


        end



        for d=WinStart:(WinEnd-1)                               %Working Arrays loading
            e=WinEnd-d;                                           

        AAshow(e)=AAAA(e)+min(BB); %Remove Background from Lorentzian and store in AA
        Dianobak(e)=AAAA(e)+(0-min(AAAA));
        NDCnobak(e)=DDDD(e)+(0-min(DDDD));  %Subtract Lorentzian from Raw store in CC

        PercentErr(e)=(100*abs((BB(e)-(AAAA(e)+DDDD(e)))/((AAAA(e)+DDDD(e)))));
        end

        EEEE=(AAAA+DDDD)-BB;



        %%Summary and  ERROR asignation
        error(b)=(gof.adjrsquare);
        error2(b)=(gofII.adjrsquare);
        error3(b)=.3*sqrt((error(b)/error2(b))/2);
        TotpercentE0(b)=mean(PercentErr);
        sigma=confint(fitresult);
        PPuncertainty0(b)=fitresult.PP-sigma(1,1);
        Wuncertainty0(b)=fitresult.w-sigma(1,4);
        IUncertainty0(b)=fitresult.a-sigma(1,2);

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

        Curve10=cat(1,Curve10,DDDD);
        Curve20=cat(1,Curve20,AAAA+DDDD);
        Curve30=cat(1,Curve30,AAAA);
        Error30=cat(1,Error30,PercentErr);
        RawData0=cat(1,RawData0,BB);

        waterfall(Curve20);

        %% Numerical Integration of fits

        Phase0(b)=trapz(AAAA)/(trapz(AAAA)+trapz(DDDD));
        PeakPosition0(b)=fitresult.PP;
        Intensity0(b)=fitresult.a;
        Width0(b)=fitresult.w;
        AreaDiamond0(b)= trapz(Dianobak);
        AreaNDC0(b)=trapz(NDCnobak);

end



      


        %% RESULTS


        %Dual Monitor

%         figure( 'Name', "Intensity vs Polarization Angle Pixel#"+showpix+" ",'position',[-1880,930,600,400] );
%         hold on;
%         errorbar(Angle,Intensity0,IUncertainty0,'MarkerEdgeColor','red');
%         errorbar(Angle,Intensity90,IUncertainty90,'MarkerEdgeColor','blue');
%         hold off;
% 
%         figure( 'Name', "Delta Omega Peak Position v Polarization Angle"+showpix+" ",'position',[-640,930,600,400] );
%         hold on;
%         errorbar(Angle,PeakPosition0,PPuncertainty0,'-s','MarkerSize',4,'MarkerEdgeColor','black','MarkerFaceColor','red')
%         errorbar(Angle,PeakPosition90,PPuncertainty90,'-s','MarkerSize',4,'MarkerEdgeColor','black','MarkerFaceColor','blue')
%         hold off;

        
        % Single Monitor
        
        figure( 'Name', "Intensity vs Polarization Angle Pixel#"+showpix+" ",'position',[1880,930,600,400] );
        hold on;
        errorbar(Angle,Intensity0,IUncertainty0,'MarkerEdgeColor','red');
        errorbar(Angle,Intensity90,IUncertainty90,'MarkerEdgeColor','blue');
        hold off;

        figure( 'Name', "Delta Omega Peak Position v Polarization Angle"+showpix+" ",'position',[640,930,600,400] );
        hold on;
        errorbar(Angle,PeakPosition0,PPuncertainty0,'-s','MarkerSize',4,'MarkerEdgeColor','black','MarkerFaceColor','red')
        errorbar(Angle,PeakPosition90,PPuncertainty90,'-s','MarkerSize',4,'MarkerEdgeColor','black','MarkerFaceColor','blue')
        hold off;

        
        
        
        
        
        
        
        
%%  Loading Export Array
pixindex=(showpix-startpixel+1)
%Even
PP0=cat(1,PP0,PeakPosition0);
PP0Err=cat(1,PP0Err,PPuncertainty0);
I0=cat(1,I0,Intensity0);
I0err=cat(1,I0Err,IUncertainty0);



end

%%  Final Plots

%Dual Monitor

% figure( 'Name', 'Peak Position v Polarization Angle 0','position',[-1260,930,600,400] );
% waterfall(PP0);
% figure( 'Name', 'Peak Position v Polarization Angle 90','position',[-1260,400,600,400] );
% waterfall(PP90);
% 
% figure( 'Name', 'Total Intensity vs Angle for 0 Degree','position',[-1880,400,600,400] );
% waterfall(I0);
% figure( 'Name', 'Total Intensity vs Anlge for 90 Degree','position',[-640,400,600,400] );
% waterfall(I90);

%Single Monitor

figure( 'Name', 'Peak Position v Polarization Angle 0','position',[1260,930,600,400] );
waterfall(PP0);
figure( 'Name', 'Peak Position v Polarization Angle 90','position',[1260,400,600,400] );
waterfall(PP90);

figure( 'Name', 'Total Intensity vs Angle for 0 Degree','position',[1880,400,600,400] );
waterfall(I0);
figure( 'Name', 'Total Intensity vs Anlge for 90 Degree','position',[640,400,600,400] );
waterfall(I90);




%% Data Save

filename=filename2
% save(filename,'Angle','PP0','PP0Err','I0','I0Err','PP90','PP90Err','I90','I90Err','Z');




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
pplow=aproxX-15;
pphigh=aproxX+15;
opts.Lower = [pplow 3 100 2];
opts.StartPoint = [aproxX 15 360 5];
opts.Upper = [pphigh Inf Inf 60];
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