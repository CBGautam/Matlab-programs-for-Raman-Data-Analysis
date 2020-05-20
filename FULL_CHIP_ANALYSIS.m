clear all;
Peaks=[];
AAAA=[];
K=[];
%% Datafiles

path='C:\Users\gauta\Desktop\11172019\5by5\11172019_ROI_5by5_600sec_14'%file path minus extension
LoadExtension='Vars';
SaveExtension='Analysis';
filename1=strcat(path,LoadExtension)
load(filename1,'H','M','MaxWave','CCDStep');
filename2=strcat(path,SaveExtension)

%% Physical Constants

% totpixels=400; %total amount of strips
% CCDDistance=35; %image size in Microns;
% pixsize= CCDDistance/totpixels;
aproxX=521;  %approximate peak location

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
startpixel=10; %starting Strip
endpixel= 390;
numbin=10;  %Number of Spectra per Bin
DatSetTot=(endpixel-startpixel)/numbin; % Total Number of Bins

range=endpixel-startpixel;
jumpstep=(range-mod(range,DatSetTot))/DatSetTot;
% HH=H(:,startpixel:endpixel);
% theResult=sepblockfun(HH, [1,numbin],'sum'); %max, min, sum, prod, mean
%dataYSEP=60
% HH=H(:,startpixel:endpixel);
% theResult=sepblockfun(HH, [1,numbin],'sum'); %max, min, sum, prod, mean

%% Main Data Analysis

figure( 'Name', 'Representative Spectra' );hold on
for b=1:DatSetTot
%  A=theResult(:,b);
  Z(b)=0.05*b;
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
% error(b)=(gof.adjrsquare);
% error2(b)=(gofII.adjrsquare);
% error3(b)=.3*sqrt((error(b)/error2(b))/2);
% TotpercentE(b)=mean(PercentErr);
% sigma=confint(fitresult);
% PPuncertainty(b)=fitresult.PP-sigma(1,1);
% Wuncertainty(b)=fitresult.w-sigma(1,4);
    errors=confint(fitresult);
    PeakPositions(b)=fitresult.PP;
    PPerr(b)=(errors(2,1)-errors(1,1))/3.92;

    Widths(b)=fitresult.w;
    Werr(b)=(errors(2,4)-errors(1,4))/3.92;
    
    Intensities(b)=fitresult.a;
    Ierr(b)=(errors(2,2)-errors(1,2))/3.92;
     for i = 493:1388

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
% figure
% waterfall(Curve1);
figure
waterfall(Curve3);
% figure( 'Name', 'Raw' );
% %waterfall(RawData);
% figure( 'Name', 'Total Fit Error' );
% waterfall(Error3);

% figure( 'Name', 'Diamond Phase' );
% errorbar(Z,Phase,error3);
% figure( 'Name', 'Diamond Phase as function o Z in microns' );
% plot(Z,Phase);
figure( 'Name', 'Peak Position v Position' );
plot(Z,PeakPositions);
figure( 'Name', 'Diamond Line Width Vs Position' );
plot(Z,Widths);

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
pplow=aproxX-5;
pphigh=aproxX+5;
%opts.Lower = [-Inf 0.1 0 1];
opts.Lower=[pplow 0.2 100 0.5];
%opts.StartPoint = [aproxX 4 200 4];
opts.StartPoint = [aproxX 0.97 150 0.8]; %0.97
opts.Upper = [pphigh Inf Inf 6];
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
function X=sepblockfun(X,blockdims,fun)
%Perform a separable operation on sub-blocks of an input array. Here, a function
%op() is said to be separable if for any array B with elements B(i,j,k,...), 
%the operation op(B(:)), resulting in a scalar, can be equivalently done by 
%applying op() first along i, then along j, then along k, etc...
%
%USAGE:
%
%      Y=sepblockfun(X,blockdims,fun)
%
%in:
%
%
%   X: A full array. If the ndSparse class defintion is on the path, then X 
%      can also be a regular sparse matrix  or ndSparse array. Performance might 
%      not be as strong as for full arrays, however.
%
%   blockdims: a vector of integers  specifying the dimensions of the
%              sub-blocks. The array X must partition evenly into blocks 
%              of this size. If blockdims(i) is set to Inf then it will be 
%              replaced with blockdims(i)=size(X,i).
%
%   fun:  function handle to an operation assumed to be separable
%         (Examples: max,min,sum,prod,mean, etc...). The function must
%         accept the input syntax fun(B,DIM) where B is an input array
%         and DIM is a dimension along which to operate.  Alternatively,
%         fun can be one of the following strings 'max','min','sum','mean',
%         'prod'.
%
%
%out:
%
%   Y: the output array. Y(i)=fun(Xi(:),1) where Xi is the i-th sub-block of
%      the input array X.
%
%
%EXAMPLE 1: Divide a 400x400x400 array into 10x10x10 blocks. Return the blockwise, 
%mean, max, and min. of each block, each organized as a 40x40x40 array.
%
%   A=rand(400,400,400);
%   Ameans=sepblockfun(A,[10,10,10],@mean);
%   Amins=sepblockfun(A,[10,10,10],'min' ); 
%   Amaxs=sepblockfun(A,[10,10,10], @(B,d) max(B,[],d)  );
%
%EXAMPLE 2: Not all operations satisfy the separability property, but
%sometimes inseparable operations can be decomposed into separable ones. As
%an example, we take the blockwise standard deviations of the  same array
%from Example 1.
%
%   Astds=sqrt( sepblockfun(A.^2,[10,10,10],'mean') - Ameans.^2 );  
%
%
% Written by Matt Jacobson, 2014

if issparse(X)&& exist('ndSparse','class') && ~isa(X,'ndSparse')

     X=ndSparse(X);
    
end


if ischar(fun)

  switch fun

    case 'max'
     
         fun=@(b,d) max(b,[],d);

    case 'min'

         fun=@(b,d) min(b,[],d);

    case 'sum'

         fun=@sum;

    case 'mean'

         fun=@mean;

    case 'prod'

        fun=@prod;


    otherwise

     error 'Unrecognized fun() selection'

  end


end


nn=max(length(blockdims),ndims(X));
blockdims(end+1:nn)=1;

[sz{1:nn}]=size(X); %M is the original array
sz=[sz{:}];

idx=~isfinite(blockdims);
blockdims(idx)=sz(idx);

newdims=sz./blockdims;

args=num2cell([blockdims;newdims]);

X=reshape(X,args{:});

for ii=1:nn


 X=fun(X,2*ii-1);
 
end

X=reshape(X,newdims);
 end