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

AproxX=567;
correction=0;
%% Physical Constants
% totpixels = 400
% CCDDistance = 40
% Pixsize=CCDDistance/totpixels


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
%%Binning setup
% startpixel=100
% endpixel=110
% numbin=1
% DatSetTot=(endpixel-startpixel)/numbin
% range=endpixel-startpixel
% jumpstep=(range-mod(range, DatSetTot))/DatSetTot

%%  Main Analysis
max=size(AproxX);   % stablishes the total number of peaks to find and fit 

figure;
hold on;
% for b=1:DatSetTot
% for num=1:numbin
%     B=H(1:1340,startpixel);
%     n=startpixel+num+(b*jumpstep)
%     A=H(1:1340,n);
%     B =A+B;
theResult = sepblockfun(H,[1,3],'sum'); %max,min,sum,prod,mean
for n=1:133        %First Loop runs through all  the pixel strips
    A=theResult(:,n);
  Z(n)=n*.1;
for peakNum=1:max(2)        %This Loop runs through all the peaks per strip and fits them 
 
 
    [fitresult, gof] = LorentzFit(XX, A, AproxX ); 
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

[fitresult,gof]=LorentzFit(XXX, BBB, AproxX);

 
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
[fitresult,gof]=LorentzFit(XXX, BBBB, AproxX); 

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
    
    PeakPositions(AproxX,n)=fitresult.PP;
    PPerr(AproxX,n)=(errors(2,1)-errors(1,1))/2;

    Widths(AproxX,n)=fitresult.w;
    Werr(AproxX,n)=(errors(2,4)-errors(1,4))/2;
    
    Intensities(AproxX,n)=fitresult.a;
    Ierr(AproxX,n)=(errors(2,2)-errors(1,2))/2;
    
       for i = 500:1400

        L(i)=(fitresult(i)-fitresult.back);
        
        end
    K{AproxX}=L;

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

%for AproxX=1:max(2)
    
   figure('name',"Peak position as function of position for Peak at "+AproxX+" wavenumber")
   errorbar(Z,PeakPositions(AproxX,:),PPerr(AproxX,:));
   
   figure('name',"Linewidth as function of position for Peak at "+AproxX+" wavenumber")
   errorbar(Z,Widths(AproxX,:),Werr(AproxX,:));
   
   figure('name',"Intensity as function of position for Peak at "+AproxX+" wavenumber")
   errorbar(Z,Intensities(AproxX,:),Ierr(AproxX,:));
%end

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