clear all;

fid = fopen('U:\HoltzGroup\CB\Caliberation\Sample_first\sammplefirst_across_secondregion_fullchip_84_1.txt', 'rt');



data = textscan(fid,'%s %f %f %f','headerlines', 1, 'delimiter', ',');
fclose(fid);
chr=data{1}{1};
chr=erase(chr,'"Wavelength"');
C = textscan(chr,'%f');
M=[];
H=M;



for n= 1:399
    
    chrr=data{1}{3+n};
    chrr=erase(chrr,'"1 "');
    D = textscan(chrr,'%f'); 
    D{1}(1)=[];

    A=cat(2,C{1},D{1});

    sensitivity=0.041;
% 
%    [pks,locs]=findpeaks(A(:,2),'MinPeakHeight',((1+sensitivity)*mean(A(:,2))));
%     findpeaks(A(:,2),'MinPeakHeight',((1+sensitivity)*mean(A(:,2))),'Annotate','extents')
%     

     M=cat(2,M,A(:,2));
   
end
H=hampel(M,5);    
figure('name',"Raw Data")
waterfall(M);

figure('name',"Hempel Outlier Filter 3 units")
waterfall(H);
filename='U:\HoltzGroup\CB\Caliberation\Sample_first\sammplefirst_across_secondregion_fullchip_84_1Vars'
save(filename,'H','M');
