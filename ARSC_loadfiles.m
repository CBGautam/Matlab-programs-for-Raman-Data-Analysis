clear all;

path='C:\Users\gauta\Desktop\12062019\Diamond_spectra_2_GOOD' %file path minus extension
LoadExtension='.txt';
SaveExtension='Vars';

filename1=strcat(path,LoadExtension)
filename2=strcat(path,SaveExtension)
A=importdata('Diamond_spectra_2_GOOD.txt');
Data=A';
Data(:,1)=[];
B=Data;
W=Data(:,1);
Data(:,1)=[];
M=Data;
% 
% fid = fopen(filename1, 'rt');
% data = textscan(fid,'%s %f %f %f','headerlines', 1, 'delimiter', ',');
% fclose(fid);
% chr=data{1}{1};
% chr=erase(chr,'"Wavelength"');
% C = textscan(chr,'%f');
% M=[];
% H=M;
% 
% 
% 
% for n= 1:399
%     
%     chrr=data{1}{3+n};
%     %chrr=erase(chrr,'"1 "');
%     D = textscan(chrr,'%f'); 
%     D{1}(1)=[];
% 
%     A=cat(2,C{1},D{1});
% 
%     sensitivity=0.041;
% 
% 
%      M=cat(2,M,A(:,2));
%    
% end

%% Wavenumber Capture
% W=str2num(chr);
 MaxWave=max(W)
 CCDStep=(W(1099)-W(1))/1099
% 
% 
% 
% %%
% 
H=hampel(M,5);    
figure('name',"Raw Data")
waterfall(M);

figure('name',"Hempel Outlier Filter 5 units")
waterfall(H);

save(filename2,'H','M','MaxWave','CCDStep');
