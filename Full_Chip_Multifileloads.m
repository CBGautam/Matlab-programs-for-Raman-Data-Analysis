numfiles = 2;
mydata = cell(1, numfiles);

for k = 1:numfiles
  myfilename = sprintf('11032019_Friday2by6l_5datas1_%d.txt', k);
  mydata{k} = importdata(myfilename);
 M=mydata{1,1}.data+mydata{1,k}.data;
end
 %M=mydata{1,1}.data+mydata{1,2}.data;
 path='C:\Users\gauta\Desktop\DESKTOP\matlabhelprequest files\11032019_Friday2by6l_5datas1' ;
 SaveExtension='Vars';
filename2=strcat(path,SaveExtension)
 chr=mydata{1,1}.data(2,:);
 MaxWave=max(chr)
CCDStep=(chr(1)-chr(1340))/1340
 M([1,2,3],:)=[];
 M=M';
 H=M';
 H=hampel(M,5);    
 figure('name',"Raw Data")
 waterfall(M')
 figure('name',"Output filter")
 waterfall(H)
save(filename2,'H','M','MaxWave','CCDStep');