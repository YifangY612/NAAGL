 
clc
clear all, close all
warning off;
addpath funs
addpath measure
%addpath clusteringMeasure
addpath(genpath('./'));


%% dataset
%ds = {'Dermatology','BBCsport','MSRC_4view','flower17','Caltech101-7','Caltech101-20','BDGP','AwA','YTF50'};        % 2386
ds = {'flower17'}; 

dsPath = 'C:\Users\yyf\Desktop\NAAGL\NAAGL-main\datasets\';

acc=[];nmi=[];purity=[];fsore=[];


for di =1:length(ds)
    dataName = ds{di}; disp(dataName);
    load(strcat(dsPath,dataName));
    k = length(unique(Y));
    n = length(Y); 

    anchor = [k 2*k 3*k];
    alpha = [1 10 10^2 10^3 10^6];

 


    %%
    
    allresult = [];
    maotu={};
    %fid = fopen('result1.txt','w');
    for ichor = 1:length(anchor)
        for id = 1:length(alpha)
            tic;
            [U,V,A,W,Z,E,iter,obj,y1] = naagl(X,Y,anchor(ichor),alpha(id));
             res = Clustering8Measure(Y,y1);
             time(ichor,id)  = toc;
            fprintf('Dataset:%-10s\t anchor:%-10d\t alpha:%-10d\t %.6f\t %.6f\t  %.6f\t  %.6f\t Time:%.6f \n',dataName, anchor(ichor), alpha(id), res(1), res(2), res(3), res(4), time(ichor,id));
            fid = fopen('result.txt','a');
            fprintf(fid,'Dataset:%-10s\t anchor:%-10d\t alpha:%-10d\t %.6f\t  %.6f\t  %.6f\t  %.6f\t Time:%.6f\n', dataName, anchor(ichor), alpha(id), res(1), res(2), res(3), res(4), time(ichor,id));
            fclose(fid);
     %fprintf(fid,'\t\n');
            %allresult = [allresult;res time(ichor,id)];
            %maotu{end+1}=Z;
            %clear A W Z E iter obj y1  
           
        end     
    end

   
    
end





