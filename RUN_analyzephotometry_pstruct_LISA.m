clear all
close all
fclose all
trate=10; %% targeting temporal resolution; [trate] data point per second
flist={
    %%% input directory of txt file here
    %%% txt file contains information of what dat files to analyze and what
    %%% experimental parameters to use
    'DR_HFD_IG_lipid_baseline_HFD.txt'
    };

for x=1:length(flist)
    [expnames,mcode_raw]=get_info_from_txtfile(flist{x}); %% extract the strucuture data from txt file, including experiment name
    %% convert mcode into array; mcode indicate how trials are organized (e.g technical repeats for the same biological repeat)
    mc1=length(mcode_raw);
    mcode_raw=str2num(mcode_raw);
    mcode=[];
    for i=1:mc1
        a=mcode_raw-floor(mcode_raw/10)*10;
        mcode=[a mcode];
        mcode_raw=floor(mcode_raw/10);
    end
    clear mc1 mcode_raw
    for i=1:length(expnames)
        expname=cell2mat(expnames(i)); % get file name
        [ mouse, mftrace,mpre,mpost,psthnomavg,mpsthnom,mpsthtime, mlowpsthnom, mtau1,mtau2,all_low_prism] = ...
            analyzephotometry_dat2mat( strcat(expname),mcode,trate); % a function reads through all the raw photometry files, preprocess them and then organize them into matrix for convenience of any further statistical analysis
        temp=flist;
        clear files
        save(fullfile(pwd,strcat(expname,'_wspace.mat')));
        flist=temp;
    end
end

function [ expname,mcode ] = get_info_from_txtfile( infile)
fid=fopen(infile);
p=fscanf(fid,'%s',[1]);
n=0;
filename=p;
fid2=fopen('t.txt','w');
expname={};
while p   
    p=fscanf(fid,'%s',[1]);
    if strcmp(p,'TREATMENT');
        n=n+1;
        p=fscanf(fid,'%s',[1]);
        outfile=strcat(filename,num2str(n),p);
        expname{n}=outfile;
        mcode=fscanf(fid,'%s',[1]);
        fid2 = fopen(strcat(filename,num2str(n),p), 'w');
    else if strcmp(p,'RETURN')
            fprintf(fid2, '\n');
        else
            fprintf(fid2, '%s\t', p);
        end
    end
end
close all;
end

