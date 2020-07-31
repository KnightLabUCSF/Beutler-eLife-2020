function [expname, mftrace,mpre,mpost,psthnomavg,mpsthnom,mpsthtime, mlowpsthnom,mtau1,mtau2,all_low_prism] = analyzephotometry_dat2mat( infile,mcode,downsampleto )
%ANALYZEPEPTIDE  read the plain text, goes to the folder, import all the .dat files and analyze them based on the stampled time.
%   These are the output parameters:
expname={}; %experiment name
mftrace=[]; %matrix containing photometry trace for each mouse (if multiple technical repeats exist for same mouse/animal, they are averaged)
mlowftrace=[]; %mftrace but low-pass filtered 
mreducedftrace=[]; %mftrace but downsampled to 1Hz
mpre=[]; %pre-stim value. each data point represents 1 biological repeat
mpost=[]; %post-stim value. each data point represents 1 biological repeat
pre_exprepeats=[]; %mpre but for all technical repeats

psthnomavg=[]; %normalized psth average
post_exprepeats=[];%mpost but for all technical repeats
psthnom_exprepeats=[]; %normalized psth for each technical repeat

mpsthnom=[]; %all normalized psth for biological repeats
mpsthtime=[]; %relative time for psth



fid = fopen(infile);
filename = fscanf(fid, '%s', [1]);
%itinerators
n=1;
i4=0;
for i2=1:length(mcode) %structuring
    %% process ftrace
    if i2>1;%structuring
        i4=sum(mcode(1:i2-1)); %structuring
    end %structuring
    
    for i3=i4+1:mcode(i2)+i4;
        fprintf(filename);
        folder=fscanf(fid, '%s', [1]); %which folder - each folder contains concatanates of the same trials data
        if isempty(dir(folder))==1
            fprintf(strcat('cannot find folder or it is empty:\t',folder)); 
            keyboard();
        end
        stimtime=fscanf(fid, '%s', [1]);%when stimulus starts
        stimtime=str2double(stimtime);
        
        srate = fscanf(fid, '%s', [1]);%what sampling rate
        srate=str2double(srate);
        
        laserp = fscanf(fid, '%s', [1]); % what laser power
        if strcmp(laserp(1),'f') %% unused functionality; keep these three lines for older txt files
            laserp=laserp(2:end);
            rigcode=1;
        end
        laserp=str2double(laserp);
        timebefore= fscanf(fid, '%s', [1]); % how long prestim time window for psth
        timebefore=str2double(timebefore);
        timeafter= fscanf(fid, '%s', [1]); % how long poststim time window for psth
        timeafter=str2double(timeafter);
        fchannel= fscanf(fid, '%s', [1]); % which data channel (we are using the same data for 2 parallel aquisition)
        fchannel=str2double(fchannel);
        %% parse data
        targetsrate=downsampleto; %temporal resolution to archieve, after downsampling
        drate=srate/targetsrate; %how much downsampling to do
        srate=targetsrate; % archieved temporal resolution
        %% saved preprocess filename
        ppfilename=strcat(filename,'_');
        ppfilename1=strcat(ppfilename,num2str(srate),'_raw.mat');
        if fchannel==3
            ppfilename2=strcat(ppfilename,num2str(srate),'_Z_raw.mat');
        else
            ppfilename2=strcat(ppfilename,num2str(srate),'_K_raw.mat');
        end
        %         keyboard();
        fprintf(num2str(fchannel));
        if isempty(which(ppfilename1)) && isempty(which(ppfilename2)) % is a preprocessed file already saved?
            [ftrace,~,time]=importdatfile_order_time( folder,[fchannel],drate );
            %         %% optional: if the trace is shorter than 1000 data points, something might be wrong and stop here to check whats going on
            %         if length(ftrace)<1000;
            %             keyboard();
            %         end
            save(ppfilename2,'ftrace','time'); % save the preprocess file
        else
            try
                load(ppfilename1);
            catch
                load(ppfilename2);
            end
        end
        %% saved low-pass filtered preprocess filename
        low_ppfilename1=strcat(ppfilename1(1:end-7),'low1d5hz_raw.mat')
        low_ppfilename2=strcat(ppfilename2(1:end-7),'low1d5hz_raw.mat')
        if isempty(which(low_ppfilename1)) && isempty(which(low_ppfilename2))
            %                 keyboard();
            lowftrace=lowpassfilter_2sec_100db(ftrace,srate);
            lowftrace=lowpassfilter_2sec_100db(ftrace,srate);
            
            save(low_ppfilename2,'lowftrace','time');
        else
            try
                load(low_ppfilename1);
            catch
                load(low_ppfilename2);
            end
        end
        expname{n}=filename; %% experiment names
        %% start analze - make psth; organize repeats into matrix/cells
        redftrace=decimate(ftrace,srate); %smaller file for illustrator
        mreducedftrace(1:length(redftrace),n)=redftrace;
        ftrace=afcorrection_firstrig(ftrace, laserp);%correct autofluorescence.
        lowftrace=afcorrection_firstrig(lowftrace, laserp);%correct autofluorescence
        
        mftrace(1:length(ftrace),n)=ftrace;
        mlowftrace(1:length(lowftrace),n)=lowftrace;
        mtime(1:length(ftrace),n)=time;
        
        halfwindow=60;
        [ ~,~,psthnom, psthtime ] = deltaps( ftrace,stimtime,timebefore, timeafter ,srate,halfwindow); %get PSTH and pre-post stim value
        [ pre,post,lowpsthnom, psthtime ] = deltaps( lowftrace,stimtime,timebefore, timeafter ,srate,halfwindow ); %get PSTH and pre-post stim values
        pre_exprepeats=[pre_exprepeats pre];
        post_exprepeats=[post_exprepeats post];
        psthnom_exprepeats=[psthnom_exprepeats; psthnom];
        filename = fscanf(fid, '%s', [1]);
    end
    
    %% organize all repeats in this experiment
    pre=mean(pre_exprepeats);
    post=mean(post_exprepeats);
    psthnom=mean(psthnom_exprepeats,1);
    pre_exprepeats=[];
    post_exprepeats=[];
    psthnom_exprepeats=[];
    mpre(1,n)=pre;
    mpost(1,n)=post;
    mpsthnom(1:length(psthnom),n)=psthnom;
    mpsthtime(1:length(psthtime),n)=psthtime;
    
    mlowpsthnom(1:length(lowpsthnom),n)=lowpsthnom;
    
    if pre>post %define the direction of change: +1 is increase
        direction=-1;
    else
        direction=1;
    end
    try
        [tau1,tau2]=calctau(psthnom,timebefore,direction,srate);
    catch
        tau1=0;
        tau2=0;
    end
    mtau1(n,1)=tau1;
    mtau2(n,1)=tau2;
    n=n+1;
end

psthnomavg=mean(mpsthnom,2);

%% prism ready
t=downsample(mlowpsthnom,10);
all_low_prism=zeros(size(t,1),3);
all_low_prism(:,1)=mean(t,2);
all_low_prism(:,2)=std(t,[],2);
all_low_prism(:,3)=t(:,1)*0+size(t,2);


fclose all

end
function y = lowpassfilter_2sec_100db(x,srate)
% low pass filter data with a filter 0.5Hz 100db
d = designfilt('lowpassfir', ...
    'PassbandFrequency',0.5,'StopbandFrequency',0.7, ...
    'PassbandRipple',1,'StopbandAttenuation',100, ...
    'DesignMethod','equiripple','SampleRate',srate);
y = filtfilt(d,x);
end


function [ cftrace ] = afcorrection_firstrig( ftrace, lp)
%AFCORRECTION corrects autofluorescence based on the data obtained from implant in non
%fluorescent brain tissue of ARC
%   Detailed explanation goes here
pw=[5:5:70];
%% IMPORTANT: this correction parameter is measured before starting experiments for fixed laser power values (5,10,15,20,25,30,35,40,45,50,55,60,65,70)
%% to do so, measure the fluorescence value of a cohort of photometry mice 5 days after surgery at each laser power. Average all reads and get this number to estimate the sum of tissue and system autofluorescence. Not meant to be accurate but does reduce internal bias.
%% another way to avoid this kind of issue is to use Z-score instead of dF/F, as Z-score is resilient to bias in baseline autofluorescence
co=[0.0918517258193450,0.178941925659473,0.256571270183852,0.340801358113509,0.416812901678658,0.489314857713830,0.565272600319742,0.642170184652280,0.710517359712228,0.792535508393284,0.861400621902477,0.927864410871303,0.994339078337328,1.07129420863309];
loc=find(pw==lp);
if isempty(loc)
    cftrace=ftrace-lp;
else
    co=co(loc);
    
    cftrace=ftrace-co;
end
end

function [ pre,post,psthnom, psthtime ] = deltaps( ftrace,stimtime,timebefore, timeafter,srate,halfavgwindow )
%CHANGEOFACTIVITY This function generate PSTH with timebefore mins and timeafter mins flanking stimtime
%pre was calculated as average of timebefore-1.5:timebefore+1.5 before the simulus and post is calculated as
%average of timeafter-1.5:timeafter+1 .5
%   Detailed explanation goes here
% keyboard();
try pre=median(ftrace(srate*(stimtime-timebefore-halfavgwindow):srate*(stimtime-timebefore+1*halfavgwindow)));
catch
    keyboard();
    pre=median(ftrace(1:srate*(stimtime-timebefore+halfavgwindow)));
end
try post=median(ftrace(srate*(stimtime+timeafter-halfavgwindow):srate*(stimtime+timeafter+halfavgwindow)));
catch
    fprintf(num2str(length(ftrace)-timeafter));
    post=median(ftrace(srate*(stimtime+timeafter-halfavgwindow):length(ftrace)));
    fprintf(num2str(length(ftrace)-timeafter));
end
try
    psthnom=ftrace(stimtime*srate-srate*timebefore:stimtime*srate+srate*timeafter)/pre;
catch
    psthnom=ftrace(stimtime*srate-srate*timebefore:end)/pre;
end
psthtime=[1:1/srate:timebefore+timeafter+1]-timebefore;
end

function [ ftrace1,ottl1,time1 ] = importdatfile_order_time( folder, ch2,decimation )
%Concatenate .dat files from same experiments based on their temporal sequence
%   Detailed explanation goes here
list=dir(folder);
ftrace1=[];
ottl1=[];
time1=[];
%% reads the date/time of the dat files. System is IN US FORMAT
for i=1:length(list)
    a=list(i).date(end-7:end);
    d=list(i).date(1:2);
    d=str2double(d);
    a=strcat(a(1:2),a(4:5),a(7:8));
    a=str2double(a);
    date(i)=timeconvert(a)+d*1000000;
end

%% Concatenate data based on time acquired
for i=1:length(list);
    b=find(date==min(date));
    date(b)=inf;
    namet=list(b).name;
    if isempty(strfind(namet,'.dat'))==0
        clear datat
        namet
        datat=importdata(namet);
        ftracet=datat.data(1:end,ch2(1));
        ftracet=decimate(ftracet,decimation);
        ftrace1=cmbarray(ftrace1,ftracet);
        if length(ch2)>1
            ttlt=datat.data(1:end,ch2(2));
            ttime=find(ttlt>3);
            ttime=round(ttime/decimation);
            newttl=zeros(1,round(length(ttlt)/decimation));
            for x=1:length(ttime)
                try newttl(ttime(x))=1;
                end
            end
            ottl1=cmbarray(ottl1,newttl);
        end
        timet=datat.data(1:end,1);
        timet=downsample(timet,decimation);
        time1=cmbarray(time1,timet);
        clear ttime ttlt ftracet newttl1 timet
    end
end
end

function [ parray ] = timeconvert( tarray )
%Covert file time (US format) to a number
%   Detailed explanation goes here
for i=1:length(tarray);
    str=sprintf('%0.0f',tarray(i));
    sec=str(end-1:end);
    sec=str2num(sec);
    tt=sec;
    if length(str)>=3
        if length(str)==3;
            minu=str(end-2);
        else
            minu=str(end-3:end-2);
        end
        minu=str2num(minu);
        tt=minu*60+sec;
        if length(str)>=5
             if length(str)==5;
             hr=str(end-4);
             else
             hr=str(end-5:end-4);
             end  
             hr=str2num(hr);
            tt=tt+60*60*hr;
        end
    end
    hr=0;
    parray(i)=tt;
    tt=0;
end
end

function [ aout ] = cmbarray( a1,a2 )
%CMBARRAY Summary of this function goes here
%   Detailed explanation goes here
long=length(a1);
aout=a1;
for  i=1:length(a2);
    aout(long+i)=a2(i);
end
end




