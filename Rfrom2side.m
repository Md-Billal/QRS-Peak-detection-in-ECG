function [y_rec,R_locs11]=Rfrom2side(X,Fs)

%Reference paper
%M. B. Hossain, S. K. Bashar, A. J. Walkey, D. D. McManus, and K. H.
% Chon, “An Accurate QRS Complex and P Wave Detection in ECG
% Signals Using Complete Ensemble Empirical Mode Decomposition with
% Adaptive Noise Approach,” IEEE Access, vol. 7, pp. 128869–128880,
% 2019.

%X is the input data
%Fs is sampling frequency
%R_locs11 is the R peak locations
%y_rec is the reconstructed signal using CEEMDAN decomposition
 


N=length(X);
X1=X(1:N);
X1=X1-mean(X1);
X1=X1/max(abs(X1));

%% remove baseline wandering
epsilon1= round(0.6*Fs);
epsilon2= round(0.2*Fs);
X2=[zeros(1,epsilon1-1) X1' zeros(1,epsilon1-1)];
out1=movmedian(X2,epsilon1);
y1=out1(epsilon1:epsilon1+N-1);
X3=[zeros(1,epsilon2-1) X1' zeros(1,epsilon2-1)];
out2=movmedian(X3,epsilon2);
y2=out2(epsilon2:epsilon2+N-1);
yy=X1(1:N)'-0.5*y1-0.5*y2;

mpeaks=find(abs(yy)>0.2*max(abs(yy)));
av=mean(yy(mpeaks));
if av<0
    yy=(-1)*yy;
end

sd=sqrt(var(yy));
[modes,its]= ceemdan1(yy,2*sd,10,100,1);
y_rec=sum(modes(2:end,:),1);
% y_rec=y_rec-mean(y_rec);
% y=X1';
yy=sum(modes(1:4,:),1);
y=yy;
mpeaks=find(abs(yy)>0.2*max(abs(yy)));
sorted_peaks=sort(y(mpeaks),'descend');
if length(sorted_peaks)>200
th=mean(sorted_peaks(1:200));
else
th=mean(sorted_peaks(1:end));
end
% plot(yy)
yy(1:5)=0;
y_rec(1:5)=0;
yy(N-5:N)=0;
y_rec(N-5:N)=0;
%% Create mask
yy1=yy;
yy2=yy;
yy1(yy<0)=0;
yy2(yy>0)=0;
yy2=yy2*(-1);



[val1 ind1]=findpeaks(abs(yy1));
[val2 ind2]=findpeaks(abs(yy2));

% figure
% histogram(val1,170);
[a b]=histcounts(val1,200);
[A B]=histcounts(val2,200);
for k=1:length(b)-1
    b1(k)=(b(k)+b(k+1))/2;
end
for k=1:length(B)-1
    B1(k)=(B(k)+B(k+1))/2;
end
sump=sum(a.*b1);
sumy=sum(a);
t=sump/sumy;

sump1=sum(A.*B1);
sumy1=sum(A);
t1=sump1/sumy1;
mask=1*(abs(yy1)>=t)+0*(abs(yy1)<t);
mask1=1*(abs(yy2)>=t1)+0*(abs(yy2)<t1);
K=find(mask==1);
K1=find(mask1==1);
for i=1:length(K)-1
    if K(i+1)-K(i)<=0.07*Fs
        mask(K(i):K(i+1))=1;
    end
end

for i=1:length(K1)-1
    if K1(i+1)-K1(i)<=0.07*Fs
        mask1(K1(i):K1(i+1))=1;
    end
end

%%find R peaks and do post processing
MaskR=mask.*yy1;
MaskR1=mask1.*yy2;
[Rval R_loc]=findpeaks(mask.*yy1,'MinPeakHeight',1.5*t);
[Rval1 R_loc_neg]=findpeaks(mask1.*yy2,'MinPeakHeight',1.5*t1);
allR=[R_loc R_loc_neg];
allR=sort(allR,'ascend');


R_locs(1)=allR(1);
temp=allR(1);
k=1;
for i=1:length(allR)-1
    if allR(i+1)<(temp+0.3*Fs)
        if abs(yy(allR(i+1)))>abs(yy(temp))
            temp=allR(i+1);
            R_locs(k)=temp;
        end
    else
            k=k+1;
            temp=allR(i+1);
            R_locs(k)=temp;
            
    end
end
R_locup=[];
for k=1:length(R_locs)-1
     if (R_locs(k+1)-R_locs(k))>3*mean(diff(R_locs))
        [~,pos]=findpeaks(abs(yy(R_locs(k)+round(0.07*Fs):R_locs(k+1)-round(0.07*Fs))),'MinPeakHeight',0.1*mean(abs(yy(R_locs))));
        if ~isempty(pos)
            [~,pos1]=max(abs(yy(R_locs(k)+round(0.07*Fs)+pos-1)));
            R_locup=[R_locup pos+R_locs(k)+round(0.07*Fs)-1];
        end
    
    elseif (R_locs(k+1)-R_locs(k))>2*mean(diff(R_locs))
        [~,pos]=findpeaks(abs(yy(R_locs(k)+round(0.07*Fs):R_locs(k+1)-round(0.07*Fs))),'MinPeakHeight',0.28*mean(abs(yy(R_locs))));
        if ~isempty(pos)
            [~,pos1]=max(abs(yy(R_locs(k)+round(0.07*Fs)+pos-1)));
            R_locup=[R_locup pos(pos1)+R_locs(k)+round(0.07*Fs)-1];
        end
    elseif (R_locs(k+1)-R_locs(k))>1.4*mean(diff(R_locs))
        [~,pos]=findpeaks(abs(yy(R_locs(k)+round(0.07*Fs):R_locs(k+1)-round(0.07*Fs))),'MinPeakHeight',0.28*mean(abs(yy(R_locs))));
        if ~isempty(pos)
            [~,pos1]=max(abs(yy(R_locs(k)+round(0.07*Fs)+pos-1)));
            R_locup=[R_locup pos(pos1)+R_locs(k)+round(0.07*Fs)-1];
        end
    end
end
R_locs=[R_locs R_locup];
R_locs=sort(R_locs,'ascend');

R_locs11=R_locs(1);
temp=R_locs(1);
k=1;
for i=1:length(R_locs)-1
    if R_locs(i+1)<(temp+round(0.27*Fs))
        if abs(yy(R_locs(i+1)))>abs(yy(temp))
            temp=R_locs(i+1);
            R_locs11(k)=temp;
        end
    elseif R_locs(i+1)<(temp+round(0.4*Fs))
        if abs(yy(R_locs(i+1)))> 0.7*mean(abs(yy(R_locs)))
            k=k+1;
            temp=R_locs(i+1);
            R_locs11(k)=temp;
        elseif abs(yy(R_locs(i+1)))>abs(yy(temp))
            temp=R_locs(i+1);
            R_locs11(k)=temp;
        end
    else
            k=k+1;
            temp=R_locs(i+1);
            R_locs11(k)=temp;
    end
end
for j=1:length(R_locs11)
    if R_locs11(j)-10>0 && R_locs11(j)+10<length(y)
        [~,pos]=max(abs(y_rec(R_locs11(j)-10:R_locs11(j)+10)));
        R_locs11(j)=pos+R_locs11(j)-10-1;
    elseif R_locs11(j)-10<0 && R_locs11(j)+10<length(y)
        [~,pos]=max(abs(y_rec(1:R_locs11(j)+10)));
        R_locs11(j)=pos;
    elseif R_locs11(j)-10>0 && R_locs11(j)+10>length(y)
        [~,pos]=max(abs(y_rec(R_locs11(j)-10:end)));
        R_locs11(j)=pos+R_locs11(j)-10-1;
    end
end
MeanR=mean(abs(y_rec(R_locs11)));
R_locs11=R_locs11(find(abs(y_rec(R_locs11))>0.35*MeanR));


