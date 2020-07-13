%% Theoretical F calculation for M-PAM
M=2;
y=[];
b_in = randi([0 M-1], 2000, 1);
x = pammod(b_in, M);
x=x/rms(x);
theoretical_F_2PAM=[];
t_range=linspace(-1.5,1.5,2008);
for t=1:length(t_range)
    temp=x<=t_range(t);
    temp_avg=mean(temp);
    theoretical_F_2PAM(t)=temp_avg;
end
figure
plot(t_range,theoretical_F_2PAM)

M=4;
y=[];
b_in = randi([0 M-1], 2000, 1);
x = pammod(b_in, M);
x=x/rms(x);
theoretical_F_4PAM=[];
t_range=linspace(-1.5,1.5,2008);
for t=1:length(t_range)
    temp=x<=t_range(t);
    temp_avg=mean(temp);
    theoretical_F_4PAM(t)=temp_avg;
end
figure
plot(t_range,theoretical_F_4PAM)
%% PAM
fc=2*10^6;
fs=4*10^6;
Ts=1/fs;
SNR_range=0:20;
P_PAM=[];
M_range=[2,4];
t_range=linspace(-1.5,1.5,2008);

for s=1:length(SNR_range)
    disp(SNR_range(s))
    SNR_lin=10^(SNR_range(s)/10);
    correct_class=0;
    for m=1:1000
        in_type=randi(length(M_range),1);
        M=M_range(in_type);
%         fprintf('The tx modulation level is %d \n',M)
        y=[];
        
        b_in = randi([0 M-1], 2000, 1);
        x = pammod(b_in, M);
        x=x/rms(x);
        xpow=rms(x);
        x_upsampled=upsample(x,4);

        filtr=rcosdesign(0.5,8,4);
        filtr=filtr/norm(filtr);
        
        x_filtered = conv(x_upsampled,filtr);
        x_filtered=x_filtered/rms(x_filtered);

        for n=1:length(x_filtered)
            noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*200))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
            y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
        end
        y=downsample(y,4);
        fun=[];
        
        for t=1:length(t_range)
            temp=y<=t_range(t);
            temp_avg=mean(temp);
            fun(t)=temp_avg;
        end

        dist=[];
        dist=[dist norm(fun-theoretical_F_2PAM)];
        dist=[dist norm(fun-theoretical_F_4PAM)];
        [dist_val, dist_index]=min(dist);
%         fprintf('The tx modulation level is %d \n',M_range(dist_index))
        if M==M_range(dist_index)
        correct_class=correct_class+1;
        end
    end
    
    acc=correct_class/m;
    P_PAM=[P_PAM acc];
    
end
%%
figure
plot(SNR_range,P_PAM,'linewidth',2)
ylim([0 1])
xlabel('SNR (dB)','FontSize',13), ylabel('Probability of Classification','FontSize',13)
title('Performance of MLC block for PAM','FontSize',15)