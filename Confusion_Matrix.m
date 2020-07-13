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
plot(t_range,theoretical_F_4PAM)
%% Theoretical F calculation for M-QAM
M=4;
b_in = randi([0 1], log2(M)*2000, 1);
x = qammod(b_in, M,'InputType','bit','UnitAveragePower', true);
theoretical_F_4QAM=[];
t_range=linspace(-1.5,1.5,2008);
for t=1:length(t_range)
    temp=x<=t_range(t);
    temp_avg=mean(temp);
    theoretical_F_4QAM(t)=temp_avg;
end
plot(t_range,theoretical_F_4QAM)

M=16;
y=[];
b_in = randi([0 1], log2(M)*2000, 1);
x = qammod(b_in, M,'InputType','bit','UnitAveragePower', true);
theoretical_F_16QAM=[];
t_range=linspace(-1.5,1.5,2008);
for t=1:length(t_range)
    temp=x<=t_range(t);
    temp_avg=mean(temp);
    theoretical_F_16QAM(t)=temp_avg;
end
plot(t_range,theoretical_F_16QAM)
%%
fc=2*10^6;
fs=4*10^6;
Ts=1/fs;
Nm=8000;

SNR_range=0:20;
MQAM_range=[4,16];
MPAM_range=[2,4];
V=[1 0 0 0 0 0; 1 1 0 1 0 1; 1 0 1 0 1 0];
t_range=linspace(-1.5,1.5,2008);
fin_prob=[];
gamma=-0.2;
for s=1:length(SNR_range)
    disp(s)
    SNR_lin=10^(SNR_range(s)/10);
    confusionmatrix=zeros(5);
    for m=1:100
        input_type = randi([1 5]); % OFDM:1, 4-QAM:2, 16-QAM:3, 4-PAM:4, 8-PAM:5
%         input_type=3;
        if input_type==1
            %OFDM
            for o=1:2:49
                for k=1:64
                    b_in(o:o+1,k) = randi([0 1], log2(4)*1, 1);
                    Xk((o+1)/2,k) = qammod(b_in(o:o+1,k), 4, 'InputType', 'bit', 'UnitAveragePower', true);
                end
            end
            for o=1:25
                Xkifft(o,:)= ifft(Xk(o,:), 64);
                xdatacp(o,:)=[Xkifft(o,end-16+1:end) Xkifft(o,:)];
            end
            xdata=[];
            for i=1:25
                xdata=[xdata xdatacp(i,:)];
            end
            xdata=upsample(xdata,4);
            y=[];
            for n=1:length(xdata)
                noise = sqrt(rms(xdata(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(xdata(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(xdata(n)*exp(1j*2*pi*fc*n*Ts))));
                y(n)=xdata(n)*exp(1j*2*pi*fc*n*Ts)+noise;
            end

        elseif input_type==2
            % 4-QAM
%             disp('TX MQAM')
            M=4;
%             fprintf('The tx modulation level is %d \n',M)
            y=[];
            b_in = randi([0 1], log2(M)*2000, 1);
            x = qammod(b_in, M,'InputType','bit','UnitAveragePower', true);

            xpow=rms(x);
            x_upsampled=upsample(x,4);

            filtr=rcosdesign(0.5,8,4);
            filtr=filtr/norm(filtr);
            filtpow=rms(filtr);

            x_filtered = conv(x_upsampled,filtr);
            x_filtered = x_filtered/rms(x_filtered);
            x_fipow=rms(x_filtered);
            % x_filteredpow=rms(x_filtered);

            for n=1:length(x_filtered)
                noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
            end
            
        elseif input_type==3
            % 16-QAM
%             disp('TX MQAM')
            M=16;
            y=[];
            b_in = randi([0 1], log2(M)*2000, 1);
            x = qammod(b_in, M,'InputType','bit','UnitAveragePower', true);

            xpow=rms(x);
            x_upsampled=upsample(x,4);

            filtr=rcosdesign(0.5,8,4);
            filtr=filtr/norm(filtr);
            filtpow=rms(filtr);

            x_filtered = conv(x_upsampled,filtr);
            x_filtered = x_filtered/rms(x_filtered);
            x_fipow=rms(x_filtered);
            % x_filteredpow=rms(x_filtered);

            for n=1:length(x_filtered)
                noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;

            end

        elseif input_type==4
            % 2-PAM
            y=[];
            M=2;
            b_in = randi([0 M-1], 2000, 1);
            x = pammod(b_in, M);
            x=x/rms(x);
            xpow=rms(x);
            x_upsampled=upsample(x,4);

            filtr=rcosdesign(0.5,8,4);
            filtr=filtr/norm(filtr);
            x_filtered = conv(x_upsampled,filtr);
            x_filtered=x_filtered/rms(x_filtered);
            x_fipow=rms(x_filtered);

            for n=1:length(x_filtered)
                noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
            end

        elseif input_type==5
            % 4-PAM
            M=4;
            b_in = randi([0 M-1], 2000, 1);
            x = pammod(b_in, M);
            x=x/rms(x);
            x_upsampled=upsample(x,4);

            filtr=rcosdesign(0.5,8,4);
            filtr=filtr/norm(filtr);
            x_filtered = conv(x_upsampled,filtr);
            x_filtered=x_filtered/rms(x_filtered);
            y=[];
            for n=1:length(x_filtered)
                noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
            end

        end
        
        C20=mean(y(1:Nm).^2);
        C21=mean(abs(y(1:Nm).^2));
        C42=mean(abs(y(1:Nm)).^4)-abs(C20)^2-2*C21^2;

        if C42>gamma
            confusionmatrix(input_type,1)=confusionmatrix(input_type,1)+1;
            
        else
            temp_sum=0;
            T=1*10^-6;
            F=[];
            N=2000;
            alpha_set=[1/T,2*fc-1/T,2*fc-1/(2*T),2*fc,2*fc+1/(2*T),2*fc+1/T];
            for alpha=1:6
                temp_sum=0;
                if alpha_set(alpha)==1/T
                    for n=0:N-1
                        temp_sum=temp_sum+(y(n+1)*conj(y(n+1))*exp(-1j*2*pi*alpha_set(alpha)*n*Ts));
                    end
                    temp_avg=temp_sum/N;
                    F=[F abs(temp_avg)];
                else
                    for n=0:N-1
                        temp_sum=temp_sum+(y(n+1)*y(n+1)*exp(-1j*2*pi*alpha_set(alpha)*n*Ts));
                    end
                    temp_avg=temp_sum/N;
                    F=[F abs(temp_avg)];
                end
            end
            Fbar=F./norm(F);
            
            dist=[];
            for i=1:3
                dist=[dist (norm(Fbar-V(i,:)./norm(V(i,:))))^2];
            end
            [value, index]=min(dist);
            
            y=downsample(y,4);
            if index==1
%                 disp('Rx MQAM')
                %QAM detected
                fun=[];
                t_range=linspace(-1.5,1.5,2008);
                for t=1:length(t_range)
                    temp=x<=t_range(t);
                    temp_avg=mean(temp);
                    fun(t)=temp_avg;
                end

                dist=[];
                dist=[dist norm(fun-theoretical_F_4QAM)];
                dist=[dist norm(fun-theoretical_F_16QAM)];
                [dist_val, dist_index]=min(dist);
        
%                 fprintf('The rx modulation level is %d \n',M_range(dist_index))
                
%                 fprintf('The rx modulation level is %d \n',MQAM_range(dist_index))
                if MQAM_range(dist_index)==4
                    confusionmatrix(input_type,2)=confusionmatrix(input_type,2)+1;
                elseif MQAM_range(dist_index)==16
                    confusionmatrix(input_type,3)=confusionmatrix(input_type,3)+1;
                end

            elseif index==2
%                 disp('Rx MPAM')
                %PAM detected
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
                
                if MPAM_range(dist_index)==2
                    confusionmatrix(input_type,4)=confusionmatrix(input_type,4)+1;
                elseif MPAM_range(dist_index)==4
                    confusionmatrix(input_type,5)=confusionmatrix(input_type,5)+1;
                end
            end
        end
    end
    prob_matrix=diag(1./sum(confusionmatrix,2))*confusionmatrix;
    if mod(SNR_range(s),10)==0
        fprintf('The confusion matrix for SNR %d dB and \x03b3=%.1f is \n',SNR_range(s),gamma);
        T = array2table(prob_matrix,'VariableNames',{'OFDM','4-QAM','16-QAM','2-PAM','4-PAM'},'RowNames',{'OFDM','4-QAM','16-QAM','2-PAM','4-PAM'});
        disp(T)
    end
    fin_prob=[fin_prob trace(prob_matrix)/5];
    
end
%%
plot(SNR_range,fin_prob,'linewidth',2)
xlabel('SNR (dB)','FontSize',13), ylabel('Probability of Classification','FontSize',13)
titl=sprintf('Performance of the entire classifier for \x03b3 = %.1f and Nm=%d',gamma,Nm);
title(titl,'FontSize',15)