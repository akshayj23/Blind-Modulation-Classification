fc=2*10^6;
fs=4*10^6;
Ts=1/fs;
Nm=8000;

SNR_range=[0,10,20];
P_detected=zeros(3,21); 
P_falsealarm=zeros(3,21);
gamma=linspace(-2,0,20);
for s=1:length(SNR_range)
    SNR_lin=10^(SNR_range(s)/10);
    g=1;
    for ga=1:length(gamma)
        OFDM_false_detected_count=0;
        OFDM_correct_detected_count=0;
        OFDM_tx_count=0;
        non_OFDM_tx_count=0;
        disp(gamma(ga))
    %     if mod(gamma,0.001)==0
    %         disp(gamma)
    %     end  
        for m=1:100
            input_type = randi([1 5]); % OFDM:1, 4-QAM:2, 16-QAM:3, 2-PAM:4, 4-PAM:5
            if input_type==1
                OFDM_tx_count=OFDM_tx_count+1;
        %         disp('The true txed signal is OFDM');
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
                non_OFDM_tx_count=non_OFDM_tx_count+1;
        %         disp('The true txed signal is QAM')
                % QAM
                M=4;
                b_in = randi([0 1], log2(M)*2000, 1);
                x = qammod(b_in, M,'InputType','bit','UnitAveragePower', true);

                xpow=rms(x);
                x_upsampled=upsample(x,4);

                filtr=rcosdesign(0.5,8,4,'sqrt');
                filtr=filtr/norm(filtr);
                filtpow=rms(filtr);

                x_filtered = conv(x_upsampled,filtr);
                x_filtered = x_filtered/rms(x_filtered);
                x_fipow=rms(x_filtered);
                
                y=[];
                for n=1:length(x_filtered)
                    noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                    y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
                end
                
                
            elseif input_type==3
                non_OFDM_tx_count=non_OFDM_tx_count+1;
        %         disp('The true txed signal is QAM')
                % QAM
                M=16;
                b_in = randi([0 1], log2(M)*2000, 1);
                x = qammod(b_in, M,'InputType','bit','UnitAveragePower', true);

                xpow=rms(x);
                x_upsampled=upsample(x,4);

                filtr=rcosdesign(0.5,8,4,'sqrt');
                filtr=filtr/norm(filtr);
                filtpow=rms(filtr);

                x_filtered = conv(x_upsampled,filtr);
                x_filtered = x_filtered/rms(x_filtered);
                x_fipow=rms(x_filtered);
                
                y=[];
                for n=1:length(x_filtered)
                    noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                    y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
                end

            elseif input_type==4
                non_OFDM_tx_count=non_OFDM_tx_count+1;
        %         disp('The true txed signal is PAM')
                % PAM
                M=2;
                b_in = randi([0 M-1], 2000, 1);
                x = pammod(b_in, M);
                x=x/rms(x);
                x=upsample(x,4);
                
                filtr=rcosdesign(0.5,8,4,'sqrt');
                filtr=filtr/norm(filtr);
                x_filtered = conv(x,filtr);
                x_filtered=x_filtered/rms(x_filtered);
                
                y=[];
                for n=1:length(x_filtered)
                    noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                    y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
                end
                
            elseif input_type==5
                non_OFDM_tx_count=non_OFDM_tx_count+1;
        %         disp('The true txed signal is PAM')
                % PAM
                M=4;
                b_in = randi([0 M-1], 2000, 1);
                x = pammod(b_in, M);
                x=x/rms(x);
                x=upsample(x,4);
                
                filtr=rcosdesign(0.5,8,4,'sqrt');
                filtr=filtr/norm(filtr);
                x_filtered = conv(x,filtr);
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

            if C42>gamma(ga)
                det=1;
            else
                det=2;
            end

            if input_type==1 && det==1
                OFDM_correct_detected_count=OFDM_correct_detected_count+1;
            elseif input_type~=1 && det==1
                OFDM_false_detected_count=OFDM_false_detected_count+1;
            end
        end

        P_detected(s,g)=OFDM_correct_detected_count/OFDM_tx_count;
        P_falsealarm(s,g)=OFDM_false_detected_count/non_OFDM_tx_count;
        g=g+1;

    end
end
%% 
figure(3)
plot(P_falsealarm(1,:), P_detected(1,:),'-x', 'linewidth',2)
hold on
plot(P_falsealarm(2,:), P_detected(2,:),'-o', 'linewidth',2)
hold on
plot(P_falsealarm(3,:), P_detected(3,:), 'linewidth',2)
legend('SNR 0dB','SNR 10dB','SNR 20dB', 'FontSize',12);
xlabel('Probability of false alarm','FontSize',13), ylabel('Probability of detection','FontSize',13)
title('Performance of MC-SC block with Nm=8000','FontSize',15)
%% MTC

T=1*10^-6;
N=8000;
SNR_full_range=0:20;
P_MTC=zeros(1,length(SNR_full_range));
V=[1 0 0 0 0 0; 1 1 0 1 0 1; 1 0 1 0 1 0];
for s=1:21
    SNR_lin=10^(SNR_full_range(s)/10);
    correctly_classified=0;
    disp(s)
    for m=1:1000
        input_type = randi([1 2]); % QAM-1, PAM-2
        if input_type==1
%             disp('The true txed signal is QAM Class 1')
            % QAM
            M=16;
            b_in = randi([0 1], log2(M)*2000, 1);
            x = qammod(b_in, M,'InputType','bit','UnitAveragePower', true);

            xpow=rms(x);
            x_upsampled=upsample(x,4);

            filtr=rcosdesign(0.5,8,4,'sqrt');
            filtr=filtr/norm(filtr);
            filtpow=rms(filtr);

            x_filtered = conv(x_upsampled,filtr);
            x_filtered = x_filtered/rms(x_filtered);
            x_fipow=rms(x_filtered);

            y=[];
            for n=1:length(x_filtered)
                noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
            end

        elseif input_type==2
            non_OFDM_tx_count=non_OFDM_tx_count+1;
%             disp('The true txed signal is PAM Class 2')
            % PAM
            M=2;
            b_in = randi([0 M-1], 2000, 1);
            x = pammod(b_in, M);
            x=x/rms(x);
            x=upsample(x,4);

            filtr=rcosdesign(0.5,8,4,'sqrt');
            filtr=filtr/norm(filtr);
            x_filtered = conv(x,filtr);
            x_filtered=x_filtered/rms(x_filtered);

            y=[];
            for n=1:length(x_filtered)
                noise = sqrt(rms(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))^2/(2*SNR_lin))*(randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))) + 1j* randn(size(x_filtered(n)*exp(1j*2*pi*fc*n*Ts))));
                y(n)=x_filtered(n)*exp(1j*2*pi*fc*n*Ts)+noise;
            end
        end
        
        F=[];
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
        
        if (index==1 && input_type==1) || (index==2 && input_type==2)
            correctly_classified=correctly_classified+1;
%             fprintf('Det class %d',index)
        end
    end
    P_MTC(s)=correctly_classified/m;
end
%%
figure
plot(SNR_full_range,P_MTC,'linewidth',2)
ylim([0 1])
xlabel('SNR (dB)','FontSize',13), ylabel('Probability of Classification','FontSize',13)
title('Performance of MTC block with N=8000','FontSize',15)