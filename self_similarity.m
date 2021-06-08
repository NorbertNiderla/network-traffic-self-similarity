close all
clear all

ss_hist = @(a,r,d,c) histcounts(d(d>a & d<(a+c*r)),c);
ss_xcov = @(s, h, k) (s^2)/2*((k+1).^(2*h) - 2*k.^(2*h) + (k-1).^(2*h));


traffic_data_type = 'TCP'; %DNS | TCP

switch traffic_data_type
    case 'DNS'
        load('dns_timestamps.mat');
        time = dns_timestamps;
    case 'TCP'
        load('tcp_timestamps.mat'); 
        time = tcp_timestamps;
    otherwise
        error('self similarity: wrong traffic data type')
end

%% HISTOGRAMS
normal_time = cumsum(abs(randn(length(time),1)*std(diff(time))));
normal_time = normal_time / max(normal_time) * max(time);
a = 30;

switch traffic_data_type
    case 'TCP'
        r = [0.01 0.1 1 10];
    case 'DNS'
        r = [0.1 1 10 60];
end

histograms_fig = figure;
for d = 1:4
    hist = ss_hist(a,r(d),time,50);
    hist_g = ss_hist(a,r(d),normal_time,50);
    t = linspace(a,a+length(hist)*r(d),length(hist));
    subplot(2,2,d)
    stairs(t', [hist', hist_g'],'LineWidth',3);
    xlabel('t [s]')
    ylabel('Packet count')
    legend('Traffic data', 'Generated data (Gauss)')
    title(strcat('Resolution: ',string(r(d)),'s'));
    grid on
end

clear hist hist_g r t d a

%% AUTOCOVARIANCE
b = 1200;
histd = (histcounts(time,b));
s = std(histd);
histd = xcov(histd,'normalized');
histd = histd(b:end);
k = 0:b-1;
h = 0.7;
ada = ss_xcov(s,h,k);
ada = ada/max(ada);

autocovariance_fig = figure;
plot(k,ada,k,histd, 'LineWidth', 3)
title('Autocovariance of packet count');
xlabel('k')
legend('Second-order self similarity, H = 0.7', 'Traffic data')
grid on

clear b histd s k h ada

%% H PARAMETER
h_fig = figure;
apr_h = [];
apr_k = [];
for x = 1:6
    switch traffic_data_type
        case 'DNS'
            b = 400;
            histd = (histcounts(time(time<(x*600)&time>(x-1)*600),b));
        case 'TCP'
            b = 1200;
            histd = (histcounts(time(time<(x*200)&time>(x-1)*200),b));
    end
    s = std(histd);
    histd = xcov(histd,'normalized');
    histd = histd(b:end);

    k = 0:b-1;
    [~,i]=find(histd<0);
    k = k(1:i(1));
    clear i
    h = 0.51:0.001:0.99;

    k_idx = 1;
    h_p = [];
    v_p = [];
    for kx = k
       ada_h = ss_xcov(s,h,kx);
       [v,i] = min(abs(abs(histd(k_idx))-ada_h/max(ada_h)));
       v_p = [v_p v];
       h_p = [h_p h(i)];
       k_idx = k_idx + 1;
    end
    
    apr_h = [apr_h h_p];
    apr_k = [apr_k k];
    plot(k,h_p,'.', 'LineWidth', 3)
    hold on
    
end

approx = polyfit(apr_k,apr_h,3);
plot(0:max(apr_k),polyval(approx,0:max(apr_k)),'-', 'LineWidth', 3)
grid on
ylim([0.5 1])
ylabel('H')
xlabel('k')
title('H parameter estimation')
legend('H_T_1','H_T_2','H_T_3','H_T_4','H_T_5','H_T_6','Approx.','Location','southwest')
clear b histd k h k_idx h_p v_p ada_h k_idx apr_h apr_v