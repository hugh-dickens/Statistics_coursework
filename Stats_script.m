% Load data
load 'CLNR_Jan_2012.mat';

%Rename matrices from extracted data
[G_4am_cn, ~] = size(clnr_X1(:,1));      %Mosaic G at 4am number of customers
G_4am_ec = clnr_X1(:,2);                 %Mosaic G at 4am electricity consumption
[G_10am_cn,~] = size(clnr_X2(:,1));      %Mosaic G at 10am number of customers
G_10am_ec = clnr_X2(:,2);                %Mosaic G at 10am electricity consumption
[I_4am_cn,~] = size(clnr_Y1(:,1));       %Mosaic I at 4am number of customers
I_4am_ec = clnr_Y1(:,2);                 %Mosaic I at 4am electricity consumption
[I_10am_cn,~] = size(clnr_Y2(:,1));      %Mosaic I at 10am number of customers
I_10am_ec = clnr_Y2(:,2);                %Mosaic I at 10am electricity consumption

%Calculate means
G_4am_mean = mean(G_4am_ec);
G_10am_mean = mean(G_10am_ec);
I_4am_mean = mean(I_4am_ec);
I_10am_mean = mean(I_10am_ec);

%Calculate standard deviations
G_4am_std = std(G_4am_ec);
G_10am_std = std(G_10am_ec);
I_4am_std = std(I_4am_ec);
I_10am_std = std(I_10am_ec);

% Find an interval that contains 95% of the values from a standard normal distribution.
y = norminv(0.975);

% Find sample mean within 95% confidence interval
G_4am_mu1 = G_4am_mean + y*((G_4am_std) /(sqrt(G_4am_cn)));
G_4am_mu2 = G_4am_mean - y*(G_4am_std /sqrt(G_4am_cn));
G_10am_mu1 = G_10am_mean + y*(G_10am_std /sqrt(G_10am_cn));
G_10am_mu2 = G_10am_mean - y*(G_10am_std /sqrt(G_10am_cn));
I_4am_mu1 = I_4am_mean + y*(I_4am_std /sqrt(I_4am_cn));
I_4am_mu2 = I_4am_mean - y*(I_4am_std /sqrt(I_4am_cn));
I_10am_mu1 = I_10am_mean + y*(I_10am_std /sqrt(I_10am_cn));
I_10am_mu2 = I_10am_mean - y*(I_10am_std /sqrt(I_10am_cn));

%Create title headings and matrices to put calculated values straight into a table 
Sample = {'Mosaic G at 4am';'Mosaic G at 10am';'Mosaic I at 4am';'Mosaic I at 10am'};
Mean = [G_4am_mean;G_10am_mean;I_4am_mean;I_10am_mean];
Std = [G_4am_std;G_10am_std;I_4am_std;I_10am_std];
NumberDataPoints = [G_4am_cn;G_10am_cn;I_4am_cn;I_10am_cn];
ConfidenceInterval1 = [G_4am_mu1;G_10am_mu1;I_4am_mu1;I_10am_mu1];
ConfidenceInterval2 = [G_4am_mu2;G_10am_mu2;I_4am_mu2;I_10am_mu2];

%Print table for user to identify results
T = table(Mean,Std,NumberDataPoints,ConfidenceInterval1,ConfidenceInterval2,'RowNames',Sample)
           
%Filter out any data above 1.0KWh
inRange1 = G_4am_ec(G_4am_ec < 1); 
inRange2 = G_10am_ec (G_10am_ec < 1);
inRange3 = I_4am_ec (I_4am_ec < 1);
inRange4 = I_10am_ec (I_10am_ec < 1);

%Create equally spaced vector to plot log normal probability density
%function against on x-axis (making smooth curve)
xt = 0:0.01:1;

%Plot for Mosaic G at 4am
subplot (2,2,1)
histogram(inRange1,100,'Normalization','pdf')
title ('Mosaic G at 4am')
xlim(0:1)
ylabel('Probability Density')
xlabel('Energy Consumption/ kWh')

%Plot for Mosaic G at 10am
subplot (2,2,2)
histogram(inRange2,100,'Normalization','pdf')
title ('Mosaic G at 10am')
xlim(0:1)
ylabel('Probability Density')
xlabel('Energy Consumption/ kWh')

%Plot for Mosaic I at 4am
subplot (2,2,3)
histogram(inRange3,100,'Normalization','pdf')
title ('Mosaic I at 4am')
xlim(0:1)
ylabel('Probability Density')
xlabel('Energy Consumption/ kWh')
hold on
ln_p3 = lognfit(I_4am_ec);
plot (xt,lognpdf(xt,ln_p3(1),ln_p3(2)),'r')
hold off

%Plot for Mosaic I at 10am
subplot (2,2,4)
histogram(inRange4,100,'Normalization','pdf')
title ('Mosaic I at 10am')
xlim(0:1)
ylabel('Probability Density')
xlabel('Energy Consumption/ kWh')
hold on
ln_p4 = lognfit(I_10am_ec);
plot (xt,lognpdf(xt,ln_p4(1),ln_p4(2)),'r')
hold off

%Mosaic I at 4am
%Calculate the cumulative probability greater than 1kWh for log normal
%graph
yt = 1-logncdf(1,ln_p3(1),ln_p3(2));
%Calculate number of customers greater than 1kWh
No_cust_4am = yt*I_4am_cn
%Calculate actual number of customers greater than 1kWh from original data
No_cust_act_4am = sum(I_4am_ec>1)

%Mosaic I at 10am
%Calculate the cumulative probability greater than 1kWh for log normal
%graph
yt2 = 1-logncdf(1,ln_p4(1),ln_p4(2));
%Calculate number of customers greater than 1kWh
No_cust_10am = yt2*I_10am_cn
%Calculate actual number of customers greater than 1kWh from original data
No_cust_act_10am = sum(I_10am_ec>1)