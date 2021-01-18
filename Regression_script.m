% HK building cost regression

% load data
load 'HK_build_cost.mat'

% Data into matrices: rc_dat, steel_dat

% Data:
% Col 1: Average Floor area
% Col 2: Total floor area
% Col 3: Story height
% Col 4: Adjusted construction cost (OUTPUT to be estimated)

% Redefine variable names to make them more user-friendly
rc_ave_fl = rc_dat(:,1);
rc_tot_fl = rc_dat(:,2);
rc_height = rc_dat(:,3);
rc_cost = rc_dat(:,4);
steel_ave_fl = steel_dat(:,1);
steel_tot_fl = steel_dat(:,2);
steel_height = steel_dat(:,3);
steel_cost = steel_dat(:,4);

% Plot data in figure
figure('Name','Raw Data Plots','units','normalized','outerposition',[0 0 1 1])

% col 1 Average Floor Area
subplot (4,4,1)
histogram(rc_ave_fl, 'FaceColor','b','EdgeColor','k','FaceAlpha',1)
xlabel('Average Floor Area/m^2')
xlim([0 4000])
hold on 
histogram(steel_ave_fl,'FaceColor','r','EdgeColor','k','FaceAlpha',0.6)
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,5)
plot(rc_ave_fl,rc_tot_fl,'x')
xlabel('Average Floor Area')
ylabel('Total Floor Area')
yticks(0:50000:200000)
hold on
plot(steel_ave_fl,steel_tot_fl,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,9)
plot(rc_ave_fl,rc_height,'x')
xlabel('Average Floor Area')
ylabel('Story Height')
ylim([3.5 5])
hold on
plot(steel_ave_fl,steel_height,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,13)
plot(rc_ave_fl,rc_cost,'x')
xlabel('Average Floor Area')
ylabel('Cost')
ylim ([0 3000])
hold on
plot(steel_ave_fl,steel_cost,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

% col 2 Total Floor Area
subplot (4,4,2)
plot(rc_tot_fl, rc_ave_fl, 'x')
xlabel('Total Floor Area')
ylabel('Average Floor Area')
yticks(0:1000:4000)
hold on
plot(steel_tot_fl,steel_ave_fl,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,6)
histogram(rc_tot_fl, 'FaceColor','b','EdgeColor','k','FaceAlpha',1)
xlabel('Total Floor Area')
ylim([0 8])
yticks(0:2:8)
hold on
histogram(steel_tot_fl,'FaceColor','r','EdgeColor','k','FaceAlpha',0.6)
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,10)
plot(rc_tot_fl, rc_height, 'x')
xlabel('Total Floor Area')
ylabel('Story Height')
ylim([3.5 5])
hold on
plot(steel_tot_fl,steel_height,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,14)
plot(rc_tot_fl, rc_cost, 'x')
xlabel('Total Floor Area')
ylabel('Cost')
ylim ([0 3000])
hold on
plot(steel_tot_fl,steel_cost,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

% col 3 story height
subplot(4,4,3)
plot(rc_height, rc_ave_fl, 'x')
xlabel('Story Height')
ylabel('Average Floor Area')
yticks(0:1000:4000)
hold on
plot(steel_height,steel_ave_fl,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,7)
plot(rc_height, rc_tot_fl, 'x')
xlabel('Story Height')
ylabel('Total Floor Area')
yticks(0:50000:200000)
hold on
plot(steel_height,steel_tot_fl,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,11)
histogram(rc_height, 'FaceColor','b','EdgeColor','k','FaceAlpha',1)
xlabel('Story Height')
xlim([3.5 5])
ylim([0 15])
hold on
histogram(steel_height,'FaceColor','r','EdgeColor','k','FaceAlpha',0.6)
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,15)
plot(rc_height, rc_cost, 'x')
xlabel('Story Height')
ylabel('Cost')
ylim ([0 3000])
hold on
plot(steel_height,steel_cost,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

% Col 4: adjusted construction cost
subplot(4,4,4)
plot(rc_cost, rc_ave_fl, 'x')
xlabel('Cost')
ylabel('Average Floor Area')
yticks(0:1000:4000)
hold on
plot(steel_cost,steel_ave_fl,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,8)
plot(rc_cost, rc_tot_fl, 'x')
xlabel('Cost')
ylabel('Total Floor Area')
yticks(0:50000:200000)
hold on
plot(steel_cost,steel_tot_fl,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,12)
plot(rc_cost, rc_height, 'x')
xlabel('Cost')
ylabel('Story Height')
ylim([3.5 5])
hold on
plot(steel_cost,steel_height,'x')
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

subplot(4,4,16)
histogram(rc_cost, 'FaceColor','b','EdgeColor','k','FaceAlpha',1)
xlabel('Cost')
xlim([0 3000])
hold on
histogram(steel_cost,'FaceColor','r','EdgeColor','k','FaceAlpha',0.6)
hold off
legend({'Reinf. Concr.','Steel'},'Location','best','FontSize',6)

%Calculate the R^2 and R^2 adjusted values for the reinforced
%concrete and steel data with the response variable being the cost and the
%other variables being the factors in the regression equation.

% Explanatory variables for regression line of reinforced concrete data
x1 = rc_ave_fl;
x2 = rc_tot_fl;
x3 = rc_height;
% Response variable
y_rc = rc_cost;

% Put all the explanatory variables (including a constant term) into one
% matrix for each of the possible combinations of factors; average floor
% area, total floor area and story height
X = [ones(size(x1)) x1 x2 x3]; %all three factors
X1= [ones(size(x1)) x1 x2]; %total and average floor area
X2= [ones(size(x1)) x1 x3]; %average floor area and story height
X3= [ones(size(x1)) x2 x3]; %total floor area and story height
X4= [ones(size(x1)) x1]; %only average floor area
X5= [ones(size(x1)) x2]; %only total floor area
X6= [ones(size(x1)) x3]; %only story height

% Calculate regression model manually
% Estimate the parameters for each model individually using regression
% function
b = regress(y_rc,X);
b1 = regress(y_rc,X1);
b2 = regress(y_rc,X2);
b3 = regress(y_rc,X3);
b4 = regress(y_rc,X4);
b5 = regress(y_rc,X5);
b6 = regress(y_rc,X6);

% Use these results to calculate the mean y values for each model 
y_hat = X * b;
y_hat1 = X1 * b1;
y_hat2 = X2 * b2;
y_hat3 = X3 * b3;
y_hat4 = X4 * b4;
y_hat5 = X5 * b5;
y_hat6 = X6 * b6;

% Calculate the residuals for each model
e = y_rc - y_hat;
e1 = y_rc - y_hat1;
e2 = y_rc - y_hat2;
e3 = y_rc - y_hat3;
e4 = y_rc - y_hat4;
e5 = y_rc - y_hat5;
e6 = y_rc - y_hat6;

% Compute residual sum of squares for each model
sse = sum(e.^2);
sse1 = sum(e1.^2);
sse2 = sum(e2.^2);
sse3 = sum(e3.^2);
sse4 = sum(e4.^2);
sse5 = sum(e5.^2);
sse6 = sum(e6.^2);

% Compute the total sum of squares for the model
v = y_rc - mean(y_rc);
sst = sum(v.*v);

%Calculate biased R^2 values
R_squared = 1-(sse/sst);
R_squared1 = 1-(sse1/sst);
R_squared2 = 1-(sse2/sst);
R_squared3 = 1-(sse3/sst);
R_squared4 = 1-(sse4/sst);
R_squared5 = 1-(sse5/sst);
R_squared6 = 1-(sse6/sst);

% Calculate the number of data points
n = size(X,1);

%Calculate number of regressors (deleting the const term)
p = size(X,2)-1;
p1 = size(X1,2)-1; 
p2 = size(X2,2)-1;
p3 = size(X3,2)-1;
p4 = size(X4,2)-1;
p5 = size(X5,2)-1;
p6 = size(X6,2)-1;

%Calculate the R^2 adj for each model
R2_adj = 1-((sse/(n-p)) / (sst/(n-1)));
R2_adj1 = 1-((sse1/(n-p1)) / (sst/(n-1)));
R2_adj2 = 1-((sse2/(n-p2)) / (sst/(n-1)));
R2_adj3 = 1-((sse3/(n-p3)) / (sst/(n-1)));
R2_adj4 = 1-((sse4/(n-p4)) / (sst/(n-1)));
R2_adj5 = 1-((sse5/(n-p5)) / (sst/(n-1)));
R2_adj6 = 1-((sse6/(n-p6)) / (sst/(n-1)));

% Explanatory variables for regression line of steel data with the same
% analysis to follow. Different nomenclature used for regression line:
z1 = steel_ave_fl;
z2 = steel_tot_fl;
z3 = steel_height;
% Response variable
y_steel = steel_cost;

% Put all the explanatory variables (including a constant term) into one matrix
Z = [ones(size(z1)) z1 z2 z3]; %all three factors
Z1= [ones(size(z1)) z1 z2]; %total and average floor area
Z2= [ones(size(z1)) z1 z3]; %average floor area and story height
Z3= [ones(size(z1)) z2 z3]; %total floor area and story height
Z4= [ones(size(z1)) z1]; %only average floor area
Z5= [ones(size(z1)) z2]; %only total floor area
Z6= [ones(size(z1)) z3]; %only story height

% Estimate the parameters 
bz = regress(y_steel,Z);
bz1 = regress(y_steel,Z1);
bz2 = regress(y_steel,Z2);
bz3 = regress(y_steel,Z3);
bz4 = regress(y_steel,Z4);
bz5 = regress(y_steel,Z5);
bz6 = regress(y_steel,Z6);

% Use these results to calculate the mean y values for each model 
y_hatz = Z * bz;
y_hatz1 = Z1 * bz1;
y_hatz2 = Z2 * bz2;
y_hatz3 = Z3 * bz3;
y_hatz4 = Z4 * bz4;
y_hatz5 = Z5 * bz5;
y_hatz6 = Z6 * bz6;

% Residuals
ez = y_steel - y_hatz;
ez1 = y_steel - y_hatz1;
ez2 = y_steel - y_hatz2;
ez3 = y_steel - y_hatz3;
ez4 = y_steel - y_hatz4;
ez5 = y_steel - y_hatz5;
ez6 = y_steel - y_hatz6;

% Sum of squares
ssez = sum(ez.^2);
ssez1 = sum(ez1.^2);
ssez2 = sum(ez2.^2);
ssez3 = sum(ez3.^2);
ssez4 = sum(ez4.^2);
ssez5 = sum(ez5.^2);
ssez6 = sum(ez6.^2);

% Compute the total sum of squares for the model
vz = y_steel - mean(y_steel);
sst1 = sum(vz.*vz);

%Calculate biased R^2 values
R_squaredz = 1-(ssez/sst1);
R_squaredz1 = 1-(ssez1/sst1);
R_squaredz2 = 1-(ssez2/sst1);
R_squaredz3 = 1-(ssez3/sst1);
R_squaredz4 = 1-(ssez4/sst1);
R_squaredz5 = 1-(ssez5/sst1);
R_squaredz6 = 1-(ssez6/sst1);

% Number of data points
nz = size(Z,1);

% Number of regressors (deleting the const term)
pz = size(Z,2)-1; 
pz1 = size(Z1,2)-1; 
pz2 = size(Z2,2)-1; 
pz3 = size(Z3,2)-1; 
pz4 = size(Z4,2)-1; 
pz5 = size(Z5,2)-1; 
pz6 = size(Z6,2)-1; 

%Calculate the R^2 adj for each model
R2_adjz = 1-((ssez/(nz-pz)) / (sst1/(nz-1)));
R2_adjz1 = 1-((ssez1/(nz-pz1)) / (sst1/(nz-1)));
R2_adjz2 = 1-((ssez2/(nz-pz2)) / (sst1/(nz-1)));
R2_adjz3 = 1-((ssez3/(nz-pz3)) / (sst1/(nz-1)));
R2_adjz4 = 1-((ssez4/(nz-pz4)) / (sst1/(nz-1)));
R2_adjz5 = 1-((ssez5/(nz-pz5)) / (sst1/(nz-1)));
R2_adjz6 = 1-((ssez6/(nz-pz6)) / (sst1/(nz-1)));

%Create matrices with R^2 and R^2 adjusted results and table headings 
RegressorVariables = {'Ave fl area, Tot fl area, Story ht';'Ave fl area, Tot fl area';'Ave fl area and story height';'Ave fl area, Story height';'Ave fl area';'Tot fl area';'Story ht'};
ReinfConcr_R2 = round([R_squared;R_squared1;R_squared2;R_squared3;R_squared4;R_squared5;R_squared6],4,'significant');
ReinfConcr_R2_adj = round([R2_adj;R2_adj1;R2_adj2;R2_adj3;R2_adj4;R2_adj5;R2_adj6],4,'significant');
Steel_R2 = round([R_squaredz;R_squaredz1;R_squaredz2;R_squaredz3;R_squaredz4;R_squaredz5;R_squaredz6],4,'significant');
Steel_R2_adj = round([R2_adjz;R2_adjz1;R2_adjz2;R2_adjz3;R2_adjz4;R2_adjz5;R2_adjz6],4,'significant');

%Print table for user to identify results
T = table(ReinfConcr_R2,ReinfConcr_R2_adj,Steel_R2,Steel_R2_adj,'RowNames',RegressorVariables)

% Plot residuals of the best and worst models against cost in order to evaluate the models
figure('Name','Residual Plots ','units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
plot(rc_cost, e, 'x')
title ('Best Reinforced Concrete Model')
xlabel('Cost/ million HK$')
ylabel({'Residuals for model';'with all variables'},'Fontsize',14)

subplot(2,2,2)
plot(rc_cost, e6, 'x')
title ('Worst Reinforced Concrete Model')
xlabel('Cost/ million HK$')
ylabel({'Residuals for model';'with story height'},'Fontsize',14)

subplot(2,2,3)
plot(steel_cost, ez, 'x')
title ('Best Steel Model')
xlabel('Cost/ million HK$')
ylabel({'Residuals for model';'with all variables'},'Fontsize',14)

subplot(2,2,4)
plot(steel_cost, ez6, 'x')
title ('Worst Steel Model')
xlabel('Cost/ million HK$')
ylabel({'Residuals for model with';'story height'},'Fontsize',14)
