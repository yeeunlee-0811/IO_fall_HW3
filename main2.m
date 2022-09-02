clear

%% data loading

data = importdata('fpa.dat');

T = size(data,1);

%% Check the min and max
minmaxmat = zeros(2,4);

for i = 1:4
    indiv = data(:,i);
    
    minmaxmat(1,i) = min(indiv); % 0~
    minmaxmat(2,i) = max(indiv); % ~150
end

%% A
grid_indiv = [0:10:150]'; %16x1
ngrid_indiv = size(grid_indiv,2);
bw = 10;
grid = zeros(ngrid_indiv^4,4);
ngrid = size(grid,2);

pdfjoint = zeros(ngrid,1);

for i = 1:ngrid
    vgrid = grid(i,:);
    xvec = (data - vgrid.*ones(T,1))/bw;
    pdfvec = normpdf(xvec,0,1); % used Normal distn 
    pdfvec1 = prod(pdfvec,2); % product Kernel
    pdfjoint(i,1) = mean(pdfvec1); 
end

%% B: Marginal CDF

v = [0:0.5:150]';
nv = size(v,1);

F = zeros(nv,4);

for j = 1:4
    indiv = data(:,j);
    for i = 1:nv
        vval = v(i,1);
        iobs = indiv < vval;
        F(i,j) = (1/T)*sum(iobs);
    end
end

plot1 =  plot(v,F(:,1));
saveas(plot1,'plot1.png')
plot2 =  plot(v,F(:,2));
saveas(plot2,'plot2.png')
plot3 =  plot(v,F(:,3));
saveas(plot3,'plot3.png')
plot4 =  plot(v,F(:,4));
saveas(plot4,'plot4.png')
%%
figure;

subplot(2,2,1)
plot(v,F(:,1));
subplot(2,2,2)
plot(v,F(:,2));
subplot(2,2,3)
plot(v,F(:,3));
subplot(2,2,4)
plot(v,F(:,4));

%% C

%% quantmat
quantmat = zeros(2,4);

for j = 1:4
    indiv = data(:,j);
    
    quantmat(1,j) = quantile(indiv,0.25);
    quantmat(2,j) = quantile(indiv,0.75);
end

%% joint CDF

jcdfval = zeros(16,4);

jcdfval(1:8,1) = quantmat(1,1);
jcdfval(9:16,1) = quantmat(2,1);
jcdfval([1:4,9:12],2) = quantmat(1,2);
jcdfval([5:8,13:16],2) = quantmat(2,2);
jcdfval([1,2,5,6,9,10,13,14],3) = quantmat(1,3);
jcdfval([3,4,7,8,11,12,15,16],3) = quantmat(2,3);
jcdfval([1,3,5,7,9,11,12,13,15],4) = quantmat(1,4);
jcdfval([2,4,6,8,10,12,14,16],4) = quantmat(2,4);

jcdf = zeros(16,1);

for i = 1:16
    id = data < jcdfval(i,:);
    iid = sum(id,2) == 4;
    jcdf(i,1) = (1/T)*sum(iid);
end

%% Checking symmetry

meanF = mean(F,2);
diffF = F - meanF.*ones(1,4);
diffFsqt = diffF.^2;
mean_diffFsqt = mean(diffFsqt,2);
maxdiffFsqt = max(mean_diffFsqt); % [0.000726749999999999]

%% F: symmetry & independence

data1 = reshape(data,[500*4,1]);

% distn 
bw = 0.5;

pdfsym = zeros(nv,1);

for i = 1:nv
    x = v(i,1);
    id = data1 > x - (1/2)*bw & data1 < x +(1/2)*bw;
    pdfsym(i,1) = (1/(4*500))*sum(id);
end


% cdf
cdfsym = zeros(nv,1);

for i = 1:nv
    x = v(i,1);
    id = data1<= x;
    cdfsym(i,1) = (1/(T*4))*sum(id);
end

%%
figure;

subplot(2,3,1)
plot(v,F(:,1));
subplot(2,3,2)
plot(v,F(:,2));
subplot(2,3,3)
plot(v,F(:,3));
subplot(2,3,4)
plot(v,F(:,4));
subplot(2,3,5)
plot(v,cdfsym);
    
    
    

        




        



