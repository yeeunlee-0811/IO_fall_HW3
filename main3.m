clear

%% data loading

data3 = importdata('ascending_data.dat');
nauction = size(data3,1);
nbidders = data3(:,1);
winprice = data3(:,2);

nbidmin = min(nbidders); %3
nbidmax = max(nbidders); %5
wpmin = min(winprice); %22
wpmax = max(winprice); %146

%% price grid
pgrid = [22:1:146]';
npgrid = size(pgrid,1);

%% gprice
gprice = zeros(npgrid,3);

for i = 1:3
    ni = i+2;
    for j = 1:npgrid
        vprice = pgrid(j,1);
        
        idni = nbidders == ni;
        idprice = winprice < vprice;
        
        gprice(j,i) = (1/sum(idni))*sum(idni.*idprice);
    end
end

%% CDF
cdf = zeros(npgrid,1);

resphi = zeros(nauction,npgrid);

for i =1:npgrid
    for t = 1:nauction
        nt = nbidders(t,1);
        
        g = gprice(i,nt-2);
        
        solve_phi = @(x) phi(x,g,2,nt);
        
        x0 = 0.5;
        options = optimoptions('fsolve','Display','iter',...
            'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-12,...
            'StepTolerance',1e-12,'FunctionTolerance',1e-12,...
            'MaxIterations',25);
        
        xres = fsolve(solve_phi,x0,options);
        resphi(t,i) = xres;
    end
    
    cdf(i,1) = (1/nauction)*sum(resphi(:,i));
end

%%
cdf3 =  plot(pgrid,cdf);
saveas(cdf3,'cdf3.png')

