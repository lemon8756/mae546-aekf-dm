% DM_AEKF_KF.m

% Adaptive Extended Kalman Filter and Kalman Filter Comparison

% Created on May 9, 2015, 12:25 by Aaron James Lemmer

% close all;

%% Initialize real deformable mirror
D = 1;
nxy = 64; 
Nact = 16;  % number of actuators across the DM
dxy = D/nxy;
x = (-(nxy-1)/2:(nxy-1)/2)*dxy;
[X,Y] = meshgrid(x,x);
var = (4*D/Nact/nxy)/(2*sqrt(2*log(2)));
DMbasmat = zeros(nxy,nxy,Nact^2);
DMbasvec = zeros(nxy^2,Nact^2);

paramtru = zeros(Nact);

actlocidx = linspace(-0.5,0.5,Nact);
for actx = 1:Nact
    for acty = 1:Nact
        p = 0 + 0.3*randn(1);
        paramtru(actx,acty) = p;
        thisact = exp(-((X-actlocidx(actx)).^2/(2*var*(1 + p)) + (Y-actlocidx(acty)).^2/(2*var*(1 + p))));
        DMbasmat(:,:,(actx-1)*Nact+acty) = thisact;
        DMbasvec(:,(actx-1)*Nact+acty) = thisact(1:nxy^2);
    end
end

DMall = zeros(nxy);
for index = 1:Nact^2
    DMall = DMall + DMbasmat(:,:,index);
end
figure(10); imagesc(x,x,DMall);
colorbar; axis square; axis xy;

figure(21);
imagesc(1:Nact,1:Nact,paramtru);
axis square; axis xy; colorbar;

% %% Decompose Kolmogorov into actuator strengths
% z = Kolm(1:nxy^2)';   %reindex Kolmogorov screen to be [nxy^2 x 1]
% 
% %z = DMbasvec*a, so do the left pseudoinverse to find 
% comvec = (DMbasvec'*DMbasvec)\(DMbasvec'*z);
% 
% A = reshape(comvec,[Nact,Nact]);  % check it by plotting
% figure(3);
% imagesc(1:Nact,1:Nact,A);
% axis square; axis xy; colorbar;

% %% Actuate DM given actuator strengths
% DMactuated = DMbasvec*comvec;
% 
% seeDMact = reshape(DMactuated,[nxy,nxy]);
% figure(7);
% imagesc(xkolm,xkolm,seeDMact,[-1,1]);
% axis square; axis xy; colorbar;

%% Initialize model deformable mirror
varxmod = (4*D/Nact/nxy)/(2*sqrt(2*log(2)));
varymod = varxmod;
DMbasmodmat = zeros(nxy,nxy,Nact^2);
DMbasmodvec = zeros(nxy^2,Nact^2);
DMbasmoddelmat = zeros(nxy,nxy,Nact^2);
DMbasmoddelvec = zeros(nxy^2,Nact^2);

for actx = 1:Nact
    for acty = 1:Nact
        thisact = exp(-((X-actlocidx(actx)).^2/(2*varxmod) + (Y-actlocidx(acty)).^2/(2*varymod)));
        DMbasmodmat(:,:,(actx-1)*Nact+acty) = thisact;
        DMbasmodvec(:,(actx-1)*Nact+acty) = thisact(1:nxy^2);
        
        thisactdelta = thisact.*((X-actlocidx(actx)).^2 + (Y-actlocidx(acty)).^2)/(2*varxmod^2);
        DMbasmoddelmat(:,:,(actx-1)*Nact+acty) = thisactdelta;
        DMbasmoddelvec(:,(actx-1)*Nact+acty) = thisactdelta(1:nxy^2);
    end
end

% DMall = zeros(nxy);
% for index = 1:Nact^2
%     DMall = DMall + DMbasmodmat(:,:,index);
% end
% figure(18); imagesc(x,x,DMbasmodmat(:,:,34));
% colorbar; axis square; axis xy;

%% Initialize KF control loop
numit = 10;  % set number of loop iterations
contvec = zeros(numit,2);

% Set DM actuator commands to zero
KFcomvec = zeros(Nact^2,1);
KFcomtot = zeros(Nact^2,1);

% Set initial DM shape (unactuated)
KFDMshapvec = DMbasvec*KFcomvec;

% Generate and plot initial Kolmogorov phase screen
r0 = D/0.008;
[Kolm,~,xkolm] = kolmogorov(D,nxy,r0);

% % Propagate initial Kolmogorov phase screen to image plane
% [imout,xout] = DM2implane_prop(Kolm,nxy);

% Initialize the state estimate
KFphsest = zeros(nxy^2,1);  % initial value is zero

% Initialize the covariance estimate
actvar = 0.015^2;  % variance in stroke for a single actuator
KFcovest = zeros(nxy^2);

% Populate R with the sensor variances for each pixel
senvar = 0.005^2;  % variance in phase measurement per pixel
KFR = senvar*eye(nxy^2);

%% Initialize AEKF control loop
% Set DM actuator commands to zero
comvec = zeros(Nact^2,1);
comtot = zeros(Nact^2,1);

% Set initial DM shape (unactuated)
DMshapvec = DMbasvec*comvec;

% Initialize the state estimate
phsest = zeros(nxy^2,1);  % initial value is zero
paramest = zeros(Nact^2,1);
phsaugest = [phsest;paramest];

% Initialize the augmented observation matrix
obsaugmat = [eye(nxy^2) zeros(nxy^2,Nact^2)];

% Initialize the covariance estimate
paramvar = 0.02^2;
covest = zeros(nxy^2+Nact^2);

% Populate R with the sensor variances for each pixel
R = senvar*eye(nxy^2);


%% Control loop
for iteration = 1:numit
    % 0. Calculate new Kolmogorov phase additive for input
    [newKolm,~,~] = kolmogorov(D,nxy,r0);
%     addKolm = newKolm*0.5*exp(-(iteration-1)/10);
    addKolm = 0.05*newKolm;
    Kolm = Kolm + addKolm;
%     Kolm = newKolm;
    Kolmvec = Kolm(1:nxy^2)';  % reindex Kolm to be [nxy^2 x 1]
    
    figure(1);
    imagesc(xkolm,xkolm,Kolm,[-0.1,0.1]);
    xlabel('[m]');
    ylabel('[m]');
    colorbar;
    axis square;
    axis xy;
    
    % 1. Calculate state (phase error between Kolomogorov input and current DM shape)
    phserr = Kolmvec - DMshapvec;
    phserrmat = reshape(phserr,[nxy,nxy]);
    
    KFphserr = Kolmvec - KFDMshapvec;
    KFphserrmat = reshape(KFphserr,[nxy,nxy]);
    
    figure(2);
    imagesc(xkolm,xkolm,phserrmat,[-0.02,0.02]);
    title('True phase error (AEKF)');
    axis xy; axis square; colorbar;
    
%     figure(52);
%     imagesc(xkolm,xkolm,KFphserrmat,[-0.02,0.02]);
%     title('True phase error (KF)');
%     axis xy; axis square; colorbar;
    
    % 2. Measure the state (phase error)
    phsmeas = phserr + sqrt(senvar)*randn(size(phserr));
    phsmeasmat = reshape(phsmeas, [nxy,nxy]);
    
    KFphsmeas = KFphserr + sqrt(senvar)*randn(size(KFphserr));
    KFphsmeasmat = reshape(KFphsmeas, [nxy,nxy]);
    
    %     figure(11);
    %     imagesc(xkolm,xkolm,phsmeasmat);
    %     title('Measured phase error');
    %     axis xy; axis square; colorbar;
    
    %     figure(15);  % with my 'sensor' this is just a picture of white noise
    %     imagesc(xkolm,xkolm,phsmeasmat-phserrmat);
    %     title('Measurement error, z-x');
    %     axis xy; axis square; colorbar;
    
    % 3. Estimate the state (phase error) with a Kalman filter
    % Extrapolate the state estimate
    phsaugest = blkdiag(eye(nxy^2),eye(Nact^2))*phsaugest + [DMbasmodvec;zeros(Nact^2)]*comvec + [DMbasmoddelvec;zeros(Nact^2)]*diag(paramest)*comvec;
    phsest = phsaugest(1:nxy^2);
    phsestmat = reshape(phsest,[nxy,nxy]);
    
    figure(12);
    imagesc(xkolm,xkolm,phsestmat);
    title('Estimated phase error, xhat(-), AEKF');
    axis xy; axis square; colorbar;
    
    figure(13);
    imagesc(xkolm,xkolm,phserrmat-phsestmat);
    title('State (phase error) residual, x-xhat(-), AEKF');
    axis xy; axis square; colorbar;
    
    figure(14);
    imagesc(xkolm,xkolm,phsmeasmat-phsestmat);
    title('Measurement residual, z-xhat(-), AEKF');
    axis xy; axis square; colorbar;
    
        %KF
    KFphsest = KFphsest + DMbasmodvec*KFcomvec;
    KFphsestmat = reshape(KFphsest,[nxy,nxy]);
    
%     figure(512);
%     imagesc(xkolm,xkolm,KFphsestmat);
%     title('Estimated phase error, xhat(-), KF');
%     axis xy; axis square; colorbar;
%     
%     figure(513);
%     imagesc(xkolm,xkolm,KFphserrmat-KFphsestmat);
%     title('State (phase error) residual, x-xhat(-), KF');
%     axis xy; axis square; colorbar;
%     
%     figure(514);
%     imagesc(xkolm,xkolm,KFphsmeasmat-KFphsestmat);
%     title('Measurement residual, z-xhat(-), KF');
%     axis xy; axis square; colorbar;
    
    % Calculate the Jacobian
    jacob = [eye(nxy^2) DMbasmoddelvec*diag(comvec); zeros(Nact^2, nxy^2), eye(Nact^2, Nact^2)];
    
    % Extrapolate the covariance estimate
    covest = jacob*covest*jacob' + blkdiag(actvar*(DMbasmodvec*DMbasmodvec'),paramvar*eye(Nact^2));
    
    KFcovest = KFcovest + actvar*(DMbasmodvec*DMbasmodvec');
    
    % Compute the filter gains
    gainK = covest*obsaugmat'/(obsaugmat*covest*obsaugmat' + R);
    
    KFgainK = KFcovest/(KFcovest + KFR);  % KF
    
    % Update the state estimate with a measurement
    phsaugest = phsaugest + gainK*(phsmeas - obsaugmat*phsaugest);
    phsest = phsaugest(1:nxy^2);
    paramest = phsaugest((nxy^2+1):(nxy^2+Nact^2));
    
    phsestmat = reshape(phsest,[nxy,nxy]);
    
    figure(15);
    imagesc(xkolm,xkolm,phsestmat);
    title('Updated estimate of phase error, xhat(+), AEKF');
    axis xy; axis square; colorbar;
    
    figure(16);
    imagesc(xkolm,xkolm,phserrmat-phsestmat);
    title('Updated state (phase error) residual, x-xhat(+), AEKF');
    axis xy; axis square; colorbar;
    
    figure(17);
    imagesc(xkolm,xkolm,phsmeasmat-phsestmat);
    title('Updated masurement residual, z-xhat(+), AEKF');
    axis xy; axis square; colorbar;
    
    paramestmat = reshape(paramest,[Nact,Nact]);
    figure(20);
    imagesc(1:Nact,1:Nact,paramestmat);
    axis square; axis xy; colorbar;
    
        % KF
    KFphsest = KFphsest + KFgainK*(KFphsmeas - KFphsest);
    KFphsestmat = reshape(KFphsest,[nxy,nxy]);
    
%     figure(515);
%     imagesc(xkolm,xkolm,KFphsestmat);
%     title('Updated estimate of phase error, xhat(+), KF');
%     axis xy; axis square; colorbar;
%     
%     figure(516);
%     imagesc(xkolm,xkolm,KFphserrmat-KFphsestmat);
%     title('Updated state (phase error) residual, x-xhat(+), KF');
%     axis xy; axis square; colorbar;
%     
%     figure(517);
%     imagesc(xkolm,xkolm,KFphsmeasmat-KFphsestmat);
%     title('Updated masurement residual, z-xhat(+), KF');
%     axis xy; axis square; colorbar;
    
    % Update the covariance estimate
%     covest = inv(inv(covest) + inv(R));
    covest = (eye(nxy^2 + Nact^2) - gainK*obsaugmat)*covest;
    
    KFcovest = (eye(nxy^2) - KFgainK)*KFcovest;
    
    % 3. Decompose estimated phase error into new differential DM commands
    comvec = (DMbasmodvec'*DMbasmodvec)\(DMbasmodvec'*phsest);  % left pseudoinverse to find comvec, since errvec = DMbasvec*comvec
    comtot = comtot + comvec;
    
    commat = reshape(comtot,[Nact,Nact]);  % check DM commands by plotting
    figure(5);
    imagesc(1:Nact,1:Nact,commat,[-0.06,0.06]);
    title('AEKF DM commands');
    axis square; axis xy; colorbar;
    
    KFcomvec = (DMbasmodvec'*DMbasmodvec)\(DMbasmodvec'*KFphsest);  % left pseudoinverse to find comvec, since errvec = DMbasvec*comvec
    KFcomtot = KFcomtot + KFcomvec;
    
%     KFcommat = reshape(KFcomtot,[Nact,Nact]);  % check DM commands by plotting
%     figure(55);
%     imagesc(1:Nact,1:Nact,KFcommat,[-0.06,0.06]);
%     title('KF DM commands');
%     axis square; axis xy; colorbar;
    
    % 4. Update the true DM shape with new actuator commands
    DMshapvec = DMshapvec + DMbasvec*comvec;
    
    DMshapmat = reshape(DMshapvec,[nxy,nxy]);  % check DM shape by plotting
    figure(6);
    imagesc(xkolm,xkolm,DMshapmat,[-0.1,0.1]);
    title('Applied DM shape, AEKF');
    axis square; axis xy; colorbar;
    
    KFDMshapvec = KFDMshapvec + DMbasvec*KFcomvec;
    
%     KFDMshapmat = reshape(KFDMshapvec,[nxy,nxy]);  % check DM shape by plotting
%     figure(56);
%     imagesc(xkolm,xkolm,KFDMshapmat,[-0.1,0.1]);
%     title('Applied DM shape, KF');
%     axis square; axis xy; colorbar;
    
    % 5. Propagate to the image plane and calculate contrast
    [imout,xout,contrast] = DM2implane_prop(phserrmat,nxy);
    contvec(iteration,1) = log10(contrast);

    figure(4)
    imagesc(xout,xout,log10(imout),[-7,-3]);
    colorbar;
    axis square;
    axis xy;
    title('AEKF Image Plane Intensity, Iout');
    xlabel('Horizontal Dimension [m]');
    ylabel('Vertical Dimension [m]');
    
    [KFimout,KFxout,KFcontrast] = DM2implane_prop(KFphserrmat,nxy);
    contvec(iteration,2) = log10(KFcontrast);

%     figure(54)
%     imagesc(linspace(-1,1,length(KFxout)),linspace(-1,1,length(KFxout)),log10(KFimout),[-7,-3]);
%     colorbar;
%     axis square;
%     axis xy;
%     title('KF Image Plane Intensity, Iout');

    pause(0.005);
end

%%
figure(8);
plot(1:numit,contvec(:,1), 1:numit,contvec(:,2));
xlabel('Iterations');
ylabel('Log10 of Contrast');
legend('AEKF','KF');

%     [Kolm, RMSKolm,~] = kolmogorov(D,nxy,r0);
% %     additive = Kolm*0.5*exp(-(iteration-1)/10);
% %     input = input + additive;
% %     
% %     figure(1);
% %     imagesc(x,x,input);
% %     xlabel('[m]');
% %     ylabel('[m]');
% %     colorbar;
% %     axis square;
% %     axis xy;
% %     
% %     pause(0.01);