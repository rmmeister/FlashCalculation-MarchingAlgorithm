%% INITIALIZE MATLAB
clear all;
clc;
close all;
format long;

%% DASHBOARD
gc = 32.174;
% 1. Entering Well Data
mixtureRate = 100000; % lb/day
pipeDiameter = 1.995/12; % ft
pipeLength = 4000; % ft
liquidViscosity = 17; % cp
gasViscosity = .013; % cp
surfaceTension = 30; % dynes/cm (or lbf/ft?)
Temperature = 130; % degrees Fahrenheit
e = 0; % ft
Area = pi*pipeDiameter^2/4; % ft2
gasGravity = .75;
API = 30;
oilGravity = 141.5/(API + 131.5);
waterDensity = 62.4; % PCF
oilDensity = waterDensity*oilGravity; % PCF
inputDeviation = input('Enter deviation angle: (0 to 90)\n-> ');
deviationAngle = degtorad(inputDeviation); % tubing inclination angle
% INCREMENTS
dL = 400; % ft
gridSize = pipeLength/dL + 1;
Pressure = zeros(gridSize, 1);
Pressure(gridSize) = 1000; % psia, pressure at depth
dP = zeros(gridSize,1);
dP(gridSize) = 10;
lengthVector = 0:dL:pipeLength;
Iteration = zeros(gridSize,1);

%% MARCHING ALGORITHM
fprintf('\ncreating pressure profile using the Marching Algorithm...\n');
for i = gridSize:-1:2
    criteria(i) = inf;
    while criteria(i) > 1e-6
        
        % 2. Guess the Pressure
        Pressure(i-1) = Pressure(i) - dP(i);
        
        % 3. Calculate Average P
        averagePressure(i) = (Pressure(i) + Pressure(i-1))/2;
        
        % 4. Determine Fluid Properties @ Pav & T 
        % the compositional approach using the data provided in assingment
        % 2 is used and composition in each step is calculated by flash
        % calculation
        [x, y, gasFraction, liquidFraction, gasMW, liquidMW, mixtureMW ] = ...
            calcFlash( averagePressure(i)/145, (Temperature-32)*5/9 + 273.15 );
        
        % Obtain phase mass flow-rates
        liquidRate = liquidFraction * liquidMW/mixtureMW  * mixtureRate; % lb/day of liquid
        gasRate = gasFraction * gasMW/mixtureMW * mixtureRate; % lb/day of gas
        
        % Only gas density is used in this function
        [ Bo, Rs, Pb, Z, gasDensity ] =  ...
            calcStanding(Temperature, averagePressure(i), gasGravity, API, NaN, 1);
        
        % 5. Calc. Pressure Gradient 
        % TWO-PHASE PROPERTIES
        % calculating average two-phase fluid properties
        qL = liquidRate/oilDensity/5.615; % bbl/day
        qG = gasRate/gasDensity/5.615; % bbl/day
        lambdaL(i) = qL/(qL+qG); % no-slip liquid hold-up
        lambdaG = 1 - lambdaL(i); % gas void fraction
        mixtureDensity = oilDensity *lambdaL(i) + gasDensity*lambdaG; % PCF
        liquidSuperficialVelocity = qL/Area*5.615/24/60/60; % ft/s
        gasSuperficialVelocity = qG/Area*5.615/24/60/60; % ft/s
        mixtureVelocity = gasSuperficialVelocity + liquidSuperficialVelocity; % ft/s
        mixtureViscosity = liquidViscosity*lambdaL(i) + gasViscosity*lambdaG;
        
        % % PRESSURE GRADIENT % % % % % % % % % % % % % % % % % % % % % % %
        % calculating RhoBar, Holdup, and Flow Regime using the predefined
        % Beggs & Brill function
        [rhoBar(i), flowRegime, holdUp(i)] = calcBB(lambdaL(i), mixtureVelocity, ...
            liquidSuperficialVelocity, pipeDiameter,...
            oilDensity, gasDensity, surfaceTension, deviationAngle, Pressure(i));
        FR{i} = flowRegime;
        
        % calculating Friction Factor using predefined function
        [f_tp(i), f_n(i)] = calcFrictionPressure(mixtureDensity, mixtureVelocity, ...
            mixtureViscosity, pipeDiameter, lambdaL(i), holdUp(i));
        
        % now pressure gradients are derived based on rho_bar and friction
        % factor
        frictionalGradient(i) = (f_tp(i) * mixtureDensity * mixtureVelocity^2)/2/gc/pipeDiameter;
        hydrostaticGradient(i) = rhoBar(i)*sin(deviationAngle);
        pressureGradient = (hydrostaticGradient(i) + frictionalGradient(i)); % PCF
        
        % 6. Determining Error Criteria
        dpGuess = Pressure(i)-Pressure(i-1); % psi
        dP(i) = pressureGradient*dL/144; % psi
        criteria(i) = abs(dP(i) - dpGuess);
        Iteration(i) = Iteration(i) + 1;
    end
end

%% RESULTS
figure
plot(Pressure, -lengthVector, '+r' ,'LineWidth', 2)
grid on
title('Tubing Pressure Profile');
ylabel('Depth, ft');
xlabel('Pressure, psi');

figure
I = imread('beggs.png');
if strcmp(FR{gridSize}, 'Segregated') == 1
    I(170:270, 170:270, 1) = I(170:270, 170:270, 1) - 100;
elseif strcmp(FR{gridSize}, 'Intermittent')
    I(105:155, 340:440, 1) = I(105:155, 340:440, 1) - 100;
elseif strcmp(FR{gridSize}, 'Transition') == 1
    I(265:315, 365:465, 1) = I(265:315, 365:465, 1) - 120;
elseif strcmp(FR{gridSize}, 'Distributed') == 1
    I(20:70, 220:320, 1) = I(20:70, 220:320, 1) - 100;
end
imshow(I);

% figure
% plot(holdUp, -lengthVector, '+k', 'LineWidth', 2);
% grid on
% xlabel('Liquid Hold-up');
% ylabel('Depth, ft');
% title('Liquid Hold-up vs. Depth');

[xData, yData] = prepareCurveData( Pressure, lengthVector' );
ft = 'linearinterp';
[h_fit, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
pipeLength = abs(h_fit(1000) - h_fit(500));
fprintf('\nThe length of pipe between points P = 1000 and P = 500 is: %g ft\n', pipeLength);