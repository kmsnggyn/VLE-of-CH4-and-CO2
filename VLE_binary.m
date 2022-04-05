clear all, clc, close all

load('datamix.mat')
set(groot,'defaultLineMarkerSize', 10, ...
    'defaultLineLineWidth', 2, ...
    'defaultAxesFontName', 'Times');

global R

R = 0.0831446261815324; % L bar / K mol
T_data = [230, 250, 270]; % K

% CH4: species 1
Tc1 = 190.564;  % K
Pc1 = 45.992;   % bar
Dc1 = 10.139;   % mol/L
Vc1 = 1/Dc1;    % L/mol
 w1 = 0.01142;

% CO2: species 2
Tc2 = 304.21;
Pc2 = 73.829955;
Dc2 = 10.6249;
Vc2 = 1/Dc2;
 w2 = 0.22394;

%% Choice of an EOS and BIP (Binary Interaction Parameter)

EOS = input('\nChoose an EOS\n(1:vdW, 2:SRK, 3:PR)\n>> ');    % Choose which EOS you want to use
EOSname = {'vdW' 'SRK' 'PR'};
EOSname = EOSname{EOS};
P_initial = [15 20 30]; % Initial values for faster calculation

fprintf('\nChoose a BIP model (0 to 4)\n');
fprintf(' 0: k12 = 0 (Classical Mixing Rule)\n')
fprintf(' 1: k12 = 0.919\n');
fprintf(' 2: k12 = 0.979\n');
fprintf(' 3: k12 = 0.5219 - 0.8254*Tr + 0.4494*Tr^3\n');

% this work
fprintf(' 4: k12 = 0.7834 - 0.9783493*Tr + 0.35693431*Tr^2\n');

BIP_model_num = input('>> ');

fprintf(['\n< ' EOSname ' EOS with BIP model no.' num2str(BIP_model_num) ' >\n'])

%% Calculate for three temperatures: 230K, 250K, 270K


for Ti = 1:3
    T = T_data(Ti);
    Tr1 = T/Tc1;
    Tr2 = T/Tc2;
    Tr = (Tr1+Tr2)/2;
    fprintf('Calculating VLE for %dK...\n',T)
    

    if EOS == 1
        [Omega, Ksi, eps, sig, alpha1, alpha2] = vdW;
    elseif EOS == 2
        [Omega, Ksi, eps, sig, alpha1, alpha2] = SRK(Tr1,Tr2,w1,w2);
    elseif EOS == 3
        [Omega, Ksi, eps, sig, alpha1, alpha2] = PR(Tr1,Tr2,w1,w2);
    end

    k12 = BIP(BIP_model_num,Tr);

    a1 = Ksi*alpha1*R^2*Tc1^2/Pc1;
    a2 = Ksi*alpha2*R^2*Tc2^2/Pc2;
    a12 = (1-k12)*a1^(1/2)*a2*(1/2);
    
    b1 = Omega*R*Tc1/Pc1;
    b2 = Omega*R*Tc2/Pc2;
    
    % x1 = 0.0001, 0.0002, 0.0003, ...
    n = 10^4;

    for k = 1:n
    
        x1 = k/n;
        x2 = 1 - x1;
        
        if k == 1
            P_new(1) = P_initial(Ti);
            y1_new(1) = 0.1;
        else
            P_new(1) = P_res(k-1);
            y1_new(1) = y1_res(k-1);
        end
            
        power = - 2;
        
        for j = 1:10000
        
            P = P_new(j);
            y1 = y1_new(j);
            
            parameters = [a1, a2, a12, b1, b2, eps, sig];
            [K(j), y1_new(j+1)] = K_converge(T, x1, P, y1, parameters);

            if j > 1
                if ( K(j) -1 )*( K(j-1) -1 ) < 0
                    power = power - 1;
                end
            end
        
            if abs( K(j) - 1 ) < 10^(-10)
                P_res(k) = P_new(j);
                y1_res(k) = y1_new(j+1);
                break
            elseif K(j) > 1
                P_new(j + 1) = P + 10^power;
            elseif K(j) < 1
                P_new(j + 1) = P - 10^power;
            end
    
            if j == 10000
                fprintf('K did not converge to 1\n')
            end
        end
        
        if k > 1        % Stop the procedure if P decreases with increasing x1
            if P_res(k) < P_res(k-1)
                break
            end
        end
    
    end

    x1 = linspace(1/n, 1, n);
    lim = k-2;
    x1 = x1(1:lim);
    
    VLE = figure('Position',[0 10000 500 500]);
    plot(dataset{Ti}(:,1),dataset{Ti}(:,3),'.'); hold on
    plot(dataset{Ti}(:,2),dataset{Ti}(:,3),'.');
    plot(x1,P_res(1:lim)); hold on
    plot(y1_res(1:lim),P_res(1:lim));
    axis([0 1 0 100])
    pbaspect([1 1 1]) 
    legend('Exp. BUBL','Exp. DEW',[EOSname ' BUBL'],[EOSname ' DEW'],'Location','southeast')
    xlabel('$x_{\mathrm{CH_4}}$, $y_{\mathrm{CH_4}}$','Interpreter','latex')
    ylabel('$P$ ${\mathrm {[bar]}}$','Interpreter','latex')
    exportgraphics(gca,['mixVLE_', EOSname, '_BIP_', num2str(BIP_model_num), '_', num2str(T), 'K.jpg'],'Resolution',300)

end

fprintf('\nComplete!\n\n')
%% Functions

function k12 = BIP(num, Tr)

    if num == 0
        k12 = 0;
    elseif num == 1
        k12 = 0.919;
    elseif num == 2
        k12 = 0.979;
    elseif num == 3
        k12 = 0.5219 - 0.8254*Tr + 0.4494*Tr^3;
    elseif num == 4
        k12 = 0.7834 - 0.9783493*Tr + 0.35693431*Tr^2;
    end

end

function [Omega, Ksi, eps, sig, alpha1, alpha2] = vdW
    eps = 0;
    sig = 0;
    alpha1 = 1;
    alpha2 = 1;
    Omega = 1/8;
    Ksi = 27/64;
end

function [Omega, Ksi, eps, sig, alpha1, alpha2] = SRK(Tr1,Tr2,w1,w2)
    eps = 0;
    sig = 1;
    alpha1 = ( 1 + (0.480 + 1.574*w1 - 0.176*w1^2) * (1 - Tr1^(0.5)) )^2;
    alpha2 = ( 1 + (0.480 + 1.574*w2 - 0.176*w2^2) * (1 - Tr2^(0.5)) )^2;
    Omega = 0.08664;
    Ksi = 0.42748;
end

function [Omega, Ksi, eps, sig, alpha1, alpha2] = PR(Tr1,Tr2,w1,w2)
    eps = 1-sqrt(2);
    sig = 1+sqrt(2);
    alpha1 = (1+(0.37464+1.57226*w1-0.26992*w1^2)*(1-Tr1^(1/2)))^2;
    alpha2 = (1+(0.37464+1.57226*w2-0.26992*w2^2)*(1-Tr2^(1/2)))^2;
    Omega = 0.07780;
    Ksi = 0.45724;
end

function I = I_calc(eps, sig, beta, Z)
    if sig == eps
        I = beta/(Z+eps*beta);
    else
        I = (1/(sig - eps)) * (log((Z + sig*beta)/(Z + eps*beta)));
    end
end

function Z_out = Zv_CEOS(eps,sig,beta,q)
    
    Z = roots([1 ...
        ( (eps+sig)*beta - 1 - beta ) ...
        ( eps*sig*beta^2 - (eps+sig)*beta - (eps+sig)*beta^2 + q*beta) ...
        ( - eps*sig*beta^2 - eps*sig*beta^3 - q*beta^2) ]);

    Zi = Z==real(Z);
    Z_real = Z(Zi);
    Z_out = max(Z_real);

end

function Z_out = Zl_CEOS(eps,sig,beta,q)
    
    Z = roots([1 ...
        ( (eps+sig)*beta - 1 - beta ) ...
        ( eps*sig*beta^2 - (eps+sig)*beta - (eps+sig)*beta^2 + q*beta) ...
        ( - eps*sig*beta^2 - eps*sig*beta^3 - q*beta^2) ]);

    Zi = Z==real(Z);
    Z_real = Z(Zi);
    Z_out = min(Z_real);

end

function [K_res, y1_res] = K_converge(T, x1, P, y1_in, parameters)
    global R

    y1_new(1) = y1_in;
    x2 = 1 - x1;

    a1  = parameters(1);
    a2  = parameters(2);
    a12 = parameters(3);
    b1  = parameters(4);
    b2  = parameters(5);
    eps = parameters(6);
    sig = parameters(7);

    for i = 1:100
    
        y1 = y1_new(i);
        y2 = 1 - y1;
                
        a_l = x1^2*a1 + 2*x1*x2*a12 + x2^2*a2;
        a_v = y1^2*a1 + 2*y1*y2*a12 + y2^2*a2;
        b_l = x1*b1 + x2*b2;
        b_v = y1*b1 + y2*b2;
        beta_l = b_l*P/(R*T);
        beta_v = b_v*P/(R*T);
        q_l = a_l/(b_l*R*T);
        q_v = a_v/(b_v*R*T);
        qbar1_l = q_l*( (2*x1*a1 + 2*x2*a12)/a_l - b1/b_l );
        qbar2_l = q_l*( (2*x2*a2 + 2*x1*a12)/a_l - b2/b_l );
        qbar1_v = q_v*( (2*y1*a1 + 2*y2*a12)/a_v - b1/b_v );
        qbar2_v = q_v*( (2*y2*a2 + 2*y1*a12)/a_v - b2/b_v );
        
        Z_l = Zl_CEOS(eps, sig, beta_l, q_l);
        Z_v = Zv_CEOS(eps, sig, beta_v, q_v);
        
        I_l = I_calc(eps, sig, beta_l, Z_l);
        I_v = I_calc(eps, sig, beta_v, Z_v);
        
        phi1_l = exp( (b1/b_l)*(Z_l - 1) - log(Z_l - beta_l) - qbar1_l*I_l);
        phi2_l = exp( (b2/b_l)*(Z_l - 1) - log(Z_l - beta_l) - qbar2_l*I_l);
        
        phi1_v = exp( (b1/b_v)*(Z_v - 1) - log(Z_v - beta_v) - qbar1_v*I_v);
        phi2_v = exp( (b2/b_v)*(Z_v - 1) - log(Z_v - beta_v) - qbar2_v*I_v);
        
        K1 = phi1_l/phi1_v;
        K2 = phi2_l/phi2_v;
        
        K(i) = K1*x1 + K2*x2;
        y1_new(i + 1) = (K1*x1) / (K1*x1 + K2*x2);

        if i > 1
            if abs ( (K(i)-K(i-1)) / K(i-1)) < 10^(-10)
                break
            end
        end
    end

    K_res = K(i);
    y1_res = y1_new(i);
   
end


