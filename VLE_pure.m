clear all, clc
close all

load('datapure.mat')
set(groot,'defaultLineMarkerSize', 10, ...
    'defaultLineLineWidth', 2, ...
    'defaultAxesFontName', 'Times');

global R
R = 0.0831446261815324; % L bar / K mol

% CH4: species 1
Tc1 = 190.564;  % K
Pc1 = 45.992;   % bar
Dc1 = 10.139;   % mol/L
Vc1 = 1/Dc1;    % L/mol
 w1 = 0.01142;
% AC1 = [6.61184 389.9278 265.99];

% CO2: species 2
Tc2 = 304.1282;
Pc2 = 73.773;
Dc2 = 10.6249;
Vc2 = 1/Dc2;
 w2 = 0.22394;
AC2 = [7.5788, 863.35, 273.15];

% FIG1 = figure('Position',[0 10000 500 500]);

% EOS = input('vdW : 1\n SRK : 2\nPR : 3\n\n');
for EOS = 1:3
    
    EOSname = {'vdW' 'SRK' 'PR'};
    EOSname = EOSname{EOS};
    fprintf(['\n< ' EOSname ' EOS >\n'])
    plotsettings = {'interpreter','latex','fontsize',14};
    
    T_lim1 = [120 200];
    T_lim2 = [220 320];
    n = 10^3;
    
    [rho1_v, rho1_l, T1] = VLE(T_lim1(1), T_lim1(2), n, 10, Pc1, Tc1, w1, EOS);
    
    VLE1 = figure('Position',[0 10000 400 400]);
    plot(CH4_vap(:,3),CH4_vap(:,1),'.'); hold on
    plot(CH4_liq(:,3),CH4_liq(:,1),'.');
    plot(rho1_v,T1,'-'); hold on
    plot(rho1_l,T1,'-');
    axis([0 25 T_lim1(1) T_lim1(2)])
    pbaspect([1 1 1])
    xlabel('$\rho \mathrm{ [mol/l]}$',plotsettings{:})
    ylabel('$T \mathrm{ [K]}$',plotsettings{:})
    legend('Exp. DEW' ,'Exp. BUBL',[EOSname ' DEW'], [EOSname ' BUBL'])
    exportgraphics(gca,['VLE_1_' EOSname '.jpg'],'Resolution',300)
    
    [rho2_v, rho2_l, T2] = VLE(T_lim2(1), T_lim2(2), n, 20, Pc2, Tc2, w2, EOS);
    
    VLE2 = figure('Position',[0 0 400 400]);
    plot(CO2_vap(:,3),CO2_vap(:,1),'.'); hold on
    plot(CO2_liq(:,3),CO2_liq(:,1),'.');
    plot(rho2_v,T2,'-'); hold on
    plot(rho2_l,T2,'-');
    axis([0 30 T_lim2(1) T_lim2(2)])
    pbaspect([1 1 1])
    xlabel('$\rho \mathrm{ [mol/l]}$',plotsettings{:})
    ylabel('$T \mathrm{ [K]}$',plotsettings{:})
    legend('Exp. DEW' ,'Exp. BUBL',[EOSname ' DEW'], [EOSname ' BUBL'])
    exportgraphics(gca,['VLE_2_' EOSname '.jpg'],'Resolution',300)

end



function [rho_v, rho_l, T] = VLE(Tmin, Tmax, n, P_in, Pc, Tc, w, EOS)

    global R
    
    for i = 1:n
        T(i) = Tmin + (Tmax - Tmin)*(i/n);
        Tr = T(i)/Tc;
        
%         fprintf('\n Calculating for T = %.4fK >',T(i))
        
        if EOS == 1
            [Omega, Ksi, eps, sig, alpha] = vdW;
        elseif EOS == 2
            [Omega, Ksi, eps, sig, alpha] = SRK(Tr,w);
        elseif EOS == 3
            [Omega, Ksi, eps, sig, alpha] = PR(Tr,w);
        end
        
        a = Ksi*alpha*R^2*Tc^2/Pc;
        b = Omega*R*Tc/Pc;
            
        [rho_v(i), rho_l(i), P_sat(i), stop] = Psat(eps, sig, a, b, T(i), P_in);
        P_in = P_sat(i);
%         fprintf(' %7f, %7f, %7f', rho_v(i), rho_l(i), P_sat(i))
        
        if stop == 1
            break
        end
    end
    
    rho_v = rho_v(1:i-1);
    rho_l = rho_l(1:i-1);
    T = T(1:i-1);

end

function [rho_v, rho_l, Psat, stop] = Psat(eps, sig, a, b, T, P_in)
    global R
    P(1) = P_in;
    power = -3;
    stop = 0;
    for i = 1:10000
        
        beta = b*P(i)/(R*T);
        q = a/(b*R*T);
        
        [Z_l, Z_v] = Z_CEOS(eps, sig, beta, q);
        
        I_l = I_calc(eps, sig, beta, Z_l);
        I_v = I_calc(eps, sig, beta, Z_v);
        
        phi_l = exp( Z_l - 1 - log(Z_l - beta) - q*I_l);
        phi_v = exp( Z_v - 1 - log(Z_v - beta) - q*I_v);
        
        K(i) = phi_l / phi_v;
    
        if abs( K(i)-1 ) < 10^(-15)
            if Z_l == Z_v
                stop = 1;
                break
            end
            break
        elseif K(i) > 1
            P(i+1) = P(i) + 10^power;
        elseif K(i) < 1
            P(i+1) = P(i) - 10^power;
        end
        
        if i > 2
            if ( K(i) - 1 )*( K(i-1) - 1 ) < 0
                power = power - 1;
            end
        end

    end

    Psat = P(i);
    rho_l = Psat/(Z_l*R*T);
    rho_v = Psat/(Z_v*R*T);

end

function [Omega, Ksi, eps, sig, alpha] = vdW
    eps = 0;
    sig = 0;
    alpha = 1;
    Omega = 1/8;
    Ksi = 27/64;
end

function [Omega, Ksi, eps, sig, alpha] = SRK(Tr,w)
    eps = 0;
    sig = 1;
    alpha = ( 1 + (0.480 + 1.574*w - 0.176*w^2) * (1 - Tr^(0.5)) )^2;
%     alpha0 = Tr^(-0.201158) * exp(0.141599*(1-Tr^2.29528));
%     alpha1 = Tr^(-0.660145) * exp(0.500315*(1-Tr^2.63165));
%     alpha = alpha0 + w*(alpha1 - alpha0);
    Omega = 0.08664;
    Ksi = 0.42748;
end

function [Omega, Ksi, eps, sig, alpha] = PR(Tr,w)
    eps = 1-sqrt(2);
    sig = 1+sqrt(2);
    alpha = (1+(0.37464+1.57226*w-0.26992*w^2)*(1-Tr^(1/2)))^2;
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

function [Zl, Zv] = Z_CEOS(eps,sig,beta,q)
    
    Z = roots([1 ...
        ( (eps+sig)*beta - 1 - beta ) ...
        ( eps*sig*beta^2 - (eps+sig)*beta - (eps+sig)*beta^2 + q*beta) ...
        ( - eps*sig*beta^2 - eps*sig*beta^3 - q*beta^2) ]);

    Zi = Z==real(Z);
    Z_real = Z(Zi);
    Zv = max(Z_real);
    Zl = min(Z_real);
end