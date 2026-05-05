%% Copyright (c) Adam Ouzeri 2022
%% Calculate the different tension contribution

ViscoelasticEnergy = zeros(length(time),1);
ActiveEnergy = zeros(length(time),1);
SpringpotEnergy = zeros(length(time),1);

for i = 1:length(time)

    gamma_ap = rho_ap(i)*xiap;
    gamma_lat_p = rho_lat_p(i)*xilat;
    gamma_lat_t = rho_lat_t(i)*xilat;
    
    %% metric components
    Gsharp_ap = [MetricTensorSharp(i,1), 0; 0, MetricTensorSharp(i,2)];
    Gsharp_lat_p = [MetricTensorSharp(i,3), 0; 0, MetricTensorSharp(i,4)];
    Gsharp_lat_t = [MetricTensorSharp(i,5), 0; 0, MetricTensorSharp(i,6)];
        
    Gflat_ap = Gsharp_ap\eye(2);
    Gflat_lat_p = Gsharp_lat_p\eye(2);
    Gflat_lat_t = Gsharp_lat_t\eye(2);
    
    %% kinematics
    J_ap = kappa(i);
    J_lat_p = 1;
    J_lat_t = sqrt(kappa(i)^2 + 3)/2/kappa(i);
    
    Cflat_ap = [kappa(i)^2, 0; 0, 1];
    Cflat_lat_p = [kappa(i)^2, 0; 0, 1/kappa(i)^2];
    Cflat_lat_t = [(kappa(i)^2 + 3)/4, 0; 0, 1/kappa(i)^2];
    
    C_inverse_sharp_ap = Cflat_ap\eye(2);
    C_inverse_sharp_lat_p = Cflat_lat_p\eye(2);
    C_inverse_sharp_lat_t = Cflat_lat_t\eye(2);
    
    Cmixed_ap = Cflat_ap*Gsharp_ap;
    Cmixed_lat_p  = Cflat_lat_p *Gsharp_lat_p;
    Cmixed_lat_t  = Cflat_lat_t *Gsharp_lat_t;
    
    I1_ap = trace(Cmixed_ap);
    I3_ap = sqrt(det(Cmixed_ap));
    I1_lat_p = trace(Cmixed_lat_p);
    I3_lat_p = sqrt(det(Cmixed_lat_p));
    I1_lat_t = trace(Cmixed_lat_t);
    I3_lat_t = sqrt(det(Cmixed_lat_t));
    
    dCflat_ap_dkappa = [2*kappa(i),0;0,0];
    dCflat_lat_p_dkappa = [2*kappa(i),0;0,-2/kappa(i)^3];
    dCflat_lat_t_dkappa = [kappa(i)/2,0;0,-2/kappa(i)^3];
    
    %% free energy
%     Psi_ap = mu*(I1_ap/I3_ap -2) + lambda*(I3_ap -1)^2;
%     Psi_lat_p = mu*(I1_lat_p/I3_lat_p - 2) + lambda*(I3_lat_p - 1)^2;
%     Psi_lat_t = mu*(I1_lat_t/I3_lat_t -2) + lambda*(I3_lat_t - 1)^2;
    
    %% variation of free energy with respect to Gsharpdot
    dPsi_ap_dI1_ap = mu/I3_ap;
    dPsi_ap_dI3_ap = 2*lambda*(I3_ap - 1) - mu*I1_ap/I3_ap^2;
    dPsi_lat_p_dI1_lat_p = mu/I3_lat_p;
    dPsi_lat_p_dI3_lat_p = 2*lambda*(I3_lat_p - 1) - mu*I1_lat_p/I3_lat_p^2;
    dPsi_lat_t_dI1_lat_t = mu/I3_lat_t;
    dPsi_lat_t_dI3_lat_t = 2*lambda*(I3_lat_t - 1) - mu*I1_lat_t/I3_lat_t^2;
    
    dI1_ap_dGsharp_ap = Cflat_ap;
    dI1_lat_p_dGsharp_lat_p = Cflat_lat_p;
    dI1_lat_t_dGsharp_lat_t = Cflat_lat_t;
    
    dI3_ap_dGsharp_ap = I3_ap*Gflat_ap/2;
    dI3_lat_p_dGsharp_lat_p = I3_lat_p*Gflat_lat_p/2;
    dI3_lat_t_dGsharp_lat_t = I3_lat_t*Gflat_lat_t/2;
    
    dPsi_ap_dGsharp_ap = dPsi_ap_dI1_ap*dI1_ap_dGsharp_ap + dPsi_ap_dI3_ap*dI3_ap_dGsharp_ap;
    dPsi_lat_p_dGsharp_lat_p = dPsi_lat_p_dI1_lat_p*dI1_lat_p_dGsharp_lat_p + dPsi_lat_p_dI3_lat_p*dI3_lat_p_dGsharp_lat_p;
    dPsi_lat_t_dGsharp_lat_t = dPsi_lat_t_dI1_lat_t*dI1_lat_t_dGsharp_lat_t + dPsi_lat_t_dI3_lat_t*dI3_lat_t_dGsharp_lat_t;
    
    
    %% variation of free energy with respect to kappadot
    dI1_ap_dCflat_ap = Gsharp_ap;
    dI1_lat_p_dCflat_lat_p = Gsharp_lat_p;
    dI1_lat_t_dCflat_lat_t = Gsharp_lat_t;
    
    dI3_ap_dCflat_ap = I3_ap*C_inverse_sharp_ap/2;
    dI3_lat_p_dCflat_lat_p = I3_lat_p*C_inverse_sharp_lat_p/2;
    dI3_lat_t_dCflat_lat_t = I3_lat_t*C_inverse_sharp_lat_t/2;
    
    dPsi_ap_dCflat_ap = dPsi_ap_dI1_ap*dI1_ap_dCflat_ap + dPsi_ap_dI3_ap*dI3_ap_dCflat_ap;
    dPsi_lat_p_dCflat_lat_p = dPsi_lat_p_dI1_lat_p*dI1_lat_p_dCflat_lat_p + dPsi_lat_p_dI3_lat_p*dI3_lat_p_dCflat_lat_p;
    dPsi_lat_t_dCflat_lat_t = dPsi_lat_t_dI1_lat_t*dI1_lat_t_dCflat_lat_t + dPsi_lat_t_dI3_lat_t*dI3_lat_t_dCflat_lat_t;
    
    E_ap = rho_ap(i)*trace(dPsi_ap_dCflat_ap'*dCflat_ap_dkappa)*J_ap*Aap0;
    E_lat_p = rho_lat_p(i)*trace(dPsi_lat_p_dCflat_lat_p'*dCflat_lat_p_dkappa)*J_lat_p*A_lat_p0;
    E_lat_t = rho_lat_t(i)*trace(dPsi_lat_t_dCflat_lat_t'*dCflat_lat_t_dkappa)*J_lat_t*A_lat_t0;
    
    ViscoelasticEnergy(i) = 2*E_ap + E_lat_p + E_lat_t;
    
    %% variation of active power input with respect to kappadot
    dJdot_ap_dkappadot = 1;
    dJdot_lat_p_dkappadot = 0;
    dJdot_lat_t_dkappadot = -3/(2*kappa(i)^2*sqrt(kappa(i)^2 + 3));
    
    dPa_ap_dkappadot = gamma_ap*dJdot_ap_dkappadot*Aap0;
    dPa_lat_p_dkappadot = gamma_lat_p*dJdot_lat_p_dkappadot*A_lat_p0;
    dPa_lat_t_dkappadot = gamma_lat_t*dJdot_lat_t_dkappadot*A_lat_t0;
    
    ActiveEnergy(i) = 2*dPa_ap_dkappadot + dPa_lat_p_dkappadot + dPa_lat_t_dkappadot;
    
end

%% Fractional element
for i = 1:length(time) - 1
    passedF = (time(i) + params.deltat - time(1:i)').^(-params.alphaPower) * (kappadot(1:i)./kappa(1:i));
    SpringpotEnergy(i + 1) = (passedF + params.deltat^(-params.alphaPower)*kappadot(i + 1)/kappa(i + 1))*2*params.materialConstant/gamma(1-params.alphaPower)/kappa(i + 1)*params.deltat*params.V0;
end

ViscoelasticTension = ViscoelasticEnergy/Aap0;
ActiveTension = ActiveEnergy/Aap0;
SpringpotTension = SpringpotEnergy/Aap0;
