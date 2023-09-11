% let , f'=g;
% f''=h;
% f'''=-f*f-B*(1-g^2);


eta_max=8; %recommended in notes
N=1000; %let it be 1000
delta_eta=eta_max/(N-1); %from notes
B= -0.2; %pressure gradient parametre beta
l = 10; %length of the plate
mu=8.90*10^-4 ; %coefficient of viscosity
u_inf=20; %free stream velocity guess value

eta_values=0:delta_eta:eta_max; %for plotting
f = zeros(size(eta_values));
g = zeros(size(eta_values));
h = zeros(size(eta_values));

%boundary conditions
f(1)=0;
g(1)=0;
h(1)=-B; %beta

%itration1
for i=1:N-1
    eta = eta_values(i);

    k1_f = g(i);
    k1_g = h(i);
    k1_h = -f(i) * h(i) - B * (1 - g(i)^2);
    
    k2_f = g(i) + 0.5 * delta_eta * k1_h;
    k2_g = h(i) + 0.5 * delta_eta * k1_h;
    k2_h = - (f(i) + 0.5 * delta_eta * k1_f) * (h(i) + 0.5 * delta_eta * k1_h) - B * (1 - (g(i) + 0.5 * delta_eta * k1_g)^2);
    
    k3_f = g(i) + 0.5 * delta_eta * k2_h;
    k3_g = h(i) + 0.5 * delta_eta * k2_h;
    k3_h = - (f(i) + 0.5 * delta_eta * k2_f) * (h(i) + 0.5 * delta_eta * k2_h) - B * (1 - (g(i) + 0.5 * delta_eta * k2_g)^2);
    
    k4_f = g(i) + delta_eta * k3_h;
    k4_g = h(i) + delta_eta * k3_h;
    k4_h = - (f(i) + delta_eta * k3_f) * (h(i) + delta_eta * k3_h) - B * (1 - (g(i) + delta_eta * k3_g)^2);
    
    % Update values
    f(i+1) = f(i) + (delta_eta / 6) * (k1_f + 2*k2_f + 2*k3_f + k4_f);
    g(i+1) = g(i) + (delta_eta / 6) * (k1_g + 2*k2_g + 2*k3_g + k4_g);
    h(i+1) = h(i) + (delta_eta / 6) * (k1_h + 2*k2_h + 2*k3_h + k4_h);
end


P=g(N)-1;
epslon=10^-6;


% Plot f, f', and f'' as functions of η
figure;
subplot(3,1,1);
plot(eta_values, f,"LineWidth",2);
title('f as a function of η');
subplot(3,1,2);
plot(eta_values, g,"LineWidth",2);
title('f'' as a function of η');
subplot(3,1,3);
plot(eta_values, h,"LineWidth",2);
title('f'''' as a function of η');
xlabel('η');


%displacement thickness
function dt = disp_thic(x,g)
    d=1-g;
    c_delta = simpsons_integral(d,1,12,13);
    dt = x*c_delta/sqrt(reynolds(x));
return
end

%momentum thickness
function mt = momt_thic(x,g)
    d=g(1-g);
    c_delta = simpsons_integral(d,1,12,13);
    mt = x*c_delta/sqrt(reynolds(x));
return
end

%shape factor
function sf = shape_fac(x,g)
    d=1/g;
    c_delta = simpsons_integral(d,1,12,13);
    sf = x*c_delta/sqrt(reynolds(x));
return
end

%skin friction
function sk = skin_fric(x,g)
    sk = 2*g(1)/sqrt(reynolds(x));
return
end

%simpsons rule of integration
function in = simpsons_integral(y,low,high,steps)
    h = (low+high)/steps;
    % Initialize the sum for Simpson's rule
    integral_sum = y(low) + y(high); % Include the endpoints
    
    % Sum for even indices
    for i = 2:2:(steps-1)
        x_i = low + i * h;
        integral_sum = integral_sum + 4 * y(x_i);
    end
    
    % Sum for odd indices
    for i = 1:2:(steps-2)
        x_i = low + i * h;
        integral_sum = integral_sum + 2 * y(x_i);
    end

    in = (h/3)*integral_sum;
return
end
%reynolds number
function re = reynolds(x)
    if(x<=l)
        re = u_inf*x/mu;
    else
        error("x is greater than lenth of plate");
    end
return
end


