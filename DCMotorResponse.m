%PARAMETER SETTING

GUI = true; %Set false to avoid generatig plots
fig_cnt = 1;

%Dc motor parameters
J   = 6.76e-6;  %Kg*m^2
K_e = 38.9e-3;  %Kg*m^2*s^-2*A^-1
K_m = 38.9e-3;  %Kg*m^2*s^-2*A^-1
R   = 1.23;     %Kg*m^2*s^-3*A^-2
b   = 2e-6;     %Kg*m^2*s^-1
L   = 1.34e-3;  %Kg*m^2*s^-2*A^-2

%Time interval
t = 0:0.0001:0.1;
%t = 0:1:100;

%Laplace variable
s = tf([1 0],[1]);

%%
%COMPUTING TF

%tf parameters
w_n = sqrt( (R*b+K_m*K_e)/(J*L) );  
z   = (R*J+b*L)/(2*sqrt( (R*b+K_m*K_e)*J*L ));

G_0 = K_m/(R*b+K_m*K_e);

%Display results
disp("z = " + z + " - w_n = " + w_n);
disp("Gain (G(0)) = " + G_0);

%tf nominator and denominator
nomG = w_n^2;
denG = [1, 2*z*w_n, w_n^2];

%compute system (tf)
G = G_0*(1/s)*tf(nomG,denG);

%COMPUTE SYSTEM RESPONSES

%%
%IMPULSE RESPONSE

%Compute impulse response
plant_y_imp = impulse(G,t);

%Plot

if(GUI)
    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    plot(t,plant_y_imp, 'LineWidth', 3);

    %Plot h-line on gain
    hold on;

    plot(xlim,[G_0,G_0], 'LineWidth', 3);

    hold off;

    %Set plot caption
    title("Impulse response for DC motor");

    xlabel("time (s)");
    ylabel("displacement (rad)");

    legend("step response", "gain (G_0)");

    %Set plot style
    setPlotStyle();
end

%%
%STEP RESPONSE

%Compute step response
plant_y_step = step(G,t);

%Plot
if(GUI)
    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    plot(t,plant_y_step, 'LineWidth', 3);

    %Set plot caption
    title("Step response for DC motor");

    xlabel("time (s)");
    ylabel("displacement (rad)");

    %Set plot style
    setPlotStyle();
end

%%
%SYSTEM ROOT LOCUS

%Compute poles of L
delta   = -z*w_n;
w_d     = w_n*sqrt(z^2-1);

poles = [delta+(w_d), delta-(w_d)];

%Compue loci

if(GUI)
    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    rlocus(G);
end

%Get data
[Roots,Gains] = rlocus(G);

if(GUI)
    %plot poles
    hold on;

    if imag(poles(1)) == 0
        plot(poles(1),0, 'square', 'Color', 'r');
    else
        plot(poles(1), 'square', 'Color', 'r');
    end

    if imag(poles(2)) == 0
        plot(poles(2),0, 'square', 'Color', 'r');
    else
        plot(poles(2), 'square', 'Color', 'r');
    end

    root0 = plot(0,0,'square', rere'Color', 'r');

    hold off;


    %Set plot captions
    legend(root0, "poles of L(s)");

    %COMOPUTE TWO TANGENTS OF LOCI

    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    rlocus(G);
end

%Compute centroid
alpha = sum(Roots(:,1))/3;

%Compute angles
ang_1 = deg2rad( 180/3 );
ang_2 = deg2rad( (180-360)/3);

%Display results
disp("Centroid: " + alpha + " with angles: " ...
    + ang_1 + ", " + ang_2 + ", and 180 degrees");

%Compute tangents
x_lim = xlim;
x = alpha:1:x_lim(2);

line_1 = tan(ang_1)*(x-alpha);
line_2 = tan(ang_2)*(x-alpha);

%Plot centroid
if(GUI)
    hold on;

    alpha_plot = plot(alpha, 0, 'square');
    legend(alpha_plot, 'centroid');

    %Plot tangents
    line_1_plot = plot(x,line_1,'DisplayName','first tangent');
    line_2_plot = plot(x,line_2,'DisplayName','second tangent');
end

%Disp stability condition
disp("For the system to be stable 0<k_p<" + (2*w_n*z)/G_0 );

%%
%STABLE SINUSOIDAL RESPONSE

%Set new time interval
t_sin = 0:0.001:100;

%Set stable kp value
k_p = 20;

%Compute controller and closed loop system
C = k_p;
G_cl = (G*C)/(1+G*C);

%Plot sinusoidal response

if(GUI)
    figure(fig_cnt);
    fig_cnt = fig_cnt+1;

    plot(t_sin,lsim(G_cl,sin(t_sin),t_sin));

    %Set plot caption
    title({"(Marginally) stable sinusoidal response for DC motor," , ...
        "for C = K_p = " + k_p});

    xlabel("time (s)");
    ylabel("displacement (rad)");

    %Set plot style
    setPlotStyle();
end

%%
%UNSTABLE SINUSOIDAL RESPONSE

%Set new time interval
t_sin = 29:0.0001:30;

%Set stable kp value
k_p = 40;

%Compute controller and closed loop system
C = k_p;
G_cl = (G*C)/(1+G*C);

%Plot sinusoidal response

if(GUI)
    figure(fig_cnt);
    fig_cnt = fig_cnt+1;
    
    plot(t_sin,lsim(G_cl,sin(t_sin),t_sin));

    %Set plot caption
    title({"Unstable sinusoidal response for DC motor,", ...
        " for C = K_p = " + k_p});

    xlabel("time (s)");
    ylabel("displacement (rad)");

    %Set plot style
    setPlotStyle();
end

%%
%PI CONTROLLER FREQUENCY and SINUSOIDAL RESPONSE

%Set frequency interval
w =logspace(-3,3,600);

%Set parameters
k_p = (1/2)*(2*w_n*z)/G_0; %Average of boundaries
k_i_vals = -400:200:400;

cl_mags = [];
cl_phases = [];

cl_imp = [];
cl_step = [];

for k_i = k_i_vals
    %Define controller
    C = k_p + (1/s)*k_i;

    %compute system (tf)
    G_cl = (G*C)/(1+G*C);

    %Compute frequency response
    [cl_mag,cl_phase]=bode(G_cl,w);
    
    %Store results
    cl_mags    = [cl_mags squeeze(cl_mag)];
    cl_phases  = [cl_phases squeeze(cl_phase)];
    
    %Store impulse response
    cl_imp = [cl_imp impulse(G_cl, t)];
    
     %Store step response
    cl_step = [cl_step step(G_cl, t)];
    
end

if(GUI)
    %PLOT IMPULSE RESPONSES
    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    plot(t, cl_imp','LineWidth', 3);

    %Set plot caption
    title("Impulse response of closed loop system, w.r.t. K_i");

    ylabel("displacement (rad)");
    xlabel("time (s)");

    leg = legend(strsplit(num2str([k_i_vals]),' '));
    title(leg, 'K_i values');

    %Set plot style
    setPlotStyle();

    %PLOT IMPULSE RESPONSES
    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    plot(t, cl_step','LineWidth', 3);

    %Set plot caption
    title("Step response of closed loop system, w.r.t. K_i");

    ylabel("displacement (rad)");
    xlabel("time (s)");

    leg = legend(strsplit(num2str([k_i_vals]),' '));
    title(leg, 'K_i values');

    %Set plot style
    setPlotStyle();

    %PLOT MAGNITUDES

    %PLot results
    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    loglog(w,cl_mags','LineWidth', 3);

    %Set plot caption
    title("Magnitudes (in function of input frquency) w.r.t. K_i");

    xlabel("normalized input frequency");
    ylabel("magnitude");

    leg = legend(strsplit(num2str([k_i_vals]),' '));
    title(leg, 'K_i values');

    %Set plot style
    setPlotStyle();

    %PLOT PHASES

    %PLot results
    figure(fig_cnt);
    fig_cnt = fig_cnt + 1;

    semilogx(w,cl_phases','LineWidth', 3);

    %Set plot caption
    title("Phases (in function of input frquency) w.r.t. K_i");

    xlabel("normalized input frequency");
    ylabel("phase");

    leg = legend(strsplit(num2str([k_i_vals]),' '));
    title(leg, 'K_i values');

    %Set plot style
    setPlotStyle();
end

