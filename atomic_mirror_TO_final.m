clear all; clc
%%
% % please place this file and the files 'qed_am_TO.m', 'qed_am.m' and 'eqed0.m' in the same folder 
%%
p = 0.5 * pi;  % tunable parameter of dimerized interactions between atoms
p1 = 0;  % phase factor between the waveguide-cavity junction and the atoms in mirror
p2 = 1.5*pi; % phase factor between the atoms in mirror
dw = 0;   % atom-cavity detuning
y0 = 1;   % spontaneous emission rate of atom
y = 5*y0;  % waveguide-induced atom decay \Gamma
k = 20*y0;  % cavity decay
g = 20*y0;  % QE-cavity coupling
J0 = 0;  % interaction strength between the atoms in mirror
N = 3; % number of atoms in mirror
x0 = zeros( N + 3, 1 ); x0(1) = 1; % initial condition

tmax = 0.5;
tlist = linspace( 0, tmax, 10001 );

% % generate string equations
s = strings( N + 3, 1 );
for j = 1:N+3
    s( j ) = strcat( 'dPdt(', num2str( j ), ') = - 1i * ( ' );
    for l = 1:N+2
        s( j ) = strcat( s( j ), 'Pm( ', num2str( j ), ', ', num2str( l ), ' ) * d( ', num2str( l ), ' ) + ' );
    end
    s( j ) = strcat( s( j ), 'Pm( ', num2str( j ), ', ', num2str( l + 1 ), ' ) * d( ', num2str( l + 1 ), ' ) );' );
end

% % numerical solutions of string equations
[t,x]=ode45( @(t,x) qed_am_TO(t,x,dw,y,y0,k,g,p1,p2,p,J0,N,s), tlist, x0 );
[tm,xm]=ode45( @(t,x) qed_am(t,x,dw,y,y0,k,g,p1,pi,N,s), tlist, x0 );
[~,x0]=ode45( @eqed0, tlist, [ 1; 0; 0 ], [ ], y0, k, g );


% % final plot
figure(1); plot( t, abs(x(:,1)).^2, 'LineWidth', 2 ); hold on;
figure(1); plot( t, abs(x(:,2)).^2, 'LineWidth', 2 ); hold on;
figure(1); plot( t, abs(x(:,3)).^2, 'LineWidth', 2 ); hold on;
% figure(1); plot( t, exp(-y0*t), 'k', 'LineWidth', 2 ); hold on;
ylim( [ 0, 1 ] );
legend( 'QE', 'CCW mode', 'CW mode', 'Free space' );
xlabel( '{\it{\gamma_0}t}' );
ylabel( 'Occupation' );
title( '{\it{N}} = 3, J_0 = 0' );
set( gca, 'Fontname', 'Arial' );
set( gca, 'Fontsize', 14 );

figure(2); plot( t, abs( x( :, 4 ) ).^2, 'LineWidth', 2 ); hold on;
figure(2); plot( t, abs( x( :, 5 ) ).^2, 'LineWidth', 2 ); hold on;
figure(2); plot( t, abs( x( :, 6 ) ).^2, 'LineWidth', 2 ); hold on;
legend( 'Atom 1', 'Atom 2', 'Atom 3' );
title( '{\it{N}} = 3, J_0 = 0' );
ylabel( 'Population of QEs' );
xlabel( '{\it{\gamma_0}t}' );
set( gca, 'Fontname', 'Arial' );
set( gca, 'Fontsize', 14 );


