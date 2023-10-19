function dPdt = eqed0(t,d,y,k,g)

dPdt = zeros( 3, 1 );
dPdt(1) = - 1i * ( 0 - 1i * y / 2 ) * d(1) - 1i * g * d(2) - 1i * g * d(3);
dPdt(2) = - 1i * ( 0 - 1i * k / 2 ) * d(2) - 1i * g * d(1);
dPdt(3) = - 1i * ( 0 - 1i * k / 2 ) * d(3) - 1i * g * d(1);
end