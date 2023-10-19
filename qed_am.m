function dPdt = qed_am(t,d,dw,y,y0,k,g,p1,p2,N,s)

Pm = zeros( N + 3, N + 3 );

Pm( 1, 1 ) = - 1i * y0 / 2; 
Pm( 2, 2 ) = - 1i * k / 2; Pm( 3, 3 ) = Pm( 2, 2 );
Pm( 1, 2:3 ) = g; Pm ( 2:3, 1 ) = g;

for j = 4:N+3
    for l = 4:N+3
        if j == l
            Pm( j, l ) = dw - 1i * ( y + y0 / 2 );
        else
            Pm( j, l ) = - 1i * y * exp( 1i * abs( p2 - p1 ) * abs( j - l ) );
        end
    end
end

for j = 4:N+3
%     Pm( 3, j ) = - 1i * sqrt( k * y ) * exp( 1i * abs( p2 - p1 ) * ( j - 2 ) );
%     Pm( j, 2 ) = Pm( 3, j );
    if j == 4
        Pm( 3, j ) = - 1i * sqrt( k * y ) * exp( 1i * p1 );
        Pm( j, 2 ) = Pm( 3, j );
    end
    if j == 5
        Pm( 3, j ) = - 1i * sqrt( k * y ) * exp( 1i * p2 );
        Pm( j, 2 ) = Pm( 3, j );
    end
    if j >= 6
        Pm( 3, j ) = - 1i * sqrt( k * y ) * exp( 1i * ( p2 + ( p2 - p1 ) * ( j - 5 ) ) );
        Pm( j, 2 ) = Pm( 3, j );
    end
end

dPdt = zeros( N + 3, 1 );

% dPdt(1) = - 1i * ( Pm( 1, 1 ) * d(1) + Pm( 1, 2 ) * d(2) + Pm( 1, 3 ) * d(3) + Pm( 1, 4 ) * d(4) + Pm( 1, 5 ) * d(5) + Pm( 1, 6 ) * d(6) );
% dPdt(2) = - 1i * ( Pm( 2, 1 ) * d(1) + Pm( 2, 2 ) * d(2) + Pm( 2, 3 ) * d(3) + Pm( 2, 4 ) * d(4) + Pm( 2, 5 ) * d(5) + Pm( 2, 6 ) * d(6) );
% dPdt(3) = - 1i * ( Pm( 3, 1 ) * d(1) + Pm( 3, 2 ) * d(2) + Pm( 3, 3 ) * d(3) + Pm( 3, 4 ) * d(4) + Pm( 3, 5 ) * d(5) + Pm( 3, 6 ) * d(6) );
% dPdt(4) = - 1i * ( Pm( 4, 1 ) * d(1) + Pm( 4, 2 ) * d(2) + Pm( 4, 3 ) * d(3) + Pm( 4, 4 ) * d(4) + Pm( 4, 5 ) * d(5) + Pm( 4, 6 ) * d(6) );
% dPdt(5) = - 1i * ( Pm( 5, 1 ) * d(1) + Pm( 5, 2 ) * d(2) + Pm( 5, 3 ) * d(3) + Pm( 5, 4 ) * d(4) + Pm( 5, 5 ) * d(5) + Pm( 5, 6 ) * d(6) );
% dPdt(6) = - 1i * ( Pm( 6, 1 ) * d(1) + Pm( 6, 2 ) * d(2) + Pm( 6, 3 ) * d(3) + Pm( 6, 4 ) * d(4) + Pm( 6, 5 ) * d(5) + Pm( 6, 6 ) * d(6) );

for j = 1:N+3
    eval( s( j ) );
end

end




