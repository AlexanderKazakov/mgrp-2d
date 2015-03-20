% Prudnikov, Brychkov, Marichev p.46 f.25, p.47 f.30, p.37 f.19

clear all

syms A B C cosAl sinAl cosAlp sinAlp ln arctg Two z
cosAl = -B/(2*sqrt(A*C));
sinAl = sqrt(1 - cosAl*cosAl);
sinAlp = sqrt( (1-cosAl)/2 );
cosAlp = sqrt( (1+cosAl)/2 );
ln = log( ( z*z + 2*z*(A/C)^(1/4)*cosAlp + sqrt(A/C) ) / ( z*z - 2*z*(A/C)^(1/4)*cosAlp + sqrt(A/C) ) );
arctg = atan( ( z*z - sqrt(A/C) ) / ( 2*z*sqrt(A/C)*sinAlp ) );
Two = 1/(3*A*sinAl) * ( sinAlp * ln + 2*cosAlp * arctg );
%pretty(Two)

% syms der err
% der = 1/( A*z^4 + B*z^2 + C );
% err = der - diff(Two, z);
% err = subs(err, A, 1);
% err = subs(err, B, 1);
% err = subs(err, C, 1);
% zArr = 0:0.5:100;
% Err = subs(err, z, zArr);
% plot(zArr, Err);

syms a b c Three
a = 1;
b = 2*(sqrt(C/4/A) - B/4/A);
c = 2*sqrt(C/4/A);
Three = 1 / ( 4*A*sqrt(sqrt(C/4/A) - B/4/A) ) * ( 1/2/a*log(abs(( a*z*z - b*z + c )/( a*z*z + b*z + c ))) + b/( a*sqrt(4*a*c - b*b) )*atan( z*sqrt(4*a*c - b*b)/(c-a*z*z) ));
pretty(Three)

% syms der err
% der = z^2/( A*z^4 + B*z^2 + C );
% err = der - diff(Three, z);
% err = subs(err, A, 1);
% err = subs(err, B, 1);
% err = subs(err, C, 1);
% zArr = 0:0.5:100;
% Err = subs(err, z, zArr);
% plot(zArr, Err);

syms b y
Two = subs(Two, A, 1);
Two = subs(Two, B, 2*b);
Two = subs(Two, C, (b*b + y*y));
Three = subs(Three, A, 1);
Three = subs(Three, B, 2*b);
Three = subs(Three, C, (b*b + y*y));

syms x z a t ksi One J Ja Jma f 

One = 2*b*z - 2/3 * z^3 - 2*(b*y*y + b*b*b) * Two - 2*(b*b - y*y) * Three;
One = subs(One, b, x - a);
One = subs(One, z, sqrt(a - x + t));
One = subs(One, t, (x - ksi));

syms der err
der = (a - ksi)^(3/2)*(x - ksi) / ( (x-ksi)^2 + y*y );
err = der - diff(One, ksi);
err = subs(err, a, 1);
err = subs(err, x, 5);
err = subs(err, y, 5);
ksiArr = -1:0.1:1;
Err = subs(err, ksi, ksiArr);
plot(ksiArr, Err);

J = -2/3 * log( (x-ksi)^2 + y^2 ) * (a-ksi)^(3/2) - 4/3 * One;


% syms arctg
% arctg = atan((ksi - a + (1/(y^2 + (a - x)^2))^(1/2))/(2*(a - ksi)^(1/2)*(1/(y^2 + (a - x)^2))^(1/2)*(1/2 - (a - x)/(2*(y^2 + (a - x)^2)^(1/2)))^(1/2)))
% Ja = subs(J, arctg, pi/2)
% Ja = subs(Ja, ksi, a)
% Jma = subs(J, ksi, -a)
% f = Ja - Jma
% fx = diff(f, x)
