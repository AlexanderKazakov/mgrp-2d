clear all
syms x y a phi rho arctg ln A1 A2

arctg = atan( 2*sqrt(2*a)*rho*sin(phi) / (rho*rho - 2*a) );
ln = log(  ( 2*a - 2*sqrt(2*a)*rho*cos(phi) + rho*rho )  /  ( 2*a + 2*sqrt(2*a)*rho*cos(phi) + rho*rho )  );
A1 = 1/rho*( 1/2*( cos(phi) - x/y*sin(phi) ) * ln + ( sin(phi) + x/y*cos(phi) ) * arctg );
A2 = rho*( 1/2*( cos(phi) + x/y*sin(phi) ) * ln + ( sin(phi) + x/y*cos(phi) ) * arctg );

A1 = subs(A1, rho, sqrt(sqrt(x*x+y*y)));
A2 = subs(A2, rho, sqrt(sqrt(x*x+y*y)));
A1 = subs(A1, phi, 1/2*atan(y/x));
A2 = subs(A2, phi, 1/2*atan(y/x));

A1y0 = subs(simplify(limit(A1, y, 0)), (x^2)^(1/2), abs(x));
A1y0 = subs(A1y0, 4*(x^2)^(1/4), 4*sqrt(abs(x))); 
pretty(A1y0)
A2y0 = subs(simplify(limit(A2, y, 0)), (x^2)^(1/2), abs(x));
A2y0 = subs(A2y0, 4*(x^2)^(1/4), 4*sqrt(abs(x))); 
pretty(A2y0)

syms Ix Iy I
Ix = A1*y;
Iy = A1*x - A2;

err = ( diff(Ix, y) - diff(Iy, x) );
deriv = subs(subs(subs(diff(Ix, y), x, 10), y, 10), a, 1)
err = subs(subs(subs(err, x, 10), y, 10), a, 1)

Ix = subs(Ix, x, a+x);
Iy = subs(Iy, x, a+x);

syms M_PI nu F2 F3 F4 F5 F6 F7 rho phi
F2 = -1/(4*M_PI*(1-nu)*sqrt(a))*Ix;
pretty(F2)
F3 = -1/(4*M_PI*(1-nu)*sqrt(a))*Iy;
F4 = diff(F2, y);
F5 = diff(F2, x);
F6 = diff(F4, y);
F7 = diff(diff(F3, y), y);

Fs = cell(1,7);
Fs{2} = char(subs(subs(F2, ((a + x)^2 + y^2), rho), atan(y/(a + x)), phi));
Fs{3} = char(subs(subs(F3, ((a + x)^2 + y^2), rho), atan(y/(a + x)), phi));
Fs{4} = char(subs(subs(F4, ((a + x)^2 + y^2), rho), atan(y/(a + x)), phi));
Fs{5} = char(subs(subs(F5, ((a + x)^2 + y^2), rho), atan(y/(a + x)), phi));
Fs{6} = char(subs(subs(F6, ((a + x)^2 + y^2), rho), atan(y/(a + x)), phi));
Fs{7} = char(subs(subs(F7, ((a + x)^2 + y^2), rho), atan(y/(a + x)), phi));

for i = 2:1:7
    Fs{i} = strrep(Fs{i}, 'rho^(1/2)', 'sqrt(rho)');
    Fs{i} = strrep(Fs{i}, 'rho^(1/4)', 'sqrt(sqrt(rho))');
    Fs{i} = strrep(Fs{i}, '2^(1/2)*a^(1/2)', 'sqrt(2*a)');
    Fs{i} = strrep(Fs{i}, 'a^(1/2)', 'sqrt(a)');
    Fs{i} = strrep(Fs{i}, 'y^2', '(y*y)');
    Fs{i} = strrep(Fs{i}, 'y^3', '(y*y*y)');
    Fs{i} = strrep(Fs{i}, '(a + x)^2', '((x+a)*(x+a))');
    Fs{i} = strrep(Fs{i}, '(a + x)^3', '((x+a)*(x+a)*(x+a))');
    Fs{i} = strrep(Fs{i}, 'rho^(3/4)', 'pow(rho, 0.75)');
    Fs{i} = strrep(Fs{i}, 'rho^(5/4)', 'pow(rho, 1.25)');
    Fs{i} = strrep(Fs{i}, 'rho^(7/4)', 'pow(rho, 1.75)');
    Fs{i} = strrep(Fs{i}, 'rho^(9/4)', 'pow(rho, 2.25)');
    Fs{i} = strrep(Fs{i}, 'rho^(3/2)', 'pow(rho, 1.5)');
    Fs{i} = strrep(Fs{i}, '(2*a - sqrt(rho))^2', '((2*a - sqrt(rho))*(2*a - sqrt(rho)))');
    Fs{i} = strrep(Fs{i}, 'sin(phi/2)^2', '(sin(phi/2)*sin(phi/2))');
    Fs{i} = strrep(Fs{i}, 'cos(phi/2)^2', '(cos(phi/2)*cos(phi/2))');
    Fs{i} = strrep(Fs{i}, '(2*a + sqrt(rho) + 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))^2', '((2*a + sqrt(rho) + 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))*(2*a + sqrt(rho) + 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2)))');
    Fs{i} = strrep(Fs{i}, '((y*y)/((x+a)*(x+a)) + 1)^2', '(((y*y)/((x+a)*(x+a)) + 1)*((y*y)/((x+a)*(x+a)) + 1))');
    Fs{i} = strrep(Fs{i}, '(2*a - sqrt(rho))^3', '((2*a - sqrt(rho))*(2*a - sqrt(rho))*(2*a - sqrt(rho)))');
    Fs{i} = strrep(Fs{i}, '(2*a + sqrt(rho) + 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))^3', '((2*a + sqrt(rho) + 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))*(2*a + sqrt(rho) + 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))*(2*a + sqrt(rho) + 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2)))');
    Fs{i} = strrep(Fs{i}, '(2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))^2', '(2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))*(2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))');
    Fs{i} = strrep(Fs{i}, '((8*a*sqrt(rho)*(sin(phi/2)*sin(phi/2)))/((2*a - sqrt(rho))*(2*a - sqrt(rho))) + 1)^2', '(((8*a*sqrt(rho)*(sin(phi/2)*sin(phi/2)))/((2*a - sqrt(rho))*(2*a - sqrt(rho))) + 1)*((8*a*sqrt(rho)*(sin(phi/2)*sin(phi/2)))/((2*a - sqrt(rho))*(2*a - sqrt(rho))) + 1))');
    Fs{i} = strrep(Fs{i}, '(y/sqrt(rho) + (sqrt(2*a)*y*cos(phi/2))/pow(rho, 0.75) - (sqrt(2*a)*sqrt(sqrt(rho))*sin(phi/2))/(((y*y)/((x+a)*(x+a)) + 1)*(a + x)))^2', '((y/sqrt(rho) + (sqrt(2*a)*y*cos(phi/2))/pow(rho, 0.75) - (sqrt(2*a)*sqrt(sqrt(rho))*sin(phi/2))/(((y*y)/((x+a)*(x+a)) + 1)*(a + x)))*(y/sqrt(rho) + (sqrt(2*a)*y*cos(phi/2))/pow(rho, 0.75) - (sqrt(2*a)*sqrt(sqrt(rho))*sin(phi/2))/(((y*y)/((x+a)*(x+a)) + 1)*(a + x))))');
    Fs{i} = strrep(Fs{i}, '((2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))*(2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2)))^2', '(((2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))*(2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2)))*((2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))*(2*a + sqrt(rho) - 2*sqrt(2*a)*sqrt(sqrt(rho))*cos(phi/2))))');
end

file = fopen('derivatives.cpp', 'w');
for i = 2:1:7
    fprintf(file, 'double Break::sqrtF%d(const double& x, const double& y) const {\n', i);
    fprintf(file, '\tdouble rho = (a + x)*(a + x) + y*y;\n');
    fprintf(file, '\tdouble phi = atan(y/(a + x));\n');
    fprintf(file, '\treturn %s;\n}\n\n', Fs{i});
end
fclose(file);


