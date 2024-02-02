function dydt = PG_eq(t,y,z) 
syms p1 p2

p2=(y(1));
p1=(y(2));

EQ2=makfun_PG(z,p1,p2,t);
t
dydt = [p1; EQ2];
end