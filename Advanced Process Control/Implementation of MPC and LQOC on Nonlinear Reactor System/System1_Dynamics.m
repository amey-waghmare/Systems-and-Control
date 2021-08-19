function fofX = System1_Dynamics(t, X)

global sys

fofX = zeros(3,1) ; 

k = sys.alfa(1) * exp(-sys.alfa(2)/X(2));

fofX(1) = (sys.Uk(1)/100)*(sys.Dk-X(1)) - 2* k* X(1)^2;

fofX(2) = (sys.Uk(1)/100)* (275-X(2)) + sys.alfa(3)* k * X(1)^2 - sys.alfa(4) * (X(2)-X(3)) ;

fofX(3) = (sys.Uk(2)/10)*(250 -X(3)) + sys.alfa(5) * (X(2)-X(3));



