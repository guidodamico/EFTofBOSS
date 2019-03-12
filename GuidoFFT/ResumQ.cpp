#include "ResumEFT.h"

double ResumQ0 (const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & X1Infinity, const double & f1) {

	double X1q = X1(q, InterpX1) ;
	double Y1q = Y1(q, InterpY1) ;

	int i = l/2, j = lp/2 ;

switch (i) { 
	case 0:
		switch (j) { 
			case 0:
				return -(pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*Pi*q*(96 + 16*(-1 + 3*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity + (4*(3 - 10*beta + \
15*pow(beta,2))*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1Infinity,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + \
beta)*beta))*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1Infinity,3))/35.)*sin(k*q))/(24.*k) + \
(2*Pi*pow(q,2)*(((96 + 16*(-1 + 3*beta)*f1*(2 + f1)*pow(k,2)*X1q + \
(4*(3 - 10*beta + 15*pow(beta,2))*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + \
beta)*beta))*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/35.)*sin(k*q))/(48.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) + (alpha*k*(96 + 16*(-1 + 3*beta)*f1*(2 + \
f1)*pow(k,2)*X1q + (4*(3 - 10*beta + 15*pow(beta,2))*pow(f1,2)*pow(2 \
+ f1,2)*pow(k,4)*pow(X1q,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + \
beta)*beta))*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/35.)*Y1q*sin(k*q))/(96.*pow(E,(beta*f1*(2 \
+ f1)*pow(k,2)*X1q)/2.)*q) - (f1*((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(48.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) + ((96 + 16*(-1 + 3*beta)*f1*(2 + \
f1)*pow(k,2)*X1q + (4*(3 - 10*beta + 15*pow(beta,2))*pow(f1,2)*pow(2 \
+ f1,2)*pow(k,4)*pow(X1q,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + \
beta)*beta))*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/35.)*Y1q*(-(k*q*(-6 + \
pow(k,2)*pow(q,2))*cos(k*q)) + 3*(-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(96.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,4)) - \
(pow(f1,2)*pow(k,2)*Y1q*((((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*sin(k*q))/(144.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) + ((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(62370.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/2.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.);
				break ;

			case 1:
				return (pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*f1*(2 + f1)*Pi*X1Infinity*(-168 + f1*(2 \
+ f1)*pow(k,2)*X1Infinity*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + f1)*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(630.*k*q) + (2*Pi*pow(q,2)*((f1*(2 + \
f1)*(-12 + pow(k,2)*pow(q,2))*X1q*(-168 + f1*(2 + \
f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(k*q*(-15 + pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 \
- 2*pow(k,2)*pow(q,2))*sin(k*q)))/(2520.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,6)) - (f1*(2 + f1)*X1q*(-168 + \
f1*(2 + f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + f1)*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(1260.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) - (alpha*f1*(2 + f1)*k*X1q*(-168 + \
f1*(2 + f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(2520.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) - \
(pow(f1,2)*pow(k,2)*Y1q*(-(((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*sin(k*q))/(360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) - ((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(79380.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (4*f1*(2 + f1)*X1q*(-1144 \
+ f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(525525.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5))))/2. - \
f1*pow(k,2)*Y1q*(-(((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(120.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (f1*(2 + f1)*X1q*(-792 + \
f1*(2 + f1)*pow(k,2)*X1q*(44*(5 - 9*beta) + (-35 + 11*(10 - \
9*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(23100.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.);
				break ;

			case 2:
				return (-2*pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(f1,2)*pow(2 + \
f1,2)*Pi*pow(X1Infinity,2)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(5*k*q*(-21 + 2*pow(k,2)*pow(q,2))*cos(k*q) \
+ (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(3465.*k*pow(q,3)) + \
(2*Pi*pow(q,2)*((pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(22 + (-5 + \
11*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(3465.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,5)) + (alpha*pow(f1,2)*pow(2 + \
f1,2)*k*pow(X1q,2)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(5*k*q*(-21 + 2*pow(k,2)*pow(q,2))*cos(k*q) + \
(105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(6930.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,5)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(-(k*q*(-22050 + 3225*pow(k,2)*pow(q,2) - \
111*pow(k,4)*pow(q,4) + pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + \
10575*pow(k,2)*pow(q,2) - 696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(6930.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,8)) - \
(pow(f1,2)*pow(k,2)*Y1q*(((11088 + 1584*(-6 + 7*beta)*f1*pow(k,2)*X1q \
+ 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(218295.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + (f1*(2 + f1)*X1q*(-1144 + \
f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(88935.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + (4*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-28 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(891891.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7))))/2. - \
f1*pow(k,2)*Y1q*((f1*(2 + f1)*X1q*(-792 + f1*(2 + \
f1)*pow(k,2)*X1q*(44*(5 - 9*beta) + (-35 + 11*(10 - \
9*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(31185.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + (5*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(26 + (-7 + 13*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(81081.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 3:
				return (-4*pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(f1,3)*pow(2 + \
f1,3)*Pi*pow(X1Infinity,3)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9009.*k*pow(q,5)) + \
(2*Pi*pow(q,2)*((2*pow(f1,3)*pow(2 + f1,3)*pow(X1q,3)*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9009.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,7)) + (alpha*pow(f1,3)*pow(2 + \
f1,3)*k*pow(X1q,3)*Y1q*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9009.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,7)) - (pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*Y1q*(k*q*(5239080 - 773955*pow(k,2)*pow(q,2) + \
27909*pow(k,4)*pow(q,4) - 342*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - 840105*pow(k,2)*pow(q,2) + \
56490*pow(k,4)*pow(q,4) - 1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(9009.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,10)) - \
(pow(f1,2)*pow(k,2)*Y1q*((-2*f1*(2 + f1)*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(429429.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) - (166*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-28 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.2297275e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) + (896*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(6.3996075e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9))))/2. - \
f1*pow(k,2)*Y1q*((-2*pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(26 + (-7 + \
13*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(k*q*(5670 - 735*pow(k,2)*pow(q,2) \
+ 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(39039.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) + (14*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(250965.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 4:
				return (2*Pi*pow(q,2)*((16*pow(f1,4)*pow(2 + \
f1,3)*pow(X1q,3)*Y1q*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(328185.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,9)) - \
(pow(f1,2)*pow(k,2)*Y1q*((16*pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(90 + \
(-28 + 45*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(4.922775e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) - (16*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(654075.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9))))/2.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

		}
		break ;

	case 1:
		switch (j) { 
			case 0:
				return -(pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*f1*(2 + f1)*k*Pi*q*X1Infinity*(-168 + \
f1*(2 + f1)*pow(k,2)*X1Infinity*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + f1)*pow(k,2)*X1Infinity))*sin(k*q))/126. + \
(10*Pi*pow(q,2)*((f1*(2 + f1)*k*X1q*(-168 + f1*(2 + \
f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*sin(k*q))/(1260.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) + (alpha*f1*(2 + f1)*pow(k,3)*X1q*(-168 + \
f1*(2 + f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*sin(k*q))/(2520.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) - (f1*(11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(41580.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) + (f1*(2 + f1)*X1q*(-168 + f1*(2 + \
f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(-(k*q*(-6 + pow(k,2)*pow(q,2))*cos(k*q)) + \
3*(-2 + pow(k,2)*pow(q,2))*sin(k*q)))/(2520.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) - (pow(f1,2)*pow(k,2)*Y1q*(((11088 + \
1584*(-6 + 7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(124740.*pow(E,(beta*f1*(2 \
+ f1)*pow(k,2)*X1q)/2.)*k*q) + ((226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(1.62162e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/2.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 1:
				return (pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*Pi*(11088 + 528*(-11 + \
21*beta)*f1*pow(k,2)*X1Infinity + 6*(-85 + 297*beta - 363*pow(beta,2) \
+ 231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-85 + \
297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
264*pow(f1,2)*pow(k,2)*X1Infinity*(-11 + 21*beta + (9 - 22*beta + \
21*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(33*(9 - 22*beta + \
21*pow(beta,2)) + (-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(11*(9 - 22*beta + \
21*pow(beta,2)) + 2*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(2772.*pow(k,3)*q) + \
(10*Pi*pow(q,2)*(((-12 + pow(k,2)*pow(q,2))*(11088 + 528*(-11 + \
21*beta)*f1*pow(k,2)*X1q + 6*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-85 + 297*beta - \
363*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
264*pow(f1,2)*pow(k,2)*X1q*(-11 + 21*beta + (9 - 22*beta + \
21*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(33*(9 - 22*beta + 21*pow(beta,2)) + \
(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(11*(9 - 22*beta + 21*pow(beta,2)) + \
2*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*Y1q*(k*q*(-15 + \
pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 - \
2*pow(k,2)*pow(q,2))*sin(k*q)))/(55440.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,6)) - ((11088 + 528*(-11 + \
21*beta)*f1*pow(k,2)*X1q + 6*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-85 + 297*beta - \
363*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
264*pow(f1,2)*pow(k,2)*X1q*(-11 + 21*beta + (9 - 22*beta + \
21*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(33*(9 - 22*beta + 21*pow(beta,2)) + \
(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(11*(9 - 22*beta + 21*pow(beta,2)) + \
2*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(27720.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (alpha*(11088 + 528*(-11 + \
21*beta)*f1*pow(k,2)*X1q + 6*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-85 + 297*beta - \
363*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
264*pow(f1,2)*pow(k,2)*X1q*(-11 + 21*beta + (9 - 22*beta + \
21*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(33*(9 - 22*beta + 21*pow(beta,2)) + \
(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(11*(9 - 22*beta + 21*pow(beta,2)) + \
2*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*Y1q*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(55440.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) - (pow(f1,2)*pow(k,2)*Y1q*(-((11088 \
+ 1584*(-6 + 7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(311850.*pow(E,(beta*f1*(2 \
+ f1)*pow(k,2)*X1q)/2.)*k*q) - ((226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(2.06388e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (2*(6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(525525.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/2. - \
f1*pow(k,2)*Y1q*(-((11088 + 1584*(-6 + 7*beta)*f1*pow(k,2)*X1q + \
6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(103950.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((61776 + 6864*(-7 + \
9*beta)*f1*pow(k,2)*X1q + 6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-805 + 2665*beta - \
3003*pow(beta,2) + 1287*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-7 + 9*beta) + (205 - 462*beta + \
297*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + (-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + 6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(600600.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 2:
				return -(pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*f1*(2 + f1)*Pi*X1Infinity*(-3432 + f1*(2 \
+ f1)*pow(k,2)*X1Infinity*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(9009.*pow(k,3)*pow(q,3)) + \
(10*Pi*pow(q,2)*((f1*(2 + f1)*X1q*(-3432 + f1*(2 + \
f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + (alpha*f1*(2 + \
f1)*X1q*(-3432 + f1*(2 + f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + \
13*(34 - 33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,5)) + (f1*(2 + f1)*X1q*(-3432 + f1*(2 + \
f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(-(k*q*(-22050 + \
3225*pow(k,2)*pow(q,2) - 111*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + 10575*pow(k,2)*pow(q,2) - \
696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,8)) - \
(pow(f1,2)*pow(k,2)*Y1q*(((226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(5.67567e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + ((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(177870.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (2*f1*(2 + f1)*X1q*(-18360 \
+ f1*(2 + f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 + 153*(23 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.5162147e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7))))/2. - \
f1*pow(k,2)*Y1q*(((61776 + 6864*(-7 + 9*beta)*f1*pow(k,2)*X1q + \
6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-805 + 2665*beta - \
3003*pow(beta,2) + 1287*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-7 + 9*beta) + (205 - 462*beta + \
297*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + (-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + 6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(810810.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (f1*(2 + f1)*X1q*(-1560 + \
f1*(2 + f1)*pow(k,2)*X1q*(500 - 780*beta + (-91 + 5*(50 - \
39*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(5670 - \
735*pow(k,2)*pow(q,2) + 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + \
2625*pow(k,2)*pow(q,2) - 135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(162162.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 3:
				return (2*pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(f1,2)*pow(2 + \
f1,2)*Pi*pow(X1Infinity,2)*(90 + (-23 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9009.*pow(k,3)*pow(q,5)) + \
(10*Pi*pow(q,2)*(-(pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(90 + (-23 + \
45*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(45045.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) - (alpha*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-23 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,7)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-23 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(k*q*(5239080 - 773955*pow(k,2)*pow(q,2) + \
27909*pow(k,4)*pow(q,4) - 342*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - 840105*pow(k,2)*pow(q,2) + \
56490*pow(k,4)*pow(q,4) - 1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,10)) - \
(pow(f1,2)*pow(k,2)*Y1q*(-((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(429429.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) - (83*f1*(2 + \
f1)*X1q*(-18360 + f1*(2 + f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 \
+ 153*(23 - 15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(3.79053675e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (448*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(342 + (-116 + 171*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.215925425e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9))))/2. - \
f1*pow(k,2)*Y1q*(-(f1*(2 + f1)*X1q*(-1560 + f1*(2 + \
f1)*pow(k,2)*X1q*(500 - 780*beta + (-91 + 5*(50 - \
39*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(5670 - \
735*pow(k,2)*pow(q,2) + 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + \
2625*pow(k,2)*pow(q,2) - 135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(195195.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (7*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(306 + (-91 + 153*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.266405e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 4:
				return (16*pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(f1,3)*pow(2 + \
f1,3)*Pi*pow(X1Infinity,3)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) \
- 770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(21879.*pow(k,3)*pow(q,7)) + \
(10*Pi*pow(q,2)*((-8*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(109395.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9)) - (4*alpha*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*Y1q*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(109395.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,9)) - (4*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*Y1q*(-(k*q*(-2006754750 + \
300405105*pow(k,2)*pow(q,2) - 11320155*pow(k,4)*pow(q,4) + \
158679*pow(k,6)*pow(q,6) - 852*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*cos(k*q)) + 3*(-668918250 + \
323107785*pow(k,2)*pow(q,2) - 22286880*pow(k,4)*pow(q,4) + \
501165*pow(k,6)*pow(q,6) - 4418*pow(k,8)*pow(q,8) + \
13*pow(k,10)*pow(q,10))*sin(k*q)))/(109395.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,12)) - \
(pow(f1,2)*pow(k,2)*Y1q*((8*f1*(2 + f1)*X1q*(-18360 + f1*(2 + \
f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 + 153*(23 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(8.3687175e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) + (8*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(342 + (-116 + 171*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.2427425e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) - (480*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(55*k*q*(11904165 - 1670760*pow(k,2)*pow(q,2) + \
51597*pow(k,4)*pow(q,4) - 468*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) + (-654729075 + \
310134825*pow(k,2)*pow(q,2) - 18918900*pow(k,4)*pow(q,4) + \
315315*pow(k,6)*pow(q,6) - 1485*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(1.04433329e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,11))))/2. - \
f1*pow(k,2)*Y1q*((8*pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(306 + (-91 + \
153*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(5.579145e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) - (72*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(k*q*(344594250 - 49324275*pow(k,2)*pow(q,2) + \
1621620*pow(k,4)*pow(q,4) - 16830*pow(k,6)*pow(q,6) + \
46*pow(k,8)*pow(q,8))*cos(k*q) + (-344594250 + \
164189025*pow(k,2)*pow(q,2) - 10405395*pow(k,4)*pow(q,4) + \
190575*pow(k,6)*pow(q,6) - 1080*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(3.926065e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,11)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;	
				break ;

		}
		break ;

	case 2:
		switch (j) { 
			case 0:
				return (-2*pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(f1,2)*pow(2 + \
f1,2)*pow(k,3)*Pi*q*pow(X1Infinity,2)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*sin(k*q))/385. + \
(18*Pi*pow(q,2)*((pow(f1,2)*pow(2 + f1,2)*pow(k,3)*pow(X1q,2)*(22 + \
(-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*sin(k*q))/(3465.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) + (alpha*pow(f1,2)*pow(2 + \
f1,2)*pow(k,5)*pow(X1q,2)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*sin(k*q))/(6930.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) - (pow(f1,2)*(2 + f1)*k*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(45045.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(k,2)*pow(X1q,2)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(-(k*q*(-6 + pow(k,2)*pow(q,2))*cos(k*q)) + \
3*(-2 + pow(k,2)*pow(q,2))*sin(k*q)))/(6930.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) - (pow(f1,2)*pow(k,2)*Y1q*((f1*(2 + \
f1)*k*X1q*(-1144 + f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + \
13*(15 - 11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*sin(k*q))/(135135.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) + ((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(135135.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/2.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 1:
				return  (pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*f1*(2 + f1)*Pi*X1Infinity*(-3432 + f1*(2 \
+ f1)*pow(k,2)*X1Infinity*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + \
(-3 + pow(k,2)*pow(q,2))*sin(k*q)))/(5005.*k*q) + \
(18*Pi*pow(q,2)*((f1*(2 + f1)*(-12 + pow(k,2)*pow(q,2))*X1q*(-3432 + \
f1*(2 + f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(k*q*(-15 + \
pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 - \
2*pow(k,2)*pow(q,2))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,6)) - (f1*(2 + f1)*X1q*(-3432 + \
f1*(2 + f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) - (alpha*f1*(2 + f1)*k*X1q*(-3432 + \
f1*(2 + f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) - (pow(f1,2)*pow(k,2)*Y1q*((-2*f1*(2 \
+ f1)*k*X1q*(-1144 + f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 \
+ 13*(15 - 11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*sin(k*q))/(675675.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) - ((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(171990.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((689520 + 816*(-643 + \
845*beta)*f1*pow(k,2)*X1q + 6*(-8295 + 27795*beta - 32793*pow(beta,2) \
+ 14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.786785e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/2. - \
f1*pow(k,2)*Y1q*((-2*f1*(2 + f1)*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(225225.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) - ((34320 + 624*(-42 + \
55*beta)*f1*pow(k,2)*X1q + 6*(-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-546 + 715*beta + (465 - 1092*beta + \
715*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(465 - 1092*beta + \
715*pow(beta,2)) + (-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(465 \
- 1092*beta + 715*pow(beta,2) + 2*(-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(450450.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 2:
				return  -(pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*Pi*(240240 + 3120*(-39 + \
77*beta)*f1*pow(k,2)*X1Infinity + 6*(-1635 + 5787*beta - \
7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-1635 + \
5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
24*pow(f1,2)*pow(k,2)*X1Infinity*(65*(-39 + 77*beta) + (1929 - \
5070*beta + 5005*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(3*(1929 - 5070*beta + \
5005*pow(beta,2)) + (-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(1929 - 5070*beta + \
5005*pow(beta,2) + 2*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1Infinity))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(60060.*pow(k,5)*pow(q,3)) + \
(18*Pi*pow(q,2)*(((240240 + 3120*(-39 + 77*beta)*f1*pow(k,2)*X1q + \
6*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-1635 + 5787*beta \
- 7605*pow(beta,2) + 5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 24*pow(f1,2)*pow(k,2)*X1q*(65*(-39 + 77*beta) + (1929 - 5070*beta + \
5005*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(1929 - 5070*beta + \
5005*pow(beta,2)) + (-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(1929 - 5070*beta + 5005*pow(beta,2) \
+ 2*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.08108e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (alpha*(240240 + 3120*(-39 \
+ 77*beta)*f1*pow(k,2)*X1q + 6*(-1635 + 5787*beta - 7605*pow(beta,2) \
+ 5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-1635 + \
5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(65*(-39 + 77*beta) + (1929 - 5070*beta + \
5005*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(1929 - 5070*beta + \
5005*pow(beta,2)) + (-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(1929 - 5070*beta + 5005*pow(beta,2) \
+ 2*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q))*Y1q*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.16216e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + ((240240 + 3120*(-39 + \
77*beta)*f1*pow(k,2)*X1q + 6*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-1635 + 5787*beta \
- 7605*pow(beta,2) + 5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 24*pow(f1,2)*pow(k,2)*X1q*(65*(-39 + 77*beta) + (1929 - 5070*beta + \
5005*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(1929 - 5070*beta + \
5005*pow(beta,2)) + (-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(1929 - 5070*beta + 5005*pow(beta,2) \
+ 2*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q))*Y1q*(-(k*q*(-22050 + \
3225*pow(k,2)*pow(q,2) - 111*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + 10575*pow(k,2)*pow(q,2) - \
696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(2.16216e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,6)*pow(q,8)) - \
(pow(f1,2)*pow(k,2)*Y1q*((2*(6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(945945.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + ((689520 + 816*(-643 + \
845*beta)*f1*pow(k,2)*X1q + 6*(-8295 + 27795*beta - 32793*pow(beta,2) \
+ 14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.209516e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + ((2713200 + 25840*(-106 + \
105*beta)*f1*pow(k,2)*X1q + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + 56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) \
+ (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.92053862e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7))))/2. - \
f1*pow(k,2)*Y1q*((2*(34320 + 624*(-42 + 55*beta)*f1*pow(k,2)*X1q + \
6*(-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-546 + 715*beta + (465 - 1092*beta + \
715*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(465 - 1092*beta + \
715*pow(beta,2)) + (-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(465 \
- 1092*beta + 715*pow(beta,2) + 2*(-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.216215e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (5*(371280 + 4080*(-69 + \
91*beta)*f1*pow(k,2)*X1q + 6*(-4389 + 14841*beta - 17595*pow(beta,2) \
+ 7735*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-4389 + \
14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(5*(-69 + 91*beta) + (291 - 690*beta + \
455*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(291 - 690*beta + \
455*pow(beta,2)) + (-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(291 - 690*beta + \
455*pow(beta,2)) + 2*(-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(3.3081048e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 3:
				return (pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*f1*(2 + f1)*Pi*X1Infinity*(-14280 + \
f1*(2 + f1)*pow(k,2)*X1Infinity*(3604 - 7140*beta + (-569 + 17*(106 - \
105*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(34034.*pow(k,5)*pow(q,5)) + \
(18*Pi*pow(q,2)*(-(f1*(2 + f1)*X1q*(-14280 + f1*(2 + \
f1)*pow(k,2)*X1q*(3604 - 7140*beta + (-569 + 17*(106 - \
105*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(612612.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (alpha*f1*(2 + \
f1)*X1q*(-14280 + f1*(2 + f1)*pow(k,2)*X1q*(3604 - 7140*beta + (-569 \
+ 17*(106 - 105*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.225224e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) + (f1*(2 + f1)*X1q*(-14280 + \
f1*(2 + f1)*pow(k,2)*X1q*(3604 - 7140*beta + (-569 + 17*(106 - \
105*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(k*q*(5239080 - \
773955*pow(k,2)*pow(q,2) + 27909*pow(k,4)*pow(q,4) - \
342*pow(k,6)*pow(q,6) + pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - \
840105*pow(k,2)*pow(q,2) + 56490*pow(k,4)*pow(q,4) - \
1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(1.225224e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,6)*pow(q,10)) - \
(pow(f1,2)*pow(k,2)*Y1q*(-((689520 + 816*(-643 + \
845*beta)*f1*pow(k,2)*X1q + 6*(-8295 + 27795*beta - 32793*pow(beta,2) \
+ 14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.9201172e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) - (83*(2713200 + 25840*(-106 \
+ 105*beta)*f1*pow(k,2)*X1q + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + 56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) \
+ (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9.6026931e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (112*f1*(2 + \
f1)*X1q*(-5320 + f1*(2 + f1)*pow(k,2)*X1q*(2010 - 2660*beta + (-394 + \
5*(201 - 133*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.05308475e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9))))/2. - \
f1*pow(k,2)*Y1q*(-((371280 + 4080*(-69 + 91*beta)*f1*pow(k,2)*X1q + \
6*(-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-4389 + 14841*beta \
- 17595*pow(beta,2) + 7735*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 408*pow(f1,2)*pow(k,2)*X1q*(5*(-69 + 91*beta) + (291 - 690*beta + \
455*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(291 - 690*beta + \
455*pow(beta,2)) + (-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(291 - 690*beta + \
455*pow(beta,2)) + 2*(-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(7.963956e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (7*f1*(2 + f1)*X1q*(-90440 \
+ f1*(2 + f1)*pow(k,2)*X1q*(380*(75 - 119*beta) + (-5033 + 95*(150 - \
119*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.0808226e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

			case 4:
				return (-4*pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(f1,2)*pow(2 + \
f1,2)*Pi*pow(X1Infinity,2)*(266 + (-67 + 133*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(46189.*pow(k,5)*pow(q,7)) + \
(18*Pi*pow(q,2)*((2*pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(266 + (-67 + \
133*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(415701.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) + (alpha*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(266 + (-67 + 133*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(415701.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(266 + (-67 + 133*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(-(k*q*(-2006754750 + \
300405105*pow(k,2)*pow(q,2) - 11320155*pow(k,4)*pow(q,4) + \
158679*pow(k,6)*pow(q,6) - 852*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*cos(k*q)) + 3*(-668918250 + \
323107785*pow(k,2)*pow(q,2) - 22286880*pow(k,4)*pow(q,4) + \
501165*pow(k,6)*pow(q,6) - 4418*pow(k,8)*pow(q,8) + \
13*pow(k,10)*pow(q,10))*sin(k*q)))/(415701.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,6)*pow(q,12)) - \
(pow(f1,2)*pow(k,2)*Y1q*((2*(2713200 + 25840*(-106 + \
105*beta)*f1*pow(k,2)*X1q + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + 56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) \
+ (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(5.30018775e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) + (2*f1*(2 + f1)*X1q*(-5320 \
+ f1*(2 + f1)*pow(k,2)*X1q*(2010 - 2660*beta + (-394 + 5*(201 - \
133*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.142475e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)) + (600*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(322 + (-108 + 161*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(55*k*q*(11904165 - 1670760*pow(k,2)*pow(q,2) + \
51597*pow(k,4)*pow(q,4) - 468*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) + (-654729075 + \
310134825*pow(k,2)*pow(q,2) - 18918900*pow(k,4)*pow(q,4) + \
315315*pow(k,6)*pow(q,6) - 1485*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(2.401966567e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,11))))/2. - \
f1*pow(k,2)*Y1q*((2*f1*(2 + f1)*X1q*(-90440 + f1*(2 + \
f1)*pow(k,2)*X1q*(380*(75 - 119*beta) + (-5033 + 95*(150 - \
119*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(3.5334585e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)) + (2*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(126 + (-37 + 63*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(344594250 - 49324275*pow(k,2)*pow(q,2) + \
1621620*pow(k,4)*pow(q,4) - 16830*pow(k,6)*pow(q,6) + \
46*pow(k,8)*pow(q,8))*cos(k*q) + (-344594250 + \
164189025*pow(k,2)*pow(q,2) - 10405395*pow(k,4)*pow(q,4) + \
190575*pow(k,6)*pow(q,6) - 1080*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(785213.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,11)))))/pow(E,(pow(k,2)*(X1q + \
alpha*Y1q))/2.) ;
				break ;

		}
		break ;

}

}

double ResumQ1 (const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & X1Infinity, const double & f1) {

	double X1q = X1(q, InterpX1) ;
	double Y1q = Y1(q, InterpY1) ;

	int i = l/2, j = lp/2 ;

switch (i) { 
	case 0:
		switch (j) { 
			case 0:
				return (-2*Pi*pow(q,2)*(((1 + (pow(k,2)*X1Infinity)/2.)*(96 + 16*(-1 + \
3*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity + (4*(3 - 10*beta + \
15*pow(beta,2))*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1Infinity,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + \
beta)*beta))*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1Infinity,3))/35.)*sin(k*q))/(48.*pow(E,(beta*f1*(\
2 + f1)*pow(k,2)*X1Infinity)/2.)*k*q) + ((4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1Infinity,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1Infinity,2)*(2 + beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity))/7. - (6*f1*(2 + f1)*pow(k,2)*X1Infinity*(8 \
+ beta*f1*(2 + f1)*pow(k,2)*X1Infinity*(4 + beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity*(24 + beta*f1*(2 + f1)*pow(k,2)*X1Infinity*(6 \
+ beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity))))/3.)*sin(k*q))/(192.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*k*q)))/pow(E,(pow(k,2)*X1Infinity)/2.) + \
(2*Pi*pow(q,2)*((((-2*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1q,3))/9. \
+ (6*pow(f1,2)*pow(2 + f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + \
f1)*pow(k,2)*X1q))/7. - (6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(4 + beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + \
beta*f1*(2 + f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 \
+ beta*f1*(2 + f1)*pow(k,2)*X1q))))/3.)*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*sin(k*q))/(96.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) + ((96 + 16*(-1 + 3*beta)*f1*(2 + \
f1)*pow(k,2)*X1q + (4*(3 - 10*beta + 15*pow(beta,2))*pow(f1,2)*pow(2 \
+ f1,2)*pow(k,4)*pow(X1q,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + \
beta)*beta))*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1q,3))/35.)*(1 + \
(pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*sin(k*q))/(48.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) - (f1*k*(4*f1*X1q + \
2*pow(f1,2)*X1q)*((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/11. + (2*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/3. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/7. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/5.)*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(192.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) - (f1*k*X1q*((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(96.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) + (X1q*(96 + 16*(-1 + 3*beta)*f1*(2 + \
f1)*pow(k,2)*X1q + (4*(3 - 10*beta + 15*pow(beta,2))*pow(f1,2)*pow(2 \
+ f1,2)*pow(k,4)*pow(X1q,2))/5. + (2*(-5 + 7*beta*(3 + 5*(-1 + \
beta)*beta))*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/35.)*Y1q*(-(k*q*(-6 + \
pow(k,2)*pow(q,2))*cos(k*q)) + 3*(-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(192.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) + ((4*f1*X1q + \
2*pow(f1,2)*X1q)*((-2*pow(f1,3)*pow(2 + f1,3)*pow(k,6)*pow(X1q,3))/9. \
+ (6*pow(f1,2)*pow(2 + f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + \
f1)*pow(k,2)*X1q))/7. - (6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(4 + beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + \
beta*f1*(2 + f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 \
+ beta*f1*(2 + f1)*pow(k,2)*X1q))))/3.)*Y1q*(-(k*q*(-6 + \
pow(k,2)*pow(q,2))*cos(k*q)) + 3*(-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(384.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*((((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*sin(k*q))/(144.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) + ((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(62370.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q*((((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/11. + (2*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/3. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/7. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/5.)*sin(k*q))/(144.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) + ((41184 + 6864*(-5 + \
6*beta)*f1*pow(k,2)*X1q + 6*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-525 + 1820*beta - \
2145*pow(beta,2) + 858*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-5 + 6*beta) + 2*(70 - 165*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(26*(70 - 165*beta + 99*pow(beta,2)) \
+ (-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q) + \
4*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(70 - 165*beta + 99*pow(beta,2)) \
+ 3*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(270270.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/8.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 1:
				return (-2*Pi*pow(q,2)*(-((4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1Infinity + 6*(-140 + 495*beta - 594*pow(beta,2) \
+ 231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-140 + \
495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
792*pow(f1,2)*pow(k,2)*X1Infinity*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(99*(5 - 12*beta + \
7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(33*(5 - 12*beta + \
7*pow(beta,2)) + 2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(166320.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,3)) - (f1*(2 + \
f1)*X1Infinity*(1 + (pow(k,2)*X1Infinity)/2.)*(-168 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + f1)*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(1260.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*k*pow(q,3))))/pow(E,(pow(k,2)*X1Infinity)\
/2.) + (2*Pi*pow(q,2)*(((-12 + pow(k,2)*pow(q,2))*(4*f1*X1q + \
2*pow(f1,2)*X1q)*(11088 + 1584*(-6 + 7*beta)*f1*pow(k,2)*X1q + \
6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*Y1q*(k*q*(-15 + \
pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 - \
2*pow(k,2)*pow(q,2))*sin(k*q)))/(332640.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,6)) + (f1*(2 + f1)*(-12 + \
pow(k,2)*pow(q,2))*pow(X1q,2)*(-168 + f1*(2 + f1)*pow(k,2)*X1q*(36 - \
84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(k*q*(-15 + pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 \
- 2*pow(k,2)*pow(q,2))*sin(k*q)))/(5040.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,6)) - ((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(83160.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (f1*(2 + f1)*X1q*(-168 + \
f1*(2 + f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + f1)*pow(k,2)*X1q))*(1 + (pow(k,2)*X1q)/2. + \
((alpha*pow(k,2))/2. + (alpha*pow(k,4)*X1q)/4.)*Y1q)*(3*k*q*cos(k*q) \
+ (-3 + pow(k,2)*pow(q,2))*sin(k*q)))/(1260.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) - (pow(f1,2)*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q*(-(((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/11. + (2*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/3. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/7. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/5.)*sin(k*q))/(360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) - ((41184 + 6864*(-5 + \
6*beta)*f1*pow(k,2)*X1q + 6*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-525 + 1820*beta - \
2145*pow(beta,2) + 858*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-5 + 6*beta) + 2*(70 - 165*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(26*(70 - 165*beta + 99*pow(beta,2)) \
+ (-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q) + \
4*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(70 - 165*beta + 99*pow(beta,2)) \
+ 3*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(343980.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (4*(6864 + 624*(-15 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-210 + 630*beta - \
585*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-195 + 143*beta + (210 - 390*beta + \
143*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(630 \
- 1170*beta + 429*pow(beta,2) + (-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(210 \
- 390*beta + 143*pow(beta,2) + 2*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.576575e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/8. - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*(-(((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*sin(k*q))/(360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) - ((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(79380.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (4*f1*(2 + f1)*X1q*(-1144 \
+ f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(525525.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5))))/4. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-(((-2*pow(f1,3)*pow(2 \
+ f1,3)*pow(k,6)*pow(X1q,3))/11. + (2*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/3. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/7. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/5.)*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(120.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((20592 + 2288*(-10 + \
9*beta)*f1*pow(k,2)*X1q + 6*(-420 + 1365*beta - 1430*pow(beta,2) + \
429*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-420 + 1365*beta - \
1430*pow(beta,2) + 429*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(-110 + 99*beta + (105 - 220*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(13*(105 - 220*beta + 99*pow(beta,2)) \
+ (-420 + 1365*beta - 1430*pow(beta,2) + \
429*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(105 - 220*beta + 99*pow(beta,2)) \
+ 6*(-420 + 1365*beta - 1430*pow(beta,2) + \
429*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(300300.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/4. - \
(f1*pow(k,4)*X1q*Y1q*(-(((-2*pow(f1,3)*pow(2 + \
f1,3)*pow(k,6)*pow(X1q,3))/9. + (6*pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,2)*(2 + beta*f1*(2 + f1)*pow(k,2)*X1q))/7. - \
(6*f1*(2 + f1)*pow(k,2)*X1q*(8 + beta*f1*(2 + f1)*pow(k,2)*X1q*(4 + \
beta*f1*(2 + f1)*pow(k,2)*X1q)))/5. + (2*(48 + beta*f1*(2 + \
f1)*pow(k,2)*X1q*(24 + beta*f1*(2 + f1)*pow(k,2)*X1q*(6 + beta*f1*(2 \
+ f1)*pow(k,2)*X1q))))/3.)*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(120.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (f1*(2 + f1)*X1q*(-792 + \
f1*(2 + f1)*pow(k,2)*X1q*(44*(5 - 9*beta) + (-35 + 11*(10 - \
9*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(23100.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5))))/2.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 2:
				return (-2*Pi*pow(q,2)*((pow(f1,2)*pow(2 + f1,2)*pow(X1Infinity,2)*(1 + \
(pow(k,2)*X1Infinity)/2.)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(5*k*q*(-21 + 2*pow(k,2)*pow(q,2))*cos(k*q) \
+ (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(3465.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*k*pow(q,5)) + (f1*(2 + \
f1)*X1Infinity*(4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,5))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (2*Pi*pow(q,2)*((f1*(2 + f1)*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(5*k*q*(-21 + 2*pow(k,2)*pow(q,2))*cos(k*q) \
+ (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(22 + (-5 + 11*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(1 + \
(pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(3465.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,5)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,3)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(-(k*q*(-22050 + 3225*pow(k,2)*pow(q,2) - \
111*pow(k,4)*pow(q,4) + pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + \
10575*pow(k,2)*pow(q,2) - 696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(13860.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,8)) + (f1*(2 + f1)*X1q*(4*f1*X1q + \
2*pow(f1,2)*X1q)*(-1144 + f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + \
(-70 + 13*(15 - 11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(-(k*q*(-22050 + 3225*pow(k,2)*pow(q,2) - \
111*pow(k,4)*pow(q,4) + pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + \
10575*pow(k,2)*pow(q,2) - 696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(360360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,8)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*(((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(218295.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + (f1*(2 + f1)*X1q*(-1144 + \
f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(88935.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + (4*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-28 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(891891.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(((41184 + \
6864*(-5 + 6*beta)*f1*pow(k,2)*X1q + 6*(-525 + 1820*beta - \
2145*pow(beta,2) + 858*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-5 + 6*beta) + 2*(70 - 165*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(26*(70 - 165*beta + 99*pow(beta,2)) \
+ (-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q) + \
4*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(70 - 165*beta + 99*pow(beta,2)) \
+ 3*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(945945.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + ((6864 + 624*(-15 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-210 + 630*beta - \
585*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-195 + 143*beta + (210 - 390*beta + \
143*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(630 \
- 1170*beta + 429*pow(beta,2) + (-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(210 \
- 390*beta + 143*pow(beta,2) + 2*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(266805.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (4*f1*(2 + f1)*X1q*(-2040 \
+ f1*(2 + f1)*pow(k,2)*X1q*(952 - 1020*beta + (-210 + 17*(28 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(5.054049e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7))))/8. - \
(f1*pow(k,4)*X1q*Y1q*((f1*(2 + f1)*X1q*(-792 + f1*(2 + \
f1)*pow(k,2)*X1q*(44*(5 - 9*beta) + (-35 + 11*(10 - \
9*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(31185.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + (5*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(26 + (-7 + 13*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(81081.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(((20592 + 2288*(-10 + \
9*beta)*f1*pow(k,2)*X1q + 6*(-420 + 1365*beta - 1430*pow(beta,2) + \
429*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-420 + 1365*beta - \
1430*pow(beta,2) + 429*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(-110 + 99*beta + (105 - 220*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(13*(105 - 220*beta + 99*pow(beta,2)) \
+ (-420 + 1365*beta - 1430*pow(beta,2) + \
429*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(105 - 220*beta + 99*pow(beta,2)) \
+ 6*(-420 + 1365*beta - 1430*pow(beta,2) + \
429*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(405405.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (f1*(2 + f1)*X1q*(-520 + \
f1*(2 + f1)*pow(k,2)*X1q*(210 - 260*beta + (-42 + 5*(21 - \
13*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(5670 - \
735*pow(k,2)*pow(q,2) + 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + \
2625*pow(k,2)*pow(q,2) - 135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(81081.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 3:
				return (-2*Pi*pow(q,2)*((2*pow(f1,3)*pow(2 + f1,3)*pow(X1Infinity,3)*(1 + \
(pow(k,2)*X1Infinity)/2.)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9009.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*k*pow(q,7)) - (pow(f1,2)*pow(2 + \
f1,2)*pow(X1Infinity,2)*(4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(90 + (-28 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(270270.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,7))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (2*Pi*pow(q,2)*(-(pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-28 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*((4*f1*pow(k,2)*X1q + 2*pow(f1,2)*pow(k,2)*X1q)/2. \
+ (alpha*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q)/4.)*(21*k*q*(495 \
- 60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(135135.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) + (2*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(1 + (pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9009.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,7)) - (pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,4)*Y1q*(k*q*(5239080 - 773955*pow(k,2)*pow(q,2) + \
27909*pow(k,4)*pow(q,4) - 342*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - 840105*pow(k,2)*pow(q,2) + \
56490*pow(k,4)*pow(q,4) - 1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(18018.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,10)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(4*f1*X1q + 2*pow(f1,2)*X1q)*(90 + (-28 + \
45*beta)*f1*(2 + f1)*pow(k,2)*X1q)*Y1q*(k*q*(5239080 - \
773955*pow(k,2)*pow(q,2) + 27909*pow(k,4)*pow(q,4) - \
342*pow(k,6)*pow(q,6) + pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - \
840105*pow(k,2)*pow(q,2) + 56490*pow(k,4)*pow(q,4) - \
1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(540540.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,10)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*((-2*f1*(2 + f1)*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(429429.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) - (166*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-28 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.2297275e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) + (896*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(6.3996075e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((-2*(6864 + \
624*(-15 + 11*beta)*f1*pow(k,2)*X1q + 6*(-210 + 630*beta - \
585*pow(beta,2) + 143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-195 + 143*beta + (210 - 390*beta + \
143*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(630 \
- 1170*beta + 429*pow(beta,2) + (-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(210 \
- 390*beta + 143*pow(beta,2) + 2*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.288287e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) - (166*f1*(2 + \
f1)*X1q*(-2040 + f1*(2 + f1)*pow(k,2)*X1q*(952 - 1020*beta + (-210 + \
17*(28 - 15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.26351225e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (896*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(38 + (-15 + 19*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.05308475e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9))))/8. - \
(f1*pow(k,4)*X1q*Y1q*((-2*pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(26 + \
(-7 + 13*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(k*q*(5670 - \
735*pow(k,2)*pow(q,2) + 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + \
2625*pow(k,2)*pow(q,2) - 135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(39039.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) + (14*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(250965.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((-2*f1*(2 + \
f1)*X1q*(-520 + f1*(2 + f1)*pow(k,2)*X1q*(210 - 260*beta + (-42 + \
5*(21 - 13*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(5670 - \
735*pow(k,2)*pow(q,2) + 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + \
2625*pow(k,2)*pow(q,2) - 135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(195195.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (14*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(34 + (-12 + 17*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.422135e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 4:
				return (8*pow(E,-(pow(k,2)*X1Infinity)/2. - (beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(f1,3)*pow(2 + \
f1,3)*Pi*pow(X1Infinity,3)*(4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(328185.*pow(k,3)*pow(q,7)) + \
(2*Pi*pow(q,2)*((-8*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*((4*f1*pow(k,2)*X1q + 2*pow(f1,2)*pow(k,2)*X1q)/2. + \
(alpha*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q)/4.)*(9*k*q*(-225225 \
+ 30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(328185.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9)) + (8*pow(f1,4)*pow(2 + \
f1,3)*k*pow(X1q,4)*Y1q*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(328185.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,9)) - (2*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-(k*q*(-2006754750 \
+ 300405105*pow(k,2)*pow(q,2) - 11320155*pow(k,4)*pow(q,4) + \
158679*pow(k,6)*pow(q,6) - 852*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*cos(k*q)) + 3*(-668918250 + \
323107785*pow(k,2)*pow(q,2) - 22286880*pow(k,4)*pow(q,4) + \
501165*pow(k,6)*pow(q,6) - 4418*pow(k,8)*pow(q,8) + \
13*pow(k,10)*pow(q,10))*sin(k*q)))/(328185.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,12)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*((16*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-28 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(4.922775e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) - (16*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(654075.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((16*f1*(2 + \
f1)*X1q*(-2040 + f1*(2 + f1)*pow(k,2)*X1q*(952 - 1020*beta + (-210 + \
17*(28 - 15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.7895725e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) + (16*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(38 + (-15 + 19*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.142475e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) - (320*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(55*k*q*(11904165 - 1670760*pow(k,2)*pow(q,2) + \
51597*pow(k,4)*pow(q,4) - 468*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) + (-654729075 + \
310134825*pow(k,2)*pow(q,2) - 18918900*pow(k,4)*pow(q,4) + \
315315*pow(k,6)*pow(q,6) - 1485*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(1.04433329e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,11))))/8. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((16*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(34 + (-12 + 17*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.859715e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) - (48*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(k*q*(344594250 - 49324275*pow(k,2)*pow(q,2) + \
1621620*pow(k,4)*pow(q,4) - 16830*pow(k,6)*pow(q,6) + \
46*pow(k,8)*pow(q,8))*cos(k*q) + (-344594250 + \
164189025*pow(k,2)*pow(q,2) - 10405395*pow(k,4)*pow(q,4) + \
190575*pow(k,6)*pow(q,6) - 1080*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(3.926065e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,11))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.);
				break ;

		}
		break ;

	case 1:
		switch (j) { 
			case 0:
				return (-10*Pi*pow(q,2)*(((4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1Infinity + 6*(-140 + 495*beta - 594*pow(beta,2) \
+ 231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-140 + \
495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
792*pow(f1,2)*pow(k,2)*X1Infinity*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(99*(5 - 12*beta + \
7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(33*(5 - 12*beta + \
7*pow(beta,2)) + 2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity))*sin(k*q))/(166320.*pow(E,(beta*\
f1*(2 + f1)*pow(k,2)*X1Infinity)/2.)*k*q) + (f1*(2 + \
f1)*k*X1Infinity*(1 + (pow(k,2)*X1Infinity)/2.)*(-168 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(36 - 84*beta - (5 + 3*beta*(-6 + \
7*beta))*f1*(2 + \
f1)*pow(k,2)*X1Infinity))*sin(k*q))/(1260.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*q)))/pow(E,(pow(k,2)*X1Infinity)/2.) + \
(10*Pi*pow(q,2)*(((11088 + 1584*(-6 + 7*beta)*f1*pow(k,2)*X1q + \
6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*sin(k*q))/(83160.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*q) + (f1*(2 + f1)*k*X1q*(-168 + f1*(2 + \
f1)*pow(k,2)*X1q*(36 - 84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*(1 + (pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*sin(k*q))/(1260.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) - (f1*k*X1q*(11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(83160.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) - (f1*k*(4*f1*X1q + \
2*pow(f1,2)*X1q)*(41184 + 6864*(-5 + 6*beta)*f1*pow(k,2)*X1q + \
6*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-525 + 1820*beta - \
2145*pow(beta,2) + 858*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-5 + 6*beta) + 2*(70 - 165*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(26*(70 - 165*beta + 99*pow(beta,2)) \
+ (-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q) + \
4*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(70 - 165*beta + 99*pow(beta,2)) \
+ 3*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q))*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(720720.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) + ((4*f1*X1q + \
2*pow(f1,2)*X1q)*(11088 + 1584*(-6 + 7*beta)*f1*pow(k,2)*X1q + \
6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*Y1q*(-(k*q*(-6 + \
pow(k,2)*pow(q,2))*cos(k*q)) + 3*(-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(332640.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) + (f1*(2 + \
f1)*pow(k,2)*pow(X1q,2)*(-168 + f1*(2 + f1)*pow(k,2)*X1q*(36 - \
84*beta - (5 + 3*beta*(-6 + 7*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(-(k*q*(-6 + pow(k,2)*pow(q,2))*cos(k*q)) + \
3*(-2 + pow(k,2)*pow(q,2))*sin(k*q)))/(5040.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) - (pow(f1,2)*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q*(((41184 + 6864*(-5 + 6*beta)*f1*pow(k,2)*X1q + \
6*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-525 + 1820*beta - \
2145*pow(beta,2) + 858*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-5 + 6*beta) + 2*(70 - 165*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(26*(70 - 165*beta + 99*pow(beta,2)) \
+ (-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q) + \
4*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(70 - 165*beta + 99*pow(beta,2)) \
+ 3*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(540540.*pow(E,(beta*f1*(2 \
+ f1)*pow(k,2)*X1q)/2.)*k*q) + ((61776 + 624*(-85 + \
99*beta)*f1*pow(k,2)*X1q + 6*(-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-861 + 2905*beta - \
3315*pow(beta,2) + 1287*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(39*(-85 + 99*beta) + (2905 - 6630*beta + \
3861*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(2905 - 6630*beta + 3861*pow(beta,2) \
+ (-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(2905 - 6630*beta + 3861*pow(beta,2) \
+ 6*(-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(540540.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/8. - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*(((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(124740.*pow(E,(beta*f1*(2 \
+ f1)*pow(k,2)*X1q)/2.)*k*q) + ((226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(1.62162e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 1:
				return (-10*Pi*pow(q,2)*(-((1 + (pow(k,2)*X1Infinity)/2.)*(11088 + 528*(-11 \
+ 21*beta)*f1*pow(k,2)*X1Infinity + 6*(-85 + 297*beta - \
363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-85 + \
297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
264*pow(f1,2)*pow(k,2)*X1Infinity*(-11 + 21*beta + (9 - 22*beta + \
21*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(33*(9 - 22*beta + \
21*pow(beta,2)) + (-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(11*(9 - 22*beta + \
21*pow(beta,2)) + 2*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(27720.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,3)) - \
((4*f1*pow(k,2)*X1Infinity + 2*pow(f1,2)*pow(k,2)*X1Infinity)*(226512 \
+ 20592*(-9 + 11*beta)*f1*pow(k,2)*X1Infinity + 6*(-2905 + 9945*beta \
- 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-2905 + \
9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
936*pow(f1,2)*pow(k,2)*X1Infinity*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(39*(85 - 198*beta + \
121*pow(beta,2)) + 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(4.32432e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,3))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (10*Pi*pow(q,2)*(((-12 + \
pow(k,2)*pow(q,2))*X1q*(11088 + 528*(-11 + 21*beta)*f1*pow(k,2)*X1q + \
6*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-85 + 297*beta - \
363*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
264*pow(f1,2)*pow(k,2)*X1q*(-11 + 21*beta + (9 - 22*beta + \
21*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(33*(9 - 22*beta + 21*pow(beta,2)) + \
(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(11*(9 - 22*beta + 21*pow(beta,2)) + \
2*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*Y1q*(k*q*(-15 + \
pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 - \
2*pow(k,2)*pow(q,2))*sin(k*q)))/(110880.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,6)) + ((-12 + \
pow(k,2)*pow(q,2))*(4*f1*X1q + 2*pow(f1,2)*X1q)*(226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*Y1q*(k*q*(-15 + \
pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 - \
2*pow(k,2)*pow(q,2))*sin(k*q)))/(8.64864e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,6)) - ((226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(2.16216e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((11088 + 528*(-11 + \
21*beta)*f1*pow(k,2)*X1q + 6*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-85 + 297*beta - \
363*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
264*pow(f1,2)*pow(k,2)*X1q*(-11 + 21*beta + (9 - 22*beta + \
21*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(33*(9 - 22*beta + 21*pow(beta,2)) + \
(-85 + 297*beta - 363*pow(beta,2) + 231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(11*(9 - 22*beta + 21*pow(beta,2)) + \
2*(-85 + 297*beta - 363*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(1 + (pow(k,2)*X1q)/2. + \
((alpha*pow(k,2))/2. + (alpha*pow(k,4)*X1q)/4.)*Y1q)*(3*k*q*cos(k*q) \
+ (-3 + pow(k,2)*pow(q,2))*sin(k*q)))/(27720.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*(-((11088 + 1584*(-6 + \
7*beta)*f1*pow(k,2)*X1q + 6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(311850.*pow(E,(beta*f1*(2 \
+ f1)*pow(k,2)*X1q)/2.)*k*q) - ((226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(2.06388e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (2*(6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(525525.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-((41184 + \
6864*(-5 + 6*beta)*f1*pow(k,2)*X1q + 6*(-525 + 1820*beta - \
2145*pow(beta,2) + 858*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-5 + 6*beta) + 2*(70 - 165*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(26*(70 - 165*beta + 99*pow(beta,2)) \
+ (-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q) + \
4*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(70 - 165*beta + 99*pow(beta,2)) \
+ 3*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(1.35135e6*pow(E,(beta*f1*(\
2 + f1)*pow(k,2)*X1q)/2.)*k*q) - ((61776 + 624*(-85 + \
99*beta)*f1*pow(k,2)*X1q + 6*(-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-861 + 2905*beta - \
3315*pow(beta,2) + 1287*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(39*(-85 + 99*beta) + (2905 - 6630*beta + \
3861*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(2905 - 6630*beta + 3861*pow(beta,2) \
+ (-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(2905 - 6630*beta + 3861*pow(beta,2) \
+ 6*(-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(687960.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (2*(360672 + 816*(-435 + \
442*beta)*f1*pow(k,2)*X1q + 6*(-6825 + 21420*beta - 22185*pow(beta,2) \
+ 7514*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-6825 + \
21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-435 + 442*beta + 2*(210 - 435*beta + \
221*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.6801775e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/8. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-((41184 + 6864*(-5 + \
6*beta)*f1*pow(k,2)*X1q + 6*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-525 + 1820*beta - \
2145*pow(beta,2) + 858*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-5 + 6*beta) + 2*(70 - 165*beta + \
99*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(26*(70 - 165*beta + 99*pow(beta,2)) \
+ (-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q) + \
4*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(70 - 165*beta + 99*pow(beta,2)) \
+ 3*(-525 + 1820*beta - 2145*pow(beta,2) + \
858*pow(beta,3))*pow(k,2)*X1q))*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(450450.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((48048 + 208*(-205 + \
231*beta)*f1*pow(k,2)*X1q + 6*(-735 + 2415*beta - 2665*pow(beta,2) + \
1001*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-735 + 2415*beta - \
2665*pow(beta,2) + 1001*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-2665 + 3003*beta + (2415 - 5330*beta + \
3003*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(2415 - 5330*beta + 3003*pow(beta,2) \
+ (-735 + 2415*beta - 2665*pow(beta,2) + \
1001*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(2415 - 5330*beta + 3003*pow(beta,2) \
+ 6*(-735 + 2415*beta - 2665*pow(beta,2) + \
1001*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(600600.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/4. - \
(f1*pow(k,4)*X1q*Y1q*(-((11088 + 1584*(-6 + 7*beta)*f1*pow(k,2)*X1q + \
6*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 495*beta - \
594*pow(beta,2) + 231*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
792*pow(f1,2)*pow(k,2)*X1q*(-6 + 7*beta + (5 - 12*beta + \
7*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(99*(5 \
- 12*beta + 7*pow(beta,2)) + (-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(33*(5 - 12*beta + 7*pow(beta,2)) + \
2*(-140 + 495*beta - 594*pow(beta,2) + \
231*pow(beta,3))*pow(k,2)*X1q))*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(103950.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((61776 + 6864*(-7 + \
9*beta)*f1*pow(k,2)*X1q + 6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-805 + 2665*beta - \
3003*pow(beta,2) + 1287*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-7 + 9*beta) + (205 - 462*beta + \
297*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + (-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + 6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(600600.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/2.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 2:
				return (-10*Pi*pow(q,2)*(((4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1Infinity + 6*(-140 + 435*beta - \
442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-140 + \
435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
8*pow(f1,2)*pow(k,2)*X1Infinity*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(435 - 884*beta + \
429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1Infinity) + \
2*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(435 - 884*beta + \
429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1Infinity))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(360360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,5)*pow(q,5)) + (f1*(2 + \
f1)*X1Infinity*(1 + (pow(k,2)*X1Infinity)/2.)*(-3432 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,5))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (10*Pi*pow(q,2)*(((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(5*k*q*(-21 + 2*pow(k,2)*pow(q,2))*cos(k*q) \
+ (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (f1*(2 + f1)*X1q*(-3432 + \
f1*(2 + f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(1 + (pow(k,2)*X1q)/2. + \
((alpha*pow(k,2))/2. + (alpha*pow(k,4)*X1q)/4.)*Y1q)*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,5)) + ((4*f1*X1q + \
2*pow(f1,2)*X1q)*(6864 + 208*(-34 + 33*beta)*f1*pow(k,2)*X1q + \
6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*Y1q*(-(k*q*(-22050 + \
3225*pow(k,2)*pow(q,2) - 111*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + 10575*pow(k,2)*pow(q,2) - \
696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(720720.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,8)) + (f1*(2 + \
f1)*pow(X1q,2)*(-3432 + f1*(2 + f1)*pow(k,2)*X1q*(884 - 1716*beta + \
(-145 + 13*(34 - 33*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(-(k*q*(-22050 + 3225*pow(k,2)*pow(q,2) - \
111*pow(k,4)*pow(q,4) + pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + \
10575*pow(k,2)*pow(q,2) - 696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(360360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,8)) - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(((61776 + \
624*(-85 + 99*beta)*f1*pow(k,2)*X1q + 6*(-861 + 2905*beta - \
3315*pow(beta,2) + 1287*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(39*(-85 + 99*beta) + (2905 - 6630*beta + \
3861*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(2905 - 6630*beta + 3861*pow(beta,2) \
+ (-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(2905 - 6630*beta + 3861*pow(beta,2) \
+ 6*(-861 + 2905*beta - 3315*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(1.89189e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + ((360672 + 816*(-435 + \
442*beta)*f1*pow(k,2)*X1q + 6*(-6825 + 21420*beta - 22185*pow(beta,2) \
+ 7514*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-6825 + \
21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-435 + 442*beta + 2*(210 - 435*beta + \
221*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(9.07137e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (2*(232560 + 15504*(-23 + \
15*beta)*f1*pow(k,2)*X1q + 6*(-9870 + 26866*beta - 22287*pow(beta,2) \
+ 4845*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-9870 + \
26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
152*pow(f1,2)*pow(k,2)*X1q*(51*(-23 + 15*beta) + (1414 - 2346*beta + \
765*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(19*(1414 - 2346*beta + \
765*pow(beta,2)) + (-9870 + 26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(19*(1414 - 2346*beta + \
765*pow(beta,2)) + 6*(-9870 + 26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9.6026931e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7))))/8. - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*(((226512 + 20592*(-9 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2905 + 9945*beta \
- 11583*pow(beta,2) + 4719*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 936*pow(f1,2)*pow(k,2)*X1q*(-99 + 121*beta + (85 - 198*beta + \
121*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(117*(85 - 198*beta + \
121*pow(beta,2)) + (-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(39*(85 - 198*beta + 121*pow(beta,2)) \
+ 2*(-2905 + 9945*beta - 11583*pow(beta,2) + \
4719*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(5.67567e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + ((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(177870.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (2*f1*(2 + f1)*X1q*(-18360 \
+ f1*(2 + f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 + 153*(23 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.5162147e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7))))/4. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(((48048 + 208*(-205 + \
231*beta)*f1*pow(k,2)*X1q + 6*(-735 + 2415*beta - 2665*pow(beta,2) + \
1001*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-735 + 2415*beta - \
2665*pow(beta,2) + 1001*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-2665 + 3003*beta + (2415 - 5330*beta + \
3003*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(2415 - 5330*beta + 3003*pow(beta,2) \
+ (-735 + 2415*beta - 2665*pow(beta,2) + \
1001*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(2415 - 5330*beta + 3003*pow(beta,2) \
+ 6*(-735 + 2415*beta - 2665*pow(beta,2) + \
1001*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(810810.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + ((53040 + 1360*(-50 + \
39*beta)*f1*pow(k,2)*X1q + 6*(-1596 + 4641*beta - 4250*pow(beta,2) + \
1105*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-1596 + 4641*beta \
- 4250*pow(beta,2) + 1105*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 136*pow(f1,2)*pow(k,2)*X1q*(5*(-50 + 39*beta) + (273 - 500*beta + \
195*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(17*(273 - 500*beta + \
195*pow(beta,2)) + (-1596 + 4641*beta - 4250*pow(beta,2) + \
1105*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(273 - 500*beta + \
195*pow(beta,2)) + 6*(-1596 + 4641*beta - 4250*pow(beta,2) + \
1105*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.756754e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7))))/4. - \
(f1*pow(k,4)*X1q*Y1q*(((61776 + 6864*(-7 + 9*beta)*f1*pow(k,2)*X1q + \
6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-805 + 2665*beta - \
3003*pow(beta,2) + 1287*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
104*pow(f1,2)*pow(k,2)*X1q*(33*(-7 + 9*beta) + (205 - 462*beta + \
297*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + (-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(13*(205 - 462*beta + \
297*pow(beta,2)) + 6*(-805 + 2665*beta - 3003*pow(beta,2) + \
1287*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(810810.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (f1*(2 + f1)*X1q*(-1560 + \
f1*(2 + f1)*pow(k,2)*X1q*(500 - 780*beta + (-91 + 5*(50 - \
39*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(5670 - \
735*pow(k,2)*pow(q,2) + 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + \
2625*pow(k,2)*pow(q,2) - 135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(162162.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7))))/2.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 3:
				return (-10*Pi*pow(q,2)*(-(pow(f1,2)*pow(2 + f1,2)*pow(X1Infinity,2)*(1 + \
(pow(k,2)*X1Infinity)/2.)*(90 + (-23 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(45045.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,7)) - (f1*(2 + \
f1)*X1Infinity*(4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(-18360 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(7038 - 9180*beta + (-1414 + 153*(23 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9.18918e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,5)*pow(q,7))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (10*Pi*pow(q,2)*(-(f1*(2 + f1)*X1q*(-18360 + f1*(2 \
+ f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 + 153*(23 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(4.59459e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(90 + (-23 + 45*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(1 + \
(pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(45045.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,7)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,3)*(90 + (-23 + 45*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(k*q*(5239080 - 773955*pow(k,2)*pow(q,2) + \
27909*pow(k,4)*pow(q,4) - 342*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - 840105*pow(k,2)*pow(q,2) + \
56490*pow(k,4)*pow(q,4) - 1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,10)) + (f1*(2 + \
f1)*X1q*(4*f1*X1q + 2*pow(f1,2)*X1q)*(-18360 + f1*(2 + \
f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 + 153*(23 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(k*q*(5239080 - \
773955*pow(k,2)*pow(q,2) + 27909*pow(k,4)*pow(q,4) - \
342*pow(k,6)*pow(q,6) + pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - \
840105*pow(k,2)*pow(q,2) + 56490*pow(k,4)*pow(q,4) - \
1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(1.837836e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,10)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*(-((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(429429.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) - (83*f1*(2 + \
f1)*X1q*(-18360 + f1*(2 + f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 \
+ 153*(23 - 15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(3.79053675e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (448*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(342 + (-116 + 171*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.215925425e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-((360672 + \
816*(-435 + 442*beta)*f1*pow(k,2)*X1q + 6*(-6825 + 21420*beta - \
22185*pow(beta,2) + 7514*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-435 + 442*beta + 2*(210 - 435*beta + \
221*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.1900879e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) - (83*(232560 + 15504*(-23 + \
15*beta)*f1*pow(k,2)*X1q + 6*(-9870 + 26866*beta - 22287*pow(beta,2) \
+ 4845*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-9870 + \
26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
152*pow(f1,2)*pow(k,2)*X1q*(51*(-23 + 15*beta) + (1414 - 2346*beta + \
765*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(19*(1414 - 2346*beta + \
765*pow(beta,2)) + (-9870 + 26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(19*(1414 - 2346*beta + \
765*pow(beta,2)) + 6*(-9870 + 26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.400673275e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (64*f1*(2 + f1)*X1q*(-3192 \
+ f1*(2 + f1)*pow(k,2)*X1q*(28*(58 - 57*beta) + (-390 + 7*(116 - \
57*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.05308475e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9))))/8. - \
(f1*pow(k,4)*X1q*Y1q*(-(f1*(2 + f1)*X1q*(-1560 + f1*(2 + \
f1)*pow(k,2)*X1q*(500 - 780*beta + (-91 + 5*(50 - \
39*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(5670 - \
735*pow(k,2)*pow(q,2) + 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + \
2625*pow(k,2)*pow(q,2) - 135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(195195.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) - (7*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(306 + (-91 + 153*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.266405e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-((53040 + 1360*(-50 + \
39*beta)*f1*pow(k,2)*X1q + 6*(-1596 + 4641*beta - 4250*pow(beta,2) + \
1105*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-1596 + 4641*beta \
- 4250*pow(beta,2) + 1105*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 136*pow(f1,2)*pow(k,2)*X1q*(5*(-50 + 39*beta) + (273 - 500*beta + \
195*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(17*(273 - 500*beta + \
195*pow(beta,2)) + (-1596 + 4641*beta - 4250*pow(beta,2) + \
1105*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(273 - 500*beta + \
195*pow(beta,2)) + 6*(-1596 + 4641*beta - 4250*pow(beta,2) + \
1105*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(3.318315e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (7*f1*(2 + f1)*X1q*(-7752 \
+ f1*(2 + f1)*pow(k,2)*X1q*(3458 - 3876*beta + (-762 + 19*(91 - \
51*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(2.7020565e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 4:
				return (-10*Pi*pow(q,2)*((-8*pow(f1,3)*pow(2 + f1,3)*pow(X1Infinity,3)*(1 + \
(pow(k,2)*X1Infinity)/2.)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(109395.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,9)) + (2*pow(f1,2)*pow(2 \
+ f1,2)*pow(X1Infinity,2)*(4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(342 + (-116 + 171*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(6.235515e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,5)*pow(q,9))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (10*Pi*pow(q,2)*((4*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(342 + (-116 + 171*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*((4*f1*pow(k,2)*X1q + 2*pow(f1,2)*pow(k,2)*X1q)/2. \
+ (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(6.235515e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) - (8*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(1 + (pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(109395.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,9)) - (2*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,4)*Y1q*(-(k*q*(-2006754750 + \
300405105*pow(k,2)*pow(q,2) - 11320155*pow(k,4)*pow(q,4) + \
158679*pow(k,6)*pow(q,6) - 852*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*cos(k*q)) + 3*(-668918250 + \
323107785*pow(k,2)*pow(q,2) - 22286880*pow(k,4)*pow(q,4) + \
501165*pow(k,6)*pow(q,6) - 4418*pow(k,8)*pow(q,8) + \
13*pow(k,10)*pow(q,10))*sin(k*q)))/(109395.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,12)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(4*f1*X1q + 2*pow(f1,2)*X1q)*(342 + (-116 + \
171*beta)*f1*(2 + f1)*pow(k,2)*X1q)*Y1q*(-(k*q*(-2006754750 + \
300405105*pow(k,2)*pow(q,2) - 11320155*pow(k,4)*pow(q,4) + \
158679*pow(k,6)*pow(q,6) - 852*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*cos(k*q)) + 3*(-668918250 + \
323107785*pow(k,2)*pow(q,2) - 22286880*pow(k,4)*pow(q,4) + \
501165*pow(k,6)*pow(q,6) - 4418*pow(k,8)*pow(q,8) + \
13*pow(k,10)*pow(q,10))*sin(k*q)))/(6.235515e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,12)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*((8*f1*(2 + f1)*X1q*(-18360 + f1*(2 + \
f1)*pow(k,2)*X1q*(7038 - 9180*beta + (-1414 + 153*(23 - \
15*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(8.3687175e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) + (8*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(342 + (-116 + 171*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.2427425e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) - (480*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(55*k*q*(11904165 - 1670760*pow(k,2)*pow(q,2) + \
51597*pow(k,4)*pow(q,4) - 468*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) + (-654729075 + \
310134825*pow(k,2)*pow(q,2) - 18918900*pow(k,4)*pow(q,4) + \
315315*pow(k,6)*pow(q,6) - 1485*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(1.04433329e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,11))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((8*(232560 + \
15504*(-23 + 15*beta)*f1*pow(k,2)*X1q + 6*(-9870 + 26866*beta - \
22287*pow(beta,2) + 4845*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-9870 + 26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
152*pow(f1,2)*pow(k,2)*X1q*(51*(-23 + 15*beta) + (1414 - 2346*beta + \
765*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(19*(1414 - 2346*beta + \
765*pow(beta,2)) + (-9870 + 26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(19*(1414 - 2346*beta + \
765*pow(beta,2)) + 6*(-9870 + 26866*beta - 22287*pow(beta,2) + \
4845*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(5.30018775e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) + (8*f1*(2 + f1)*X1q*(-3192 \
+ f1*(2 + f1)*pow(k,2)*X1q*(28*(58 - 57*beta) + (-390 + 7*(116 - \
57*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(2.8997325e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)) + (160*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(414 + (-175 + 207*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(55*k*q*(11904165 - 1670760*pow(k,2)*pow(q,2) + \
51597*pow(k,4)*pow(q,4) - 468*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) + (-654729075 + \
310134825*pow(k,2)*pow(q,2) - 18918900*pow(k,4)*pow(q,4) + \
315315*pow(k,6)*pow(q,6) - 1485*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(2.401966567e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,11))))/8. - \
(f1*pow(k,4)*X1q*Y1q*((8*pow(f1,2)*pow(2 + f1,2)*pow(X1q,2)*(306 + \
(-91 + 153*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(5.579145e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) - (72*pow(f1,3)*pow(2 + \
f1,3)*pow(X1q,3)*(k*q*(344594250 - 49324275*pow(k,2)*pow(q,2) + \
1621620*pow(k,4)*pow(q,4) - 16830*pow(k,6)*pow(q,6) + \
46*pow(k,8)*pow(q,8))*cos(k*q) + (-344594250 + \
164189025*pow(k,2)*pow(q,2) - 10405395*pow(k,4)*pow(q,4) + \
190575*pow(k,6)*pow(q,6) - 1080*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(3.926065e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,11))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((8*f1*(2 + \
f1)*X1q*(-7752 + f1*(2 + f1)*pow(k,2)*X1q*(3458 - 3876*beta + (-762 + \
19*(91 - 51*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(3.5334585e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)) + (72*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(42 + (-16 + 21*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(344594250 - 49324275*pow(k,2)*pow(q,2) + \
1621620*pow(k,4)*pow(q,4) - 16830*pow(k,6)*pow(q,6) + \
46*pow(k,8)*pow(q,8))*cos(k*q) + (-344594250 + \
164189025*pow(k,2)*pow(q,2) - 10405395*pow(k,4)*pow(q,4) + \
190575*pow(k,6)*pow(q,6) - 1080*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(2.7482455e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,11))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

		}
		break ;

	case 2:
		switch (j) { 
			case 0:
				return (-18*Pi*pow(q,2)*((pow(f1,2)*pow(2 + \
f1,2)*pow(k,3)*pow(X1Infinity,2)*(1 + (pow(k,2)*X1Infinity)/2.)*(22 + \
(-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*sin(k*q))/(3465.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*q) + (f1*(2 + \
f1)*k*X1Infinity*(4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity))*sin(k*q))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*q)))/pow(E,(pow(k,2)*X1Infinity)/2.) + \
(18*Pi*pow(q,2)*((f1*(2 + f1)*k*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*sin(k*q))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) + (pow(f1,2)*pow(2 + \
f1,2)*pow(k,3)*pow(X1q,2)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(1 + (pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*sin(k*q))/(3465.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) - (f1*k*(4*f1*X1q + 2*pow(f1,2)*X1q)*(6864 + \
624*(-15 + 11*beta)*f1*pow(k,2)*X1q + 6*(-210 + 630*beta - \
585*pow(beta,2) + 143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-195 + 143*beta + (210 - 390*beta + \
143*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(630 \
- 1170*beta + 429*pow(beta,2) + (-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(210 \
- 390*beta + 143*pow(beta,2) + 2*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(540540.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) - (pow(f1,2)*(2 + \
f1)*pow(k,3)*pow(X1q,2)*(-1144 + f1*(2 + f1)*pow(k,2)*X1q*(390 - \
572*beta + (-70 + 13*(15 - 11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,3)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(k,4)*pow(X1q,3)*(22 + (-5 + 11*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(-(k*q*(-6 + pow(k,2)*pow(q,2))*cos(k*q)) + \
3*(-2 + pow(k,2)*pow(q,2))*sin(k*q)))/(13860.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) + (f1*(2 + f1)*pow(k,2)*X1q*(4*f1*X1q \
+ 2*pow(f1,2)*X1q)*(-1144 + f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta \
+ (-70 + 13*(15 - 11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(-(k*q*(-6 + pow(k,2)*pow(q,2))*cos(k*q)) + \
3*(-2 + pow(k,2)*pow(q,2))*sin(k*q)))/(360360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,4)) - (pow(f1,2)*pow(k,4)*X1q*Y1q*((f1*(2 \
+ f1)*k*X1q*(-1144 + f1*(2 + f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 \
+ 13*(15 - 11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*sin(k*q))/(135135.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) + ((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(135135.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(((6864 + \
624*(-15 + 11*beta)*f1*pow(k,2)*X1q + 6*(-210 + 630*beta - \
585*pow(beta,2) + 143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-195 + 143*beta + (210 - 390*beta + \
143*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(630 \
- 1170*beta + 429*pow(beta,2) + (-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(210 \
- 390*beta + 143*pow(beta,2) + 2*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(405405.*pow(E,(beta*f1*(2 \
+ f1)*pow(k,2)*X1q)/2.)*k*q) + ((360672 + 816*(-435 + \
442*beta)*f1*pow(k,2)*X1q + 6*(-6825 + 21420*beta - 22185*pow(beta,2) \
+ 7514*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-6825 + \
21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-435 + 442*beta + 2*(210 - 435*beta + \
221*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(6.891885e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3))))/8.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 1:
				return (-18*Pi*pow(q,2)*(-((4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1Infinity + 6*(-140 + 435*beta - \
442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-140 + \
435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
8*pow(f1,2)*pow(k,2)*X1Infinity*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(435 - 884*beta + \
429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1Infinity) + \
2*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(435 - 884*beta + \
429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(360360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,3)*pow(q,3)) - (f1*(2 + \
f1)*X1Infinity*(1 + (pow(k,2)*X1Infinity)/2.)*(-3432 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(3*k*q*cos(k*q) + \
(-3 + pow(k,2)*pow(q,2))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*k*pow(q,3))))/pow(E,(pow(k,2)*X1Infinity)\
/2.) + (18*Pi*pow(q,2)*(((-12 + pow(k,2)*pow(q,2))*(4*f1*X1q + \
2*pow(f1,2)*X1q)*(6864 + 208*(-34 + 33*beta)*f1*pow(k,2)*X1q + \
6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*Y1q*(k*q*(-15 + \
pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 - \
2*pow(k,2)*pow(q,2))*sin(k*q)))/(720720.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,2)*pow(q,6)) + (f1*(2 + f1)*(-12 + \
pow(k,2)*pow(q,2))*pow(X1q,2)*(-3432 + f1*(2 + f1)*pow(k,2)*X1q*(884 \
- 1716*beta + (-145 + 13*(34 - 33*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(k*q*(-15 + pow(k,2)*pow(q,2))*cos(k*q) + 3*(5 \
- 2*pow(k,2)*pow(q,2))*sin(k*q)))/(360360.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(q,6)) - ((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(180180.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - (f1*(2 + f1)*X1q*(-3432 + \
f1*(2 + f1)*pow(k,2)*X1q*(884 - 1716*beta + (-145 + 13*(34 - \
33*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(1 + (pow(k,2)*X1q)/2. + \
((alpha*pow(k,2))/2. + (alpha*pow(k,4)*X1q)/4.)*Y1q)*(3*k*q*cos(k*q) \
+ (-3 + pow(k,2)*pow(q,2))*sin(k*q)))/(90090.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*((-2*f1*(2 + f1)*k*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*sin(k*q))/(675675.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*q) - ((6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(171990.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((689520 + 816*(-643 + \
845*beta)*f1*pow(k,2)*X1q + 6*(-8295 + 27795*beta - 32793*pow(beta,2) \
+ 14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.786785e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((-2*(6864 + \
624*(-15 + 11*beta)*f1*pow(k,2)*X1q + 6*(-210 + 630*beta - \
585*pow(beta,2) + 143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-195 + 143*beta + (210 - 390*beta + \
143*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(630 \
- 1170*beta + 429*pow(beta,2) + (-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(210 \
- 390*beta + 143*pow(beta,2) + 2*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*sin(k*q))/(2.027025e6*pow(E,(beta*f1*(\
2 + f1)*pow(k,2)*X1q)/2.)*k*q) - ((360672 + 816*(-435 + \
442*beta)*f1*pow(k,2)*X1q + 6*(-6825 + 21420*beta - 22185*pow(beta,2) \
+ 7514*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-6825 + \
21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-435 + 442*beta + 2*(210 - 435*beta + \
221*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(8.77149e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((9969072 + 15504*(-545 + \
643*beta)*f1*pow(k,2)*X1q + 6*(-145215 + 472815*beta - \
528105*pow(beta,2) + \
207689*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-145215 + \
472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
456*pow(f1,2)*pow(k,2)*X1q*(17*(-545 + 643*beta) + (8295 - 18530*beta \
+ 10931*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(57*(8295 - 18530*beta + \
10931*pow(beta,2)) + (-145215 + 472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(19*(8295 - 18530*beta + \
10931*pow(beta,2)) + 2*(-145215 + 472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(3.3948915e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/8. - \
(f1*pow(k,4)*X1q*Y1q*((-2*f1*(2 + f1)*X1q*(-1144 + f1*(2 + \
f1)*pow(k,2)*X1q*(390 - 572*beta + (-70 + 13*(15 - \
11*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(225225.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*k*pow(q,3)) - ((34320 + 624*(-42 + \
55*beta)*f1*pow(k,2)*X1q + 6*(-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-546 + 715*beta + (465 - 1092*beta + \
715*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(465 - 1092*beta + \
715*pow(beta,2)) + (-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(465 \
- 1092*beta + 715*pow(beta,2) + 2*(-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(450450.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((-2*(6864 + 624*(-15 + \
11*beta)*f1*pow(k,2)*X1q + 6*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-210 + 630*beta - \
585*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-195 + 143*beta + (210 - 390*beta + \
143*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(630 \
- 1170*beta + 429*pow(beta,2) + (-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(210 \
- 390*beta + 143*pow(beta,2) + 2*(-210 + 630*beta - 585*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(2*k*q*cos(k*q) + (-2 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(675675.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) - ((148512 + 816*(-155 + \
182*beta)*f1*pow(k,2)*X1q + 6*(-2205 + 7140*beta - 7905*pow(beta,2) + \
3094*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2205 + 7140*beta \
- 7905*pow(beta,2) + 3094*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 408*pow(f1,2)*pow(k,2)*X1q*(-155 + 182*beta + 2*(70 - 155*beta + \
91*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(70 - 155*beta + 91*pow(beta,2)) \
+ (-2205 + 7140*beta - 7905*pow(beta,2) + \
3094*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(70 - 155*beta + 91*pow(beta,2)) \
+ (-2205 + 7140*beta - 7905*pow(beta,2) + \
3094*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.55255e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 2:
				return (-18*Pi*pow(q,2)*(((1 + (pow(k,2)*X1Infinity)/2.)*(240240 + \
3120*(-39 + 77*beta)*f1*pow(k,2)*X1Infinity + 6*(-1635 + 5787*beta - \
7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-1635 + \
5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
24*pow(f1,2)*pow(k,2)*X1Infinity*(65*(-39 + 77*beta) + (1929 - \
5070*beta + 5005*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(3*(1929 - 5070*beta + \
5005*pow(beta,2)) + (-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(1929 - 5070*beta + \
5005*pow(beta,2) + 2*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1Infinity))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.08108e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,5)*pow(q,5)) + \
((4*f1*pow(k,2)*X1Infinity + 2*pow(f1,2)*pow(k,2)*X1Infinity)*(689520 \
+ 816*(-643 + 845*beta)*f1*pow(k,2)*X1Infinity + 6*(-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
408*pow(f1,2)*pow(k,2)*X1Infinity*(-643 + 845*beta + (545 - 1286*beta \
+ 845*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1Infinity) + \
6*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1Infinity))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.450448e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,5)*pow(q,5))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (18*Pi*pow(q,2)*(((689520 + 816*(-643 + \
845*beta)*f1*pow(k,2)*X1q + 6*(-8295 + 27795*beta - 32793*pow(beta,2) \
+ 14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(5*k*q*(-21 + 2*pow(k,2)*pow(q,2))*cos(k*q) \
+ (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.225224e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + ((240240 + 3120*(-39 + \
77*beta)*f1*pow(k,2)*X1q + 6*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-1635 + 5787*beta \
- 7605*pow(beta,2) + 5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 24*pow(f1,2)*pow(k,2)*X1q*(65*(-39 + 77*beta) + (1929 - 5070*beta + \
5005*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(1929 - 5070*beta + \
5005*pow(beta,2)) + (-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(1929 - 5070*beta + 5005*pow(beta,2) \
+ 2*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q))*(1 + (pow(k,2)*X1q)/2. + \
((alpha*pow(k,2))/2. + (alpha*pow(k,4)*X1q)/4.)*Y1q)*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.08108e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (X1q*(240240 + 3120*(-39 + \
77*beta)*f1*pow(k,2)*X1q + 6*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-1635 + 5787*beta \
- 7605*pow(beta,2) + 5005*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) \
+ 24*pow(f1,2)*pow(k,2)*X1q*(65*(-39 + 77*beta) + (1929 - 5070*beta + \
5005*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(1929 - 5070*beta + \
5005*pow(beta,2)) + (-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(1929 - 5070*beta + 5005*pow(beta,2) \
+ 2*(-1635 + 5787*beta - 7605*pow(beta,2) + \
5005*pow(beta,3))*pow(k,2)*X1q))*Y1q*(-(k*q*(-22050 + \
3225*pow(k,2)*pow(q,2) - 111*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + 10575*pow(k,2)*pow(q,2) - \
696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(4.32432e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,8)) + ((4*f1*X1q + \
2*pow(f1,2)*X1q)*(689520 + 816*(-643 + 845*beta)*f1*pow(k,2)*X1q + \
6*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*Y1q*(-(k*q*(-22050 + \
3225*pow(k,2)*pow(q,2) - 111*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*cos(k*q)) + (-22050 + 10575*pow(k,2)*pow(q,2) - \
696*pow(k,4)*pow(q,4) + \
13*pow(k,6)*pow(q,6))*sin(k*q)))/(4.900896e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,8)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*((2*(6864 + 208*(-34 + \
33*beta)*f1*pow(k,2)*X1q + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-140 + 435*beta - \
442*pow(beta,2) + 143*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
8*pow(f1,2)*pow(k,2)*X1q*(-442 + 429*beta + (435 - 884*beta + \
429*pow(beta,2))*pow(k,2)*X1q) + 8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + (-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q) + 2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(435 \
- 884*beta + 429*pow(beta,2) + 6*(-140 + 435*beta - 442*pow(beta,2) + \
143*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(945945.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + ((689520 + 816*(-643 + \
845*beta)*f1*pow(k,2)*X1q + 6*(-8295 + 27795*beta - 32793*pow(beta,2) \
+ 14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.209516e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + ((2713200 + 25840*(-106 + \
105*beta)*f1*pow(k,2)*X1q + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + 56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) \
+ (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.92053862e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((2*(360672 + \
816*(-435 + 442*beta)*f1*pow(k,2)*X1q + 6*(-6825 + 21420*beta - \
22185*pow(beta,2) + 7514*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + \
(-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-435 + 442*beta + 2*(210 - 435*beta + \
221*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(210 - 435*beta + \
221*pow(beta,2)) + (-6825 + 21420*beta - 22185*pow(beta,2) + \
7514*pow(beta,3))*pow(k,2)*X1q))*(3*k*q*cos(k*q) + (-3 + \
pow(k,2)*pow(q,2))*sin(k*q)))/(4.8243195e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,3)*pow(q,3)) + ((9969072 + 15504*(-545 + \
643*beta)*f1*pow(k,2)*X1q + 6*(-145215 + 472815*beta - \
528105*pow(beta,2) + \
207689*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-145215 + \
472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
456*pow(f1,2)*pow(k,2)*X1q*(17*(-545 + 643*beta) + (8295 - 18530*beta \
+ 10931*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(57*(8295 - 18530*beta + \
10931*pow(beta,2)) + (-145215 + 472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(19*(8295 - 18530*beta + \
10931*pow(beta,2)) + 2*(-145215 + 472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.2980804e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + ((8217120 + 4560*(-1707 + \
1802*beta)*f1*pow(k,2)*X1q + 6*(-145635 + 459396*beta - \
486495*pow(beta,2) + \
171190*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-145635 + \
459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(95*(-1707 + 1802*beta) + 2*(76566 - \
162165*beta + 85595*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(76566 - 162165*beta + \
85595*pow(beta,2) + (-145635 + 459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(6*(76566 - 162165*beta + \
85595*pow(beta,2)) + (-145635 + 459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(5.76161586e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7))))/8. - \
(f1*pow(k,4)*X1q*Y1q*((2*(34320 + 624*(-42 + 55*beta)*f1*pow(k,2)*X1q \
+ 6*(-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(-546 + 715*beta + (465 - 1092*beta + \
715*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(3*(465 - 1092*beta + \
715*pow(beta,2)) + (-420 + 1395*beta - 1638*pow(beta,2) + \
715*pow(beta,3))*pow(k,2)*X1q) + 6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(465 \
- 1092*beta + 715*pow(beta,2) + 2*(-420 + 1395*beta - \
1638*pow(beta,2) + 715*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(1.216215e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + (5*(371280 + 4080*(-69 + \
91*beta)*f1*pow(k,2)*X1q + 6*(-4389 + 14841*beta - 17595*pow(beta,2) \
+ 7735*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-4389 + \
14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(5*(-69 + 91*beta) + (291 - 690*beta + \
455*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(291 - 690*beta + \
455*pow(beta,2)) + (-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(291 - 690*beta + \
455*pow(beta,2)) + 2*(-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(3.3081048e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((2*(148512 + 816*(-155 \
+ 182*beta)*f1*pow(k,2)*X1q + 6*(-2205 + 7140*beta - 7905*pow(beta,2) \
+ 3094*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-2205 + \
7140*beta - 7905*pow(beta,2) + \
3094*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-155 + 182*beta + 2*(70 - 155*beta + \
91*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(70 - 155*beta + 91*pow(beta,2)) \
+ (-2205 + 7140*beta - 7905*pow(beta,2) + \
3094*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(102*(70 - 155*beta + 91*pow(beta,2)) \
+ (-2205 + 7140*beta - 7905*pow(beta,2) + \
3094*pow(beta,3))*pow(k,2)*X1q))*(k*q*(-60 + \
7*pow(k,2)*pow(q,2))*cos(k*q) + (60 - 27*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(6.891885e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) + ((8914800 + 77520*(-97 + \
115*beta)*f1*pow(k,2)*X1q + 6*(-127323 + 416955*beta - \
469965*pow(beta,2) + \
185725*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-127323 + \
416955*beta - 469965*pow(beta,2) + \
185725*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
2280*pow(f1,2)*pow(k,2)*X1q*(17*(-97 + 115*beta) + (1463 - 3298*beta \
+ 1955*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(285*(1463 - 3298*beta + \
1955*pow(beta,2)) + (-127323 + 416955*beta - 469965*pow(beta,2) + \
185725*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1463 - 3298*beta + \
1955*pow(beta,2)) + 2*(-127323 + 416955*beta - 469965*pow(beta,2) + \
185725*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) \
+ 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.09513304e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 3:
				return (-18*Pi*pow(q,2)*(-((4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(2713200 + 25840*(-106 + \
105*beta)*f1*pow(k,2)*X1Infinity + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1Infinity,3) + (-51044 + \
162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1Infinity,3) + \
760*pow(f1,2)*pow(k,2)*X1Infinity*(-1802 + 1785*beta + (1707 - \
3604*beta + 1785*pow(beta,2))*pow(k,2)*X1Infinity) + \
8*pow(f1,3)*pow(k,4)*pow(X1Infinity,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1Infinity) + \
2*pow(f1,4)*pow(k,4)*pow(X1Infinity,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1Infinity))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.3279256e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,7)*pow(q,7)) - (f1*(2 + \
f1)*X1Infinity*(1 + (pow(k,2)*X1Infinity)/2.)*(-14280 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(3604 - 7140*beta + (-569 + 17*(106 - \
105*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(612612.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,5)*pow(q,7))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (18*Pi*pow(q,2)*(-((2713200 + 25840*(-106 + \
105*beta)*f1*pow(k,2)*X1q + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + 56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) \
+ (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.1639628e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (f1*(2 + f1)*X1q*(-14280 + \
f1*(2 + f1)*pow(k,2)*X1q*(3604 - 7140*beta + (-569 + 17*(106 - \
105*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(1 + (pow(k,2)*X1q)/2. + \
((alpha*pow(k,2))/2. + (alpha*pow(k,4)*X1q)/4.)*Y1q)*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(612612.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,7)) + ((4*f1*X1q + \
2*pow(f1,2)*X1q)*(2713200 + 25840*(-106 + 105*beta)*f1*pow(k,2)*X1q + \
6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-51044 + \
162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*Y1q*(k*q*(5239080 - \
773955*pow(k,2)*pow(q,2) + 27909*pow(k,4)*pow(q,4) - \
342*pow(k,6)*pow(q,6) + pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - \
840105*pow(k,2)*pow(q,2) + 56490*pow(k,4)*pow(q,4) - \
1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(4.6558512e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,6)*pow(q,10)) + (f1*(2 + \
f1)*pow(X1q,2)*(-14280 + f1*(2 + f1)*pow(k,2)*X1q*(3604 - 7140*beta + \
(-569 + 17*(106 - 105*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*Y1q*(k*q*(5239080 - 773955*pow(k,2)*pow(q,2) + \
27909*pow(k,4)*pow(q,4) - 342*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) - 3*(1746360 - 840105*pow(k,2)*pow(q,2) + \
56490*pow(k,4)*pow(q,4) - 1178*pow(k,6)*pow(q,6) + \
8*pow(k,8)*pow(q,8))*sin(k*q)))/(2.450448e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,10)) - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-((9969072 + \
15504*(-545 + 643*beta)*f1*pow(k,2)*X1q + 6*(-145215 + 472815*beta - \
528105*pow(beta,2) + \
207689*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-145215 + \
472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
456*pow(f1,2)*pow(k,2)*X1q*(17*(-545 + 643*beta) + (8295 - 18530*beta \
+ 10931*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(57*(8295 - 18530*beta + \
10931*pow(beta,2)) + (-145215 + 472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(19*(8295 - 18530*beta + \
10931*pow(beta,2)) + 2*(-145215 + 472815*beta - 528105*pow(beta,2) + \
207689*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(5.54822268e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) - (83*(8217120 + 4560*(-1707 \
+ 1802*beta)*f1*pow(k,2)*X1q + 6*(-145635 + 459396*beta - \
486495*pow(beta,2) + \
171190*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-145635 + \
459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(95*(-1707 + 1802*beta) + 2*(76566 - \
162165*beta + 85595*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(76566 - 162165*beta + \
85595*pow(beta,2) + (-145635 + 459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(6*(76566 - 162165*beta + \
85595*pow(beta,2)) + (-145635 + 459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.88080793e10*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (16*(5139120 + 38640*(-201 \
+ 133*beta)*f1*pow(k,2)*X1q + 6*(-205710 + 570906*beta - \
485415*pow(beta,2) + \
107065*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-205710 + \
570906*beta - 485415*pow(beta,2) + \
107065*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
3864*pow(f1,2)*pow(k,2)*X1q*(-1005 + 665*beta + (1182 - 2010*beta + \
665*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(483*(1182 - 2010*beta + \
665*pow(beta,2)) + (-205710 + 570906*beta - 485415*pow(beta,2) + \
107065*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(161*(1182 - 2010*beta + \
665*pow(beta,2)) + 2*(-205710 + 570906*beta - 485415*pow(beta,2) + \
107065*pow(beta,3))*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(2.7966284775e10*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,9)*pow(q,9))))/8. - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*(-((689520 + 816*(-643 + \
845*beta)*f1*pow(k,2)*X1q + 6*(-8295 + 27795*beta - 32793*pow(beta,2) \
+ 14365*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-8295 + \
27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(-643 + 845*beta + (545 - 1286*beta + \
845*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(545 - 1286*beta + \
845*pow(beta,2)) + (-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(545 - 1286*beta + \
845*pow(beta,2)) + 2*(-8295 + 27795*beta - 32793*pow(beta,2) + \
14365*pow(beta,3))*pow(k,2)*X1q))*(5*k*q*(-21 + \
2*pow(k,2)*pow(q,2))*cos(k*q) + (105 - 45*pow(k,2)*pow(q,2) + \
pow(k,4)*pow(q,4))*sin(k*q)))/(2.9201172e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,5)) - (83*(2713200 + 25840*(-106 \
+ 105*beta)*f1*pow(k,2)*X1q + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + 56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) \
+ (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(9.6026931e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (112*f1*(2 + \
f1)*X1q*(-5320 + f1*(2 + f1)*pow(k,2)*X1q*(2010 - 2660*beta + (-394 + \
5*(201 - 133*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.05308475e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9))))/4. - \
(f1*pow(k,4)*X1q*Y1q*(-((371280 + 4080*(-69 + \
91*beta)*f1*pow(k,2)*X1q + 6*(-4389 + 14841*beta - 17595*pow(beta,2) \
+ 7735*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-4389 + \
14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
408*pow(f1,2)*pow(k,2)*X1q*(5*(-69 + 91*beta) + (291 - 690*beta + \
455*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(51*(291 - 690*beta + \
455*pow(beta,2)) + (-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(17*(291 - 690*beta + \
455*pow(beta,2)) + 2*(-4389 + 14841*beta - 17595*pow(beta,2) + \
7735*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) + \
16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(7.963956e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (7*f1*(2 + f1)*X1q*(-90440 \
+ f1*(2 + f1)*pow(k,2)*X1q*(380*(75 - 119*beta) + (-5033 + 95*(150 - \
119*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.0808226e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*(-((8914800 + \
77520*(-97 + 115*beta)*f1*pow(k,2)*X1q + 6*(-127323 + 416955*beta - \
469965*pow(beta,2) + \
185725*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-127323 + \
416955*beta - 469965*pow(beta,2) + \
185725*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
2280*pow(f1,2)*pow(k,2)*X1q*(17*(-97 + 115*beta) + (1463 - 3298*beta \
+ 1955*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(285*(1463 - 3298*beta + \
1955*pow(beta,2)) + (-127323 + 416955*beta - 469965*pow(beta,2) + \
185725*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1463 - 3298*beta + \
1955*pow(beta,2)) + 2*(-127323 + 416955*beta - 469965*pow(beta,2) + \
185725*pow(beta,3))*pow(k,2)*X1q))*(k*q*(5670 - 735*pow(k,2)*pow(q,2) \
+ 16*pow(k,4)*pow(q,4))*cos(k*q) + (-5670 + 2625*pow(k,2)*pow(q,2) - \
135*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(2.5219194e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) - (7*(542640 + f1*(2 + \
f1)*pow(k,2)*X1q*(2280*(-150 + 119*beta) + f1*(2 + \
f1)*pow(k,2)*X1q*(90594 - 15252*f1*(2 + f1)*pow(k,2)*X1q + \
11305*pow(beta,3)*f1*(2 + f1)*pow(k,2)*X1q - 570*pow(beta,2)*(-119 + \
75*f1*(2 + f1)*pow(k,2)*X1q) + 9*beta*(-19000 + 5033*f1*(2 + \
f1)*pow(k,2)*X1q))))*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(3.2424678e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,9)*pow(q,9))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

			case 4:
				return (-18*Pi*pow(q,2)*((2*pow(f1,2)*pow(2 + f1,2)*pow(X1Infinity,2)*(1 + \
(pow(k,2)*X1Infinity)/2.)*(266 + (-67 + 133*beta)*f1*(2 + \
f1)*pow(k,2)*X1Infinity)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(415701.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,5)*pow(q,9)) + (f1*(2 + \
f1)*X1Infinity*(4*f1*pow(k,2)*X1Infinity + \
2*pow(f1,2)*pow(k,2)*X1Infinity)*(-5320 + f1*(2 + \
f1)*pow(k,2)*X1Infinity*(2010 - 2660*beta + (-394 + 5*(201 - \
133*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1Infinity))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.15701e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1Infinity)/2.)*pow(k,7)*pow(q,9))))/pow(E,(pow(k,2)*\
X1Infinity)/2.) + (18*Pi*pow(q,2)*((f1*(2 + f1)*X1q*(-5320 + f1*(2 + \
f1)*pow(k,2)*X1q*(2010 - 2660*beta + (-394 + 5*(201 - \
133*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*((4*f1*pow(k,2)*X1q + \
2*pow(f1,2)*pow(k,2)*X1q)/2. + (alpha*pow(k,4)*(4*f1*X1q + \
2*pow(f1,2)*X1q)*Y1q)/4.)*(9*k*q*(-225225 + 30030*pow(k,2)*pow(q,2) - \
770*pow(k,4)*pow(q,4) + 4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - \
945945*pow(k,2)*pow(q,2) + 51975*pow(k,4)*pow(q,4) - \
630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(2.078505e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)) + (2*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(266 + (-67 + 133*beta)*f1*(2 + f1)*pow(k,2)*X1q)*(1 \
+ (pow(k,2)*X1q)/2. + ((alpha*pow(k,2))/2. + \
(alpha*pow(k,4)*X1q)/4.)*Y1q)*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(415701.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,5)*pow(q,9)) + (pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,3)*(266 + (-67 + 133*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*Y1q*(-(k*q*(-2006754750 + \
300405105*pow(k,2)*pow(q,2) - 11320155*pow(k,4)*pow(q,4) + \
158679*pow(k,6)*pow(q,6) - 852*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*cos(k*q)) + 3*(-668918250 + \
323107785*pow(k,2)*pow(q,2) - 22286880*pow(k,4)*pow(q,4) + \
501165*pow(k,6)*pow(q,6) - 4418*pow(k,8)*pow(q,8) + \
13*pow(k,10)*pow(q,10))*sin(k*q)))/(831402.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,4)*pow(q,12)) + (f1*(2 + \
f1)*X1q*(4*f1*X1q + 2*pow(f1,2)*X1q)*(-5320 + f1*(2 + \
f1)*pow(k,2)*X1q*(2010 - 2660*beta + (-394 + 5*(201 - \
133*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*Y1q*(-(k*q*(-2006754750 + \
300405105*pow(k,2)*pow(q,2) - 11320155*pow(k,4)*pow(q,4) + \
158679*pow(k,6)*pow(q,6) - 852*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*cos(k*q)) + 3*(-668918250 + \
323107785*pow(k,2)*pow(q,2) - 22286880*pow(k,4)*pow(q,4) + \
501165*pow(k,6)*pow(q,6) - 4418*pow(k,8)*pow(q,8) + \
13*pow(k,10)*pow(q,10))*sin(k*q)))/(8.31402e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,6)*pow(q,12)) - \
(pow(f1,2)*pow(k,4)*X1q*Y1q*((2*(2713200 + 25840*(-106 + \
105*beta)*f1*pow(k,2)*X1q + 6*(-51044 + 162165*beta - \
171190*pow(beta,2) + 56525*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) \
+ (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
760*pow(f1,2)*pow(k,2)*X1q*(-1802 + 1785*beta + (1707 - 3604*beta + \
1785*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + (-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q) + \
2*pow(f1,4)*pow(k,4)*pow(X1q,2)*(95*(1707 - 3604*beta + \
1785*pow(beta,2)) + 6*(-51044 + 162165*beta - 171190*pow(beta,2) + \
56525*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - 60*pow(k,2)*pow(q,2) \
+ pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + 4725*pow(k,2)*pow(q,2) - \
210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(5.30018775e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) + (2*f1*(2 + f1)*X1q*(-5320 \
+ f1*(2 + f1)*pow(k,2)*X1q*(2010 - 2660*beta + (-394 + 5*(201 - \
133*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(4.142475e6*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)) + (600*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(322 + (-108 + 161*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(55*k*q*(11904165 - 1670760*pow(k,2)*pow(q,2) + \
51597*pow(k,4)*pow(q,4) - 468*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) + (-654729075 + \
310134825*pow(k,2)*pow(q,2) - 18918900*pow(k,4)*pow(q,4) + \
315315*pow(k,6)*pow(q,6) - 1485*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(2.401966567e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,11))))/4. - \
(pow(f1,2)*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((2*(8217120 + \
4560*(-1707 + 1802*beta)*f1*pow(k,2)*X1q + 6*(-145635 + 459396*beta - \
486495*pow(beta,2) + \
171190*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-145635 + \
459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
24*pow(f1,2)*pow(k,2)*X1q*(95*(-1707 + 1802*beta) + 2*(76566 - \
162165*beta + 85595*pow(beta,2))*pow(k,2)*X1q) + \
12*pow(f1,4)*pow(k,4)*pow(X1q,2)*(76566 - 162165*beta + \
85595*pow(beta,2) + (-145635 + 459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(6*(76566 - 162165*beta + \
85595*pow(beta,2)) + (-145635 + 459396*beta - 486495*pow(beta,2) + \
171190*pow(beta,3))*pow(k,2)*X1q))*(21*k*q*(495 - \
60*pow(k,2)*pow(q,2) + pow(k,4)*pow(q,4))*cos(k*q) + (-10395 + \
4725*pow(k,2)*pow(q,2) - 210*pow(k,4)*pow(q,4) + \
pow(k,6)*pow(q,6))*sin(k*q)))/(1.590056325e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,7)) + (2*(5139120 + 38640*(-201 \
+ 133*beta)*f1*pow(k,2)*X1q + 6*(-205710 + 570906*beta - \
485415*pow(beta,2) + \
107065*pow(beta,3))*pow(f1,5)*pow(k,6)*pow(X1q,3) + (-205710 + \
570906*beta - 485415*pow(beta,2) + \
107065*pow(beta,3))*pow(f1,6)*pow(k,6)*pow(X1q,3) + \
3864*pow(f1,2)*pow(k,2)*X1q*(-1005 + 665*beta + (1182 - 2010*beta + \
665*pow(beta,2))*pow(k,2)*X1q) + \
8*pow(f1,3)*pow(k,4)*pow(X1q,2)*(483*(1182 - 2010*beta + \
665*pow(beta,2)) + (-205710 + 570906*beta - 485415*pow(beta,2) + \
107065*pow(beta,3))*pow(k,2)*X1q) + \
6*pow(f1,4)*pow(k,4)*pow(X1q,2)*(161*(1182 - 2010*beta + \
665*pow(beta,2)) + 2*(-205710 + 570906*beta - 485415*pow(beta,2) + \
107065*pow(beta,3))*pow(k,2)*X1q))*(9*k*q*(-225225 + \
30030*pow(k,2)*pow(q,2) - 770*pow(k,4)*pow(q,4) + \
4*pow(k,6)*pow(q,6))*cos(k*q) + (2027025 - 945945*pow(k,2)*pow(q,2) + \
51975*pow(k,4)*pow(q,4) - 630*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(2.000815425e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,9)*pow(q,9)) + (24*f1*(2 + \
f1)*X1q*(-32200 + f1*(2 + f1)*pow(k,2)*X1q*(16200 - 16100*beta + \
(-3822 + 25*(324 - 161*beta)*beta)*f1*(2 + \
f1)*pow(k,2)*X1q))*(55*k*q*(11904165 - 1670760*pow(k,2)*pow(q,2) + \
51597*pow(k,4)*pow(q,4) - 468*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*cos(k*q) + (-654729075 + \
310134825*pow(k,2)*pow(q,2) - 18918900*pow(k,4)*pow(q,4) + \
315315*pow(k,6)*pow(q,6) - 1485*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(2.401966567e9*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,9)*pow(q,11))))/8. - \
(f1*pow(k,4)*X1q*Y1q*((2*f1*(2 + f1)*X1q*(-90440 + f1*(2 + \
f1)*pow(k,2)*X1q*(380*(75 - 119*beta) + (-5033 + 95*(150 - \
119*beta)*beta)*f1*(2 + f1)*pow(k,2)*X1q))*(k*q*(-1081080 + \
148995*pow(k,2)*pow(q,2) - 4284*pow(k,4)*pow(q,4) + \
29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - 509355*pow(k,2)*pow(q,2) \
+ 29925*pow(k,4)*pow(q,4) - 434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(3.5334585e7*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,9)) + (2*pow(f1,2)*pow(2 + \
f1,2)*pow(X1q,2)*(126 + (-37 + 63*beta)*f1*(2 + \
f1)*pow(k,2)*X1q)*(k*q*(344594250 - 49324275*pow(k,2)*pow(q,2) + \
1621620*pow(k,4)*pow(q,4) - 16830*pow(k,6)*pow(q,6) + \
46*pow(k,8)*pow(q,8))*cos(k*q) + (-344594250 + \
164189025*pow(k,2)*pow(q,2) - 10405395*pow(k,4)*pow(q,4) + \
190575*pow(k,6)*pow(q,6) - 1080*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(785213.*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,7)*pow(q,11))))/2. - \
(f1*pow(k,4)*(4*f1*X1q + 2*pow(f1,2)*X1q)*Y1q*((2*(542640 + f1*(2 + \
f1)*pow(k,2)*X1q*(2280*(-150 + 119*beta) + f1*(2 + \
f1)*pow(k,2)*X1q*(90594 - 15252*f1*(2 + f1)*pow(k,2)*X1q + \
11305*pow(beta,3)*f1*(2 + f1)*pow(k,2)*X1q - 570*pow(beta,2)*(-119 + \
75*f1*(2 + f1)*pow(k,2)*X1q) + 9*beta*(-19000 + 5033*f1*(2 + \
f1)*pow(k,2)*X1q))))*(k*q*(-1081080 + 148995*pow(k,2)*pow(q,2) - \
4284*pow(k,4)*pow(q,4) + 29*pow(k,6)*pow(q,6))*cos(k*q) + (1081080 - \
509355*pow(k,2)*pow(q,2) + 29925*pow(k,4)*pow(q,4) - \
434*pow(k,6)*pow(q,6) + \
pow(k,8)*pow(q,8))*sin(k*q)))/(1.06003755e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,9)*pow(q,9)) + (6*f1*(2 + \
f1)*X1q*(-135240 + f1*(2 + f1)*pow(k,2)*X1q*(1610*(37 - 42*beta) + \
(-12858 - 805*beta*(-37 + 21*beta))*f1*(2 + \
f1)*pow(k,2)*X1q))*(k*q*(344594250 - 49324275*pow(k,2)*pow(q,2) + \
1621620*pow(k,4)*pow(q,4) - 16830*pow(k,6)*pow(q,6) + \
46*pow(k,8)*pow(q,8))*cos(k*q) + (-344594250 + \
164189025*pow(k,2)*pow(q,2) - 10405395*pow(k,4)*pow(q,4) + \
190575*pow(k,6)*pow(q,6) - 1080*pow(k,8)*pow(q,8) + \
pow(k,10)*pow(q,10))*sin(k*q)))/(6.32096465e8*pow(E,(beta*f1*(2 + \
f1)*pow(k,2)*X1q)/2.)*pow(k,9)*pow(q,11))))/4.))/pow(E,(pow(k,2)*(X1q \
+ alpha*Y1q))/2.) ;
				break ;

		}
		break ;

}


}