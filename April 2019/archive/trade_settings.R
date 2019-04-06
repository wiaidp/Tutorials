# Day one has harpe of 1.83

K<-240
# Spectrum
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
#weight_func[1,]<-1000
periodicity<-10
Lag<-1
L<-10*periodicity
L<-100
L<-min(2*K,L)
plot_T<-F
asset<-"USDEUR"
lag_fx<-1
sign_day_of_week<-rep(1,5)
sign_day_of_week[5]<--1
x<-FX_diff[,asset]
# MSE estimation and trading: sign rule and signal weighting
#--------------------------------------------------------------------
# Same as above but (slightly better sharpe)
Lag<-0
Lag<--1
Lag<--2
Lag<--3
# This is very good for EURUSD with gain-picking
Lag<--6
#------------------------------------------------------------
K<-240
# Spectrum
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
#weight_func[1,]<-1000
periodicity<-20
Lag<-1
L<-10*periodicity
L<-100
L<-min(2*K,L)
plot_T<-F
asset<-"USDEUR"
lag_fx<-1
sign_day_of_week<-rep(1,5)
sign_day_of_week[5]<-sign_day_of_week[3]<--1
#-----------------------------------------------------------
K<-240
# Spectrum
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
#weight_func[1,]<-1000
periodicity<-2.5
Lag<-1
L<-10*periodicity
L<-100
L<-min(2*K,L)
plot_T<-F
asset<-"USDEUR"
lag_fx<-1
sign_day_of_week<-rep(1,5)
sign_day_of_week[5]<--1#sign_day_of_week[3]<--1
#sign_day_of_week[5]<--1
#sign_day_of_week<-c(1,1,-1,1,-1)



#---------------------------------------------------------------

# Settings: equally-weighted filter 
#   weight_func large in frequency zero: non-zero level
#   one-step ahead forecast: periodicity=1 (allpass) and Lag=-1 (one-step ahead)
#   L determines length of filter i.e. coefficients=1/L
K<-240
# Spectrum
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
weight_func[1,]<-1000
periodicity<-1
Lag<--1
L<-10*periodicity
L<-20
L<-min(2*K,L)
plot_T<-T
asset<-"USDEUR"
lag_fx<-1
sign_day_of_week<-rep(1,5)
sign_day_of_week<-c(1,1,-1,1,-1)
