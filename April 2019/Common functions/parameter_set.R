K<-nrow(weight_func)-1
#weight_func[1,]<-1000
periodicity<-10
# Cutoff frequency
cutoff<-pi/periodicity
# Target 
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Real-time estimate
Lag<-0
# Filter length
L<-100
L<-min(2*K,L)
lambda<-eta<-0
lambda_smooth<-lambda_cross<-0
lambda_decay<-c(0,0)





