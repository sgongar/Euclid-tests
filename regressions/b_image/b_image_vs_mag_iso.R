fit <- lm(catalogue$B_IMAGE ~ catalogue$MAG_ISO)
print(summary(fit))


# require(splines)
# require(ISLR)
# fit<-lm(B_IMAGE ~ bs(MAG_ISO,knots = c(22,26)),data = catalogue )

# print(summary(fit))

# plot(catalogue$MAG_ISO,catalogue$B_IMAGE,col="grey",xlab="Age",ylab="Wages")
# points(MAG_ISO.grid,predict(fit,newdata = list(MAG_ISO=MAG_ISO.grid)),col="darkgreen",lwd=2,type="l")
# adding cutpoints
# abline(v=c(22, 25),lty=2,col="darkgreen")