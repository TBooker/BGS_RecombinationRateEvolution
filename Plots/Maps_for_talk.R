seq(1.75, 3.75, 0.25)
r1 <- c(1.25e-6, 1.5e-6, 2e-6, 2.25e-6, 2.5e-6,2.5e-6, 2.75e-6, 3e-6, 3.5e-6, 3.75e-6)
r1 <- rep(r1, each = 2)


r2 <- rev(r1)

pos_raw <- seq(1,10e6,1e6)
pos <- sort(c(pos_raw, pos_raw-1+1e6))

plot(  pos/1e6, r1*4000 ,
       xlab = "Position (Mbp)",
       ylab = expression(italic(rho)* " = 4N"[e]*"r"),
       type = 'l')

par(mfrow = c(2,1))
plot(  pos/1e6, r1*4000 ,
       xlab = "Position (Mbp)",
       ylab = expression(italic(rho)* " = 4N"[e]*"r"),
       type = 'l')

plot(  pos/1e6, r2*4000 ,
       xlab = "Position (Mbp)",
       ylab = expression(italic(rho)* " = 4N"[e]*"r"),
       type = 'l')



### Hotspot map



r1_raw <- c(27375,37374,74804,84803,747157,757156,763312,773311,838553,848552,1198118,1208117,1225481,1235480,1477996,1487995,1562255,1572254,1780528,1790527,1794761,1804760,1842583,1852582)
r1 <- sort(c(1,r1_raw, r1_raw-1, 2e6))
rates_1 <- c(rep(c(0,0,2.08e-5,2.08e-5), 12),0,0)*4000

r2_raw <- c(163497,            173496,            589640,            599639,            629429,            639428,            660523,            670522,            872705,            882704,            913841,            923840,            1143231,            1148585,            1153230,            1158584,            1265531,            1275530,            1336571,            1346570,            1411392,            1421391)
r2 <- sort(c(1,r2_raw, r2_raw-1, 2e6))
rates_2 <- c(rep(c(0,0,2.08e-5,2.08e-5), 11),0,0)*4000

plot( rates_1 ~ r1, type =  "l",
      ylab = expression(italic(rho)* " = 4N"[e]*"r"),
      xlab = "Position (Mbp)"
)


plot( rates_2 ~ r2, type =  "l",
      ylab = expression(italic(rho)* " = 4N"[e]*"r"),
      xlab = "Position (Mbp)"
)
