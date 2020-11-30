#!/bin/bash
# add gaussian noise to data as dictated by the errors
# now using cpu's random generator as a seed
Eerr=0.003;
Nerr=0.003;
Verr=0.006;


awk -v seed=$RANDOM ' 
function gaussrand(sigma){if(! sigma) sigma=1; rsq=0;
 while ((rsq >= 1)||(rsq == 0))
   { x=2.0*rand()-1.0
     y=2.0*rand()-1.0
     rsq=x*x+y*y
   }
 fac = sqrt(-2.0*log(rsq)/rsq);
 return sigma*x*fac 
}
BEGIN {srand(seed)}           
$1!~"point"{print $0}; $1~"point" {GausE=Ee*gaussrand();GausN=Ne*gaussrand();GausV=Ve*gaussrand(); 
#print $1,$2, $3,$4,$5,$6,$7,$8,$9,Ee,Ne,Ve,$13
print $1,$2, $3,$4,$5,$6,$7+GausE/2,$8+GausN/2,$9+GausV/2,Ee,Ne,Ve,$13
}' Ee=$Eerr Ne=$Nerr Ve=$Verr $1 
