# Steps

1. curve file.txt -> fn,a,b,n,P.x,P.y,Q.x,Q.y
1. generate test points -> python ECCurve.py p curve_file.txt
1. generate q point -> python ECCurve.py q curve_file.txt (have to add manually to curve_file.txt)
1. generate large rho iterations -> python ECCurve.py l curve_file.txt