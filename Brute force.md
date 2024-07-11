BruteForceClosest(P)
n = P.length
dmin = infinity
for i = 1 to n-1
	for j = i+1 to n
		d = ((xi - xj)^2 + (yi-yj)^2)
		if d < dmin
		dmin = d;
		index1 = i;
		index2 = j;
return index1, index2;