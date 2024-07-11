
Closest-Pair(P)
construct Px and Py  //O(nlogn)
(p0*, p1*) = Closest-Pair-Rec(Px, Py)

Closest-Pair-Rec(Px, Py)
if |P| <= 3
	find closest pair by measuring all pairwise distances
construct Qx, Qy, Rx, Ry // O(n)
(q0*, q1*) = Closest-Pair-Rec(Qx, Qy)
(r0*, r1*) = Closest-Pair-Rec(Rx, Ry)
delta = min(distance(q0*, q1*), d(r0*, r1*))
x* = maximum x-coordinate of a point in set Q
L = {(x, y): x = x*}
S = points in P within distance delta of L
construct Sy // O(n) time
for each point s in Sy // O(n)
	compute the distance from s to each of the next 15 points in Sy
let s, s' be the pair with the minimum distance
	if distance(s, s') < delta
		return (s, s')
	else if distance(q0*, q1*) < distance(r0*, r1*)
		return (q0*, q1*)
	else
		return (r0*, r1*)