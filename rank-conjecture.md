# Conjecture 1: rank(M_LWE-CC_n) = Θ(q^n)

## Evidence
| (n,q) | q^n | rank (computed) |
|--------|-----|-----------------|
| (2,5)  | 25  | 31              |
| (2,7)  | 49  | 53              |
| (3,5)  | 125 | 126             |

## Proof Strategy (incomplete)
Lower bound: show the q^n rows indexed by distinct secrets are linearly
independent over Q using the DFT structure of (s -> k*s_j mod q) functions.
The current submatrix approach gives rank n(q-1) < q^n -- needs full Bob input space.

Upper bound: follows from output space size q^m.
