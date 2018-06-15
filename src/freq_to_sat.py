#!/usr/bin/python
import sys

from pyeda.inter import *
from pyeda.boolalg.expr import expr2dimacscnf
from pyeda.logic.addition import ripple_carry_add
import math

def parse_F(filename):
    with open(filename) as f:
        f.readline()
        m = int(f.readline().split()[0])
        n = int(f.readline().split()[0])

        F_lb = [ [ 0. for c in range(n)] for p in range(m) ]
        F_ub = [ [ 0. for c in range(n)] for p in range(m) ]

        f.readline()
        for line in f:
            s = line.rstrip("\n").split("\t")
            if len(s) < 8:
                continue

            p = int(s[0])
            c = int(s[4])

            assert(0 <= p < m)
            assert(0 <= c < n)

            F_lb[p][c] = float(s[-2])
            F_ub[p][c] = float(s[-1])

        return F_lb, F_ub, m, n

def generate_G(F_lb, F_ub, m, n):
    G = []
    for c in range(n):
        G.append((n, c))

    for c in range(n):
        for d in range(n):
            if c == d: continue
            ok = True
            for p in range(m):
                if F_ub[p][c] < F_lb[p][d]:
                    ok = False
                    break
            if ok:
                G.append((c, d))

    return G

def to_DOT(G, f):
    f.write("digraph G {\n")
    for (c, d) in G:
        f.write("  %d -> %d\n" % (c,d))
    f.write("}\n")

def discretize(F, m, n, bits):
    nrBins = 1 << bits
    delta = 1. / nrBins
    FF = [ [ int(math.floor(max(0., (F[p][c] - sys.float_info.epsilon) / delta))) for c in range(n)] for p in range(m) ]

    return FF

def to_bit_array(val, bits):
    assert(0 <= val <= (1 << bits))
    B = [ 0 if (1 << (bits - i - 1)) & val == 0 else 1 for i in range(bits) ]

    return B[::-1]

def to_SAT(G, F, m, n, bits):
    # first spanning trees

    sys.stderr.write("Generating variables x[(c,d)] for each edge (c,d)...")
    # X[(c,d)] = 1 iff (c,d) is an edge in T
    X = {}
    for (c, d) in G:
        X[(c, d)] = exprvar("x_%d_%d" % (c, d))
    sys.stderr.write(" Done!\n")

    # H[p, d, c, i] = F[p, d, i] & X[c, d]
    sys.stderr.write("Generating expressions H[p, d, c, i] = F[p, d, i] & X[c, d]...")
    H = [[{} for d in range(n)] for p in range(m)]
    for p in range(m):
        for d in range(n):
            for c in [cc for (cc, dd) in G if dd == d]:
                H[p][d][c] = [And(F[p][d][i], X[(c, d)]) for i in range(bits)]
    sys.stderr.write(" Done!\n")

    # FF[p, c] = sum of frequencies of children of c in sample p
    sys.stderr.write("Generating expressions FF[p, c] = sum of frequencies of children of c in sample p...")
    #FF = exprvars('ff', m, n + 1, bits)
    FF = [ [ [ None for i in range(bits) ] for c in range(n + 1)] for p in range(m) ]
    phi = expr(1)
    for p in range(m):
        for c in range(n + 1):
            S = [ 0 for i in range(bits) ]
            carry = 0
            for d in [ dd for (cc, dd) in G if cc == c ]:
                assert((c,d) in G)
                S, C = ripple_carry_add(H[p][d][c], S, carry)
                carry = C[-1]
                phi = phi & Not(carry, simplify=False)

            for i in range(bits):
                FF[p][c][i] = S[i]
#                 phi = phi & Equal(FF[p][c][i], S[i])

            if c != n:
                # INEQUALITY CONSTRAINTS
                ineq = expr(0)
                for i in range(bits)[::-1]:
                    ineq_i = Not(F[p][c][i]) & FF[p][c][i]
                    for j in range(i + 1, bits):
                        ineq_i = ineq_i & Equal(F[p][c][j], FF[p][c][j], simplify=False)
                    ineq = ineq | ineq_i
                phi = phi & ~ineq
    sys.stderr.write(" Done\n")

    # ADD SPANNING TREE CONSTRAINTS
    for d in range(n):
        parents = [ c for (c, dd) in G if dd == d ]
        phi = phi & (Or(*map(lambda c: X[(c, d)], parents)))
        parent_pairs = [(c1, c2) for (c1, d1) in G for (c2, d2) in G if c1 < c2 and d1 == d2 == d]
        phi = phi & And(*[ Not(And(X[(c1, d)], X[(c2, d)])) for (c1, c2) in parent_pairs])

    return X, phi

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <frequencies> <bits>\n" % sys.argv[0])
        sys.exit(1)

    filename = sys.argv[1]
    bits = int(sys.argv[2])

    F_lb, F_ub, m, n = parse_F(sys.argv[1])
    G = generate_G(F_lb, F_ub, m, n)

    #to_DOT(G, sys.stderr)
    FF = discretize(F_lb, m, n, bits)
    FFF = [ [ to_bit_array(FF[p][c], bits) for c in range(n) ] for p in range(m) ]
    X, phi = to_SAT(G, FFF, m, n, bits)

    sys.stderr.write("Convert to CNF...\n")
    mapping, CNF = expr2dimacscnf(phi.to_cnf())
    #print(phi)
    #sys.stderr.write("Tseitin transform...\n")
    #mapping, CNF = expr2dimacscnf(phi.tseitin())
    for var in mapping:
        print("c", var, mapping[var])
    print(CNF)
