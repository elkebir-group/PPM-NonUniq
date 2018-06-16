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
#    FF = exprvars('ff', m, n + 1, bits)
    SS = [ [ {} for c in range(n+1) ] for p in range(m) ]
    CC = [ [ {} for c in range(n+1) ] for p in range(m) ]

    #FF = [ [ [ None for i in range(bits) ] for c in range(n + 1)] for p in range(m) ]
    phi = expr(1)
    for p in range(m):
        for c in range(n + 1):
            sys.stderr.write("Generating expressions FF[%d, %d] = sum of frequencies of children of %d in sample %d..." % (p, c, c, p))
            sys.stderr.flush()

            S = [ 0 for i in range(bits) ]
            carry = 0
            children = [ dd for (cc, dd) in G if cc == c ]
            if len(children) == 1:
                d = children[0]
                SS[p][c][d] = exprvars('ss_%d_%d_%d' % (p,c,d), bits)
                for i in range(bits):
                    if F[p][d][i]:
                        phi = phi & SS[p][c][d][i]
                    else:
                        phi = phi & Not(SS[p][c][d][i])
            else:
                for (idx, d) in enumerate(children):
                    assert((c,d) in G)

                    if idx == 0: continue
                    if idx == 1:
                        SS[p][c][d] = exprvars('ss_%d_%d_%d' % (p,c,d), bits)
                        CC[p][c][d] = exprvars('cc_%d_%d_%d' % (p,c,d), bits)
                        S, C = ripple_carry_add(H[p][d][c], H[p][children[0]][c], 0)
                        phi = phi & And(*[Equal(SS[p][c][d][i], S[i]) for i in range(bits)])
                        phi = phi & And(*[Equal(CC[p][c][d][i], C[i]) for i in range(bits)])
                        phi = phi & Not(CC[p][c][d][-1])
                    else:
                        SS[p][c][d] = exprvars('ss_%d_%d_%d' % (p,c,d), bits)
                        CC[p][c][d] = exprvars('cc_%d_%d_%d' % (p,c,d), bits)
                        S, C = ripple_carry_add(H[p][d][c], SS[p][c][children[idx-1]], CC[p][c][children[idx-1]][-1])
                        phi = phi & And(*[Equal(SS[p][c][d][i], S[i]) for i in range(bits)])
                        phi = phi & And(*[Equal(CC[p][c][d][i], C[i]) for i in range(bits)])
                        phi = phi & Not(CC[p][c][d][-1])

            sys.stderr.write(" Done\n")

            if c != n and len(children) > 0:
                # INEQUALITY CONSTRAINTS
                ineq = expr(0)
                for i in range(bits)[::-1]:
                    ineq_i = Not(F[p][c][i]) & SS[p][c][children[-1]][i]
                    for j in range(i + 1, bits):
                        ineq_i = ineq_i & Equal(F[p][c][j], SS[p][c][children[-1]][j], simplify=True)
                    ineq = ineq | ineq_i
                phi = phi & ~ineq

    # ADD SPANNING TREE CONSTRAINTS
    for d in range(n):
        parents = [ c for (c, dd) in G if dd == d ]
        phi = phi & (Or(*map(lambda c: X[(c, d)], parents)))
        parent_pairs = [(c1, c2) for (c1, d1) in G for (c2, d2) in G if c1 < c2 and d1 == d2 == d]
        phi = phi & And(*[ Not(And(X[(c1, d)], X[(c2, d)])) for (c1, c2) in parent_pairs])

    # ADD MONOCLONALITY CONSTRAINTS
    parent_pairs = [(d1, d2) for (c1, d1) in G for (c2, d2) in G if d1 < d2 and c1 == c2 == n]
    phi = phi & And(*[ Not(And(X[(c, d1)], X[(c, d2)])) for (d1, d2) in parent_pairs])

    return X, SS, phi

def solve(G, F, X, SS, phi, m, n, bits):
    l = list(phi.satisfy_all())
    for idx,sol in enumerate(l):
        with open("%d.dot" % idx, "w") as f:
            f.write("digraph G {\n")
            for c in range(n+1):
                children = [ dd for (cc, dd) in G if cc == c ]
                if len(children) > 0:
                    f.write('  %d [label="%d\\n' % (c, c))
                    for p in range(m):
                        if c == n:
                            f.write('1' * bits + ' - ')
                        else:
                            f.write("".join(map(str, F[p][c])) + ' - ')
                        for i in range(bits):
                            f.write(str(sol[SS[p][c][children[-1]][i]]))
                        f.write("\\n")
                    f.write('"]\n')
                else:
                    f.write('  %d [label="%d\\n' % (c, c))
                    for p in range(m):
                        f.write("".join(map(str, F[p][c])) + ' - ')
                        f.write("\\n")
                    f.write('"]\n')


            for (c, d) in G:
                f.write("  %d -> %d%s\n" % (c, d, " [color=red]" if sol[X[(c, d)]] == 1 else ""))
            f.write("}\n")

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
    X, SS, phi = to_SAT(G, FFF, m, n, bits)

    sys.stderr.write("Convert to CNF...\n")
    phi = phi.to_cnf()
    mapping, CNF = expr2dimacscnf(phi)
    for (idx, (c,d)) in enumerate(G):
        if idx % 10 == 0:
            if idx > 0:
                sys.stdout.write(" 0\n")
            sys.stdout.write("c ind")
        sys.stdout.write(" " + str(mapping[X[(c,d)]]))
    sys.stdout.write(" 0\n")

    for var in mapping:
        if str(var).isdigit():
            print("c", var, mapping[var])
    #print(CNF)
    #solve(G, FFF, X, SS, phi, m, n, bits)
