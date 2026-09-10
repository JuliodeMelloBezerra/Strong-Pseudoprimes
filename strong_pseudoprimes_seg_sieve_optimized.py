from sage.all import *
import time

start_time = time.time()

# We implement Arnault's method with h=3, t=64 (i.e. n is prod of 3 primes)
# as described in Prime and Prejudice paper.


def generate_basis(n):
    return list(prime_range(1, n + 1))


def miller_rabin(n, basis=None):
    """
    Miller-Rabin test over the chosen prime basis T.
    """
    if basis is None:
        basis = T
    # if n == 2 or n == 3:
    #     return True
    # if n % 2 == 0:
    #     return False

    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2

    for b in basis:
        if b >= n:
            break
        if gcd(b, n) != 1:
            return False
        x = pow(b, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def random_base_miller_rabin(n, k, seed=None):
    """
    Run k independent Miller-Rabin checks with uniformly random bases.
    """

    if seed is not None:
        set_random_seed(seed)

    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2

    for _ in range(k):
        a = randint(2, n - 2)
        if gcd(a, n) != 1:
            return False
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


t = 1200
T = generate_basis(t)
h = 3


def Sb(T):
    S = {}
    for b in T:
        S[b] = [p for p in range(3, 4 * b, 2) if kronecker_symbol(b, p) == -1]
    return S


def find_congruent(B):
    for r in (1, 3):
        solution = []
        for L in B.values():
            found = next((x for x in L if x % 4 == r), None)
            if found is None:
                break
            solution.append(found)
        if len(solution) == len(B):
            return solution
    return None


def Rt(T, k, S):
    B = {}
    for b in T:
        mod = 4 * b
        residue_sets = []
        for i in range(1, h + 1):
            inv = inverse_mod(k[i], mod)
            residue_sets.append({
                (inv * (s + k[i] - 1)) % mod
                for s in S[b]
            })
        B[b] = sorted(residue_sets[0] & residue_sets[1] & residue_sets[2])
    return B


def supplemental_k(n, basis=None):
    if basis is None:
        basis_set = set(generate_basis(n))
    else:
        basis_set = set(basis)
    return [q for q in prime_range(n + 1, 3 * n + 1) if q not in basis_set]


supplemental = supplemental_k(t, T)
residue_sets = Sb(T)

candidate_sol = None
for i in range(len(supplemental) // 2):
    k = {1: 1, 2: supplemental[i], 3: supplemental[-(i + 1)]}
    B = Rt(T, k, residue_sets)
    candidate_sol = find_congruent(B)
    if candidate_sol is not None:
        break

if candidate_sol is None:
    print('No candidate solution found')
    raise SystemExit

x2 = inverse_mod(-k[3], k[2])
x3 = inverse_mod(-k[2], k[3])

crt1 = candidate_sol + [x2, x3]
crt2 = [i * 4 for i in T] + [k[2], k[3]]
try:
    r = crt(crt1, crt2)
except Exception:
    print('candidate solution is not a solution...')
    raise

mod = 4 * prod(T) * k[2] * k[3]

print('r = %d mod %d , k2=%d, k3=%d' % (r, mod, k[2], k[3]))
current_time = time.time()
total_time = current_time - start_time
print(f'Program took {total_time:.4f} seconds to run.')


# Segmented sieve

def prepare_sieve(M, r, k2, k3, t, B):
    data = []
    for q in prime_range(t+1, B + 1):
        if q == k2 or q == k3:
            continue
        invM = inverse_mod(M, q)
        invk2 = inverse_mod(k2, q)
        invk3 = inverse_mod(k3, q)

        b1 = (-r * invM) % q
        b2 = (b1 + (1 - invk2) * invM) % q
        b3 = (b1 + (1 - invk3) * invM) % q

        data.append((q, tuple({b1, b2, b3})))

    return data


def sieve_block(k0, L, data):
    good = bytearray(b'\x01') * L
    zero = b'\x00' * L

    for q, residues in data:
        for b in residues:
            i = (b - k0) % q
            if i < L:
                n = (L - 1 - i) // q + 1
                good[i::q] = zero[:n]

    return good


rmod8 = candidate_sol[0]
k2mod8 = k[2] % 8
k3mod8 = k[3] % 8

if rmod8 in (1, 5):
    ourcase = [rmod8] * 3
else:
    ourcase = [rmod8, (k2mod8 * (rmod8 - 1) + 1) % 8, (k3mod8 * (rmod8 - 1) + 1) % 8]


def MR_1(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s

    x = pow(2, d, p)
    if x == 1 or x == p - 1:
        return True

    for _ in range(s - 3):
        x = (x * x) % p
        if x == p - 1:
            return True
    return False


def MR_5(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s
    x = pow(2, d, p)
    if (x * x) % p == p - 1:
        return True
    return False


def MR_3(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s
    if pow(2, d, p) == p - 1:
        return True
    return False


def MR_7(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s
    if pow(2, d, p) == 1:
        return True
    return False


def MR(i, p):
    if i == 1:
        return (MR_1(p) and (not p.is_square()))
    if i == 3:
        return MR_3(p)
    if i == 5:
        return MR_5(p)
    if i == 7:
        return MR_7(p)
    raise ValueError(f'unsupported MR base case: {i}')


# Lucas part

SELFRIDGE_D = tuple((-1) ** i * (2 * i + 5) for i in range(20))


def selfridge_D(n):
    n = Integer(n)
    for D in SELFRIDGE_D:
        if kronecker_symbol(D, n) == -1:
            return Integer(D)
    return 45


def lucas_parameters(n):
    D = selfridge_D(n)
    if D == 5:
        P = Integer(5)
        Q = Integer(5)
    else:
        P = Integer(1)
        Q = (1 - D) // 4
    return D, P, Q


def lucas_uv(n, k, P, Q):
    n = Integer(n)
    k = Integer(k)
    P = Integer(P)
    Q = Integer(Q)

    if k == 0:
        return Integer(0), Integer(2 % n), Integer(1 % n)

    D = P * P - 4 * Q
    inv2 = inverse_mod(2, n)

    U = Integer(1)
    V = P % n
    Qk = Q % n

    for bit in bin(k)[3:]:
        U2 = (U * V) % n
        V2 = (V * V - 2 * Qk) % n
        Q2 = (Qk * Qk) % n

        if bit == '0':
            U, V, Qk = U2, V2, Q2
        else:
            U = ((P * U2 + V2) * inv2) % n
            V = ((D * U2 + P * V2) * inv2) % n
            Qk = (Q2 * Q) % n

    return U, V, Qk


def strong_lucas(n):
    D, P, Q = lucas_parameters(n)
    m = n + 1
    s = m.valuation(2)
    d = m >> s

    U, V, Qd = lucas_uv(n, d, P, Q)

    if U == 0:
        return True
    if V == 0:
        return True

    for _ in range(1, s):
        V = (V * V - 2 * Qd) % n
        Qd = (Qd * Qd) % n
        if V == 0:
            return True
    return False


def custom_pseudoprimes(p1, k2, k3, cases):
    if not MR(cases[0], p1):
        return False
    p2 = (p1 - 1) * k2 + 1
    if not MR(cases[1], p2):
        return False
    p3 = (p1 - 1) * k3 + 1
    if not MR(cases[2], p3):
        return False
    # if not random_base_miller_rabin(p1, 2):
    #     return False
    # if not random_base_miller_rabin(p2, 2):
    #     return False
    # if not random_base_miller_rabin(p3, 2):
    #     return False
    return p1, p2, p3


def find_triple(M, r, k2, k3, t, B=10 * (t ** 2), block_size=10 * (t ** 2)):
    data = prepare_sieve(M, r, k2, k3, t, B)
    k0 = 0

    while True:
        good = sieve_block(k0, block_size, data)
        for i, x in enumerate(good):
            if not x:
                continue

            j = k0 + i
            p1 = r + j * M
            candidate = custom_pseudoprimes(p1, k2, k3, ourcase)
            if candidate is False:
                continue
            p1, p2, p3 = candidate

            pseudoprime = p1 * p2 * p3
            if not miller_rabin(pseudoprime):
                continue
            return p1, p2, p3, pseudoprime

        k0 += block_size


p1, p2, p3, pseudoprime = find_triple(mod, r, k[2], k[3], t)

print('p1 is ', p1)
print('p2 is ', p2)
print('p3 is ', p3)
print('pseudoprime: ', pseudoprime)
print('bitlength: ', pseudoprime.nbits(), 'digits ', pseudoprime.ndigits())

current_time = time.time()
total_time = current_time - start_time
print(f'Program took {total_time:.4f} seconds to run.')
