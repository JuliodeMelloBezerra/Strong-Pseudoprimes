# Constructing Strong Pseudoprimes

Our goal is to construct strong pseudoprimes with respect to the base consisting of all primes smaller than a given $t>0$. In the next chapter we will briefly explain what is a strong pseudoprime to a given base.

## The Miller-Rabin primality test

Any prime number $p$ satisfies the following conditions:

- (Fermat's little theorem) for any base $a \in \mathbb{Z}$ not divisible by $p$, we have that $a^{p-1} \equiv 1 \mathrm{mod} p$.
- The only roots of unity over $\mathbb{F}_p$ are $\pm 1$. That is, $x^2 \equiv 1 \mathrm{mod} p \Leftrightarrow x \equiv \pm 1$.

These two facts are always true for prime numbers, but often not true if we replace the prime $p$ with a composite $n$.

The idea behind the Miller-Rabin probabilistic primality test is to check, for a given odd $n \in \mathbb{Z}$ and coprime base $a \in \mathbb{Z}$, whether it satisfies the two properties above. If not, then $n$ is definitely a composite number. However, if $n$ satisfies the two conditions for a base $a$, we cannot conclude whether it is a prime or not. In this case $n$ is either a prime or a **strong pseudoprime** to the base $a$.

The test goes as follows: since $n$ is odd, we can write $n-1 = 2^sd$ for $d$ odd. Then, if $n$ is prime, we would have by Fermat's little theorem $(a^d)^{2^s} \equiv 1 \mathrm{mod} p $, and by taking $s$ successive square roots, all of the values $(a^d)^{2^r}$ for $0 \leq r \lt s $ would have to be congruent to $\pm 1$, by the second condition. The test then consists of simply checking whether this is true for all such $r$. It is clear that, starting from $r=0$, if any of the $(a^d)^{2^r}$ is congruent to $-1$, then all the subsequent values will be congruent to $1$, so it is enough to check if $a^d \equiv 1 \mathrm{mod} p$ or $(a^d)^{2^r} \equiv -1 \mathrm{mod} p$ for all $1 \leq r$ only up to $r=s-1$. 

We now remark that no composite number is a strong pseudoprime for all bases. In fact, at most $1/4$ of the bases between $1\lt a\lt n-1$ can make $n$ a strong pseudoprime. This means that, if $n$ is composite and $a$ a base chosen at random, there is less than $1/4$ chance that $n$ will pass the test. By choosing $k$ different bases at random, the chance is $(1/4)^k$, which can be made arbitrarily small. Since the test is also computationally cheap (more on this later), it is very good for determining whether a number is prime or not with very high accuracy.

However, if the bases are chosen from a fixed set, say the set $T$ of all primes smaller than a given $t>0$, then it is possible to create a strong pseudoprime to all bases in this set. A method for doing this was described by Arnault [^arnault], and we dedicate the rest of this article to implementing it.



## Arnault's method

We will not explain the method in detail for now, nor provide proofs, as that is done very well in Arnault's short original paper[^arnault], as well as in the Appendix of the Prime and Prejudice paper[^prime] . What we can say is that the method consists of finding a composite number which is the product of an odd number of primes (or even $\geq 4$) which will be a strong pseudoprime with respect to the given basis. The simplest case is when our composite $n = p_1p_2p_3$ is the product of three distinct primes, which is the only case we cover here. The first prime $p_1$ is constructed via an appropriate residue condition, and then the two other primes $p_2,p_3$ are constructed from $p_1$.

We begin by generating the list $T$ of primes smaller than a given $t>0$ by using the sage function```prime_range```, though a simple sieve of Eratosthenes is sufficient since $t$ is relatively small (see the section on asymptotic heuristics). 

<details>
<summary>Click to show code</summary>

```python
def generate_basis(n):
    return list(prime_range(1, n + 1))
```
</details> 

We then generate the set $S_b$ of odd non-quadratic residues for each prime $b \in T$ in the basis, by using the following dictionary

<details>
<summary>Click to show code</summary>

```python
def Sb(T):
    S = {}
    for b in T:
        S[b] = [p for p in range(3, 4 * b, 2) if kronecker_symbol(b, p) == -1]
    return S
```
</details> 

Now, as Prime and Prejudice[^prime] suggest, we should find two small primes $k_2,k_3$ not in $T$ such that the following intersection is non-empty in $\mathbb{Z}/4b\mathbb{Z}$ for each $b \in T$.

$$
k_2^{-1}(S_b + k_2 - 1) \cap k_3^{-1}(S_b + k_3 - 1) \neq \emptyset \hspace{10pt} \forall b \in T.   
$$

Finding these two values is not difficult, but no method other than trial-and-error is given in any of the papers. Based on the examples given, we simply generate a list of primes $t\lt q\lt 3t$ and try pairs by combining the first and last primes in the list, then second and penultimate, and so on. The code to generate this list and the new dictionary from the non-empty intersections is given below.

<details>
<summary>Click to show code</summary>

```python
def supplemental_k(n, basis=None):
    if basis is None:
        basis_set = set(generate_basis(n))
    else:
        basis_set = set(basis)
    return [q for q in prime_range(n + 1, 3 * n + 1) if q not in basis_set]
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
```
</details> 

Before we loop over the different choices of $k_2,k_3$, there is one last remark: we must now choose one member of each $S_b$ for every $b \in T$, such that they form a residue modulo $\mathrm{lcm}(4,b_1,\dots,b_m)$ for all $b_i \in T$. That is, they must be compatible in the sense that we can successfully apply the Chinese Remainder Theorem (CRT) on them to get a common solution.

Both the Prime and Prejudice article[^prime], as well as the other implementations known to the author, proceed to perform multiple CRT's and "back-tracking" to find such a common solution, but a much easier approach is possible. Indeed, simply note that all of the moduli are coprime, except for the number $4$ which is multiplying all of them. That means that the residue modulo $4$ is the only obstruction to finding a common solution, and it suffices to select elements of each $S_b$ such that they are all congruent to the same value modulo $4$. Since they are all odd, that means simply either all congruent to $1$ or $3$ modulo $4$, for all primes $b \in T$. No CRT's are needed in this process. If no such compatible list exists, we simply try different $k_i$'s until we have it.

<details>
<summary>Click to show code</summary>

```python
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

supplemental = supplemental_k(t, T)
residue_sets = Sb(T)

candidate_sol = None
for i in range(len(supplemental) // 2):
    k = {1: 1, 2: supplemental[i], 3: supplemental[-(i + 1)]}
    B = Rt(T, k, residue_sets)
    candidate_sol = find_congruent(B)
    if candidate_sol is not None:
        break
```
</details> 

We are now almost done. We only need to include two more simple relations that our residue $r$ must satisfy, namely $r \equiv k_3^{-1} \mathrm{mod} k_2$ and $r \equiv k_2^{-1} \mathrm{mod} k_3$. These are well defined by the definition of the $k_i$'s, and so is the solution to the final CRT when we include these to our previous compatible solution. 

<details>
<summary>Click to show code</summary>

```python
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
```
</details> 

We now have the congruence that our first prime $p_1$ must satisfy, namely $p_1 \equiv r \mathrm{mod} M$, for $M := \mathrm{lcm}(4,b_1,\dots,b_m,k_2,k_3)$ where again the $b_i$'s are all elements of $T$.

The next step is simply finding a prime $p_1$ satisfying this residue condition, and such that we can build two more primes from it, as explained in the next section. This next step is also the most expensive one.

## Finding the three primes
This is by far the most expensive part of our algorithm, both in terms of time and memory.

We want to build a prime $p_1 \equiv r \mathrm{mod} M$, and such that $p_2 = k_2(p_1-1) + 1$  and $p_3 = k_3(p_1-1) + 1$ are also primes. Naively, we could try $r + j \cdot M$ for each $j\in \mathbb{N}$ and check whether this defines our three primes, but we can do better.
It's clear that our three candidate numbers $p_1,p_2,p_3$ will never be divisible by any prime in $T$. However, we can also use knowledge from the structure of the three candidates to skip over values of $j \in \mathbb{N}$ for which one of the three $p$'s will be divisible by a small prime. This is best accomplished by constructing an appropriate sieve, which is done in the next section. 

Afterwards, once we have candidates that passed the sieve, we also need to check if they really are prime numbers. We don't need absolute certainty that they are prime, because if they are not, they will likely fail the final Miller-Rabin test with respect to $T$, and then we can try again. This means that we need to check if they are highly-likely prime numbers, so in other words, we use a "pseudoprimality" test. Once again, we may use the structure of these candidates to improve on known primality tests, constructing our own customized primality test. This will be done in the second section of this chapter. 
### The sieve

Let $S$ be a set of "small" primes such that $S \cap (T \cup \{k_2,k_3 \}) = \emptyset$. 

Let $A_i := k_i(r-1) + 1$, for $i=2,3$, so that $p_i = A_i + k_i  M j$. Then for any $q \in S$, we have that

$$
\begin{align}
q|p_1 \Leftrightarrow  j  \equiv -rM^{-1} \mathrm{mod} q \\

q|p_2 \Leftrightarrow   j  \equiv -A_2(k_2M)^{-1} \mathrm{mod} q \\

q|p_3 \Leftrightarrow  j  \equiv -A_3(k_3M)^{-1} \mathrm{mod} q
\end{align}
$$

Note that these inverses exist by the definition of $S$. Our sieve should therefore skip over the $j$'s satisfying any of the congruences above, for any prime $q \in S$, in which case one of the $p$'s won't be prime.

We can simplify the notation and save a couple of redundant modular reductions by setting 

$$
b := -rM^{-1} \text{ and }d_i := (1-k_i^{-1})M^{-1} \in \mathbb{Z}/q\mathbb{Z}
$$

 so that for each interval $[aq, (a+1)q)$ with $a \in \mathbb{N}$, we eliminate $b, b+d_1, b+d_2$ in that interval. This means that for each such interval, we eliminate up to 3 values (almost always exactly 3), and we do this for each $q \in S$. 

The basic math for the sieve is done, what is left is the implementation. We make the sieve **segmented**, which does not improve the time-cost, but quadratically improves its memory-cost. Later, in the third section, we will use some asymptotic heuristics to better choose our set $S$ used for the sieve.

<details>
<summary>Click to show code</summary>

```python
def prepare_sieve(M, r, k2, k3, t, B):
    data = []

    for q in prime_range(t + 1, B + 1):
        if q == k2 or q == k3:
            continue

        invM = inverse_mod(M, q)
        invk2 = inverse_mod(k2, q)
        invk3 = inverse_mod(k3, q)

        b1 = (-r * invM) % q
        b2 = (b1 + ( 1 - invk2) * invM) % q
        b3 = (b1 + ( 1 - invk3) * invM) % q

        data.append((q, tuple(set((b1, b2, b3)))))

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
```
</details> 

### Primality testing

Once the sieve is done, we test each number that passed the sieve for primality, at least in a probabilistic way. In an older version of our program, we used the Baillie-PSW primality testing.

We will not explain the history and math behind the development of Baillie-PSW, as that can be found in the original and subsequent papers. For now we will just remark that the idea is to use two different probabilistic primality tests, which are in some sense "orthogonal" to each other, in order to lower the odds of a composite number passing the test. Indeed, no known examples of composites which pass the test are known, though infinitely many are expected according to some heuristics (they must however be bigger than $2^{64}$).

The ``` is_pseudoprime() ``` function in Sage uses the similar named function from PARI/GP, which consists of three steps:
- In the first step, when applied to an integer $n$, it tests whether $n$ is divisible by small primes first, where small should be less than 101. 
- The second step is a single Miller-Rabin test with base 2.
- The third and final step, which is also the most expensive, consists of a certain almost-extra-strong-Lucas primality test. We will skip discussion of this third step for now, as it is the most complicated one. 

It is clear that for our candidate numbers $p_1,p_2,p_3$, the first step is not necessary. It turns out we can also make the second step faster, first by better implementing Miller-Rabin to the specific base 2, and then also by using knowledge about the structure of our candidates.

First, note that for an odd prime $p = 2^sd + 1$, quadratic reciprocity with respect to the prime $2$ implies the following congruences:
- $p \equiv 1 \mathrm{mod} 8 \; $   $ \; \Rightarrow \; $  $ ( \; 2^d \equiv 1 \mathrm{mod} p$ or $2^{2^rd} \equiv -1 \mathrm{mod} p$ for some $0 \leq r \leq s-2 \; )$
- $p \equiv 3 \mathrm{mod} 8 \;$   $\; \Rightarrow \;$  $ \;\;\, 2^d \equiv -1 \mathrm{mod} p$
- $p \equiv 5 \mathrm{mod} 8 \;$   $\; \Rightarrow \;$  $ \;\;\, 2^{2d} \equiv -1 \mathrm{mod} p$
- $p \equiv 7 \mathrm{mod} 8 \;$   $\; \Rightarrow \;$  $ \;\;\, 2^d \equiv 1 \mathrm{mod} p$

Therefore the "full" Miller-Rabin test with base 2 is only necessary in the first case, when $p \equiv 1 \mathrm{mod} 8$, and it's enough to test only until $r \leq s-2$, not $s-1$. In the other cases, only one congruence needs to be tested. 

The above remark also works well in our case, since we already know the value of $p_1 \mathrm{mod} 8$ from the construction, which is $p_1 \equiv r \mathrm{mod} 8$. In fact we don't even need to compute this value, since it was already computed in the construction of $r$, and it's none other than the first value in the list of candidate solutions ``` candidate_sol[0].```It is also easy to see from the definitions that $p_1 \equiv 1 \mathrm{mod} 8 \Leftrightarrow p_2,p_3 \equiv 1 \mathrm{mod} 8$ and hence $p_1 \equiv 5 \mathrm{mod} 8 \Leftrightarrow p_2,p_3 \equiv 5 \mathrm{mod} 8$. The cases where $p_1 \equiv 3 \mathrm{mod} 4$ are not so direct, and we need to compute $k_2$ and $k_3 \mathrm{mod} 8$ in order to get the residues of $p_2,p_3 \mathrm{mod} 8$. 

All of this can be done before the sieve loop, so we can already determine in which of the above cases we are, and which Miller-Rabin test with base 2 we need to perform in the loop. In fact, even the value of $s$ in $p-1 = 2^sd$ can be computed only once, for $p_1$, and it will be the same for $p_2$ and $p_3$, as per Arnault's lemma[^arnault]. This $s$ will necessarily be equal to $ 1$ if $p_1 \equiv 3$ or $7 \mathrm{mod} 8$, equal to $2$ if $p_1 \equiv 5 \mathrm{mod} 8$ and necessarily strictly greater than $2$ otherwise.

While the above discussion reveals certain small optimizations, it should be remarked that they are tiny compared to more computation heavy operations. Still, we want to optimize as much as possible everything that happens inside the big loop. 

Experimentation has shown that our optimized primality testing indeed yields a faster program than Sage's ``` is_pseudoprime() ```. In fact, we get a faster program by ignoring the Lucas' test entirely, since it is more expensive than a regular Miller-Rabin test (though asymptotically equivalent), and high confidence of primality for our three candidates is not so important, as we test their pseudoprimality in the end anyhow. As a last remark, we also experimented with random-base Miller-Rabin tests, and the single base 2 test still remained faster or equivalent.

 Here is the implementation:

<details>
<summary>Click to show code</summary>

```python
rmod8 = candidate_sol[0] 

#p1 = rmod8 mod 8, and =1,5 iff p2,p3 =1,5. For 3 and 7 depends on k2,k3 mod 8.

k2mod8 = k[2] % 8
k3mod8 = k[3] % 8

if rmod8 == 1 or rmod8 == 5:
    ourcase = [rmod8] * 3
else:
    ourcase = [rmod8,(k2mod8*(rmod8 -1) +1)%8, (k3mod8*(rmod8 -1) +1)%8  ]


#with this knowledge, Miller-Rabin for base 2 can be replaced with the following, using quadratic reciprocity and basic arithmetic:



def MR_1(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s

    x = pow(2, d, p)

    if x == 1 or x == p - 1:
        return True

    for _ in range(s - 3): #enough to check until r = (s-2). We already checked r=0 in previous step.
        x = (x * x) % p
        if x == p - 1:
            return True
    return False

def MR_5(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s
    x = pow(2,d,p)
    if (x * x)%p == p-1:
        return True
    return False

def MR_3(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s
    if pow(2,d,p) == p-1:
        return True
    return False

def MR_7(p):
    d = p - 1
    s = d.valuation(2)
    d >>= s
    if pow(2,d,p) == 1:
        return True
    return False

def MR(i,p):
    if i == 1:
        return (MR_1(p) and (not p.is_square())) #already eliminate squares, they can't be prime and necessary for lucas test. Only possible if =1 mod 8.
    if i == 3:
        return MR_3(p)
    if i == 5:
        return MR_5(p)
    if i == 7:
        return MR_7(p)

```
</details> 
 
The code for the Lucas test, which isn't used in the end, but can be added for experimentation (not the same Lucas variant as other Baillie-PSW versions!).

<details>
<summary>Click to show code</summary>

```python
SELFRIDGE_D = tuple((-1)**i * (2*i + 5) for i in range(20))

def selfridge_D(n):
    """
    Find the first D in

        5, -7, 9, -11, 13, -15, ...

    for which Jacobi(D/n) = -1.

    Returns D.
    """
    n = Integer(n)

    i = 0
    while i < 20 :
        D = SELFRIDGE_D[i]

        j = kronecker_symbol(D, n)

        if j == -1:
            return Integer(D)

        i += 1
    return 45


def lucas_parameters(n):
    """
    Return (D, P, Q) according to the requested method:

        D = first element of
            5, -7, 9, -11, 13, -15, ...

            with Jacobi(D/n) = -1.

        Normally:
            P = 1
            Q = (1-D)/4

        But for D = 5:
            P = Q = 5
    """
    D = selfridge_D(n)

    if D == 5:
        P = Integer(5)
        Q = Integer(5)
    else:
        P = Integer(1)
        Q = (1 - D) // 4

    return D, P, Q


def lucas_uv(n, k, P, Q):
    """
    Compute (U_k, V_k, Q^k) modulo n using binary Lucas arithmetic.

    Lucas sequences are defined by

        U_0 = 0
        U_1 = 1
        U_{j+2} = P U_{j+1} - Q U_j

        V_0 = 2
        V_1 = P
        V_{j+2} = P V_{j+1} - Q V_j

    Returns:
        (U_k mod n, V_k mod n, Q^k mod n)
    """
    n = Integer(n)
    k = Integer(k)
    P = Integer(P)
    Q = Integer(Q)

    if k == 0:
        return Integer(0), Integer(2 % n), Integer(1 % n)

    # Discriminant:
    D = P * P - 4 * Q

    # n is odd in the Lucas test, so 2 has an inverse mod n.
    inv2 = inverse_mod(2, n)

    # State for index 1:
    #
    # U_1 = 1
    # V_1 = P
    # Q^1 = Q
    U = Integer(1)
    V = P % n
    Qk = Q % n

    # Process binary representation of k from the second bit onward.
    for bit in bin(k)[3:]:
        # Doubling:
        #
        # U_{2j} = U_j V_j
        # V_{2j} = V_j^2 - 2 Q^j
        # Q^{2j} = (Q^j)^2

        U2 = (U * V) % n
        V2 = (V * V - 2 * Qk) % n
        Q2 = (Qk * Qk) % n

        if bit == '0':
            U, V, Qk = U2, V2, Q2

        else:
            # Addition by one:
            #
            # U_{2j+1}
            #     = (P U_{2j} + V_{2j}) / 2
            #
            # V_{2j+1}
            #     = (D U_{2j} + P V_{2j}) / 2
            #
            # Q^{2j+1} = Q^{2j} Q

            U = ((P * U2 + V2) * inv2) % n
            V = ((D * U2 + P * V2) * inv2) % n
            Qk = (Q2 * Q) % n

    return U, V, Qk


def strong_lucas(n):
    """
    Strong Lucas probable-prime test with Selfridge parameters.

    Parameter selection:

        D = first value in
            5, -7, 9, -11, 13, -15, ...

        for which Jacobi(D/n) = -1.

        If D == 5:
            P = Q = 5

        otherwise:
            P = 1
            Q = (1-D)/4

    Strong Lucas test:

        n + 1 = 2^s * d, d odd.

    Compute U_d, V_d.

    Accept iff:

        U_d == 0 mod n

    or

        V_{2^r d} == 0 mod n

    for some 0 <= r < s.
    """

    D, P, Q = lucas_parameters(n)

    # n + 1 = 2^s * d, with d odd.
    m = n + 1
    s = m.valuation(2)
    d = m >> s

    U, V, Qd = lucas_uv(n, d, P, Q)

    # First strong-Lucas condition.
    if U == 0:
        return True

    # V_d = 0 is the r=0 case.
    if V == 0:
        return True

    # Repeated doubling:
    #
    # V_{2j} = V_j^2 - 2 Q^j
    # Q^{2j} = (Q^j)^2
    #
    for _ in range(1, s):
        V = (V * V - 2 * Qd) % n
        Qd = (Qd * Qd) % n

        if V == 0:
            return True

    return False

```
</details> 
The random-base Miller-Rabin test, also for experimentation:

<details>
<summary>Click to show code</summary>

```python
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

```
</details> 



Our customized primality test then becomes:
<details>
<summary>Click to show code</summary>

```python
def custom_pseudoprimes(p1):
    #if anything fails at any point, return false. Otherwise, return True

    if not MR(ourcase[0],p1):
        return False
    p2 = (p1 - 1)*k[2] + 1
    if not MR(ourcase[1],p2):
        return False
    p3 = (p1 - 1)*k[3] + 1
    if not MR(ourcase[2],p3):
        return False
    # if not strong_lucas(p1):
    #     return False
    # if not strong_lucas(p2):
    #     return False
    # if not strong_lucas(p3):
    #     return False
    # if not random_base_miller_rabin(p1, 2):
    #     return False
    # if not random_base_miller_rabin(p2, 2):
    #     return False
    # if not random_base_miller_rabin(p3, 2):
    #     return False
    return p1, p2, p3

```
</details> 


Before we implement our code for finding the three primes, we must find the optimal parameters to be used on the sieve.
### Asymptotic heuristics and parameter optimization

First let's try to estimate the size of $M$ and of our final pseudoprime, with respect to the original input $t$.

Remember that by definition 

$$
M = 4 \cdot k_2 \cdot k_3 \cdot \prod \limits_{p \lt t} p
$$

where the product is over all primes smaller than $t$. This means that $\mathrm{ln}(M) = \vartheta(t) + \mathrm{ln}(4) + \mathrm{ln}(k_2) + \mathrm{ln}(k_3)$, where $\vartheta$ is Chebyshev's first function. Ignoring the extra terms and using the known asymptotic behavior of this function $\vartheta(t) \sim t$, we asymptotically approximate $M \sim e^{(1+o(1))t}$, which means the number of digits of $M$ should be roughly $\mathrm{log}(e) \cdot t$, that is, proportional to $t$. The same is true for our final pseudoprime, though with a different constant of proportionality, since it should be "slightly" bigger than $M^3$.

Now, a single Miller-Rabin test on a number with $t$ digits should have time complexity $\tilde O(t^2)$ using a FFT based multiplication. Assuming a similar time complexity for the Lucas test, this should be the total cost of our primality test. 

Let's now have a look at the cost of sieving. Suppose our sieve set $S$ is defined as all primes smaller than $B$, but greater than $t$. The two primes $k_2,k_3$ don't matter for these rough asymptotic estimates. Suppose further that for every prime $q \in S$, we eliminate exactly $3$ numbers from the interval of length $q$, as remarked in the section on Sieves. That means that the surviving numbers after sieving should be given by

$$
s(t,B):=\prod \limits_{t\lt q\leq B} \left(1-\frac{3}{q}\right).
$$

Now let's look at this function more closely. First note that for each $q$, we have 

$$
0 \lt  \left( 1- \frac{3}{q}\right) \lt  \left( 1- \frac{1}{q}\right)^3
$$

since the difference is exactly $3/q^2 - 1/q^3$ which is positive. Now, taking limits and products, we get

$$
0 \leq \lim_{n\rightarrow \infty} \mathrm{ln}(n)^3 \prod_{q\leq n}\left( 1- \frac{3}{q}\right) \leq \lim_{n\rightarrow \infty} \mathrm{ln}(n)^3 \prod_{q\leq n} \left( 1- \frac{1}{q}\right)^3 = e^{-3\gamma}
$$

where the last equality is Merten's third theorem. In particular, this means that asymptotically $\prod_{q\leq n}\left( 1- \frac{3}{q}\right) \sim \mathrm{ln}(n)^{-3}$. So we may conclude

$$
s(t,B) = \frac{\prod_{q\leq B}\left( 1- \frac{3}{q}\right)}{\prod_{q\leq t}\left( 1- \frac{3}{q}\right)} \sim \left(\frac{\mathrm{ln}(t)}{\mathrm{ln}(B)}\right)^3.
$$

That means that if we take $B = t^a$ for some real number $a$, the number of survivors after the sieve will be roughly $s(t,t^a) \sim 1/a^3$. Based on the Sieve of Eratosthenes, the time complexity for sieving should be $\tilde O(B) = \tilde O(t^a)$ using our choice for $B$. So as $a$ increases, so does the time for sieving, but the number of survivors on which we need to perform primality tests decreases. This means that, in order to find $a$ such that the total time is minimal, we need to solve the (approximate) equation 

$$
t^a \approx \frac{t^2}{a^3}
$$

where the left hand side corresponds to the cost of sieving and the right hand side the cost of performing primality tests on the survivors. Since $a$ is much smaller than $t$, this means that the optimal choice should be $a=2$, so we should sieve over all primes between $t$ and  $t^2$. Experimentation shows that this is indeed a good choice. We also remark that with a segmented sieve the space complexity of our sieve algorithm should be roughly $\tilde O(t)$ for $B=t^2$.

Now we are ready to implement the final part of our code:

<details>
<summary>Click to show code</summary>

```python
def find_triple(M, r, k2, k3, t,
                B=10*(t**2),
                block_size=10*(t**2)):

    data = prepare_sieve(M, r, k2, k3, t, B)

    k0 = 0

    while True:
        good = sieve_block(k0, block_size, data)

        for i, x in enumerate(good):
            if not x:
                continue

            j = k0 + i

            p1 = r + j*M
            if not custom_pseudoprimes(p1):
                continue
            p1,p2,p3 = custom_pseudoprimes(p1)
            #Commented-out version using Sage's built-in function for experimentation
            
            # if not is_pseudoprime(p1):
            #     continue

            # p2 = k[2] * p1 - (k[2]-1)
            # if not is_pseudoprime(p2): 
            #     continue

            # p3 = k[3] * p1 - (k[3]-1)
            # if not is_pseudoprime(p3): 
                # continue            




            pseudoprime = p1*p2*p3
            if not miller_rabin(pseudoprime):
                continue
            return p1,p2,p3,pseudoprime

        k0 += block_size

p1,p2,p3,pseudoprime = find_triple(mod, r, k[2], k[3], t)



print("p1 is ", p1)
print("p2 is ", p2)
print("p3 is ", p3)
print("pseudoprime: ", pseudoprime)
print("bitlength: ", pseudoprime.nbits(), "digits ", pseudoprime.ndigits())

```
</details> 

## Final remarks

Python is much slower than C, so we can optimize this program by changing the language in which it is written. We also remark that different implementations of Baillie-PSW use slightly different versions of the Lucas test.

To add:
- a better description of Arnault's method (maybe)
- explain how we can use quadratic-reciprocity to improve MR with base 2.
- explain Lucas test.
- Pritchard's wheel sieve does not make program consistently faster.
- In the future we should implement the above code in C (or even GMP) to fully optimize. (for this check https://pari.math.u-bordeaux.fr/lcov-report/basemath/prime.c.func.html but be aware of the different versions!)
- Apparently GMP is even faster. And then it can be optimized for multithreading.


[^arnault]: Arnault, François - *Constructing Carmichael Numbers which are Strong Pseudoprimes to Several Bases* (1995) -  [Link](https://doi.org/10.1006/jsco.1995.1042).

[^prime]:  Albrecht, M., Massimo, J., Paterson, K.Somorovsky, J. - *Prime and Prejudice: Primality Testing Under Adversarial Conditions* (2018) - [Link](https://doi.org/10.1145/3243734.3243787)

