#import "template.typ": minimal-document
#import "eq.typ": eq, nonumber
#import "@preview/ouset:0.2.0": ouset, overset

#show: minimal-document.with(
  title: "Documentation for Matrixfuncs",
  authors: (
    (name: "Nextdorf", subtitle: link("https://github.com/nextdorf/matrixfuncs")[matrixfuncs\@github.com],),
  ),
  // Insert your abstract after the colon, wrapped in brackets.
  // Example: `abstract: [This is my abstract...]`
  // abstract: lorem(55),
  keywords: ("First keyword", "Second keyword", "etc."),
  date: "January 28, 2025",
)
#set cite(style: "chicago-author-date")



= Mathematical background
== 2 dimensional case
The Cayley-Hamilton theorem states for $A in CC^(n times n): p_A (A)=0$\
For $n=2$:

#eq[
$p_A (lambda)=&det(lambda bb(1)-A)=(lambda-lambda_1)(lambda-lambda_2)=lambda^2-lambda tr A+det A$

$=> A^2=&A tr A-bb(1)det A$
]
So $A^2$ is a linear combination of $A$ and $bb(1)$. By complete induction it follows that $A^n$ is also such a linear combination of $A$ and $bb(1)$ for any $n in NN_0$:

$
  A^n =: a_n A+b_n bb(1), quad A^(n+1)=(a_n tr A+b_n)A-(a_n det A)bb(1)
$
$
  a_(n+1) = a_n tr A+b_n, quad b_(n+1)=-a_n det A, quad a_0=0, quad a_1=1, quad b_0=1, quad b_1=0 
$ <recurrent-2dim>

In order to solve the recurrence equation we use a shifting operator $Delta$ for $a_n$, s.t. $a_(n+1) = Delta a_n$

$
  0=(Delta^2-Delta tr A+det A)a_n=p_A (Delta)a_n=(Delta-lambda_1)(Delta-lambda_2)a_n
$

As an Ansatz we choose $a_n$ as a linear combination of the solutions of $0=(Delta-lambda_1)a_n$ and $0=(Delta-lambda_2)a_n$. In general this ansatz only works if $lambda_1 != lambda_2$. We will consider the case $lambda_1 = lambda_2$ again in the end of this section and in @ch-math-general-case:

$
  a_n = c_1 lambda_1^n +c_2 lambda_2^n
$
As the base case for induction, we consider the trivial cases $n=0$ and $n=1$:

#eq(numberlast: true)[
  $  a_0 =& c_1+c_2=0, quad a_1 = c_1 lambda_1 + c_2 lambda c_2 = 1$
  $=> c_2 =& -c_1, quad c_1(lambda_1-lambda_2) = 1$
  $=> c_2 =& -c_1, quad c_1(lambda_1-lambda_2) = 1$
]

We found the solutions for $a_n$:

$
  a_n = (lambda_1^n-lambda_2^n)/(lambda_1-lambda_2), quad b_n ouset(=,#ref(<recurrent-2dim>)) -(lambda_2lambda_1^n - lambda_1lambda_2^n)/(lambda_1 - lambda_2)
$

and thereby the linear combination of $A^n$:

$
  A^n = (lambda_1^n-lambda_2^n)/(lambda_1-lambda_2)A -(lambda_2lambda_1^n - lambda_1lambda_2^n)/(lambda_1 - lambda_2)bb(1)
$

Since $A^n$ is linear in $lambda_i^n$ we know how to evaluate arbitrary polynomials of $A$ in the Krylov space. Let $q(x)$ be a polynomial:

$
  q(A) = 1/(lambda_1-lambda_2)[(q(lambda_1) - q(lambda_2))A - (lambda_2 q(lambda_1) - lambda_1 q(lambda_2))bb(1)]
$

Since the coefficients only depend on $p(lambda_i)$ the formula can be generalized to all functions which are analytical in $lambda_1$ and $lambda_2$. Let $f(x)$ be a analytical in the eigenvalues of $A$:

$
  f(A) = 1/(lambda_1-lambda_2)[(f(lambda_1) - f(lambda_2))A - (lambda_2 f(lambda_1) - lambda_1 f(lambda_2))bb(1)]
$ <func-2dim>
The generalization of @func-2dim is the core of this package. The important properties are that it yields the result to arbitrary precision in constant time and that the application of $f$ commutes with any linear function applied to the Krylov space. With other words this means there exists a rank-3 tensor $F$, which only depends on $A$, s.t. $f(A)_(i j) = F_(i j k) f(lambda_k)$.

Since $f(A)$ is continuous in the eigenvalues and the set of square matrices with unequal eigenvalues is dense in the set of all square matrices, we can just consider the limit of $lambda_1 -> lambda_2$ to calculate the case $lambda_1 = lambda_2$:

#eq(numberlast: true)[
  $f(A | lambda_1& = lambda + epsilon, lambda_2 = lambda) = 1/epsilon [epsilon f'(lambda)A - epsilon(lambda f'(lambda) - f(lambda))bb(1)] + o(epsilon)$
  $ouset(-->, epsilon -> 0, &) f'(lambda)(A - lambda bb(1)) + f(lambda)bb(1)$
  $=> f(A) =& f'(lambda)(A - lambda bb(1)) + f(lambda)bb(1) wide #[if $lambda$ is the only eigenvalue of $A$]$ <func-2dim-1eigval>
]

So when $A$ is not diagonal $f(A)$ depends on $f'(lambda)$. In general we will find that $f(A)$ depends on $f(lambda_i), f'(lambda_i), ..., f^((m))(lambda_i)$ where $m$ is the difference of the algebraic and geometric multiplicity of $lambda_i$.


== Examples for 2 dimensions
=== Generating $sin x$ from a difference equation

There are many applications and possible use cases for matrix functions, but the main application this package is written for are difference equations. Assuming you want to consider a sequence of vectors $(v_i)_(i=0)^inf subset FF^n$ which satisfies $v_(i+1) = M v_i$ for some matrix $M in FF^(n times n)$. An obvious conclusion is $v_k = M^k v_0$. With equations @func-2dim or @func-2dim-1eigval we can express $M^k$ without any matrix multiplications. Instead with just apply $f(lambda) = lambda^k$ to all eigenvalues. We even can apply $v_0$ or any other tensor contraction to $f(M)$ before we specify $k$. Calculating these things beforehand can speedup the total calculation if we need to calculate $v_k$ for many different $k$.

Another neat application is the analytical continuation of $v_k$ in $k$ by setting $f(lambda) = e^(k ln(lambda))$. We will use this in the example here by solving a difference equation for a sampled $sin$ function and then evaluating the numerically solved $sin$ between the sampled points.

Consider the sum identity of $sin$:

$
  sin(x+a) = sin(x)cos(a) + cos(x)sin(a)
$ <sin-ab>

For fixed $a$ $sin(x+a)$ can be expressed as a linear combination of the basis $sin x$ and $cos b$. Therefore $sin x$ can be written as a linear combination of $sin(x+a)$ and $sin(x+2a)$ for almost all $a$:

#eq(numberlast: none)[
  $sin(x) =:& alpha sin(x+a) + beta sin(x+2a)$
  $=& alpha(sin(x)cos(a) + sin(a)cos(x)) + beta(sin(x)cos(2a) + sin(2a)cos(x))$
  $=& sin x(alpha cos a + beta cos(2a)) + cos x(alpha sin(a) + beta sin(2a))$
  $=& vec(sin x, cos x) mat(cos a, cos(2a); sin a, sin(2a)) vec(alpha, beta)$
  $=> vec(alpha, beta) =& mat(cos a, cos(2a); sin a, sin(2a))^(-1) vec(1, 0)$
  $ouset(=, #ref(<sin-ab>), &) vec(2cos(a), -1)$
]

From this we can construct a difference equation which shifts $sin$ by $a$.

#eq(numberlast: true)[
  $vec(sin x, sin(x-a)) =& underbrace(mat(2cos a, -1; 1, 0), =: M(a)) vec(sin(x-a), sin(x-2a))$
  $=> vec(sin(x+n a), sin(x+(n-1)a)) =& M(a)^n vec(sin x, sin(x-a))$
  $=> sin(x+n a) =& hat(e)_1 dot M(a)^n vec(sin x, sin(x-a))$
  $=> sin(x+y) =& hat(e)_1 dot exp(y/a ln(M(a))) vec(sin x, sin(x-a))$
]

We can now use the formula for two dimensions @func-2dim in order to solve for $sin(x+y)$:


$
  sin(x+y) ouset(=, #ref(<func-2dim>)) hat(e)_1 dot 1/(lambda_1-lambda_2)&[(f(lambda_1) - f(lambda_2))M(a) - (lambda_2 f(lambda_1) - lambda_1 f(lambda_2))bb(1)] vec(sin x, sin(x-a))
$ <sin-ab-numerical>

#eq(numberlast: none)[
  $#[with ] f(lambda) =& exp(y/a ln(lambda))$
  $#[and ] lambda_(1,2) =& 1/2 tr M(a) plus.minus sqrt((1/2 tr M(a))^2 - det M(a))$
  $=& cos a plus.minus sqrt(cos^2 a - 1)$
  $=& e^(plus.minus i a)$
  $#[s.t. ] f(lambda_(1,2)) =&e^(plus.minus i y)$
]
Simplifying @sin-ab-numerical will yield @sin-ab.


== General case <ch-math-general-case>

Similarly to the 2d-case we use the Cayley-Hamilton-theorem: $A in CC^(n times n), p_A (A) = 0$

#eq[
  #nonumber($p_A (lambda) =& det(lambda bb(1) - A) = product_(k=1)^n (lambda - lambda_k) =: lambda^n - sum_(k=0)^(n-1) Lambda_k lambda^k$)
  $=> Lambda_k =& sum_(j_1...j_(n-k)=1,\ j_1<...<j_(n-k))^n (-1)^(n-k+1) product_(i=1)^(n-k) lambda_(j_i)$
  #nonumber($A^n =& sum_(k=0)^(n-1) Lambda_k A^k$)
  $A^m =:& sum_(k=0)^(n-1) alpha_(m k) A^k$
  #nonumber($A^(m+1) =& sum_(k=1)^(n-1)(alpha_(m,k-1) + alpha_(m, n-1) Lambda_k)A^k + alpha_(m, n-1) Lambda_0 A^0$)
  $=> alpha_(m+1, 0) =& alpha_(m,n-1) Lambda_0, quad alpha_(m+1, k) = alpha_(m, k-1) + alpha_(m, n-1) Lambda_k$ <recurrent-ndim>
]

The recurrence equation @recurrent-ndim can be solved by noticing that all $alpha_(m, n-1)$ can be reduced to a function of $alpha_(m-1, n-1), alpha_(m-2, n-1), ...$

#eq(numberlast: true)[
  $alpha_(m,n-1) =& sum_(k=1)^(n-1)alpha_(m-k,n-1)Lambda_(n-k) + alpha_(m+1-n, 0)$
  $=& sum_(k=1)^n alpha_(m-k,n-1)Lambda_(n-k)$
  $=& sum_(k=0)^(n-1) alpha_(m-n+k,n-1)Lambda_k$
  $=> 0=& alpha_(m+n, n-1) - sum_(k=0)^(n-1) alpha_(m-k,n-1)Lambda_k$
  $=& (Delta^n - sum_(k=0)^(n-1) Delta^k Lambda_k) alpha_(m, n-1) wide #[with ] Delta alpha_(m k) = alpha_(m+1, k)$
  $=& p_A (Delta) alpha_(m, n-1)$
  $=& product_(k=1)^r (Delta - lambda_k)^(mu_k) alpha_(m, n-1)\ &#[with $mu_k$ being the algebraic multiplicity of $lambda_k$]$
]

The general solution is $alpha_(m, n-1) = sum_(k=1)^r lambda_k^m p_k (m)$ with $p_k$ being a polynomial of degree $mu_(k-1)$.

Proof:

#eq(numberlast: none)[
  $#[Induction Start: ] 0 =& (Delta - lambda)c_n => lambda c_(n-1) = lambda^n c_0$
  $#[Assume ] 0 =& (Delta - lambda)^m lambda^n sum_(k=0)^(m-1) c_k n^k quad forall c_0, ..., c_(m-1)$
  $=> (Delta - lambda)^(m-1) lambda^n sum_(k=0)^m c_k n^k =& (Delta - lambda)^m (Delta - lambda) lambda^n c_m n^m + (Delta - lambda)(Delta - lambda)^m lambda^n sum_(k=0)^(m-1) c_k n^k$
  $=& (Delta - lambda)^m ((n+1)^m - n^m) lambda^(n+1) c_m$
  $=& (Delta - lambda)^m lambda^n sum_(k=0)^(m-1) binom(m, k) lambda c_m n^k$
  $=& 0$
  $=> 0 =& (Delta - lambda)^m lambda^n sum_(k=0)^m' c_k n^k quad forall c_0, ..., c_(m'-1), m > m'$
]

Now consider $overline(c) = (Delta - lambda) c_n => c_(n+1) = overline(c) + lambda c_n = sum_(k=0)^n lambda^k overline(c)_(n-k) + lambda^(n+1) c_0$

Since the solution of $c_n$ is linear in $overline(c)_n$ we can consider the dimension of the solution space if $overline(c)_n$ is itself a solution of a similar equation. The solution space of $0 = (Delta - lambda_1) c_n^((1))$ is 1-dimensional. Therefore the solution space of $c_n^((m)) = (Delta - lambda_1) c_n^((m+1))$ is either of the same dimension as the solution space of $c_n^((m))$ or the dimension increases by 1. So the dimension of the solution space of $p_A (Delta) alpha_(m, n-1)$ is at most $n$. Since $sum_{k=1}^r lambda_k^m p_k(m)$ is a $sum_(k=1)^r dim(p_k) = n$ dimensional solution it is the general solution.

#place(right, $qed$)
\ 
#eq[
  $alpha_(m, n-1) =& sum_(k=1)^r sum_(l=0)^(min(mu_k-1, m)) overline(beta)_(k l) lambda_k^(m-l) m!/((m-l)!) = sum_(k=1)^r sum_(l=0)^(min(mu_k-1, m)) overline(beta)_(k l) partial_lambda^l lambda_k^m bar_(lambda=lambda_k) =: sum_(k=1)^n beta_k lambda_k^((m))$
  #nonumber($beta_k =& (overline(beta)_10, ..., overline(beta)_(1,mu_1-1), overline(beta)_20, ..., overline(beta)_(2,mu_2-1), ..., overline(beta)_(r 0), ..., overline(beta)_(r,mu_r-1))_k$)
  #nonumber($lambda_k^((m)) =& (lambda_1^m, ..., partial_(lambda_1)^(mu_1-1) lambda_1^m, lambda_2^m, ..., partial_(lambda_2)^(mu_2-1) lambda_2^m, ..., lambda_r^m, ..., partial_(lambda_r)^(mu_r-1) lambda_r^m)_k$)
]

For $m<n: alpha_(m k) = delta_(m k) => delta_(m,n-1) = beta_k lambda_k^((m))$

#eq(numberlast: true)[
  $arrow(beta) perp& sum_(k=1)^n hat(e)_k lambda_k^m =: arrow(lambda)^((m)) #[ for ] m = 0 ... n-2$
  $tilde(beta)_k :=& det(hat(e)_k | arrow(lambda)^((n-2)) | arrow(lambda)^((n-3)) | ... | arrow(lambda)^((0)))$
  $beta =& tilde(beta) / (tilde(beta) dot arrow(lambda)^((n-1)))$
]

@recurrent-ndim yields for $alpha_(m k)$:

$
  alpha_(m k) = sum_(j=1)^(k+1) alpha_(m-j, n-1) Lambda_(k+1-j) = sum_(j=1)^(k+1) arrow(beta) dot arrow(lambda)^((n-1)) Lambda_(k+1-j)
$

As in the two dimensional case this defines $f(A)$ if $f$ is analytical in the eigenvalues of $A$. If none of the eigenvalues is $0$:

#eq[
  $A^m =& sum_(k=0) A^k sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(l=1)^r sum_(p=0)^(min(mu_l-1, m-j)) overline(beta)_(l p) lambda_l^(m-j-p) (m-j)!/((m-j-p)!)$
  #nonumber($f(A) =& sum_(k=0) A^k sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(l=1)^r sum_(p=0)^(mu_l-1) overline(beta)_(l p) partial_(lambda_l)^p lambda_l^(-j) f(lambda_l)$)
  $f(A) =& sum_(k=0) A^k sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(l=1)^r sum_(p=0)^(mu_l-1) overline(beta)_(l p) sum_(q=0)^p binom(p, q) (-1)^(p-q) (j-1+p-q)!/((j-1)!) lambda_l^(-j-p+q) f^((q))(lambda_l)$
]

Since $f(A)$ is linear in $A^k$ and $f^((q))(lambda_l)$ we can define the tensors $phi$ and $phi.alt$, s.t.

$
  f(A) = sum_(i=0)^(n-1) sum_(j=1)^r sum_(k=0)^(mu_j-1) phi_(i j)^((k)) A^i f^((k))(lambda_j), quad phi.alt_(i j k)^((l)) = sum_(m=0)^(n-1) phi_(m k)^((l))(A^m)_(i j)
$
#eq(numberlast: true)[
  $phi_(k l)^((q)) =& partial_(x_q) partial_(y_l) sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(l'=1)^r sum_(p=0)^(mu_(l')-1) overline(beta)_(l' p) sum_(q'=0)^p binom(p, q') (-1)^(p-q') (j-1+p-q')!/((j-1)!) lambda_(l')^(-j-p+q') x_(q') y_(l')$
  $=& sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(p=q)^(mu_l-1) overline(beta)_(l p) binom(p, q) (-1)^(p-q) (j-1+p-q)!/((j-1)!) lambda_l^(-j-p+q)$
]

=== Sanity check for $n=2$

#list(
[
  $lambda_1 != lambda_2 => mu_1 = mu_2 = 1, r = 2:$

  #math.equation(block: true, numbering: _ => "", [
    // $lambda_k^((m)) =& lambda_k^m, quad beta_(dot 0) = binom(0, 1), quad Lambda_0 = -lambda^2, quad Lambda_1 = 2lambda$\
    $lambda_k^((m)) =& lambda_k^m, quad beta_(dot 0) = 1/(lambda_1 - lambda_2)binom(1, -1), quad Lambda_0 = -lambda_1 lambda_2_, quad Lambda_1 = lambda_1 + lambda_2$\
    $phi^((0)) =& mat(Lambda_0 overline(beta)_10 slash lambda_1, Lambda_0 overline(beta)_20 slash lambda_2; Lambda_0 overline(beta)_10 slash lambda_1^2 + Lambda_1 overline(beta)_10 slash lambda_1, Lambda_0 overline(beta)_20 slash lambda_2^2 + Lambda_1 overline(beta)_20 slash lambda_2) = 1/(lambda_1 - lambda_2) mat(-lambda_2, lambda_1; 1, -1)$\
    $=> f(A) =& 1/(lambda_1 - lambda_2)[(f(lambda_1) - f(lambda_2))A - (lambda_2 f(lambda_1) - lambda_1 f(lambda_2))bb(1)] wide checkmark$
  ])
], [
  $lambda_1 = lambda_2 = lambda => mu_1 = 2, r = 1:$

  #math.equation(block: true, numbering: _ => "", [
    $lambda_k^((m)) =& (lambda^m, m lambda^(m-1))_k, quad beta_(dot 1) = binom(0, 1), quad Lambda_0 = -lambda^2, quad Lambda_1 = 2lambda$\
    $phi^((0)) =& mat(Lambda_0(overline(beta)_10 slash lambda - overline(beta)_11 slash lambda^2); Lambda_1(overline(beta)_10 slash lambda - overline(beta)_11 slash lambda^2) + Lambda_0(overline(beta)_10 slash lambda^2 - 2overline(beta)_11 slash lambda^3)) = mat(1; 0)$\
    $phi^((1)) =& mat(Lambda_0 overline(beta)_11 slash lambda; Lambda_1 overline(beta)_11 slash lambda + Lambda_0 overline(beta)_11 slash lambda^2) = mat(-lambda; 1)$\
    $=> f(A) =& f(lambda)bb(1) + f'(lambda)(A - lambda bb(1)) wide checkmark$
  ])
])


