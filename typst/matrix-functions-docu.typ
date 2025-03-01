#import "template.typ": minimal-document, appendices
#import "eq.typ": eq, nonumber
#import "@preview/ouset:0.2.0": ouset, overset
#import "@preview/wrap-it:0.1.1": wrap-content
#import "@preview/codelst:2.0.2": sourcecode, sourcefile

#show: minimal-document.with(
  title: "Documentation for Matrixfuncs",
  authors: (
    (
      name: "Nextdorf",
      // subtitle: link("https://github.com/nextdorf/matrix-functions")[matrixfuncs\@github.com],  
      subtitle: link("https://github.com/nextdorf/matrix-functions")[
        #grid(columns: (auto, auto, auto), rows: 0pt, align: alignment.bottom, [matrixfuncs\@ #h(.05em)], [], [#figure(image("github-logo.svg", height: .9em))])
      ],
    ),
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
The Cayley-Hamilton theorem states that any square matrix $A in CC^(n times n)$ satisfies its own characteristic equation. If $p_A (lambda)=det(lambda bb(1) - A)$ is the characteristic polynomial of $A$, then substituting $A$ into it results in $p_A (A)=0$.

For $n=2$:

#eq[
$p_A (lambda)=det(lambda bb(1)-A)=(lambda-lambda_1)(lambda-lambda_2)=lambda^2-lambda tr A+det A$

$=> A^2=A tr A-bb(1)det A$
]
So $A^2$ is a linear combination of $A$ and $bb(1)$. By multiplication with $A$ we find, that $A^3$ is a linear combination of $A^2$ and $A$ which reduces according to the above relation to a linear combination of $A$ and $bb(1)$. By repeatedly  multiplying with $A$, we find equivalent results for $A^4$, $A^5$, ..., as well. By complete induction it follows that $A^n$ is also such a linear combination of $A$ and $bb(1)$ for any $n in NN_0$ (the base case of $n=0, 1$ is trivial):

$
  A^n =: a_n A+b_n bb(1), quad A^(n+1)=(a_n tr A+b_n)A-(a_n det A)bb(1)
$
$
  a_(n+1) = a_n tr A+b_n, quad b_(n+1)=-a_n det A,\ quad a_0=0, quad a_1=1, quad b_0=1, quad b_1=0 
$ <recurrent-2dim>

In order to solve the recurrence equation we use a shift operator $Delta$ for $a_n$, s.t. $a_(n+1) = Delta a_n$

$
  0=(Delta^2-Delta tr A+det A)a_n=p_A (Delta)a_n=(Delta-lambda_1)(Delta-lambda_2)a_n
$

As an Ansatz we choose $a_n$ as a linear combination of the solutions of $0=(Delta-lambda_1)a_n$ and $0=(Delta-lambda_2)a_n$. We are motivated to pick this Ansatz because clearly any solution to $0=(Delta-lambda_(1 #[ or ] 2))a_n$ solves the above equation and the $Delta^2$ hints at a 2-dimensional solution-space. In addition the problem is linear. This approach assumes distinct eigenvalues $lambda_1, lambda_2$. When the eigenvalues coincide, additional techniques are required, which will be discussed at the end of this section and in @ch-math-general-case:

$
  a_n = c_1 lambda_1^n +c_2 lambda_2^n
$
As the base case for induction, we consider the trivial cases $n=0$ and $n=1$:

#eq(numberlast: true)[
  $  a_0 =& c_1+c_2=0, quad a_1 = c_1 lambda_1 + c_2 lambda c_2 = 1$
  $=> c_2 =& -c_1, quad c_1(lambda_1-lambda_2) = 1$
  $=> c_1 =& 1/(lambda_1-lambda_2), quad c_2 = (-1)/(lambda_1-lambda_2)$
]

Thus, the explicit closed-form solutions for $a_n$ and $b_n$ are:

$
  a_n = (lambda_1^n-lambda_2^n)/(lambda_1-lambda_2), quad b_n ouset(=,#ref(<recurrent-2dim>)) -(lambda_2lambda_1^n - lambda_1lambda_2^n)/(lambda_1 - lambda_2)
$

Thus, we obtain the explicit linear combination representation of $A^n$:

$
  A^n = (lambda_1^n-lambda_2^n)/(lambda_1-lambda_2)A -(lambda_2lambda_1^n - lambda_1lambda_2^n)/(lambda_1 - lambda_2)bb(1)
$

Since $A^n$ is linear in $lambda_i^n$ we know how to evaluate arbitrary polynomials of $A$. Let $q(x)$ be a polynomial:

$
  q(A) = 1/(lambda_1-lambda_2)[(q(lambda_1) - q(lambda_2))A - (lambda_2 q(lambda_1) - lambda_1 q(lambda_2))bb(1)]
$

Since the coefficients only depend on $q(lambda_i)$ the formula can be generalized to all functions which are analytic in $lambda_1$ and $lambda_2$. Let $f(x)$ be an analytic function in the eigenvalues of $A$:

$
  f(A) = 1/(lambda_1-lambda_2)[(f(lambda_1) - f(lambda_2))A - (lambda_2 f(lambda_1) - lambda_1 f(lambda_2))bb(1)]
$ <func-2dim>
The generalization of @func-2dim is the core of this package. The important properties are that it allows efficient computation of $f(A)$ with arbitrary precision and that the application of $f$ commutes with any linear function applied on the vector space. With other words this means there exists a rank-3 tensor $F$, which only depends on $A$, s.t. $f(A)_(i j) = F_(i j k) f(lambda_k)$.

Since $f(A)$ is continuous in the eigenvalues and the set of square matrices with distinct eigenvalues is dense, we evaluate $f(A)$ for equal eigenvalues using a limiting process:

#eq(numberlast: true)[
  $f(A | lambda_1& = lambda + epsilon, lambda_2 = lambda) = 1/epsilon [epsilon f'(lambda)A - epsilon(lambda f'(lambda) - f(lambda))bb(1)] + o(epsilon)$
  $ouset(-->, epsilon -> 0, &) f'(lambda)(A - lambda bb(1)) + f(lambda)bb(1)$
  $=> f(A) =& f'(lambda)(A - lambda bb(1)) + f(lambda)bb(1) wide #[if $lambda$ is the only eigenvalue of $A$]$ <func-2dim-1eigval>
]

So when $A$ is not diagonal $f(A)$ depends on $f'(lambda)$. In general we will find that $f(A)$ depends on $f(lambda_i), f'(lambda_i), ..., f^((m))(lambda_i)$ where $m$ is the difference of the algebraic and geometric multiplicity of $lambda_i$.


== Applications and Examples in 2 dimensions
=== Comparison to Jordan form approach
...
=== Relation to Krylov subspaces <ch-krylov>
...
=== Generating $sin x$ from a difference equation <ch-sin-generation>

Matrix functions have many applications, but the primary focus of this package is solving difference equations. Assuming you want to consider a multi-dimensional first-order linear recurrence, i.e. a sequence of vectors $(v_i)_(i=0)^infinity subset FF^n$ satisfying the recurrence relation $v_(i+1) = M v_i$ for some matrix $M in FF^(n times n)$. An obvious conclusion is $v_k = M^k v_0$. Using equations @func-2dim or @func-2dim-1eigval, we can express $M^k$ without any matrix multiplications. Instead we just need to compute $f(lambda) = lambda^k$ for all eigenvalues. As discussed in @ch-krylov, we can precompute $f(M)v_0$ or any other tensor contraction to $f(M)$ before we specify $k$, allowing for efficient evaluation for different or delayed values of $k$.

Another neat application is the analytic continuation of $v_k$ in $k$ by setting $f(lambda) = e^(k ln(lambda))$, allowing for non-integer shifts. We apply this method in the example here by solving a difference equation for a sampled $sin$ function and then evaluating the numerically solved $sin$ between the sampled points.

For the setup of the difference equation, consider the sum identity of $sin$:

$
  sin(x+a) = sin(x)cos(a) + cos(x)sin(a)
$ <sin-ab>

For fixed $a$ $sin(x+a)$ can be expressed as a linear combination of the basis $sin x$ and $cos x$. As shown in the following derivation, $sin x$ can be written as a linear combination of $sin(x+a)$ and $sin(x+2a)$ for almost all $a$:

#eq(numberlast: none)[
  $sin(x) =:& alpha sin(x+a) + beta sin(x+2a)$
  $=& alpha(sin(x)cos(a) + sin(a)cos(x)) + beta(sin(x)cos(2a) + sin(2a)cos(x))$
  $=& sin x(alpha cos a + beta cos(2a)) + cos x(alpha sin(a) + beta sin(2a))$
  $=& vec(sin x, cos x) mat(cos a, cos(2a); sin a, sin(2a)) vec(alpha, beta)$
  $=> vec(alpha, beta) =& mat(cos a, cos(2a); sin a, sin(2a))^(-1) vec(1, 0)$
  $ouset(=, #ref(<sin-ab>), &) vec(2cos(a), -1)$
]

From this we can construct a difference equation which shifts $sin x$ by $a$.

#eq(numberlast: true)[
  $vec(sin x, sin(x-a)) =& underbrace(mat(2cos a, -1; 1, 0), =: M(a)) vec(sin(x-a), sin(x-2a))$
  $=> vec(sin(x+n a), sin(x+(n-1)a)) =& M(a)^n vec(sin x, sin(x-a))$
  $=> sin(x+n a) =& hat(e)_1 dot M(a)^n vec(sin x, sin(x-a))$
  $=> sin(x+y) =& hat(e)_1 dot exp(y/a ln(M(a))) vec(sin x, sin(x-a))$
]

Applying the 2D recurrence formula from @func-2dim, we compute $sin(x+y)$ as follows:


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
Simplifying @sin-ab-numerical yields @sin-ab.


#wrap-content(
  figure(image("../sinFromMFunc.svg")),
  align: right,
)[
#box(height: 1em)

The example code at #link("https://github.com/nextdorf/matrix-functions/blob/main/example.py")["/example.py"] illustrates how to do the above computation using the library. It generates the figure on the right:

For my examples check out the snippets at #link("https://github.com/nextdorf/matrix-functions/tree/main/examples")["/examples/"] or at @ch-code-examples.
]


== General case <ch-math-general-case>

Similarly to the 2d-case we use the Cayley-Hamilton-theorem: $A in CC^(n times n), p_A (A) = 0$

#eq[
  #nonumber($p_A (lambda) =& det(lambda bb(1) - A) = product_(k=1)^n (lambda - lambda_k) =: lambda^n - sum_(k=0)^(n-1) Lambda_k lambda^k$)
  $=> Lambda_k =& sum_(j_1...j_(n-k)=1,\ j_1<...<j_(n-k))^n (-1)^(n-k+1) product_(i=1)^(n-k) lambda_(j_i)$
  $A^n =& sum_(k=0)^(n-1) Lambda_k A^k$ <matrix-power-lin-n>
  $A^m =:& sum_(k=0)^(n-1) alpha_(m k) A^k$ <matrix-power-lin-ansatz>
]

As before we multiply equation @matrix-power-lin-ansatz with $A$ which will generate an $A^n$ on the rhs. That generated term is substituted using @matrix-power-lin-n, resulting in a recurrence relation for $alpha_(m k)$:

#eq[
  #nonumber($A^(m+1) =& sum_(k=1)^(n-1)(alpha_(m,k-1) + alpha_(m, n-1) Lambda_k)A^k + alpha_(m, n-1) Lambda_0 A^0$)
  $=> alpha_(m+1, 0) =& alpha_(m,n-1) Lambda_0, quad alpha_(m+1, k) = alpha_(m, k-1) + alpha_(m, n-1) Lambda_k$ <recurrent-ndim>
]

In order to solve the recurrence equation @recurrent-ndim, note that $alpha_(m k)$ can be expressed in terms of $alpha_(m-1, k-1)$ and $alpha_(m-1, n-1)$. $alpha_(m-1, k-1)$ can be further expressed in terms of $alpha_(m-2, k-2)$ and $alpha_(m-2, n-1)$, and so forth, until eventually $alpha_(m-k, 0)$ can be expressed in terms of $alpha_(m-k-1, n-1)$. Hence, the recurrence equation @recurrent-ndim can be solved by noticing that all $alpha_(m, n-1)$ can be reduced to a function of $alpha_(m-1, n-1), alpha_(m-2, n-1), ...$. Thus, we can simplify the recurrence relation to a recurrence relation in the first argument only and fixing the second argument to $n-1$:

#eq(numberlast: true)[
  $alpha_(m,n-1) =& sum_(k=1)^(n-1)alpha_(m-k,n-1)Lambda_(n-k) + alpha_(m+1-n, 0)$
  $=& sum_(k=1)^n alpha_(m-k,n-1)Lambda_(n-k)$
  $=& sum_(k=0)^(n-1) alpha_(m-n+k,n-1)Lambda_k$
  $=> 0=& alpha_(m+n, n-1) - sum_(k=0)^(n-1) alpha_(m+k,n-1)Lambda_k$
  $=& (Delta^n - sum_(k=0)^(n-1) Delta^k Lambda_k) alpha_(m, n-1) wide #grid(rows: (1.25em, 0em), align: left, $#[with ] Delta alpha_(m k) = alpha_(m+1, k)$, [i.e. $Delta$ acts on the first index])$
  //#[with $Delta alpha_(m k) = alpha_(m+1, k)$ i.e. $Delta$ acts on the first index]$
  $=& p_A (Delta) alpha_(m, n-1)$
  $=& product_(k=1)^r (Delta - lambda_k)^(mu_k) alpha_(m, n-1)\ &#[with $mu_k$ being the algebraic multiplicity of $lambda_k$]$
]

The general solution is $alpha_(m, n-1) = sum_(k=1)^r lambda_k^m p_k (m)$ with $p_k$ being an arbitrary polynomial of degree $mu_k-1$.

Proof:

#eq(numberlast: none)[
  $#[Induction Start: ] 0 =& (Delta - lambda)c_n => c_n =lambda c_(n-1) = lambda^n c_0$
  $#[Assume ] 0 =& (Delta - lambda)^m lambda^n sum_(k=0)^(m-1) c_k n^k quad forall c_0, ..., c_(m-1)$
  $=> (Delta - lambda)^(m+1) lambda^n sum_(k=0)^m c_k n^k =& (Delta - lambda)^m (Delta - lambda) lambda^n c_m n^m + (Delta - lambda)(Delta - lambda)^m lambda^n sum_(k=0)^(m-1) c_k n^k$
  $=& (Delta - lambda)^m ((n+1)^m - n^m) lambda^(n+1) c_m$
  $=& (Delta - lambda)^m lambda^n sum_(k=0)^(m-1) binom(m, k) lambda c_m n^k$
  $=& 0$
  $=> 0 =& (Delta - lambda)^m lambda^n sum_(k=0)^m' c_k n^k quad forall c_0, ..., c_(m'), m > m'$
]

What is left to show is that the above solution is general. Considering $(Delta - lambda)^n$ as a linear operator, the solution space we are interested in is effectively the operator's kernel. Thus, the above $n$-dimensional solution is general if the dimension of the operator's kernel is $n$ as well.

Now consider $overline(c)_n = (Delta - lambda) c_n => c_(n+1) = overline(c)_n + lambda c_n = sum_(k=0)^n lambda^k overline(c)_(n-k) + lambda^(n+1) c_0$

Since the solution of $c_(n+1)$ is linear in $overline(c)_n$ we can consider the dimension of the solution space if $overline(c)_n$ is itself a solution of a similar equation. The solution space of $0 = (Delta - lambda_1) c_n^((1))$ is 1-dimensional. Therefore the solution space of $c_n^((m)) = (Delta - lambda_1) c_n^((m+1))$ is either of the same dimension as the solution space of $c_n^((m))$ or the dimension increases by 1. So the dimension of the solution space of $p_A (Delta) alpha_(m, n-1)$ is at most $n$. Since $sum_(k=1)^r lambda_k^m p_k(m)$ is a $sum_(k=1)^r dim(p_k) = n$ dimensional solution it is the general solution.

#place(right, $qed$)
\ 
#eq[
  $alpha_(m, n-1) =& sum_(k=1)^r sum_(l=0)^(min(mu_k-1, m)) overline(beta)_(k l) lambda_k^(m-l) m!/((m-l)!) = sum_(k=1)^r sum_(l=0)^(min(mu_k-1, m)) overline(beta)_(k l) partial_lambda^l lambda_k^m bar_(lambda=lambda_k) =: sum_(k=1)^n beta_k lambda_k^((m))$
  #nonumber($beta_k =& (overline(beta)_10, ..., overline(beta)_(1,mu_1-1), overline(beta)_20, ..., overline(beta)_(2,mu_2-1), ..., overline(beta)_(r 0), ..., overline(beta)_(r,mu_r-1))_k$)
  #nonumber($lambda_k^((m)) =& (lambda_1^m, ..., partial_(lambda_1)^(mu_1-1) lambda_1^m, lambda_2^m, ..., partial_(lambda_2)^(mu_2-1) lambda_2^m, ..., lambda_r^m, ..., partial_(lambda_r)^(mu_r-1) lambda_r^m)_k$)
]

For $m<n: alpha_(m k) = delta_(m k) => delta_(m,n-1) = beta_k lambda_k^((m))$

$
  => arrow(beta) perp& sum_(k=1)^n hat(e)_k lambda_k^m =: arrow(lambda)^((m)) #[ for ] m = 0 ... n-2, quad #[and ] arrow(beta) dot arrow(lambda)^(n-1) = 1
$

// There are several ways to construct $beta$. Here we will list two ways. The first way relies on some properties of the determinant. It is the approach used in an older version of this package and is often more useful for analytical calculations. It introduces however precision-errors for numerical calculations and becomes unstable for large dimensions and/or eigenvalues.

// This is why in a more recent version, the change was made to directly finding the orthogonal complement of the space spanned by $arrow(lambda)^((0)), ..., arrow(lambda)^((n-2))$. We expect the orthogonal complement to be of dimension 1, s.t. $beta$'s relation to $arrow(lambda)^((n-1))$ uniquely defines $beta$. In both cases we rely on the fact that $arrow(lambda)^((0)), ..., arrow(lambda)^((n-1))$ are linear independent.

// #list([
//   *Construction via the determinant*
  
  If the rows (or columns) of a matrix are linearly dependent then the matrix is not invertable and its determinant will vanish. On the flip side if the entries are linear independent then the matrix's determinant will be non-zero. Thus, for a set of independent vectors ${w_(k+1), ..., w_n}$ an orthogonal tensor $T$ with
  
  #nonumber($T_(m_1...m_k) = det(hat(e)_(m_1) | ... | hat(e)_(m_k) | w_(k+1) | ... | w_n)$)

  can be constructed. If contracted with arbitrary vectors $(w_1, ..., w_k)$ then the orthogonal tensor will yield $det(w_1 | ... | w_n)$ which will be 0 iff the vectors $w_1, ..., w_n$ are linear dependent. This means that $T$ is non-zero, and means in particular for $k=1$ that $T$ becomes the tensor perpendicular to the vectors $w_2, ..., w_n$

  Since $arrow(lambda)^((0)), ..., arrow(lambda)^((n-1))$ are linear independent, we find:

  #math.equation(block: true, /*numbering: n => "("+[#n]+"a)"*/)[
    $beta_k =& det(hat(e)_k | arrow(lambda)^((n-2)) | arrow(lambda)^((n-3)) | ... | arrow(lambda)^((0))) / det(arrow(lambda)^((n-1)) | ... | arrow(lambda)^((0)))$
  ]
// ], [
//   *Subtracting the span*

//   As elegant the above approach might be, for numerical problems it turns out to be more stable to just spanning the embedding $n$-dimensional vector space by some vectors and then subtracting the spanned space from the embedding space. The subtracted vectors will then only span the orthogonal complement: 

//   Consider some arbitrary vectors $(v_1, v_2, ...)$ with $"span"(v_1, v_2, ...) = FF^n$, and define:

//   #nonumber($v'_i =& v_i - sum_(m=0)^(n-2) arrow(lambda)^((m)) (arrow(lambda)^((m)) dot v_i)/(arrow(lambda)^((m)) dot arrow(lambda)^((m)))$)

//   #math.equation(block: true, numbering: n => "("+[#(n - 1)]+"b)")[
//     $=> "span"(arrow(beta)) =& "span"(v'_1, v'_2, ...)$
//   ]

//   The condition $arrow(beta) dot arrow(lambda)^((n-1)) ouset(=,!) 1$ fixes the scale of $beta$.
// ])
// #counter(math.equation).update(n => n - 1)

@recurrent-ndim yields for $alpha_(m k)$:

$
  alpha_(m k) = sum_(j=1)^(k+1) alpha_(m-j, n-1) Lambda_(k+1-j) = sum_(j=1)^(k+1) arrow(beta) dot arrow(lambda)^((m-j)) Lambda_(k+1-j)
$

As in the two dimensional case this defines $f(A)$ if $f$ is analytic in the eigenvalues of $A$:

#eq[
  $A^m =& sum_(k=0) A^k sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(l=1)^r sum_(p=0)^(min(mu_l-1, m-j)) overline(beta)_(l p) lambda_l^(m-j-p) (m-j)!/((m-j-p)!)$
  #nonumber($f(A) =& sum_(k=0) A^k sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(l=1)^r sum_(p=0)^(mu_l-1) overline(beta)_(l p) partial_(lambda_l)^p lambda_l^(-j) f(lambda_l)$)
  $=& sum_(k=0) A^k sum_(j=1)^(k+1) Lambda_(k+1-j) sum_(l=1)^r sum_(p=0)^(mu_l-1) overline(beta)_(l p) sum_(q=0)^p binom(p, q) (-1)^(p-q) (j-1+p-q)!/((j-1)!) lambda_l^(-j-p+q) f^((q))(lambda_l)$
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

// == Vector space properties
// #let sdot(dual, vec) = $angle.l dual, vec angle.r$
// Scalar product:
// $
//   sdot(A, B) = 1/n tr(A^dagger B), quad norm(A) = sqrt(sdot(A, A))
// $
// ONS:
// $
//   N_0 = bb(1), quad N_(k+1) = (A N_k - sum_(j=0)^k N_j sdot(N_j, A N_k))/norm(A N_k - sum_(j=0)^k N_j sdot(N_j, A N_k))
// $

= Code Examples <ch-code-examples>
== Example from @ch-sin-generation

#sourcefile(read("../example.py"), file: "example.py",lang: "python")
#figure(image("../sinFromMFunc.svg")),

== More advanced Fcuntion Continution from random Samples

#sourcefile(read("../examples/many_frequencies.py"), file: "example.py",lang: "python")
#figure(image("../examples/many_frequencies.svg"))

// #show: appendices

