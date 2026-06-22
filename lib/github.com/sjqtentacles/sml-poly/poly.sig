(* poly.sig

   Univariate polynomials over an arbitrary commutative ring.

   The library is built around a `Poly` functor: given any structure matching
   `RING`, `Poly(R)` produces polynomials whose coefficients are drawn from
   `R`. Crucially, `Poly(R)` itself matches `RING`, so the construction nests:
   `Poly(Poly(R))` gives bivariate polynomials.

   For coefficient *fields* (where division is total except by zero), the
   `PolyField` functor additionally provides Euclidean division (`divMod`) and
   the polynomial GCD. *)

(* A commutative ring with unit. `equal` is required because polynomial
   normalization (dropping trailing zero coefficients) needs to detect the
   ring's zero, and the built-in `=` cannot be used on abstract types. *)
signature RING =
sig
  type t
  val zero  : t
  val one   : t
  val add   : t * t -> t
  val neg   : t -> t
  val mul   : t * t -> t
  val equal : t * t -> bool
  val toString : t -> string
end

(* A field: a ring in which every nonzero element has a multiplicative
   inverse. `inv` may raise on zero. *)
signature FIELD =
sig
  include RING
  val inv : t -> t          (* multiplicative inverse; may raise on zero *)
end

signature POLY =
sig
  (* Poly(R) is itself a ring, so these are exactly the RING members... *)
  type t
  val zero  : t
  val one   : t
  val add   : t * t -> t
  val neg   : t -> t
  val mul   : t * t -> t
  val equal : t * t -> bool
  val toString : t -> string

  (* ...plus polynomial-specific operations. *)
  type coeff

  val sub      : t * t -> t
  val scale    : coeff * t -> t          (* multiply by a scalar coefficient *)

  (* Build x^n, the monomial with coefficient `c`. `monomial (c, 0)` is the
     constant polynomial `c`. *)
  val monomial : coeff * int -> t

  (* The polynomial `x`. *)
  val x : t

  (* Construct from a low-to-high coefficient list: `fromList [a0, a1, a2]`
     is a0 + a1 x + a2 x^2. *)
  val fromList : coeff list -> t
  (* The normalized low-to-high coefficient list (no trailing zeros). *)
  val toList   : t -> coeff list

  (* `degree zero` is ~1 by convention; otherwise the highest power present. *)
  val degree   : t -> int
  (* The coefficient of x^n (zero if out of range). *)
  val coeff    : t * int -> coeff
  (* Leading coefficient (zero for the zero polynomial). *)
  val leading  : t -> coeff

  (* Horner evaluation at a point. *)
  val eval     : t -> coeff -> coeff

  (* Formal derivative. *)
  val derivative : t -> t
end

(* Euclidean-domain operations available when coefficients form a field. *)
signature POLY_FIELD =
sig
  include POLY

  (* `divMod (a, b)` returns (q, r) with a = q*b + r and degree r < degree b.
     Raises Div if b is the zero polynomial. *)
  val divMod : t * t -> t * t
  val div    : t * t -> t
  val rem    : t * t -> t

  (* Monic GCD of two polynomials (zero for gcd(0,0)). *)
  val gcd    : t * t -> t

  (* Divide through by the leading coefficient to make the polynomial monic.
     `monic zero` is `zero`. *)
  val monic  : t -> t

  (* ---- Interpolation ----

     Given a list of (x_i, y_i) sample points with *distinct* abscissae, build
     the unique interpolating polynomial of degree < (number of points) passing
     through them.  `lagrange` uses the Lagrange basis; `newton` uses Newton's
     divided differences; both return the same polynomial.  The empty list maps
     to `zero`.  Raises `Div` if two points share an x-coordinate.

     Because coefficients form a field, interpolation is *exact* (e.g. over the
     rationals it reproduces each sample point with no rounding). *)
  val lagrange : (coeff * coeff) list -> t
  val newton   : (coeff * coeff) list -> t

  (* ---- Root refinement (Newton's method) ----

     `findRoot (p, seed, iters)` performs `iters` Newton steps,
     x' = x - p(x)/p'(x), starting from `seed`, and returns the refined
     estimate.  It works over any field and stays exact over the rationals; if
     the derivative evaluates to the field's zero the iteration stops early and
     returns the current estimate.

     A general real-root *solver* (`roots`/`realRoots`) is intentionally not
     provided: the library is parameterized over an arbitrary field with no
     ordering or notion of "real", so a tolerance-based solver has no meaning at
     this level.  `findRoot` is the design-compatible primitive; callers over an
     ordered/`real` field can iterate it from chosen seeds. *)
  val findRoot : t * coeff * int -> coeff
end
