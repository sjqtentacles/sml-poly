# sml-poly

[![CI](https://github.com/sjqtentacles/sml-poly/actions/workflows/ci.yml/badge.svg)](https://github.com/sjqtentacles/sml-poly/actions/workflows/ci.yml)

Univariate polynomials over an arbitrary commutative ring, for Standard ML.

`sml-poly` is built around a single functor:

```sml
functor Poly (R : RING) :> POLY where type coeff = R.t
```

Given any structure matching `RING`, `Poly(R)` gives you polynomials with
coefficients drawn from `R` -- with `add`, `sub`, `mul`, `scale`, `eval`
(Horner), `derivative`, `degree`, and normalization built in.

Because `Poly(R)` *itself* matches `RING`, the construction nests:
`Poly(Poly(R))` gives bivariate polynomials, `Poly(Poly(Poly(R)))` trivariate,
and so on -- all from the same code.

When coefficients form a `FIELD`, the `PolyField` functor adds Euclidean
division (`divMod`/`div`/`rem`), the monic polynomial `gcd`, and `monic`.

## Portability

Pure Standard ML using only the Basis library. Verified on:

- **MLton**
- **Poly/ML**

The sources are shared via an [ML Basis](http://mlton.org/MLBasis) (`.mlb`)
file. MLton consumes it natively; for Poly/ML the test target simply `use`s
the sources in order.

## Building and testing

```sh
make test        # build + run the suite under MLton (default)
make test-poly   # run the suite under Poly/ML
make all-tests   # run under both
make clean
```

## Installing with smlpkg

`sml-poly` follows the conventions of the
[`smlpkg`](https://github.com/diku-dk/smlpkg) package manager. There is no
registry or account to sign up for -- packages are referenced directly by
their git URL. In your own project's directory:

```sh
smlpkg add github.com/sjqtentacles/sml-poly
smlpkg sync
```

This downloads the library into `lib/github.com/sjqtentacles/sml-poly/`.
Reference it from your own `.mlb` with a relative path to `poly.mlb`:

```
lib/github.com/sjqtentacles/sml-poly/poly.mlb
```

For Poly/ML, `use` the sources in order:

```sml
use "lib/github.com/sjqtentacles/sml-poly/poly.sig";
use "lib/github.com/sjqtentacles/sml-poly/poly.sml";
```

## Usage

Two instances are provided out of the box: `IntPoly` (polynomials over the
integers) and `RatPoly` (polynomials over the rationals, a `POLY_FIELD`).

```sml
structure P = IntPoly

(* (x + 1)(x - 1) = x^2 - 1 *)
val xp1 = P.add (P.x, P.one)
val xm1 = P.sub (P.x, P.one)
val sq  = P.mul (xp1, xm1)

val () = print (P.toString sq ^ "\n")          (* ~1 + 1*x^2     *)
val () = print (Int.toString (P.eval sq 3) ^ "\n")   (* 8           *)
val () = print (P.toString (P.derivative sq) ^ "\n") (* 2*x         *)
```

Defining your own coefficient ring is just a matter of matching `RING`:

```sml
structure Mod7 : RING =
struct
  type t = int
  val zero = 0 and one = 1
  fun add (a, b) = (a + b) mod 7
  fun neg a = (7 - a) mod 7
  fun mul (a, b) = (a * b) mod 7
  fun equal (a, b) = a = b
  val toString = Int.toString
end
structure Mod7Poly = Poly (Mod7)
```

Bivariate polynomials come for free by nesting:

```sml
structure BiPoly = Poly (IntPoly)   (* coefficients are themselves polynomials *)
```

## Project layout

```
sml.pkg                                         smlpkg manifest
Makefile                                        build + test
lib/github.com/sjqtentacles/sml-poly/
  poly.sig                                      RING/FIELD/POLY/POLY_FIELD signatures
  poly.sml                                      Poly + PolyField functors, instances
  poly.mlb                                      MLB for consumers
test/
  test.mlb                                      test basis (MLton)
  test.sml                                      assertion suite
.github/workflows/ci.yml                        CI (MLton + Poly/ML)
```

## License

MIT. See [LICENSE](LICENSE).
