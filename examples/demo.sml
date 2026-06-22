(* demo.sml - univariate polynomials over a ring (integers) and over a field
   (rationals). Everything is exact: integer coefficients print as integers,
   rational coefficients as "num/den" via the library's own toString, so the
   output is identical on every run and on both compilers. No reals involved. *)

(* ---- Integer polynomials (Poly over IntRing) ---- *)
structure P = IntPoly

val f = P.fromList [1, 2, 3]            (* 1 + 2x + 3x^2 *)
val () = print "Integer polynomials:\n"
val () = print ("  f(x)        = " ^ P.toString f ^ "\n")
val () = print ("  f(2)        = " ^ Int.toString (P.eval f 2) ^ "\n")
val () = print ("  f'(x)       = " ^ P.toString (P.derivative f) ^ "\n")
val () = print ("  (x+1)(x-1)  = "
                ^ P.toString (P.mul (P.add (P.x, P.one), P.sub (P.x, P.one)))
                ^ "\n")

(* ---- Rational polynomials (PolyField over RatField) ---- *)
structure F = RatPoly

fun rat n = (IntInf.fromInt n, IntInf.fromInt 1) : RatField.t
fun fp ns = F.fromList (map rat ns)

val a = fp [~1, 0, 1]                   (* x^2 - 1 *)
val b = fp [~1, 1]                      (* x - 1   *)
val (q, r) = F.divMod (a, b)
val () = print "\nRational polynomials:\n"
val () = print ("  (x^2-1) / (x-1) = " ^ F.toString q
                ^ "  rem " ^ F.toString r ^ "\n")
val () = print ("  gcd(x^2-1, (x-1)^2) = "
                ^ F.toString (F.gcd (a, F.mul (b, b))) ^ "\n")

(* recover q(x) = 2 + 3x - x^2 from three exact sample points *)
val target = fp [2, 3, ~1]
val pts = List.map (fn xi => (rat xi, F.eval target (rat xi))) [0, 1, 2]
val () = print ("  lagrange of {(0,?),(1,?),(2,?)} = "
                ^ F.toString (F.lagrange pts) ^ "\n")
