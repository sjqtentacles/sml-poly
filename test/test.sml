(* Dependency-free test runner for the polynomial library.
 * Prints one line per assertion and exits non-zero if any assertion fails. *)

val passed = ref 0
val failed = ref 0

fun check (name : string) (cond : bool) : unit =
    if cond
    then (passed := !passed + 1; print ("ok   - " ^ name ^ "\n"))
    else (failed := !failed + 1; print ("FAIL - " ^ name ^ "\n"))

fun raisesDiv (thunk : unit -> 'a) : bool =
    (ignore (thunk ()); false) handle General.Div => true | _ => false

(* ---- Integer polynomials ----------------------------------------------- *)
structure P = IntPoly

(* helpers to build small polynomials from low-to-high int lists *)
fun p cs = P.fromList cs

fun runInt () =
  let
    val () = check "zero degree is ~1" (P.degree P.zero = ~1)
    val () = check "one is constant 1" (P.toList P.one = [1])
    val () = check "x is [0,1]" (P.toList P.x = [0, 1])
    val () = check "fromList drops trailing zeros"
                   (P.toList (p [1, 2, 0, 0]) = [1, 2])
    val () = check "fromList all zeros is zero"
                   (P.equal (p [0, 0, 0], P.zero))

    val () = check "monomial 3 x^2" (P.toList (P.monomial (3, 2)) = [0, 0, 3])
    val () = check "monomial coeff 0 is zero"
                   (P.equal (P.monomial (0, 5), P.zero))

    (* (1 + 2x) + (3 + 4x + 5x^2) = 4 + 6x + 5x^2 *)
    val () = check "add"
                   (P.toList (P.add (p [1, 2], p [3, 4, 5])) = [4, 6, 5])
    (* cancellation normalizes away the high term *)
    val () = check "add cancels high terms"
                   (P.toList (P.add (p [1, 2, 3], p [0, 0, ~3])) = [1, 2])
    val () = check "sub"
                   (P.toList (P.sub (p [5, 5, 5], p [1, 2, 3])) = [4, 3, 2])
    val () = check "neg" (P.toList (P.neg (p [1, ~2, 3])) = [~1, 2, ~3])
    val () = check "scale" (P.toList (P.scale (2, p [1, 2, 3])) = [2, 4, 6])
    val () = check "scale by zero is zero"
                   (P.equal (P.scale (0, p [1, 2, 3]), P.zero))

    (* (x + 1)(x - 1) = x^2 - 1 *)
    val xp1 = P.add (P.x, P.one)
    val xm1 = P.sub (P.x, P.one)
    val () = check "(x+1)(x-1) = x^2 - 1"
                   (P.toList (P.mul (xp1, xm1)) = [~1, 0, 1])
    (* (x+1)^2 = x^2 + 2x + 1 *)
    val () = check "(x+1)^2"
                   (P.toList (P.mul (xp1, xp1)) = [1, 2, 1])
    val () = check "mul by zero is zero"
                   (P.equal (P.mul (xp1, P.zero), P.zero))
    val () = check "mul by one is identity"
                   (P.equal (P.mul (xp1, P.one), xp1))

    val () = check "degree of x^2-1 is 2"
                   (P.degree (P.mul (xp1, xm1)) = 2)
    val () = check "leading coeff" (P.leading (p [1, 2, 7]) = 7)
    val () = check "coeff in range" (P.coeff (p [1, 2, 3], 1) = 2)
    val () = check "coeff out of range is zero" (P.coeff (p [1, 2, 3], 9) = 0)
    val () = check "coeff negative index is zero" (P.coeff (p [1, 2, 3], ~1) = 0)

    (* eval: (x^2 - 1) at 3 = 8 ; at 1 = 0 *)
    val sq = P.mul (xp1, xm1)
    val () = check "eval (x^2-1) at 3 = 8" (P.eval sq 3 = 8)
    val () = check "eval (x^2-1) at 1 = 0" (P.eval sq 1 = 0)
    val () = check "eval constant" (P.eval (p [7]) 100 = 7)

    (* derivative: d/dx (1 + 2x + 3x^2) = 2 + 6x *)
    val () = check "derivative"
                   (P.toList (P.derivative (p [1, 2, 3])) = [2, 6])
    val () = check "derivative of constant is zero"
                   (P.equal (P.derivative (p [42]), P.zero))
    val () = check "derivative of x^3 is 3x^2"
                   (P.toList (P.derivative (P.monomial (1, 3))) = [0, 0, 3])

    val () = check "equal ignores trailing zeros"
                   (P.equal (p [1, 2], p [1, 2, 0, 0]))
    val () = check "not equal" (not (P.equal (p [1, 2], p [1, 3])))

    val () = check "toString x^2-1" (P.toString sq = "~1 + 1*x^2")
  in () end

(* ---- Bivariate: Poly(Poly(IntRing)) ------------------------------------ *)
structure PP = Poly (IntPoly)

fun runBivariate () =
  let
    (* Treat outer variable as y, inner as x.
       p(x,y) = (x) + (1)*y  i.e. coefficients in x: [x, 1]
       Represent x as IntPoly.x, the constant 1 as IntPoly.one. *)
    val cX  = IntPoly.x                 (* the polynomial "x" as a coeff *)
    val c1  = IntPoly.one
    (* poly in y:  x + 1*y  *)
    val q = PP.fromList [cX, c1]
    val () = check "bivariate degree in y is 1" (PP.degree q = 1)
    (* (x + y)^2 = x^2 + 2xy + y^2 ; coefficients in y:
       [x^2, 2x, 1] *)
    val sq = PP.mul (q, q)
    val coeffs = PP.toList sq
    val () = check "bivariate (x+y)^2 has 3 y-terms" (length coeffs = 3)
    val () = check "  y^0 coeff is x^2"
                   (IntPoly.toList (List.nth (coeffs, 0)) = [0, 0, 1])
    val () = check "  y^1 coeff is 2x"
                   (IntPoly.toList (List.nth (coeffs, 1)) = [0, 2])
    val () = check "  y^2 coeff is 1"
                   (IntPoly.toList (List.nth (coeffs, 2)) = [1])
    (* evaluate the outer poly at y = (the constant poly 3): result in x is
       x^2 + 6x + 9 = (x+3)^2 *)
    val atY3 = PP.eval sq (IntPoly.fromList [3])
    val () = check "  (x+y)^2 at y=3 is (x+3)^2"
                   (IntPoly.toList atY3 = [9, 6, 1])
  in () end

(* ---- Field polynomials: divMod / gcd ----------------------------------- *)
structure FP = RatPoly

fun rat n = (IntInf.fromInt n, IntInf.fromInt 1) : RatField.t
fun fp ns = FP.fromList (map rat ns)

fun runField () =
  let
    (* divide x^2 - 1 by x - 1 -> quotient x + 1, remainder 0 *)
    val a = fp [~1, 0, 1]
    val b = fp [~1, 1]
    val (q, r) = FP.divMod (a, b)
    val () = check "divMod quotient (x^2-1)/(x-1) = x+1"
                   (FP.equal (q, fp [1, 1]))
    val () = check "divMod remainder is zero" (FP.equal (r, FP.zero))

    (* divide with nonzero remainder: (x^2 + 1) / (x - 1) = x + 1 rem 2 *)
    val a2 = fp [1, 0, 1]
    val (q2, r2) = FP.divMod (a2, b)
    val () = check "divMod quotient (x^2+1)/(x-1) = x+1"
                   (FP.equal (q2, fp [1, 1]))
    val () = check "divMod remainder = 2" (FP.equal (r2, fp [2]))

    (* a = q*b + r reconstruction *)
    val () = check "a = q*b + r"
                   (FP.equal (FP.add (FP.mul (q2, b), r2), a2))

    val () = check "divMod by zero raises Div"
                   (raisesDiv (fn () => FP.divMod (a, FP.zero)))

    (* gcd((x^2-1),(x-1)) = x - 1 made monic = x - 1 *)
    val g = FP.gcd (a, b)
    val () = check "gcd(x^2-1, x-1) = x-1" (FP.equal (g, fp [~1, 1]))

    (* gcd of (x^2-1) and (x^2 - 2x + 1)=(x-1)^2 is (x-1) (monic) *)
    val c = FP.mul (fp [~1, 1], fp [~1, 1])      (* (x-1)^2 *)
    val g2 = FP.gcd (a, c)
    val () = check "gcd(x^2-1,(x-1)^2) = x-1" (FP.equal (g2, fp [~1, 1]))

    (* monic: 2x + 4 -> x + 2 *)
    val () = check "monic scales leading to 1"
                   (FP.equal (FP.monic (fp [4, 2]), fp [2, 1]))
  in () end

(* ---- Interpolation and Newton root refinement (exact, over rationals) --- *)

(* a rational literal a/b, normalized through the field constructor *)
fun ratFrac (a, b) = RatField.mul ((IntInf.fromInt a, IntInf.fromInt 1),
                                   RatField.inv (IntInf.fromInt b, IntInf.fromInt 1))

(* |a| < b, assuming b > 0 (RatField keeps the denominator positive) *)
fun ratLessAbs (a : RatField.t, b : RatField.t) =
    let
      val (an, ad) = a
      val (bn, bd) = b
      val an = if an < 0 then ~an else an
    in
      an * bd < bn * ad
    end

fun runInterp () =
  let
    val eps = ratFrac (1, 1000000)   (* 1e-6 *)

    (* recover a known degree-2 polynomial Q = 2 + 3x - x^2 from samples *)
    val q = fp [2, 3, ~1]
    fun sampleAt xi = (rat xi, FP.eval q (rat xi))
    val pts3 = [sampleAt 0, sampleAt 1, sampleAt 2]
    val lag = FP.lagrange pts3
    val newt = FP.newton pts3
    val () = check "lagrange recovers 2+3x-x^2" (FP.equal (lag, q))
    val () = check "newton recovers 2+3x-x^2"   (FP.equal (newt, q))
    val () = check "lagrange = newton (deg 2)"  (FP.equal (lag, newt))

    (* an extra collinear-degree sample must not change the result *)
    val pts4 = pts3 @ [sampleAt 5]
    val () = check "lagrange stable with extra point" (FP.equal (FP.lagrange pts4, q))
    val () = check "newton stable with extra point"   (FP.equal (FP.newton pts4, q))

    (* arbitrary points: interpolant must reproduce each sample exactly *)
    val pts = [(rat 0, rat 1), (rat 1, rat 3), (rat 2, rat 2), (rat 3, rat 5)]
    val li = FP.lagrange pts
    val ni = FP.newton pts
    val () = check "lagrange = newton (arbitrary)" (FP.equal (li, ni))
    val () = List.app
               (fn (xi, yi) =>
                  check ("interp reproduces sample x=" ^ RatField.toString xi)
                        (RatField.equal (FP.eval li xi, yi)))
               pts
    val () = check "interpolant degree < #points" (FP.degree li <= 3)

    (* Newton root refinement on (x-1)(x-2)(x-3) *)
    val cubic = FP.mul (FP.mul (fp [~1, 1], fp [~2, 1]), fp [~3, 1])
    fun nearRoot (r, seedOffsetNum) =
        let
          val seed = RatField.add (rat r, ratFrac (seedOffsetNum, 10))
          val approx = FP.findRoot (cubic, seed, 8)
          val distToR = RatField.add (approx, RatField.neg (rat r))
          val pval = FP.eval cubic approx
        in
          check ("findRoot -> " ^ Int.toString r ^ " (within eps)")
                (ratLessAbs (distToR, eps));
          check ("p(root~=" ^ Int.toString r ^ ") ~= 0")
                (ratLessAbs (pval, eps))
        end
    val () = nearRoot (1, ~1)   (* seed 0.9 *)
    val () = nearRoot (2, 1)    (* seed 2.1 *)
    val () = nearRoot (3, 1)    (* seed 3.1 *)
    (* seeding exactly at a root is a fixed point *)
    val () = check "findRoot fixed at exact root"
                   (RatField.equal (FP.findRoot (cubic, rat 2, 5), rat 2))
    val () = check "findRoot constant-deriv guard (deriv 0 stops)"
                   (RatField.equal (FP.findRoot (fp [5], rat 3, 4), rat 3))
  in () end

fun run () =
  let
    val () = runInt ()
    val () = runBivariate ()
    val () = runField ()
    val () = runInterp ()
  in
    print ("\n" ^ Int.toString (!passed) ^ " passed, "
           ^ Int.toString (!failed) ^ " failed\n");
    OS.Process.exit (if !failed = 0 then OS.Process.success else OS.Process.failure)
  end

val () = run ()
