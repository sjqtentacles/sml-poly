(* poly.sml

   The Poly and PolyField functors plus a few ready-made instances.

   Representation: a polynomial is a dense, low-to-high coefficient list that
   is always kept normalized (no trailing coefficients equal to R.zero). The
   zero polynomial is the empty list. This invariant is established by `norm`
   and maintained by every operation that can produce trailing zeros. *)

functor Poly (R : RING) :> POLY where type coeff = R.t =
struct
  type coeff = R.t

  (* low-to-high; invariant: no trailing R.zero (so [] is the only zero). *)
  type t = R.t list

  fun isZero c = R.equal (c, R.zero)

  (* Drop trailing zeros. Implemented by recursing to the tail and dropping a
     leading zero only when everything after it is already empty. *)
  fun norm xs =
      let
        fun strip [] = []
          | strip (y :: ys) =
              (case strip ys of
                   [] => if isZero y then [] else [y]
                 | r  => y :: r)
      in
        strip xs
      end

  val zero = []
  val one  = norm [R.one]

  fun degree p = length (norm p) - 1

  fun coeff (p, n) =
      if n < 0 then R.zero
      else (case List.drop (p, n) of c :: _ => c | [] => R.zero)
           handle Subscript => R.zero

  fun leading p =
      case norm p of
          [] => R.zero
        | xs => List.last xs

  fun add (a, b) =
      let
        fun go ([], ys) = ys
          | go (xs, []) = xs
          | go (x :: xs, y :: ys) = R.add (x, y) :: go (xs, ys)
      in
        norm (go (a, b))
      end

  fun neg a = map R.neg a

  fun sub (a, b) = add (a, neg b)

  fun scale (c, a) =
      if isZero c then zero else norm (map (fn x => R.mul (c, x)) a)

  (* multiply by x *)
  fun shift a = case norm a of [] => [] | xs => R.zero :: xs

  fun mul (a, b) =
      let
        fun go [] = zero
          | go (x :: xs) = add (scale (x, b), shift (go xs))
      in
        norm (go (norm a))
      end

  fun eqList ([], []) = true
    | eqList (x :: xs, y :: ys) = R.equal (x, y) andalso eqList (xs, ys)
    | eqList _ = false

  fun equal (a, b) = eqList (norm a, norm b)

  fun monomial (c, n) =
      if n < 0 then raise Domain
      else if isZero c then zero
      else norm (List.tabulate (n + 1, fn i => if i = n then c else R.zero))

  val x = monomial (R.one, 1)

  fun fromList cs = norm cs
  fun toList p = norm p

  (* Horner: a0 + x*(a1 + x*(a2 + ...)) *)
  fun eval p v =
      List.foldr (fn (c, acc) => R.add (c, R.mul (v, acc))) R.zero p

  fun derivative p =
      case p of
          [] => []
        | (_ :: rest) =>
            let
              (* term a_k x^k -> k*a_k x^(k-1); k counts from 1 over `rest`. *)
              fun timesN (c, n) =
                  let
                    fun loop (0, acc) = acc
                      | loop (k, acc) = loop (k - 1, R.add (acc, c))
                  in
                    loop (n, R.zero)
                  end
              fun go (_, []) = []
                | go (k, c :: cs) = timesN (c, k) :: go (k + 1, cs)
            in
              norm (go (1, rest))
            end

  fun toString p =
      case norm p of
          [] => R.toString R.zero
        | xs =>
            let
              fun term (c, 0) = R.toString c
                | term (c, 1) = R.toString c ^ "*x"
                | term (c, n) = R.toString c ^ "*x^" ^ Int.toString n
              fun go (_, []) = []
                | go (i, c :: cs) =
                    if isZero c then go (i + 1, cs)
                    else term (c, i) :: go (i + 1, cs)
            in
              case go (0, xs) of
                  [] => R.toString R.zero
                | parts => String.concatWith " + " parts
            end
end

functor PolyField (F : FIELD) :> POLY_FIELD where type coeff = F.t =
struct
  structure P = Poly (F)
  open P

  fun monic p =
      case toList p of
          [] => zero
        | _  =>
            let val lc = leading p
            in if F.equal (lc, F.one) then p
               else scale (F.inv lc, p)
            end

  (* Long division over a field. *)
  fun divMod (a, b) =
      let
        val bl = toList b
      in
        if null bl then raise Div
        else
          let
            val db = degree b
            val lcb = leading b
            val invB = F.inv lcb
            fun loop (r, q) =
                if equal (r, zero) orelse degree r < db
                then (q, r)
                else
                  let
                    val dr = degree r
                    val factorCoeff = F.mul (leading r, invB)
                    val factor = monomial (factorCoeff, dr - db)
                    val r' = sub (r, mul (factor, b))
                    val q' = add (q, factor)
                  in
                    loop (r', q')
                  end
          in
            loop (a, zero)
          end
      end

  fun op div (a, b) = #1 (divMod (a, b))
  fun rem (a, b) = #2 (divMod (a, b))

  fun gcd (a, b) =
      let
        fun loop (a, b) =
            if equal (b, zero) then a
            else loop (b, rem (a, b))
        val g = loop (a, b)
      in
        monic g
      end
end

(* ---- Ready-made instances ---------------------------------------------- *)

structure IntRing : RING =
struct
  type t = int
  val zero = 0
  val one  = 1
  fun add (a, b) = a + b
  fun neg a = ~a
  fun mul (a, b) = a * b
  fun equal (a, b) = a = b
  fun toString a = Int.toString a
end

(* Rationals over IntInf, used as the field instance for divMod/gcd tests. *)
structure RatField : FIELD =
struct
  (* normalized: denominator > 0, gcd(|num|,den) = 1, zero is (0,1) *)
  type t = IntInf.int * IntInf.int

  val izero : IntInf.int = 0
  val ione  : IntInf.int = 1

  fun igcd (a : IntInf.int, b : IntInf.int) : IntInf.int =
      if b = izero then a else igcd (b, a mod b)

  fun normalize (n : IntInf.int, d : IntInf.int) : t =
      if d = izero then raise Div
      else
        let
          val s = if d < izero then ~ione else ione
          val n = n * s
          val d = d * s
          val g = igcd (if n < izero then ~n else n, d)
          val g = if g = izero then ione else g
        in
          (n div g, d div g)
        end

  val zero = (izero, ione) : t
  val one  = (ione, ione) : t
  fun add ((a, b) : t, (c, d) : t) = normalize (a * d + c * b, b * d)
  fun neg ((a, b) : t) : t = (~a, b)
  fun mul ((a, b) : t, (c, d) : t) = normalize (a * c, b * d)
  fun equal ((a, b) : t, (c, d) : t) = a * d = c * b
  fun inv ((a, b) : t) = if a = izero then raise Div else normalize (b, a)
  fun toString ((a, b) : t) =
      if b = ione then IntInf.toString a
      else IntInf.toString a ^ "/" ^ IntInf.toString b
end

structure IntPoly  = Poly (IntRing)
structure RatPoly  = PolyField (RatField)
