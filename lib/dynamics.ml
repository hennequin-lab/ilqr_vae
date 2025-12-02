open Base
open Misc

module type T = Dynamics_intf.T

module Integrate (D : T) = struct
  open AD.Maths

  let integrate ~prms =
    let dyn_k = D.dyn ~theta:prms in
    fun ~n ~u ->
      (* assume u is n_samples x n_steps x m *)
      assert (Array.length (AD.shape u) = 3);
      let u = transpose ~axis:[| 1; 0; 2 |] u in
      (* now u is T x K x M *)
      let n_steps = AD.(shape u).(0) in
      let n_samples = AD.(shape u).(1) in
      let x0 = AD.Mat.zeros n_samples n in
      let us =
        let u = reshape u [| n_steps; -1 |] in
        split ~axis:0 (Array.init n_steps ~f:(fun _ -> 1)) u
        |> Array.map ~f:(fun v -> reshape v [| n_samples; -1 |])
        |> Array.to_list
      in
      let rec dyn k x xs us =
        match us with
        | [] -> List.rev xs
        | u :: unexts ->
          let new_x = dyn_k ~k ~x ~u in
          dyn Int.(k + 1) new_x (new_x :: xs) unexts
      in
      dyn 0 x0 [] us
      |> Array.of_list
      |> Array.map ~f:(fun v -> reshape v [| 1; n_samples; n |])
      |> concatenate ~axis:0 (* T x K x N *)
      |> transpose ~axis:[| 1; 0; 2 |]

  (* result KxTxN *)
end

let b_rescaled b =
  let open AD.Maths in
  Option.map b ~f:(function b -> b / sqrt (sum ~axis:0 (sqr b)))


module Linear = struct
  module P = Dynamics_intf.Linear_P
  open P
  open AD.Maths

  let requires_linesearch = false

  (* alpha is the spectral abscissa of the equivalent continuous-time system
     beta is the spectral radius of the random S *)
  let init ~dt_over_tau ~alpha ~beta ~n ~m () =
    (* exp (dt_over_tau * (W-I))
       where W = alpha*I + S *)
    let d =
      let tmp = Float.(exp (-2. * dt_over_tau * (1.0 - alpha))) in
      AA.init [| 1; n |] (fun _ -> Float.(tmp / (1. - tmp)))
    in
    let u = AD.Mat.eye n in
    let q =
      let s =
        AA.(Float.(beta * dt_over_tau / sqrt (2. * of_int n)) $* gaussian [| n; n |])
      in
      AA.Linalg.expm AA.(s - transpose s)
    in
    let b = if n = m then None else Some (Prms.create (AD.Mat.gaussian m n)) in
    { d = Prms.create ~above:(F 1E-5) (AD.pack_arr d)
    ; u = Prms.create u
    ; q = Prms.create (AD.pack_arr q)
    ; b
    }


  let unpack_a ~prms =
    let q =
      let q, r = AD.Linalg.qr prms.q in
      q * signum (diag r)
    in
    let u =
      let q, r = AD.Linalg.qr prms.u in
      q * signum (diag r)
    in
    let dp1_sqrt_inv = F 1. / sqrt (F 1. + prms.d) in
    let d_sqrt = sqrt prms.d in
    u * d_sqrt *@ (q * dp1_sqrt_inv) *@ transpose u


  let extract_b ~theta ~n =
    match b_rescaled theta.b with
    | None -> AD.Mat.(eye n)
    | Some b -> b


  let dyn ~theta =
    let a = unpack_a ~prms:theta in
    let n = AD.Mat.row_num a in
    let b = extract_b ~theta ~n in
    fun ~k ~x ~u -> if k = 0 then (x *@ a) + u else (x *@ a) + (u *@ b)


  let dyn_x =
    let dyn_x ~theta =
      let a = unpack_a ~prms:theta in
      fun ~k:_ ~x:_ ~u:_ -> a
    in
    Some dyn_x


  let dyn_u =
    let dyn_u ~theta =
      let n = AD.Mat.row_num theta.q in
      let id_n = AD.Mat.eye n in
      let b = extract_b ~theta ~n in
      fun ~k ~x:_ ~u:_ -> if k = 0 then id_n else b
    in
    Some dyn_u
end

module Nonlinear (X : sig
    val phi : [ `linear | `nonlinear of (AD.t -> AD.t) * (AD.t -> AD.t) ]
  end) =
struct
  module P = Dynamics_intf.Nonlinear_P
  open P
  open X
  open AD.Maths

  let phi, d_phi, requires_linesearch =
    match phi with
    | `linear -> (fun x -> x), (fun x -> AD.Arr.(ones (shape x))), false
    | `nonlinear (f, df) -> f, df, true


  let init ?(radius = 0.1) ~n ~m () =
    let sigma = Float.(radius / sqrt (of_int n)) in
    { a = Prms.create (AD.Mat.gaussian ~sigma n n)
    ; bias = Prms.create (AD.Mat.zeros 1 n)
    ; b = Some (Prms.create (AD.Mat.gaussian ~sigma:Float.(1. / sqrt (of_int m)) m n))
    }


  let u_eff ~prms =
    match b_rescaled prms.b with
    | None -> fun u -> u
    | Some b -> fun u -> u *@ b


  let dyn ~theta =
    let u_eff = u_eff ~prms:theta in
    let passive x = (phi x *@ theta.a) + theta.bias in
    fun ~k ~x ~u -> if k = 0 then passive x + u else passive x + u_eff u


  let dyn_x =
    let dyn_x ~theta = fun ~k:_ ~x ~u:_ -> transpose (d_phi x) * theta.a in
    Some dyn_x


  let dyn_u =
    let dyn_u ~theta =
      let n = AD.Mat.row_num theta.a in
      let id_n = AD.Mat.eye n in
      let b =
        match b_rescaled theta.b with
        | None -> id_n
        | Some b -> b
      in
      fun ~k ~x:_ ~u:_ -> if k = 0 then id_n else b
    in
    Some dyn_u
end

module Linear_nonlinear (X : sig
    val phi : (AD.t -> AD.t) * (AD.t -> AD.t)
    val dt : float
    val tau : float
  end) =
struct
  module P = Dynamics_intf.Linear_nonlinear_P
  open P
  open X
  open AD.Maths

  let phi, d_phi = phi
  let requires_linesearch = true

  let init ?(radius = 0.1) ~n ~nh ~m () =
    let sigma1 = Float.(1. / of_int n) in
    let sigma2 = Float.(radius / sqrt (of_int nh)) in
    { a0 = Prms.create (F Float.(1. - (dt / tau)) * AD.Mat.eye n)
    ; a1 = Prms.create (AD.Mat.gaussian ~sigma:sigma1 n nh)
    ; a2 = Prms.create (AD.Mat.gaussian ~sigma:sigma2 nh n)
    ; bias1 = Prms.create (AD.Mat.gaussian 1 nh)
    ; bias2 = Prms.create (AD.Mat.zeros 1 n)
    ; b = Some (Prms.create (AD.Mat.gaussian ~sigma:Float.(1. / sqrt (of_int m)) m n))
    }


  let u_eff ~prms =
    match b_rescaled prms.b with
    | None -> fun u -> u
    | Some b -> fun u -> AD.Maths.(u *@ b)


  let dyn ~theta =
    let u_eff = u_eff ~prms:theta in
    let passive x =
      (x *@ theta.a0)
      + (F Float.(dt / tau)
         * ((phi ((x *@ theta.a1) + theta.bias1) *@ theta.a2) + theta.bias2))
    in
    fun ~k ~x ~u ->
      if k = 0 then passive x + u else passive x + (F Float.(dt / tau) * u_eff u)


  let dyn_x =
    let dyn_x ~theta =
      fun ~k:_ ~x ~u:_ ->
      let h = (x *@ theta.a1) + theta.bias1 in
      (F Float.(dt / tau) * (theta.a1 * d_phi h *@ theta.a2)) + theta.a0
    in
    Some dyn_x


  let dyn_u =
    let dyn_u ~theta =
      let n = AD.Mat.row_num theta.a1 in
      let id_n = AD.Mat.eye n in
      let b =
        match b_rescaled theta.b with
        | None -> id_n
        | Some b -> b
      in
      let b = F Float.(dt / tau) * b in
      fun ~k ~x:_ ~u:_ -> if k = 0 then id_n else b
    in
    Some dyn_u
end
