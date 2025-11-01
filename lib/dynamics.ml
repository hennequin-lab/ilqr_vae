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


module Linear (X : sig
    val n_beg : int Option.t
  end) =
struct
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


  let generate_bs ~n ~m =
    Option.map X.n_beg ~f:(fun nb ->
      let nr = Int.(n / nb) in
      assert (nr = m);
      ( nb
      , Array.init nb ~f:(fun i ->
          let inr = Int.(i * nr) in
          let rnr = Int.(n - ((i + 1) * nr)) in
          transpose
            (concatenate
               ~axis:0
               [| AD.Mat.zeros inr m; AD.Mat.eye nr; AD.Mat.zeros rnr m |])) ))


  let extract_b ~theta ~n =
    match b_rescaled theta.b with
    | None -> AD.Mat.(eye n)
    | Some b -> b


  let dyn ~theta =
    let a = unpack_a ~prms:theta in
    let n = AD.Mat.row_num a in
    let b = extract_b ~theta ~n in
    let m = AD.Mat.row_num b in
    let beg_bs = generate_bs ~n ~m in
    let default x u = (x *@ a) + (u *@ b) in
    fun ~k ~x ~u ->
      match beg_bs with
      | None -> default x u
      | Some (i, beg_b) -> if k < i then x + (u *@ beg_b.(k)) else default x u


  let dyn_x =
    (* Marine to check this *)
    let dyn_x ~theta =
      let a = unpack_a ~prms:theta in
      let n = AD.Mat.row_num a in
      let id_n = AD.Mat.eye n in
      fun ~k ~x:_ ~u:_ ->
        match X.n_beg with
        | None -> a
        | Some i -> if k < i then id_n else a
    in
    Some dyn_x


  let dyn_u =
    let dyn_u ~theta =
      let n = AD.Mat.row_num theta.q in
      let b = extract_b ~theta ~n in
      let m = AD.Mat.row_num b in
      let beg_bs = generate_bs ~n ~m in
      fun ~k ~x:_ ~u:_ ->
        match beg_bs with
        | None -> b
        | Some (i, beg_b) -> if k < i then beg_b.(k) else b
    in
    Some dyn_u
end

module Nonlinear (X : sig
    val phi : [ `linear | `nonlinear of (AD.t -> AD.t) * (AD.t -> AD.t) ]
    val n_beg : int Option.t
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


  let generate_bs ~n ~m =
    Option.map X.n_beg ~f:(fun nb ->
      let nr = Int.(n / nb) in
      assert (nr = m);
      ( nb
      , Array.init nb ~f:(fun i ->
          let inr = Int.(i * nr) in
          let rnr = Int.(n - ((i + 1) * nr)) in
          AD.Maths.(
            transpose
              (concatenate
                 ~axis:0
                 [| AD.Mat.zeros inr m; AD.Mat.eye nr; AD.Mat.zeros rnr m |]))) ))


  let u_eff ~prms =
    match b_rescaled prms.b with
    | None -> fun u -> u
    | Some b -> fun u -> AD.Maths.(u *@ b)


  let dyn ~theta =
    let n = AD.Mat.row_num theta.a in
    let m =
      match theta.b with
      | None -> n
      | Some b -> AD.Mat.row_num b
    in
    let u_eff = u_eff ~prms:theta in
    let beg_bs = generate_bs ~n ~m in
    let default x u = (phi x *@ theta.a) + u_eff u + theta.bias in
    fun ~k ~x ~u ->
      match beg_bs with
      | None -> default x u
      | Some (nb, beg_b) -> if Int.(k < nb) then x + (u *@ beg_b.(k)) else default x u


  let dyn_x =
    let dyn_x ~theta =
      let n = AD.Mat.row_num theta.a in
      let id_n = AD.Mat.eye n in
      let default x = AD.Maths.(transpose (d_phi x) * theta.a) in
      fun ~k ~x ~u:_ ->
        match X.n_beg with
        | None -> default x
        | Some nb -> if Int.(k < nb) then id_n else default x
    in
    Some dyn_x


  let dyn_u =
    let dyn_u ~theta =
      let n = AD.Mat.row_num theta.a in
      let b =
        match b_rescaled theta.b with
        | None -> AD.Mat.eye n
        | Some b -> b
      in
      let m = AD.Mat.row_num b in
      let beg_bs = generate_bs ~n ~m in
      fun ~k ~x:_ ~u:_ ->
        match beg_bs with
        | None -> b
        | Some (nb, beg_b) -> if Int.(k < nb) then beg_b.(k) else b
    in
    Some dyn_u
end

module InvertedBottleneck (X : sig
    val phi : (AD.t -> AD.t) * (AD.t -> AD.t)
    val n_beg : int Option.t
  end) =
struct
  module P = Dynamics_intf.InvertedBottleneck_P
  open P
  open X
  open AD.Maths

  let phi, d_phi = phi
  let requires_linesearch = true

  let init ?(radius = 0.1) ?(decay = 0.5) ~n ~nh ~m () =
    let sigma = Float.(1. / sqrt (of_int n)) in
    let sigma_h = Float.(radius / sqrt (of_int nh)) in
    { a0 = Prms.create (F decay * AD.Mat.eye n)
    ; a1 = Prms.create (AD.Mat.gaussian ~sigma n nh)
    ; a2 = Prms.create (AD.Mat.gaussian ~sigma:sigma_h nh n)
    ; bias1 = Prms.create (AD.Mat.gaussian 1 nh)
    ; bias2 = Prms.create (AD.Mat.zeros 1 n)
    ; b = Some (Prms.create (AD.Mat.gaussian ~sigma:Float.(1. / sqrt (of_int m)) m n))
    }


  let generate_bs ~n ~m =
    Option.map X.n_beg ~f:(fun nb ->
      let nr = Int.(n / nb) in
      assert (nr = m);
      ( nb
      , Array.init nb ~f:(fun i ->
          let inr = Int.(i * nr) in
          let rnr = Int.(n - ((i + 1) * nr)) in
          AD.Maths.(
            transpose
              (concatenate
                 ~axis:0
                 [| AD.Mat.zeros inr m; AD.Mat.eye nr; AD.Mat.zeros rnr m |]))) ))


  let u_eff ~prms =
    match b_rescaled prms.b with
    | None -> fun u -> u
    | Some b -> fun u -> AD.Maths.(u *@ b)


  let dyn ~theta =
    let n = AD.Mat.row_num theta.a1 in
    let m =
      match theta.b with
      | None -> n
      | Some b -> AD.Mat.row_num b
    in
    let u_eff = u_eff ~prms:theta in
    let beg_bs = generate_bs ~n ~m in
    let default x u =
      (x *@ theta.a0)
      + (phi ((x *@ theta.a1) + theta.bias1) *@ theta.a2)
      + theta.bias2
      + u_eff u
    in
    fun ~k ~x ~u ->
      match beg_bs with
      | None -> default x u
      | Some (nb, beg_b) -> if Int.(k < nb) then x + (u *@ beg_b.(k)) else default x u


  let dyn_x =
    let dyn_x ~theta =
      let n = AD.Mat.row_num theta.a1 in
      let id_n = AD.Mat.eye n in
      let default x =
        let h = (x *@ theta.a1) + theta.bias1 in
        AD.Maths.((theta.a1 * d_phi h *@ theta.a2) + theta.a0)
      in
      fun ~k ~x ~u:_ ->
        match X.n_beg with
        | None -> default x
        | Some nb -> if Int.(k < nb) then id_n else default x
    in
    Some dyn_x


  let dyn_u =
    let dyn_u ~theta =
      let n = AD.Mat.row_num theta.a1 in
      let b =
        match b_rescaled theta.b with
        | None -> AD.Mat.eye n
        | Some b -> b
      in
      let m = AD.Mat.row_num b in
      let beg_bs = generate_bs ~n ~m in
      fun ~k ~x:_ ~u:_ ->
        match beg_bs with
        | None -> b
        | Some (nb, beg_b) -> if Int.(k < nb) then beg_b.(k) else b
    in
    Some dyn_u
end

module MGU (X : sig
    val phi : AD.t -> AD.t
    val d_phi : AD.t -> AD.t
    val sigma : AD.t -> AD.t
    val d_sigma : AD.t -> AD.t
    val n_beg : int Option.t
  end) =
struct
  module P = Dynamics_intf.MGU_P
  open P
  open X
  open AD.Maths

  let requires_linesearch = true

  let init ~n ~m () =
    (*h = size 1xN
     x = size 1xN (x = Bu)
     h = size 1xK
     f = size of h so 1xN
     Wf = NxN
     B = MxN *)
    { wf = Prms.create (AD.Mat.zeros m n)
    ; wh = Prms.create (AD.Mat.gaussian m n)
    ; bh = Prms.create (AD.Mat.zeros 1 n)
    ; bf = Prms.create (AD.Mat.zeros 1 n)
    ; uh = Prms.create (AD.Mat.zeros n n)
    ; uf = Prms.create (AD.Mat.zeros n n)
    }


  let with_wh_rescaled theta =
    let m, n = AD.Mat.shape theta.wh in
    let z = AD.F Float.(sqrt (of_int n / of_int m)) in
    { theta with wh = z * theta.wh / sqrt (sum ~axis:1 (sqr theta.wh)) }


  let generate_bs ~n ~m =
    Option.map X.n_beg ~f:(fun nb ->
      let nr = Int.(n / nb) in
      assert (nr = m);
      ( nb
      , Array.init nb ~f:(fun i ->
          let inr = Int.(i * nr) in
          let rnr = Int.(n - ((i + 1) * nr)) in
          transpose
            (concatenate
               ~axis:0
               [| AD.Mat.zeros inr m; AD.Mat.eye nr; AD.Mat.zeros rnr m |])) ))


  let dyn ~theta =
    let theta = with_wh_rescaled theta in
    let n = AD.Mat.col_num theta.bh in
    let m = AD.Mat.row_num theta.wh in
    let beg_bs = generate_bs ~n ~m in
    let default x u =
      let h_pred = x in
      let f = sigma ((u *@ theta.wf) + theta.bf + (h_pred *@ theta.uf)) in
      let h_hat =
        let hf = h_pred * f in
        phi ((u *@ theta.wh) + theta.bh + (hf *@ theta.uh))
      in
      ((F 1. - f) * h_pred) + (f * h_hat)
    in
    fun ~k ~x ~u ->
      match beg_bs with
      | None -> default x u
      | Some (nb, beg_b) -> if k < nb then x + (u *@ beg_b.(k)) else default x u


  let dyn_x =
    let _dyn_x ~theta =
      let theta = with_wh_rescaled theta in
      let n = AD.Mat.col_num theta.bh in
      let id_n = AD.Mat.eye n in
      let default x u =
        let h_pred = x in
        let f_pre = (u *@ theta.wf) + theta.bf + (h_pred *@ theta.uf) in
        let f = sigma f_pre in
        let h_hat_pre =
          let hf = h_pred * f in
          (u *@ theta.wh) + theta.bh + (hf *@ theta.uh)
        in
        let h_hat = phi h_hat_pre in
        diagm (F 1. - f)
        - (theta.uf * ((h_pred - h_hat) * d_sigma f_pre))
        + (((transpose f * theta.uh)
            + (theta.uf *@ (transpose (h_pred * d_sigma f_pre) * theta.uh)))
           * (f * d_phi h_hat_pre))
      in
      fun ~k ~x ~u ->
        match X.n_beg with
        | None -> default x u
        | Some i -> if k < i then id_n else default x u
    in
    Some _dyn_x


  let dyn_u =
    let _dyn_u ~theta =
      let theta = with_wh_rescaled theta in
      let m = AD.Mat.row_num theta.wh in
      let n = AD.Mat.col_num theta.bh in
      let beg_bs = generate_bs ~n ~m in
      let default x u =
        let h_pred = x in
        let f_pre = (u *@ theta.wf) + theta.bf + (h_pred *@ theta.uf) in
        let f = sigma f_pre in
        let h_hat_pre =
          let hf = h_pred * f in
          (u *@ theta.wh) + theta.bh + (hf *@ theta.uh)
        in
        let h_hat = phi h_hat_pre in
        (theta.wf * (d_sigma f_pre * (h_hat - h_pred)))
        + ((theta.wh + (theta.wf *@ (transpose (h_pred * d_sigma f_pre) * theta.uh)))
           * (f * d_phi h_hat_pre))
      in
      fun ~k ~x ~u ->
        match beg_bs with
        | None -> default x u
        | Some (nb, beg_b) -> if k < nb then beg_b.(k) else default x u
    in
    Some _dyn_u
end

module MGU2 (X : sig
    val phi : AD.t -> AD.t
    val d_phi : AD.t -> AD.t
    val sigma : AD.t -> AD.t
    val d_sigma : AD.t -> AD.t
    val n_beg : int Option.t
  end) =
struct
  module P = Dynamics_intf.MGU2_P
  open P
  open X
  open AD.Maths

  let requires_linesearch = true

  let init ~n ~m () =
    (* h : size 1xN
       x : size 1xN (x = Bu)
       h : size 1xK
       f : size of h so 1xN *)
    { wh = Prms.create (AD.Mat.gaussian m n)
    ; bh = Prms.create (AD.Mat.zeros 1 n)
    ; bf = Prms.create (AD.Mat.zeros 1 n)
    ; uh = Prms.create (AD.Mat.zeros n n)
    ; uf = Prms.create (AD.Mat.zeros n n)
    }


  let with_wh_rescaled theta =
    let m, n = AD.Mat.shape theta.wh in
    let z = AD.F Float.(sqrt (of_int n / of_int m)) in
    { theta with wh = z * theta.wh / sqrt (sum ~axis:1 (sqr theta.wh)) }


  let generate_bs ~n ~m =
    Option.map X.n_beg ~f:(fun nb ->
      let nr = Int.(n / nb) in
      assert (nr = m);
      ( nb
      , Array.init nb ~f:(fun i ->
          let inr = Int.(i * nr) in
          let rnr = Int.(n - ((i + 1) * nr)) in
          AD.Maths.(
            transpose
              (concatenate
                 ~axis:0
                 [| AD.Mat.zeros inr m; AD.Mat.eye nr; AD.Mat.zeros rnr m |]))) ))


  let dyn ~theta =
    let theta = with_wh_rescaled theta in
    let n = AD.Mat.col_num theta.bh in
    let m = AD.Mat.row_num theta.wh in
    let beg_bs = generate_bs ~n ~m in
    let default x u =
      let h_pred = x in
      let f = sigma (theta.bf + (h_pred *@ theta.uf)) in
      let h_hat =
        let hf = h_pred * f in
        phi (theta.bh + (hf *@ theta.uh)) + (u *@ theta.wh)
      in
      ((F 1. - f) * h_pred) + (f * h_hat)
    in
    fun ~k ~x ~u ->
      match beg_bs with
      | None -> default x u
      | Some (nb, beg_b) -> if k < nb then x + (u *@ beg_b.(k)) else default x u


  let dyn_x =
    let _dyn_x ~theta =
      let theta = with_wh_rescaled theta in
      let n = AD.Mat.col_num theta.bh in
      let id_n = AD.Mat.eye n in
      let default x u =
        let h_pred = x in
        let f_pre = theta.bf + (h_pred *@ theta.uf) in
        let f = sigma f_pre in
        let h_hat_pre =
          let hf = h_pred * f in
          theta.bh + (hf *@ theta.uh)
        in
        let h_hat = phi h_hat_pre + (u *@ theta.wh) in
        diagm (F 1. - f)
        - (theta.uf * ((h_pred - h_hat) * d_sigma f_pre))
        + (((transpose f * theta.uh)
            + (theta.uf *@ (transpose (h_pred * d_sigma f_pre) * theta.uh)))
           * (f * d_phi h_hat_pre))
      in
      fun ~k ~x ~u ->
        match X.n_beg with
        | None -> default x u
        | Some i -> if k < i then id_n else default x u
    in
    Some _dyn_x


  let dyn_u =
    let _dyn_u ~theta =
      let theta = with_wh_rescaled theta in
      let m = AD.Mat.row_num theta.wh in
      let n = AD.Mat.col_num theta.bh in
      let beg_bs = generate_bs ~n ~m in
      let default x =
        let h_pred = x in
        let f_pre = theta.bf + (h_pred *@ theta.uf) in
        let f = sigma f_pre in
        theta.wh * f
      in
      fun ~k ~x ~u:_ ->
        match beg_bs with
        | None -> default x
        | Some (nb, beg_b) -> if k < nb then beg_b.(k) else default x
    in
    Some _dyn_u
end
