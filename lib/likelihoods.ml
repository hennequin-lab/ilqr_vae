open Base
open Misc
open Owl

module type T = Likelihoods_intf.T

let with_prefix ?prefix s =
  match prefix with
  | None -> s
  | Some p -> p ^ "." ^ s


module Gaussian (X : sig
    val label : string
    val normalize_c : bool
  end) =
struct
  module P = Likelihoods_intf.Gaussian_P
  open P
  open AD.Maths

  type datum = AD.t
  type data = AD.t

  let requires_linesearch = false
  let label = X.label

  let init ?(sigma2 = 1.) ?(bias = 0.) ~n ~n_output () =
    { c = Prms.create (AD.Mat.gaussian ~sigma:Float.(1. / sqrt (of_int n)) n_output n)
    ; c_mask = None
    ; bias = Prms.create (AD.Mat.create 1 n_output bias)
    ; variances =
        Prms.create ~above:(F 0.001) (AD.pack_arr (Mat.create 1 n_output sigma2))
    }


  let save_data ?prefix data =
    Mat.save_txt ~out:(with_prefix ?prefix label) (AD.unpack_arr data)


  let data_slice ~k data = get_slice [ [ k ] ] data
  let to_mat_list x = [ label, x ]
  let size ~prms = AD.Mat.row_num prms.c

  let unpack_c ~prms =
    let c =
      match prms.c_mask with
      | None -> prms.c
      | Some cm -> prms.c * cm
    in
    if X.normalize_c then c / sqrt (sum ~axis:1 (sqr c)) else c


  let pre_sample ~prms ~z =
    let c = unpack_c ~prms in
    (* z is T x M *)
    prms.bias + (z *@ transpose c)


  let sample ~prms ~z =
    let mu = pre_sample ~prms ~z in
    let res =
      let xi = AD.Arr.(gaussian (shape mu)) in
      xi * sqrt prms.variances
    in
    AD.Maths.(mu + res)


  let neg_logp_t ~prms ~data_t =
    let c = unpack_c ~prms in
    let c_t = transpose c in
    let n = AD.Mat.row_num c |> Float.of_int in
    let cst = AD.F Float.(n * log Const.pi2) in
    let sum_log_var = sum' (log prms.variances) in
    fun ~k:_ ~z_t ->
      let mu_t = prms.bias + (z_t *@ c_t) in
      assert (Poly.(AD.shape mu_t = AD.shape data_t));
      F 0.5 * (cst + sum_log_var + sum' (sqr (data_t - mu_t) / prms.variances))


  let neg_jac_t =
    let neg_jac_t ~prms ~data_t =
      let c = unpack_c ~prms in
      let c_t = transpose c in
      let c_inv_variances = F 1. / transpose prms.variances * c in
      fun ~k:_ ~z_t ->
        (*z = 1*N, c = 0*N, data = 1*0*)
        let mu_t = prms.bias + (z_t *@ c_t) in
        (mu_t - data_t) *@ c_inv_variances
    in
    Some neg_jac_t


  let neg_hess_t =
    let neg_hess_t ~prms ~data_t:_ =
      let c = unpack_c ~prms in
      let c_inv_variances = AD.Maths.(F 1. / transpose prms.variances * c) in
      let tmp = transpose c *@ c_inv_variances in
      (*z = 1*N, c = 0*N, data = 1*0*)
      fun ~k:_ ~z_t:_ -> tmp
    in
    Some neg_hess_t


  let logp ~prms ~data =
    let c = unpack_c ~prms in
    let c_t = transpose c in
    let n_out = AD.Mat.row_num c in
    let n = AD.Mat.col_num c in
    let sum_log_var = sum' (log prms.variances) in
    let data = AD.expand0 data in
    fun ~z ->
      let z_s = AD.shape z in
      assert (Array.length z_s = 3);
      let n_samples = z_s.(0) in
      let n_steps = z_s.(1) in
      let z = reshape z [| -1; n |] in
      let mu = prms.bias + (z *@ c_t) in
      let diff =
        let mu = reshape mu [| n_samples; n_steps; n_out |] in
        reshape (data - mu) [| -1; n_out |]
      in
      let cst = AD.F Float.(of_int Int.(n_samples * n_steps * n_out) * log Const.pi2) in
      F (-0.5)
      * (cst
         + (F Float.(of_int Int.(n_steps * n_samples)) * sum_log_var)
         + sum' (sqr diff / prms.variances))
end

module Poisson (X : sig
    val label : string
    val dt : AD.t
    val link_function : AD.t -> AD.t
    val d_link_function : AD.t -> AD.t
    val d2_link_function : AD.t -> AD.t
  end) =
struct
  module P = Likelihoods_intf.Poisson_P
  open P
  open X
  open AD.Maths

  type datum = AD.t
  type data = AD.t

  let requires_linesearch = true
  let label = X.label

  let init ~n ~n_output () =
    { c = Prms.create (AD.Mat.gaussian ~sigma:Float.(1. / sqrt (of_int n)) n_output n)
    ; c_mask = None
    ; bias = Prms.create (AD.Mat.zeros 1 n_output)
    ; gain = Prms.create (AD.Mat.ones 1 n_output)
    }


  let save_data ?prefix data =
    Mat.save_txt ~out:(with_prefix ?prefix label) (AD.unpack_arr data)


  let data_slice ~k data = get_slice [ [ k ] ] data
  let to_mat_list x = [ label, x ]
  let size ~prms = AD.Mat.row_num prms.c

  let unpack_c ~prms =
    match prms.c_mask with
    | None -> prms.c
    | Some cm -> AD.Maths.(prms.c * cm)


  let pre_sample_before_link_function ~prms ~z =
    let c = unpack_c ~prms in
    (* z is T x M *)
    (z *@ transpose c) + prms.bias


  let pre_sample ~prms ~z = link_function (pre_sample_before_link_function ~prms ~z)

  let sample ~prms ~z =
    let t = AD.Mat.row_num z in
    (* z is T x M *)
    let mu = pre_sample ~prms ~z in
    let spikes =
      Owl_distribution_generic.poisson_rvs ~mu:(AD.unpack_arr (mu * dt)) ~n:1
    in
    reshape (AD.pack_arr spikes) [| t; -1 |]


  let logfact k =
    let rec iter k accu =
      if k <= 1 then accu else iter Int.(k - 1) Float.(accu + log (of_int k))
    in
    iter k 0.


  (* redefine the link_function to include a safe floor *)
  let link_function x = AD.Maths.(AD.F 1E-3 + link_function x)

  let neg_logp_t ~prms =
    let c = unpack_c ~prms in
    let c_t = transpose c in
    fun ~data_t ->
      let logfact =
        AD.pack_arr (Mat.map (fun x -> logfact Int.(of_float x)) (AD.unpack_arr data_t))
      in
      fun ~k:_ ~z_t ->
        let rate_t = dt * prms.gain * link_function ((z_t *@ c_t) + prms.bias) in
        let log_rate_t = log rate_t in
        assert (Poly.(AD.shape log_rate_t = AD.shape data_t));
        sum' (rate_t + logfact - (data_t * log_rate_t))


  let d_log_link_function x = d_link_function x / link_function x

  let d2_log_link_function x =
    let lx = link_function x in
    let dlx = d_link_function x in
    let ddlx = d2_link_function x in
    ((ddlx * lx) - sqr dlx) / sqr lx


  let neg_jac_t =
    (* dlogp/dz *)
    let neg_jac_t ~prms =
      let c = unpack_c ~prms in
      let c_t = transpose c in
      fun ~data_t ~k:_ ~z_t ->
        (* 1 x M *)
        let a = prms.bias + (z_t *@ c_t) in
        let tmp1 = dt * prms.gain * d_link_function a in
        let tmp2 = data_t * d_log_link_function a in
        (tmp1 - tmp2) *@ c
    in
    Some neg_jac_t


  let neg_hess_t =
    let neg_hess_t ~prms =
      let c = unpack_c ~prms in
      let c_t = transpose c in
      fun ~data_t ~k:_ ~z_t ->
        (* 1 x M *)
        let a = prms.bias + (z_t *@ c_t) in
        let tmp1 = dt * prms.gain * d2_link_function a in
        let tmp2 = data_t * d2_log_link_function a in
        transpose c * (tmp1 - tmp2) *@ c
    in
    Some neg_hess_t


  let logp ~prms =
    let c = unpack_c ~prms in
    let c_t = transpose c in
    let n_out = AD.Mat.row_num c in
    let n = AD.Mat.col_num c in
    fun ~data ->
      let logfact =
        AD.pack_arr (Mat.map (fun x -> logfact Int.(of_float x)) (AD.unpack_arr data))
      in
      let data = AD.expand0 data in
      fun ~z ->
        let z_s = AD.shape z in
        assert (Array.length z_s = 3);
        let n_samples = z_s.(0) in
        let z = reshape z [| -1; n |] in
        let rates = dt * prms.gain * link_function ((z *@ c_t) + prms.bias) in
        let rates = reshape rates [| n_samples; -1; n_out |] in
        let log_rates = log rates in
        AD.Maths.(
          sum' ((data * log_rates) - rates) - (F Float.(of_int n_samples) * sum' logfact))
end

module Pair (L1 : T) (L2 : T) = struct
  module P = Prms.Pair (L1.P) (L2.P)
  open P
  open AD.Maths

  type datum = L1.datum * L2.datum
  type data = L1.data * L2.data

  let requires_linesearch = L1.requires_linesearch || L2.requires_linesearch
  let label = Printf.sprintf "pair(%s-%s)" L1.label L2.label

  let save_data ?prefix (data1, data2) =
    L1.save_data ?prefix data1;
    L2.save_data ?prefix data2


  let data_slice ~k (data1, data2) = L1.data_slice ~k data1, L2.data_slice ~k data2

  let to_mat_list (data1, data2) =
    List.concat [ L1.to_mat_list data1; L2.to_mat_list data2 ]


  let size ~prms:(prms1, prms2) = Int.(L1.size ~prms:prms1 + L2.size ~prms:prms2)

  let pre_sample ~prms:(prms1, prms2) ~z =
    L1.pre_sample ~prms:prms1 ~z, L2.pre_sample ~prms:prms2 ~z


  let sample ~prms:(prms1, prms2) ~z = L1.sample ~prms:prms1 ~z, L2.sample ~prms:prms2 ~z

  let add f1 f2 ~prms:(prms1, prms2) =
    let f1 = f1 ~prms:prms1 in
    let f2 = f2 ~prms:prms2 in
    fun ~data_t:(d1, d2) ->
      let f1 = f1 ~data_t:d1 in
      let f2 = f2 ~data_t:d2 in
      fun ~k ~z_t -> f1 ~k ~z_t + f2 ~k ~z_t


  let neg_logp_t = add L1.neg_logp_t L2.neg_logp_t

  let neg_jac_t =
    match L1.neg_jac_t, L2.neg_jac_t with
    | Some f1, Some f2 -> Some (add f1 f2)
    | _ -> None


  let neg_hess_t =
    match L1.neg_hess_t, L2.neg_hess_t with
    | Some f1, Some f2 -> Some (add f1 f2)
    | _ -> None


  let logp ~prms:(prms1, prms2) =
    let f1 = L1.logp ~prms:prms1 in
    let f2 = L2.logp ~prms:prms2 in
    fun ~data:(d1, d2) ->
      let f1 = f1 ~data:d1 in
      let f2 = f2 ~data:d2 in
      fun ~z -> f1 ~z + f2 ~z
end
