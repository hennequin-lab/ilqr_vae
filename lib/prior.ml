open Base
open Misc
open Owl

module type T = Prior_intf.T

module Gaussian = struct
  module P = Prior_intf.Gaussian_P
  open P
  open AD.Maths

  let requires_linesearch = false

  let init ?(spatial_std = 1.) ?(first_bin = 1.) ~m () =
    let spatial_stds = AD.Mat.create 1 m spatial_std in
    { spatial_stds = Prms.create ~above:(F 1E-3) spatial_stds
    ; first_bin = Prms.create ~above:(F 1E-5) (AD.F first_bin)
    }


  (* returns a column vector *)
  let temporal_stds ~prms ~n_steps =
    let t1 = prms.first_bin * AD.Mat.ones 1 1 in
    let t2 = AD.Mat.ones Int.(n_steps - 1) 1 in
    concat ~axis:0 t1 t2


  let kl_to_gaussian =
    let kl ~prms ~mu ~space ~time =
      let ell_q_s = Covariance.to_chol_factor space in
      let ell_q_t = Covariance.to_chol_factor time in
      let m = AD.Mat.row_num ell_q_s in
      let t = AD.Mat.row_num ell_q_t in
      let mm = Float.of_int m in
      let tt = Float.of_int t in
      let dim = Float.(mm * tt) in
      let ell_p_s = prms.spatial_stds in
      let ell_p_t = temporal_stds ~prms ~n_steps:t in
      let quadratic_term = l2norm_sqr' (F 1. / ell_p_t * mu / ell_p_s) in
      let trace_term =
        let s = l2norm_sqr' (ell_q_s / ell_p_s) in
        let t = l2norm_sqr' (ell_q_t / transpose ell_p_t) in
        s * t
      in
      let logdet_term =
        let log_det_q =
          F 2. * ((F mm * sum' (log time.d)) + (F tt * sum' (log space.d)))
        in
        let log_det_p =
          let d_space = ell_p_s in
          let d_time = ell_p_t in
          F 2. * ((F mm * sum' (log d_time)) + (F tt * sum' (log d_space)))
        in
        log_det_p - log_det_q
      in
      F 0.5 * (trace_term + logdet_term + quadratic_term - F dim)
    in
    `direct kl


  let sample ~prms ~n_steps ~m =
    let ell_t = temporal_stds ~prms ~n_steps in
    let ell_s = prms.spatial_stds in
    let xi = AD.Mat.gaussian n_steps m in
    ell_t * xi * ell_s


  let neg_logp_t ~prms =
    let ell_s = prms.spatial_stds in
    let m = AD.Mat.numel ell_s in
    let mlpi2 = Float.(of_int m * log Const.pi2) in
    let ell_s_0 = prms.first_bin * ell_s in
    let cst_0 = F mlpi2 + (F 2. * sum' (log ell_s_0)) in
    let cst_k = F mlpi2 + (F 2. * sum' (log ell_s)) in
    fun ~k ~x:_ ~u ->
      let sigma = if k = 0 then ell_s_0 else ell_s in
      let cst = if k = 0 then cst_0 else cst_k in
      F 0.5 * (cst + l2norm_sqr' (u / sigma))


  let neg_jac_t =
    let jac_t ~prms =
      let ell_s = prms.spatial_stds in
      let ell_s_0 = prms.first_bin * ell_s in
      fun ~k ~x:_ ~u ->
        let sigma = if k = 0 then ell_s_0 else ell_s in
        u / sqr sigma
    in
    Some jac_t


  let neg_hess_t =
    let hess_t ~prms =
      let sigma2 = sqr prms.spatial_stds in
      let h_0 = diagm (F 1. / (sqr prms.first_bin * sigma2)) in
      let h_k = diagm (F 1. / sigma2) in
      fun ~k ~x:_ ~u:_ -> if k = 0 then h_0 else h_k
    in
    Some hess_t


  (* shouldn't ever need this *)
  let logp ~prms:_ = assert false
end

module Student = struct
  module P = Prior_intf.Student_P
  open P
  open AD.Maths

  let requires_linesearch = true

  let init ?(pin_std = false) ?(spatial_std = 1.) ?(nu = 10.) ~m () =
    let spatial_stds = Prms.create ~above:(F 1E-3) (AD.Mat.create 1 m spatial_std) in
    { spatial_stds = (if pin_std then Prms.pin spatial_stds else spatial_stds)
    ; nu = Prms.create ~above:(F 2.0) (AD.F nu)
    ; first_step = spatial_stds
    }


  let kl_to_gaussian = `sampling_based

  let get_eff_prms ~prms =
    let nu = prms.nu in
    let sigma = sqrt ((prms.nu - F 2.) / prms.nu) * prms.spatial_stds in
    nu, sigma


  (* non-differentiable *)
  let sample ~prms ~n_steps ~m =
    let nu, sigma = get_eff_prms ~prms in
    let xi = AA.(gaussian [| Int.(n_steps - 1); m |] * AD.unpack_arr sigma) in
    let u = Stats.chi2_rvs ~df:(AD.unpack_flt nu) in
    let z = Float.(sqrt (AD.unpack_flt nu / u)) in
    let z = AA.(z $* xi) in
    let z0 = AA.(gaussian [| 1; m |] * AD.unpack_arr prms.first_step) in
    AD.pack_arr (AA.concatenate ~axis:0 [| z0; z |])


  let neg_logp_t ~prms =
    let nu, sigma = get_eff_prms ~prms in
    let m = AD.Mat.numel sigma in
    let m_half = AD.F Float.(of_int m / 2.) in
    let nu_half = F 0.5 * nu in
    let nu_plus_m_half = F 0.5 * (nu + F Float.(of_int m)) in
    let sigma0 = prms.first_step in
    let cst0 = Float.(of_int m * log Const.pi2) in
    let cst =
      let cst1 = AD.loggamma nu_half - AD.loggamma nu_plus_m_half in
      let cst2 = m_half * log (F Const.pi * nu) in
      let cst3 = sum' (log sigma) in
      cst1 + cst2 + cst3
    in
    fun ~k ~x:_ ~u ->
      if k = 0
      then F 0.5 * (F cst0 + (F 2. * sum' (log sigma0)) + l2norm_sqr' (u / sigma0))
      else (
        let utilde = u / sigma in
        cst + (nu_plus_m_half * log (F 1. + (l2norm_sqr' utilde / nu))))


  let neg_jac_t =
    let jac_t ~prms =
      let nu, sigma = get_eff_prms ~prms in
      let m = AD.Mat.numel sigma in
      let nu_plus_m_half = F 0.5 * (nu + F Float.(of_int m)) in
      let sigma2 = sqr sigma in
      fun ~k ~x:_ ~u ->
        let stu =
          let tmp =
            let utilde = u / sigma in
            F 1. + (l2norm_sqr' utilde / nu)
          in
          let tmp' = F 2. * u / sigma2 / nu in
          nu_plus_m_half * tmp' / tmp
        in
        if k = 0 then u / sqr prms.first_step else stu
    in
    Some jac_t


  let neg_hess_t =
    let hess_t ~prms =
      let nu, sigma = get_eff_prms ~prms in
      let m = AD.Mat.numel sigma in
      let nu_plus_m_half = F 0.5 * (nu + F Float.(of_int m)) in
      let sigma2 = sqr sigma in
      fun ~k ~x:_ ~u ->
        let stu =
          let u_over_s2 = u / sigma2 in
          let tau = F 1. + l2norm_sqr' (u / sigma / nu) in
          let cst = F 2. * nu_plus_m_half / nu / sqr tau in
          let term1 = diagm (tau / sigma2) in
          let term2 = F 2. * (transpose u_over_s2 *@ u_over_s2) / nu in
          cst * (term1 - term2)
        in
        if k = 0 then diagm (F 1. / sqr prms.first_step) else stu
    in
    Some hess_t


  let logp ~prms ~n_steps =
    let nu, sigma = get_eff_prms ~prms in
    let sigma0 = prms.first_step in
    let m = AD.Mat.numel sigma in
    let m_half = AD.F Float.(of_int m / 2.) in
    let nu_half = F 0.5 * nu in
    let nu_plus_m_half = m_half + nu_half in
    let cst0 = F Float.(of_int m * log Const.pi2) + (F 2. * sum' (log sigma0)) in
    let cst =
      let cst1 = AD.loggamma nu_half - AD.loggamma nu_plus_m_half in
      let cst2 = m_half * log (F Const.pi * nu) in
      let cst3 = sum' (log sigma) in
      F Float.(of_int n_steps) * (cst1 + cst2 + cst3)
    in
    fun u ->
      (* u is K x T x M *)
      let u_s = AD.shape u in
      let n_samples = u_s.(0) in
      let u0 = get_slice [ []; [ 0 ]; [] ] u |> fun v -> reshape v [| -1; m |] in
      let u = get_slice [ []; [ 1; -1 ]; [] ] u in
      assert (Array.length u_s = 3);
      let cst0 = F Float.(of_int n_samples) * cst0 in
      let cst = F Float.(of_int n_samples) * cst in
      let u = reshape u [| -1; m |] in
      let first_term = F 0.5 * (cst0 + l2norm_sqr' (u0 / sigma0)) in
      let rest =
        cst + (nu_plus_m_half * sum' (log (F 1. + (sum ~axis:1 (sqr (u / sigma)) / nu))))
      in
      neg (first_term + rest)
end
