open Base

let backtrack ?(alpha = 1.) ?(alpha_min = 1e-8) ?(tau = 0.5) ?(beta = 0.1) f0 f =
  let rec backtrack alpha =
    let fv, df, prms = f alpha in
    if
      match df with
      | Some df -> Float.(f0 <= fv + (beta *. df))
      | None -> Float.(f0 < fv)
    then if Float.(alpha < alpha_min) then fv, None else backtrack Float.(tau * alpha)
    else fv, Some prms
  in
  backtrack (tau *. alpha)
