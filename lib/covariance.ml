open Base
open Misc
module P = Covariance_intf.P
open P

let init ?(floor = 1e-4) ?(no_triangle = false) ?(pin_diag = false) ?(sigma2 = 1.) n =
  let d = Prms.create ~above:(F floor) (AD.Mat.create 1 n Float.(sqrt sigma2)) in
  let t = Prms.create (AD.Mat.zeros n n) in
  { d = (if pin_diag then Prms.pin d else d)
  ; t = (if no_triangle then Prms.pin t else t)
  }


let to_chol_factor c = AD.Maths.(triu ~k:1 c.t + diagm c.d)
