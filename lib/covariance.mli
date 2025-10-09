open Misc
module P = Covariance_intf.P

val init
  :  ?floor:float
  -> ?no_triangle:bool
  -> ?pin_diag:bool
  -> ?sigma2:float
  -> int
  -> P.t

val to_chol_factor : P.t' -> AD.t
