open Misc

module type T = Dynamics_intf.T

module Integrate (D : T) : sig
  val integrate : prms:D.P.t' -> n:int -> u:AD.t -> AD.t
end

module Linear : sig
  include T with module P = Dynamics_intf.Linear_P

  val init
    :  dt_over_tau:float
    -> alpha:float
    -> beta:float
    -> n:int
    -> m:int
    -> unit
    -> P.t

  val unpack_a : prms:P.t' -> AD.t
end

module Nonlinear (X : sig
    val phi : [ `linear | `nonlinear of (AD.t -> AD.t) * (AD.t -> AD.t) ]
  end) : sig
  include T with module P = Dynamics_intf.Nonlinear_P

  val init : ?radius:float -> n:int -> m:int -> unit -> P.t
end

module Linear_nonlinear (X : sig
    val phi : (AD.t -> AD.t) * (AD.t -> AD.t)
  end) : sig
  include T with module P = Dynamics_intf.Linear_nonlinear_P

  val init : ?radius:float -> ?decay:float -> n:int -> nh:int -> m:int -> unit -> P.t
end
