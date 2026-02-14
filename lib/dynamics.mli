open Misc

module type T = Dynamics_intf.T

module Integrate (D : T) : sig
  val integrate : prms:D.P.t' -> n:int -> u:AD.t -> AD.t
end

module Linear (X : sig
    val n_beg : int Option.t
  end) : sig
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

module Linear_unconstrained (X : sig
    val n_beg : int Option.t
  end) : sig
  include T with module P = Dynamics_intf.Linear_unconstrained_P

  val init
    :  dt_over_tau:float
    -> alpha:float
    -> beta:float
    -> n:int
    -> m:int
    -> unit
    -> P.t
end

module Nonlinear (X : sig
    val phi : [ `linear | `nonlinear of (AD.t -> AD.t) * (AD.t -> AD.t) ]
    val n_beg : int Option.t
  end) : sig
  include T with module P = Dynamics_intf.Nonlinear_P

  val init : ?radius:float -> n:int -> m:int -> unit -> P.t
end

module InvertedBottleneck (X : sig
    val phi : (AD.t -> AD.t) * (AD.t -> AD.t)
    val n_beg : int Option.t
  end) : sig
  include T with module P = Dynamics_intf.InvertedBottleneck_P

  val init : ?radius:float -> ?decay:float -> n:int -> nh:int -> m:int -> unit -> P.t
end

(** Mini-GRU (Heck, 2017) *)
module MGU (X : sig
    val phi : AD.t -> AD.t
    val d_phi : AD.t -> AD.t
    val sigma : AD.t -> AD.t
    val d_sigma : AD.t -> AD.t
    val n_beg : int Option.t
  end) : sig
  include T with module P = Dynamics_intf.MGU_P

  val init : n:int -> m:int -> unit -> P.t
end

(** Mini-GRU, 2rd simplification (Heck, 2017) *)
module MGU2 (X : sig
    val phi : AD.t -> AD.t
    val d_phi : AD.t -> AD.t
    val sigma : AD.t -> AD.t
    val d_sigma : AD.t -> AD.t
    val n_beg : int Option.t
  end) : sig
  include T with module P = Dynamics_intf.MGU2_P

  val init : n:int -> m:int -> unit -> P.t
end

module GNODE (X : sig
    val n_beg : int Option.t
    val phi : (AD.t -> AD.t) * (AD.t -> AD.t)
  end) : sig
  include T with module P = Dynamics_intf.GNODE_P

  val init
    :  ?radius:float
    -> dt:float
    -> tau:float
    -> n:int
    -> m:int
    -> nh:int
    -> unit
    -> P.t
end
