open Base
open Misc

module type T = sig
  module P : Prms.T

  val requires_linesearch : bool
  val dyn : theta:P.t' -> k:int -> x:AD.t -> u:AD.t -> AD.t
  val dyn_x : (theta:P.t' -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
  val dyn_u : (theta:P.t' -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
end

module Linear_P = struct
  type ('a, 'opt) prm =
    { dt_over_tau : float
    ; d : 'a
    ; u : 'a
    ; q : 'a
    ; b : 'opt
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single) (Prms.Option (Prms.Single))
end

module Linear_nonlinear_P = struct
  type ('a, 'opt) prm =
    { dt_over_tau : 'a
    ; a0 : 'a
    ; a1 : 'a
    ; a2 : 'a
    ; bias1 : 'a
    ; bias2 : 'a
    ; b : 'opt
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single) (Prms.Option (Prms.Single))
end

module GNODE_P = struct
  type ('a, 'opt) prm =
    { dt_over_tau : float
    ; g_a1 : 'a
    ; g_bias1 : 'a
    ; h_a1 : 'a
    ; h_a2 : 'a
    ; h_bias1 : 'a
    ; h_bias2 : 'a
    ; b : 'opt
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single) (Prms.Option (Prms.Single))
end
