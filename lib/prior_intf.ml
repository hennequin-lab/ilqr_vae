open Base
open Misc

module type T = sig
  module P : Prms.T

  val requires_linesearch : bool
  (* val spatial_stds : prms:P.t' -> AD.t *)

  val kl_to_gaussian
    : [ `direct of
          prms:P.t' -> mu:AD.t -> space:Covariance.P.t' -> time:Covariance.P.t' -> AD.t
      | `sampling_based
      ]

  val sample : prms:P.t' -> n_steps:int -> m:int -> AD.t
  val neg_logp_t : prms:P.t' -> k:int -> x:AD.t -> u:AD.t -> AD.t
  val neg_jac_t : (prms:P.t' -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
  val neg_hess_t : (prms:P.t' -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
  val logp : prms:P.t' -> n_steps:int -> AD.t -> AD.t
end

module Gaussian_P = struct
  type 'a prm =
    { spatial_stds : 'a
    ; first_bin : 'a
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single)
end

module Student_P = struct
  type 'a prm =
    { spatial_stds : 'a
    ; nu : 'a
    ; first_step : 'a
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single)
end
