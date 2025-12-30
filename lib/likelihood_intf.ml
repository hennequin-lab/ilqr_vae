open Base
open Misc

module type T = sig
  module P : Prms.T

  type datum
  type data

  val requires_linesearch : bool
  val label : string
  val save_data : ?prefix:string -> data -> unit
  val data_slice : k:int -> data -> datum
  val to_mat_list : data -> (string * AD.t) list
  val size : prms:P.t' -> int
  val pre_sample : prms:P.t' -> z:AD.t -> data
  val sample_noise : mu:data -> prms:P.t' -> data
  val neg_logp_t : prms:P.t' -> data_t:datum -> k:int -> z_t:AD.t -> AD.t
  val neg_jac_t : (prms:P.t' -> data_t:datum -> k:int -> z_t:AD.t -> AD.t) option
  val neg_hess_t : (prms:P.t' -> data_t:datum -> k:int -> z_t:AD.t -> AD.t) option
  val logp : prms:P.t' -> data:data -> z:AD.t -> AD.t
end

module Gaussian_P = struct
  type 'a prm =
    { c : 'a
    ; c_mask : AD.t option
    ; bias : 'a
    ; variances : 'a (* 1 x space *)
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single)
end

module Poisson_P = struct
  type 'a prm =
    { c : 'a
    ; c_mask : AD.t option
    ; bias : 'a
    ; gain : 'a
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single)
end
