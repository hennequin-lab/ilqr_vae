open Misc
include module type of Vae_intf

module ILQR (U : Prior.T) (D : Dynamics.T) (L : Likelihood.T) : sig
  module P : module type of ILQR_P.Make (U.P) (D.P) (L.P)

  val solve
    :  ?conv_threshold:float
    -> ?n_beg:int
    -> ?saving_iter:string
    -> u_init:AA.arr option
    -> primal':(P.t' -> P.t')
    -> n:int
    -> m:int
    -> n_steps:int
    -> prms:P.t'
    -> L.data data
    -> AD.t
end

module Make
    (U : Prior.T)
    (UR : Prior.T)
    (D : Dynamics.T)
    (L : Likelihood.T)
    (X : sig
       val n : int
       val m : int
       val n_steps : int
       val n_beg : int Option.t
       val diag_time_cov : bool
     end) : sig
  module P : module type of VAE_P.Make (U.P) (UR.P) (D.P) (L.P) (Covariance.P)
  module Integrate : module type of Dynamics.Integrate (D)
  module Ilqr : module type of ILQR (UR) (D) (L)

  val broadcast_prms : P.t -> P.t

  val init
    :  ?sigma:float
    -> prior:U.P.t
    -> prior_recog:UR.P.t
    -> dynamics:D.P.t
    -> likelihood:L.P.t
    -> unit
    -> P.t

  val sample_generative : prms:P.t' -> L.data data
  val sample_generative_autonomous : sigma:float -> prms:P.t' -> L.data data

  val posterior_mean
    :  ?saving_iter:string
    -> ?conv_threshold:float
    -> ?u_init:AA.arr
    -> prms:P.t'
    -> L.data data
    -> AD.t

  val sample_recognition : prms:P.t' -> mu_u:AD.t -> int -> AD.t

  val predictions
    :  ?pre:bool
    -> n_samples:int
    -> prms:P.t'
    -> AD.t
    -> AD.t * AD.t * (string * AD.t) Array.t

  val elbo
    :  ?conv_threshold:float
    -> mu_u:[ `known of AD.t | `guess of AA.arr option ]
    -> n_samples:int
    -> ?beta:Float.t
    -> prms:P.t'
    -> L.data data
    -> AD.t * AA.arr

  (* returns [Some gradient] at the root note, [None] elsewhere *)
  val elbo_gradient
    :  ?n_samples:int
    -> ?mini_batch:int
    -> ?conv_threshold:float
    -> ?reg:(prms:P.t' -> AD.t)
    -> P.t
    -> L.data data array
    -> float * P.t' option

  (*
     val recalibrate_uncertainty
    :  ?n_samples:(int -> int)
    -> ?max_iter:int
    -> ?save_progress_to:int * int * string
    -> ?in_each_iteration:(u_init:Owl.Mat.mat option Array.t -> prms:P.t' -> int -> unit)
    -> ?eta:[ `constant of float | `of_iter of int -> float ]
    -> prms:P.t'
    -> L.data data Array.t
    -> P.t'

  val check_grad
    :  prms:P.t'
    -> L.data data Array.t
    -> [ `all | `random of int ]
    -> string
    -> unit
  *)

  val save_data : ?prefix:string -> L.data data -> unit
end
