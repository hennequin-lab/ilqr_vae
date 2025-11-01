open Misc

module type T = Likelihood_intf.T

module Gaussian (X : sig
    val label : string
    val normalize_c : bool
  end) : sig
  include
    T
    with module P = Likelihood_intf.Gaussian_P
     and type datum = AD.t
     and type data = AD.t

  val init
    :  ?scale:float
    -> ?sigma2:float
    -> ?bias:float
    -> n:int
    -> n_output:int
    -> unit
    -> P.t
end

module Poisson (X : sig
    val label : string
    val dt : AD.t
    val link_function : AD.t -> AD.t
    val d_link_function : AD.t -> AD.t
    val d2_link_function : AD.t -> AD.t
  end) : sig
  include
    T with module P = Likelihood_intf.Poisson_P and type datum = AD.t and type data = AD.t

  val init : n:int -> n_output:int -> unit -> P.t
  val pre_sample_before_link_function : prms:P.t' -> z:AD.t -> AD.t
end

module Pair (L1 : T) (L2 : T) :
  T
  with module P = Prms.Pair(L1.P)(L2.P)
   and type datum = L1.datum * L2.datum
   and type data = L1.data * L2.data
