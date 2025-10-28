module type T = Prior_intf.T

module Gaussian (X : sig
    val n_beg : int Option.t
  end) : sig
  include T with module P = Prior_intf.Gaussian_P

  val init : ?spatial_std:float -> ?first_bin:float -> m:int -> unit -> P.t
end

module Student (X : sig
    val n_beg : int Option.t
  end) : sig
  include T with module P = Prior_intf.Student_P

  val init : ?pin_std:bool -> ?spatial_std:float -> ?nu:float -> m:int -> unit -> P.t
end
