module Make (A : Prms.Intf.A) : sig
  val print_dim : A.arr -> unit
  val bmm : A.arr -> A.arr -> A.arr
  val bchol : ?upper:bool -> A.arr -> A.arr

  module AD : sig
    include module type of Owl_algodiff_generic.Make (A)

    val print_dim : t -> unit
    val bmm : t -> t -> t
  end
end
