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
    { d : 'a
    ; u : 'a
    ; q : 'a
    ; b : 'opt
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single) (Prms.Option (Prms.Single))
end

module Nonlinear_P = struct
  type ('a, 'opt) prm =
    { a : 'a
    ; bias : 'a
    ; b : 'opt
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single) (Prms.Option (Prms.Single))
end

module InvertedBottleneck_P = struct
  type ('a, 'opt) prm =
    { a1 : 'a
    ; a2 : 'a
    ; bias1 : 'a
    ; bias2 : 'a
    ; b : 'opt
    ; decay : 'a
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single) (Prms.Option (Prms.Single))
end

module MGU_P = struct
  type 'a prm =
    { wf : 'a
    ; uf : 'a
    ; wh : 'a
    ; uh : 'a
    ; bf : 'a
    ; bh : 'a
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single)
end

module MGU2_P = struct
  type 'a prm =
    { uf : 'a
    ; wh : 'a
    ; uh : 'a
    ; bf : 'a
    ; bh : 'a
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single)
end
