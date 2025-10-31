open Base
open Misc

type 'a data =
  { u : AD.t option
  ; z : AD.t option
  ; o : 'a
  }

module ILQR_P = struct
  type ('u, 'd, 'l) p =
    { prior : 'u
    ; dynamics : 'd
    ; likelihood : 'l
    }
  [@@deriving prms, accessors ~submodule:A]
end

module VAE_P = struct
  type ('u, 'ur, 'd, 'l, 'cov) p =
    { prior : 'u
    ; prior_recog : 'ur
    ; dynamics : 'd
    ; likelihood : 'l
    ; space_cov : 'cov
    ; time_cov : 'cov
    }
  [@@deriving prms, accessors ~submodule:A]
end
