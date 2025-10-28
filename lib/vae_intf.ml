open Base
open Misc

type 'a data =
  { u : AD.t option
  ; z : AD.t option
  ; o : 'a
  }

module Generative_P = struct
  type ('u, 'd, 'l) prm =
    { prior : 'u
    ; dynamics : 'd
    ; likelihood : 'l
    }
  [@@deriving prms, accessors ~submodule:A]
end

module Recognition_P = struct
  type ('g_opt, 'cov) prm =
    { generative : 'g_opt
    ; space_cov : 'cov
    ; time_cov : 'cov
    }
  [@@deriving prms, accessors ~submodule:A]
end

module VAE_P = struct
  type ('gen, 'recog) p =
    { generative : 'gen
    ; recognition : 'recog
    }
  [@@deriving prms, accessors ~submodule:A]
end
