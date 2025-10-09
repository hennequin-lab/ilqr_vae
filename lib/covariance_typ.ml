open Misc

module Covariance_P = struct
  type 'a p =
    { d : 'a
    ; t : 'a
    }
  [@@deriving prms, accessors ~submodule:A]
end
