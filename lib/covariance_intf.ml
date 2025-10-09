open Misc

module P = struct
  type 'a prm =
    { d : 'a
    ; t : 'a
    }
  [@@deriving prms, accessors ~submodule:A]

  include Make (Prms.Single)
end
