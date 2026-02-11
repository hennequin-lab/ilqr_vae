open Base

(* change this AA module to change the precision globally *)
module AA = Prms.Intf.S

module DILQR = Dilqr.Default.Make (struct
    module A = AA

    let verbose = false
  end)

module Prms = Prms.Make (AA)
module Opt = Opt.Make (AA)

module AD = struct
  include Owl_algodiff_generic.Make (AA)

  module Mat = struct
    include Mat

    let create n m x = pack_arr (AA.create [| n; m |] x)
  end

  let print_shape ~label x =
    x
    |> shape
    |> Array.to_list
    |> List.map ~f:Int.to_string
    |> String.concat ~sep:" x "
    |> Stdio.printf "[%s] shape = %s\n%!" label


  let expand0 x = Maths.reshape x (Array.append [| 1 |] (shape x))

  let squeeze0 x =
    let s = shape x in
    assert (s.(0) = 1);
    assert (Array.length s = 3);
    Maths.reshape x [| s.(1); s.(2) |]


  let rec requad x = Stdlib.Lazy.force _requad x
  and d_requad x = Stdlib.Lazy.force _d_requad x
  and d2_requad x = Stdlib.Lazy.force _d2_requad x
  and rq x = Float.(0.5 * (x + sqrt (4. + square x)))
  and rqv x = AA.((x + sqrt (4. $+ sqr x)) *$ 0.5)

  and drq x =
    let tmp = Float.(sqrt (4. + square x)) in
    Float.(0.5 * (1. + (x / tmp)))


  and drqv x =
    let tmp = Maths.(sqrt (F 4. + sqr x)) in
    Maths.(F 0.5 * (F 1. + (x / tmp)))


  and drq_arr x =
    let tmp = AA.(sqrt (4. $+ sqr x)) in
    AA.(0.5 $* (1. $+ x / tmp))


  and ddrq x =
    let tmp = Float.(4. + square x) in
    let tmp = Float.(tmp * sqrt tmp) in
    Float.(2. / tmp)


  and ddrq_arr x =
    let tmp = AA.(4. $+ sqr x) in
    let tmp = AA.(tmp * sqrt tmp) in
    AA.(2. $/ tmp)


  and ddrqv x =
    let tmp = Maths.(F 4. + sqr x) in
    let tmp = Maths.(tmp * sqrt tmp) in
    Maths.(F 2. / tmp)


  and dddrqv x =
    let tmp = Maths.(F 4. + sqr x) in
    let tmp = Maths.(sqr tmp * sqrt tmp) in
    Maths.(F (-6.) * x / tmp)


  and _requad =
    let open Builder in
    lazy
      (build_siso
         (module struct
           let label = "requad"
           let ff_f a = F (rq a)
           let ff_arr a = Arr (rqv a)
           let df _cp ap at = Maths.(at * drqv ap)
           let dr a _cp ca = Maths.(!ca * drqv a)
         end : Siso))


  and _d_requad =
    let open Builder in
    lazy
      (build_siso
         (module struct
           let label = "d_requad"
           let ff_f a = F (drq a)
           let ff_arr a = Arr (drq_arr a)
           let df _cp ap at = Maths.(at * ddrqv ap)
           let dr a _cp ca = Maths.(!ca * ddrqv a)
         end : Siso))


  and _d2_requad =
    let open Builder in
    lazy
      (build_siso
         (module struct
           let label = "d2_requad"
           let ff_f a = F (ddrq a)
           let ff_arr a = Arr (ddrq_arr a)
           let df _cp ap at = Maths.(at * dddrqv ap)
           let dr a _cp ca = Maths.(!ca * dddrqv a)
         end : Siso))


  (* log-gamma, reverse-mode only ! *)
  let rec loggamma x = Stdlib.Lazy.force _loggamma x

  and _loggamma =
    let open Builder in
    lazy
      (build_siso
         (module struct
           let label = "loggamma"
           let ff_f a = F (Owl_maths_special.loggamma a)
           let ff_arr _ = assert false
           let df _cp ap at = Maths.(at * F (Owl_maths_special.psi (unpack_flt ap)))
           let dr a _cp ca = Maths.(!ca * F (Owl_maths_special.psi (unpack_flt a)))
         end : Siso))
end

let save_bin ~out:filename m =
  let output = Stdio.Out_channel.create filename in
  Stdlib.Marshal.to_channel output m [ Stdlib.Marshal.No_sharing ];
  Stdio.Out_channel.close output


let read_bin filename =
  let input = Stdio.In_channel.create filename in
  let m = Stdlib.Marshal.from_channel input in
  Stdio.In_channel.close input;
  m


let time_this ~label f =
  let t0 = Unix.gettimeofday () in
  let res = f () in
  let t = Unix.gettimeofday () -. t0 in
  C.print_endline (Printf.sprintf "[%s] %f seconds\n%!" label t);
  res


let print s = if C.first then Stdio.print_endline (Sexp.to_string_hum s)

let with_prefix ?prefix s =
  match prefix with
  | None -> s
  | Some p -> p ^ "." ^ s
