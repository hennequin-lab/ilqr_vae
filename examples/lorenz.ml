open Base
open Ilqr_vae
open Misc
open Lorenz_common
open Vae

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let in_dir = Printf.sprintf "%s/%s" dir
let reuse_data = Cmdargs.check "-reuse_data"
let max_iter = Cmdargs.(get_int "-max_iter" |> default 40_000)
let setup = { n = 3; nh = 64; m = 3; n_trials = 112; n_steps = 100 }
let n_output = 3
let noise_std = 0.1

module M = Make_model (struct
    let setup = setup
    let n_beg = Some (setup.n / setup.m)
  end)

open M

let reg ~(prms : Model.P.t') =
  let z = Float.(1e-5 / of_int Int.(setup.n * setup.n)) in
  let part1 = AD.Maths.(F z * l2norm_sqr' prms.dynamics.a1) in
  let part2 = AD.Maths.(F z * l2norm_sqr' prms.dynamics.a2) in
  AD.Maths.(part1 + part2)


(* -----------------------------------------
   -- Generate Lorenz data
   ----------------------------------------- *)

let _ = C.print_endline "Data generation..."

(* generate training and test data right away *)
let data =
  if reuse_data
  then Misc.read_bin (in_dir "train_data.bin")
  else
    C.broadcast' (fun () ->
      let data =
        Lorenz_common.generate_from_long ~n_steps:setup.n_steps (2 * setup.n_trials)
        |> (fun v -> AA.reshape v [| -1; 3 |])
        |> (fun v -> AA.((v - mean ~axis:0 v) / sqrt (var ~axis:0 v)))
        |> (fun v -> AA.reshape v [| -1; setup.n_steps; 3 |])
        |> (fun v ->
        Array.init (2 * setup.n_trials) ~f:(fun k -> AA.(squeeze (get_slice [ [ k ] ] v))))
        |> Array.map ~f:(fun z ->
          let o = AA.(z + gaussian ~sigma:noise_std (shape z)) in
          (* here I'm hijacking z to store the Lorenz traj *)
          { u = None; z = Some (AD.pack_arr z); o = AD.pack_arr o })
      in
      let train_data = Array.sub data ~pos:0 ~len:setup.n_trials in
      let test_data = Array.sub data ~pos:setup.n_trials ~len:setup.n_trials in
      Misc.save_bin ~out:(in_dir "train_data.bin") train_data;
      Misc.save_bin ~out:(in_dir "test_data.bin") test_data;
      let save_data label data =
        Array.iteri data ~f:(fun i data ->
          let file label' = in_dir (Printf.sprintf "%s_data_%s_%i" label label' i) in
          Option.iter data.z ~f:(fun z ->
            AA.save_txt ~out:(file "latent") (AD.unpack_arr z));
          L.save_data ~prefix:(file "o") data.o)
      in
      save_data "train" train_data;
      save_data "test" test_data;
      train_data)


let _ = C.print_endline "Data generated and broadcast."

(* -----------------------------------------
   -- Initialise parameters and train
   ----------------------------------------- *)

let init_prms =
  C.broadcast' (fun () ->
    match Cmdargs.get_string "-reuse" with
    | Some file -> Misc.read_bin file
    | None ->
      let n = setup.n
      and nh = setup.nh
      and m = setup.m in
      (* let prior = U.init ~spatial_std:1.0 ~nu:20. ~m () in *)
      let prior = U.init ~spatial_std:1.0 ~m () in
      let prior_recog = UR.init ~spatial_std:1.0 ~m () in
      let dynamics = D.init ~radius:0.01 ~decay:0.2 ~n ~nh ~m () in
      let likelihood = L.init ~sigma2:Float.(square noise_std) ~n ~n_output () in
      Model.init ~prior ~prior_recog ~dynamics ~likelihood ())


let save_results prefix prms data =
  let prms = C.broadcast prms in
  let file s = prefix ^ "." ^ s in
  C.root_perform (fun () ->
    Misc.save_bin ~out:(file "params.bin") prms;
    Model.P.save_txt ~prefix prms);
  let prms = Model.P.value prms in
  Array.iteri data ~f:(fun i dat_trial ->
    if Int.(i % C.n_nodes = C.rank)
    then (
      let mu = Model.posterior_mean ~prms dat_trial in
      AA.save_txt ~out:(file (Printf.sprintf "posterior_u_%i" i)) (AD.unpack_arr mu);
      let us, zs, os = Model.predictions ~n_samples:100 ~prms mu in
      let process label a =
        let a = AD.unpack_arr a in
        AA.(mean ~axis:2 a @|| var ~axis:2 a)
        |> (fun z -> AA.reshape z [| setup.n_steps; -1 |])
        |> AA.save_txt ~out:(file (Printf.sprintf "predicted_%s_%i" label i))
      in
      process "u" us;
      process "z" zs;
      Array.iter ~f:(fun (label, x) -> process label x) os))


let _ = save_results (in_dir "init") init_prms data

module Optimizer = Opt.Adam.Make (Model.P)

let config k =
  Opt.Adam.
    { default_config with
      learning_rate = Some Float.(0.04 / sqrt (1. + of_int k))
    ; beta2 = 0.99
    ; epsilon = 1e-7
    }


let rec iter ~k state =
  let prms = C.broadcast (Optimizer.v state) in
  if Int.(k % 200 = 0) then save_results (in_dir "final") prms data;
  let loss, g =
    Model.elbo_gradient ~n_samples:10 ~mini_batch:8 ~conv_threshold:1E-4 ~reg prms data
  in
  (if C.first
   then AA.(save_txt ~append:true ~out:(in_dir "loss") (of_array [| loss |] [| 1; 1 |])));
  print [%message (k : int) (loss : float)];
  let state =
    match g with
    | None -> state
    | Some g -> Optimizer.step ~config:(config k) ~info:g state
  in
  if k < max_iter then iter ~k:(k + 1) state else Optimizer.v state


let final_prms = iter ~k:0 (Optimizer.init init_prms)
let _ = save_results (in_dir "final") final_prms data
