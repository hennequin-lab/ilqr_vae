open Base

(* let print s = Stdio.print_endline (Sexp.to_string_hum s) *)

module Make (A : Prms.Intf.A) = struct
  module Lqr = Lqr.Make (A)
  module Bmo = Bmo.Make (A)
  module AD = Bmo.AD
  open AD.Builder
  open AD.Maths

  type 'a t = theta:'a -> k:int -> x:AD.t -> u:AD.t -> AD.t
  type 'a s = theta:'a -> k:int -> x:AD.t -> AD.t
  type 'a final_loss = theta:'a -> k:int -> x:AD.t -> AD.t
  type 'a running_loss = theta:'a -> k:int -> x:AD.t -> u:AD.t -> AD.t

  let forward_for_backward
        ~theta
        ~dyn_x
        ~dyn_u
        ~rl_uu
        ~rl_xx
        ~rl_ux
        ~rl_u
        ~rl_x
        ~fl_xx
        ~fl_x
        ~dyn
    =
    let dyn = dyn ~theta
    and dyn_x = dyn_x ~theta
    and dyn_u = dyn_u ~theta
    and rl_x = rl_x ~theta
    and rl_u = rl_u ~theta
    and rl_xx = rl_xx ~theta
    and rl_uu = rl_uu ~theta
    and rl_ux = rl_ux ~theta in
    let fl_xx = fl_xx ~theta in
    let fl_x = fl_x ~theta in
    fun () x0 us ->
      let kf, xf, tape =
        List.fold us ~init:(0, x0, []) ~f:(fun (k, x, tape) u ->
          let a = dyn_x ~k ~x ~u
          and b = dyn_u ~k ~x ~u
          and rlx = rl_x ~k ~x ~u
          and rlu = rl_u ~k ~x ~u in
          let rlxx = rl_xx ~k ~x ~u in
          let rluu = rl_uu ~k ~x ~u
          and rlux = rl_ux ~k ~x ~u in
          let f = dyn ~k ~x ~u - (x *@ a) - (u *@ b) in
          let s = Lqr.{ x; u; a; b; rlx; rlu; rlxx; rluu; rlux; f } in
          let x = dyn ~k ~x ~u in
          Int.(k + 1), x, s :: tape)
      in
      let flxx = fl_xx ~k:kf ~x:xf in
      let flx = fl_x ~k:kf ~x:xf in
      flxx, flx, tape, xf


  module type P = sig
    type theta

    val primal' : theta -> theta
    val n : int
    val m : int
    val n_steps : int
    val dyn : theta t
    val final_loss : theta final_loss
    val running_loss : theta running_loss
    val dyn_x : theta t option
    val dyn_u : theta t option
    val rl_uu : theta t option
    val rl_xx : theta t option
    val rl_ux : theta t option
    val rl_u : theta t option
    val rl_x : theta t option
    val fl_xx : theta s option
    val fl_x : theta s option
  end

  module Make (P : P) = struct
    include P

    let zeros_n = AD.Mat.zeros 1 n
    let zeros_m = AD.Mat.zeros 1 m
    let zeros_nm = AD.Mat.zeros n m
    let zeros_1nmn = AD.Arr.zeros [| 1; Int.(n + m); n |]
    let zeros_nmn = AD.Mat.zeros Int.(n + m) n
    let zeros_mnm = AD.Mat.zeros m Int.(n + m)

    let dyn_u =
      let default ~theta =
        let dyn = dyn ~theta in
        fun ~k ~x ~u -> AD.jacobian (fun u -> dyn ~k ~x ~u) u
      in
      Option.value dyn_u ~default


    let dyn_x =
      let default ~theta =
        let dyn = dyn ~theta in
        fun ~k ~x ~u -> AD.jacobian (fun x -> dyn ~k ~x ~u) x
      in
      Option.value dyn_x ~default


    let rl_u =
      let default ~theta =
        let running_loss = running_loss ~theta in
        fun ~k ~x ~u -> AD.grad (fun u -> running_loss ~k ~x ~u) u
      in
      Option.value rl_u ~default


    let rl_x =
      let default ~theta =
        let running_loss = running_loss ~theta in
        fun ~k ~x ~u -> AD.grad (fun x -> running_loss ~k ~x ~u) x
      in
      Option.value rl_x ~default


    let rl_uu =
      let default ~theta =
        let rl_u = rl_u ~theta in
        fun ~k ~x ~u -> AD.jacobian (fun u -> rl_u ~k ~x ~u) u
      in
      Option.value rl_uu ~default


    let rl_xx =
      let default ~theta =
        let rl_x = rl_x ~theta in
        fun ~k ~x ~u -> AD.jacobian (fun x -> rl_x ~k ~x ~u) x
      in
      Option.value rl_xx ~default


    let rl_ux =
      let default ~theta =
        let rl_x = rl_x ~theta in
        fun ~k ~x ~u -> AD.jacobian (fun u -> rl_x ~k ~x ~u) u
      in
      Option.value rl_ux ~default


    let fl_x =
      let default ~theta =
        let final_loss = final_loss ~theta in
        fun ~k ~x -> AD.grad (fun x -> final_loss ~k ~x) x
      in
      Option.value fl_x ~default


    let fl_xx =
      let default ~theta =
        let fl_x = fl_x ~theta in
        fun ~k ~x -> AD.jacobian (fun x -> fl_x ~k ~x) x
      in
      Option.value fl_xx ~default


    let forward ~theta =
      let dyn = dyn ~theta in
      fun x0 us ->
        List.fold us ~init:(0, x0, [], []) ~f:(fun (k, x, xs, us) u ->
          let xs = x :: xs in
          let us = u :: us in
          let x = dyn ~k ~x ~u in
          Int.(k + 1), x, xs, us)


    let ffb ~theta =
      forward_for_backward
        ~theta
        ~dyn_x
        ~dyn_u
        ~rl_uu
        ~rl_xx
        ~rl_ux
        ~rl_u
        ~rl_x
        ~fl_xx
        ~fl_x
        ~dyn
        ()


    let update ~theta =
      let ffb = ffb ~theta in
      let dyn = dyn ~theta in
      fun x0 us ->
        (* xf, xs, us are in reverse *)
        let vxxf, vxf, tape, _ = ffb x0 us in
        let acc, (df1, df2) = Lqr.backward vxxf vxf tape in
        fun alpha ->
          let _, _, uhats =
            List.fold
              acc
              ~init:(0, x0, [])
              ~f:(fun (k, xhat, uhats) ((s : Lqr.t), (_K, _k)) ->
                let dx = xhat - s.x in
                let du = (dx *@ _K) + (F alpha * _k) in
                let uhat = s.u + du in
                let uhats = uhat :: uhats in
                let xhat = dyn ~k ~x:xhat ~u:uhat in
                Int.(k + 1), xhat, uhats)
          in
          let df = Float.((alpha * df1) + (0.5 * (alpha * alpha * df2))) in
          List.rev uhats, df


    let trajectory ~theta =
      let forward = forward ~theta in
      fun x0 us ->
        let _, xf, xs, _ = forward x0 us in
        let xs = List.rev xs |> Array.of_list |> concatenate ~axis:0 in
        concatenate ~axis:0 [| xs; xf |]


    let loss ~theta =
      let forward = forward ~theta in
      let running_loss = running_loss ~theta in
      let final_loss = final_loss ~theta in
      fun x0 us ->
        let kf, xf, xs, us = forward x0 us in
        let fl = final_loss ~k:kf ~x:xf in
        let _, rl =
          List.fold2_exn
            xs
            us
            ~init:(Int.(kf - 1), AD.F 0.)
            ~f:(fun (k, rl) x u -> Int.(k - 1), rl + running_loss ~k ~x ~u)
        in
        fl + rl |> AD.unpack_flt


    let differentiable_loss ~theta =
      let final_loss = final_loss ~theta in
      let running_loss = running_loss ~theta in
      fun taus_f ->
        let array_taus =
          split ~axis:0 (Array.create ~len:(AD.Arr.shape taus_f).(0) 1) taus_f
        in
        let tf = Array.length array_taus in
        Array.foldi array_taus ~init:(AD.F 0.) ~f:(fun i accu tau ->
          let tau = AD.Maths.reshape tau Int.[| 1; n + m |] in
          let x, u =
            get_slice [ []; [ 0; Int.pred n ] ] tau, get_slice [ []; [ n; -1 ] ] tau
          in
          if Int.(i = pred tf)
          then accu + final_loss ~k:(Int.succ i) ~x
          else accu + running_loss ~k:(Int.succ i) ~x ~u)


    let learn ?(max_iter = 10) ?(conv_threshold = 1e-4) ?(linesearch = true) ~theta =
      let loss = loss ~theta in
      let update = update ~theta in
      fun x0 us ->
        let rec loop iter (c : float) us =
          let has_converged c' = Float.(abs ((c' - c) / c) < conv_threshold) in
          if iter > max_iter
          then us
          else (
            let f0 = loss x0 us in
            let update = update x0 us in
            let f alpha =
              let us, df = update alpha in
              let fv = loss x0 us in
              fv, Some df, us
            in
            if not linesearch
            then (
              let c', _, us = f 1. in
              if has_converged c' then us else loop (Int.succ iter) c' us)
            else (
              match Linesearch.backtrack f0 f with
              | c', Some us when has_converged c' -> us
              | c', Some us -> loop (Int.succ iter) c' us
              | _, None -> failwith "linesearch did not converge"))
        in
        loop 0 1E8 us


    let g2 =
      let swap_out_tape tape tau_bar =
        (* swapping out the tape *)
        let _, tape =
          (* the tape is backward in time hence we fold_right *)
          List.fold_right tape ~init:(0, []) ~f:(fun (s : Lqr.t) (k, tape) ->
            let rlx =
              reshape (get_slice [ [ k ]; []; [ 0; Int.pred n ] ] tau_bar) [| 1; n |]
            in
            let rlu = reshape (get_slice [ [ k ]; []; [ n; -1 ] ] tau_bar) [| 1; m |] in
            Int.succ k, Lqr.{ s with rlu; rlx } :: tape)
        in
        let flx =
          reshape (get_slice [ [ -1 ]; []; [ 0; Int.pred n ] ] tau_bar) [| 1; n |]
        in
        flx, tape
      in
      fun ~theta ->
        let ffb = ffb ~theta in
        fun ~taus ~ustars ~lambdas ->
          let ds ~x0 ~tau_bar =
            (* recreating tape, pass as argument in the future *)
            let flxx, _, tape, _ = ffb x0 ustars in
            let flx, tape = swap_out_tape tape tau_bar in
            let acc, _ = Lqr.backward flxx flx tape in
            let ctbars_xf, ctbars_tape = Lqr.forward acc zeros_n in
            let dlambda0, dlambdas = Lqr.adjoint_back ctbars_xf flxx flx ctbars_tape in
            let ctbars =
              List.map ctbars_tape ~f:(fun (s : Lqr.t) ->
                concatenate ~axis:1 [| s.x; s.u |])
              |> List.cons (concatenate ~axis:1 [| ctbars_xf; zeros_m |])
              |> List.rev
            in
            ( stack ~axis:0 (Array.of_list ctbars)
            , stack ~axis:0 (Array.of_list (dlambda0 :: dlambdas)) )
          in
          let big_ft_bar ~taus ~lambdas ~dlambdas ~ctbars () =
            let tdl =
              Bmo.AD.bmm
                (transpose ~axis:[| 0; 2; 1 |] (get_slice [ [ 0; -2 ]; []; [] ] taus))
                (get_slice [ [ 1; -1 ]; []; [] ] dlambdas)
            in
            let dtl =
              Bmo.AD.bmm
                (transpose ~axis:[| 0; 2; 1 |] (get_slice [ [ 0; -2 ]; []; [] ] ctbars))
                (get_slice [ [ 1; -1 ]; []; [] ] lambdas)
            in
            let output = tdl + dtl in
            concatenate ~axis:0 [| output; zeros_1nmn |]
          in
          let big_ct_bar ~taus ~ctbars () =
            let tdt = Bmo.AD.bmm (transpose ~axis:[| 0; 2; 1 |] ctbars) taus in
            F 0.5 * (tdt + transpose ~axis:[| 0; 2; 1 |] tdt)
          in
          build_aiso
            (module struct
              let label = "g2"
              let ff _ = AD.primal' taus
              let df _ _ _ _ = raise (Owl_exception.NOT_IMPLEMENTED "g2 forward mode")

              let dr idxs x _ ybar =
                let x0 = x.(4) in
                let ctbars, dlambdas = ds ~x0 ~tau_bar:!ybar in
                List.map idxs ~f:(fun idx ->
                  if idx = 0
                  then big_ft_bar ~taus ~lambdas ~dlambdas ~ctbars ()
                  else if idx = 1
                  then big_ct_bar ~taus ~ctbars ()
                  else if idx = 2
                  then ctbars
                  else if idx = 3
                  then get_slice [ [ 1; -1 ] ] dlambdas
                  else get_slice [ [ 0 ] ] dlambdas |> fun x -> reshape x [| 1; -1 |])
            end : Aiso)


    let g1 ~theta =
      let ffb = ffb ~theta in
      fun ~x0 ~ustars ->
        let flxx, flx, tape, xf = ffb x0 ustars in
        let lambda0, lambdas = Lqr.adjoint flx tape in
        let lambdas = stack ~axis:0 (Array.of_list (lambda0 :: lambdas)) in
        let big_taus = [ concatenate ~axis:1 [| xf; zeros_m |] ] in
        let big_fs = [ zeros_nmn ] in
        let big_cs =
          let row1 = concatenate ~axis:1 [| flxx; zeros_nm |] in
          let row2 = zeros_mnm in
          [ concatenate ~axis:0 [| row1; row2 |] ]
        in
        let cs = [ concatenate ~axis:1 [| flx - (xf *@ flxx); zeros_m |] ] in
        let fs = [] in
        let big_taus, big_fs, big_cs, cs, fs, _ =
          List.fold
            tape
            ~init:(big_taus, big_fs, big_cs, cs, fs, xf)
            ~f:(fun (taus, big_fs, big_cs, cs, fs, next_x) (s : Lqr.t) ->
              ignore next_x;
              let taus =
                let tau = concatenate ~axis:1 [| s.x; s.u |] in
                tau :: taus
              in
              let big_f = concatenate ~axis:0 [| s.a; s.b |] in
              let big_c =
                let row1 = concatenate ~axis:1 [| s.rlxx; transpose s.rlux |] in
                let row2 = concatenate ~axis:1 [| s.rlux; s.rluu |] in
                concatenate ~axis:0 [| row1; row2 |]
              in
              let c =
                concatenate
                  ~axis:1
                  [| s.rlx - (s.x *@ s.rlxx) - (s.u *@ s.rlux)
                   ; s.rlu - (s.u *@ s.rluu) - (s.x *@ transpose s.rlux)
                  |]
              in
              taus, big_f :: big_fs, big_c :: big_cs, c :: cs, s.f :: fs, s.x)
        in
        let taus = stack ~axis:0 Array.(of_list big_taus) in
        let big_fs = stack ~axis:0 Array.(of_list big_fs) in
        let big_cs = stack ~axis:0 Array.(of_list big_cs) in
        let cs = stack ~axis:0 Array.(of_list cs) in
        let fs = stack ~axis:0 Array.(of_list fs) in
        taus, big_fs, big_cs, cs, lambdas, fs


    let ilqr ?max_iter ?conv_threshold ?(linesearch = true) ~theta =
      let theta' = primal' theta in
      let g1 = g1 ~theta in
      fun ~us ~x0 () ->
        let ustars =
          learn ?max_iter ?conv_threshold ~linesearch ~theta:theta' AD.(primal' x0) us
          |> List.map ~f:AD.primal'
        in
        let taus, big_fs, big_cs, cs, lambdas, fs = g1 ~x0:(AD.primal' x0) ~ustars in
        let inp = [| big_fs; big_cs; cs; fs; x0 |] in
        g2 ~lambdas:(AD.primal' lambdas) ~taus:(AD.primal' taus) ~ustars ~theta:theta' inp
  end
end
