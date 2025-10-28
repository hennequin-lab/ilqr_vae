module Make (A : Prms.Intf.A) = struct
  module AD = struct
    include Owl_algodiff_generic.Make (A)

    module Primal = struct
      let print_dim x =
        let shp = A.shape x in
        Array.iter (fun s -> Printf.printf "%i " s) shp;
        print_newline ();
        flush_all ()


      (* batch matrix multiplication *)
      let bmm =
        let check_dims_match ndim shpx shpy =
          if shpx.(ndim - 1) <> shpy.(ndim - 2)
          then failwith "bmm: last two dimensions mismatch"
          else
            for i = 0 to ndim - 3 do
              if shpx.(i) <> shpy.(i) then failwith "bmm: do not support broadcast"
            done
        in
        fun x y ->
          let shpx = A.shape x in
          let shpy = A.shape y in
          let ndimx = Array.length shpx in
          let ndimy = Array.length shpy in
          if ndimx <> ndimy
          then failwith "bmm: dimensions of [x] and [y] must be the same";
          if ndimx < 2
          then failwith "bmm: dimensions must be greater than 2"
          else if ndimx = 2
          then A.dot x y
          else (
            let ndim = ndimx in
            check_dims_match ndim shpx shpy;
            let shp = Array.copy shpx in
            let m = shpx.(ndim - 2) in
            let k = shpx.(ndim - 1) in
            let l = shpy.(ndim - 2) in
            assert (k == l);
            let n = shpy.(ndim - 1) in
            shp.(ndim - 1) <- n;
            shp.(ndim - 2) <- m;
            let batch_size = Array.fold_left ( * ) 1 shp / m / n in
            let z = A.empty shp in
            (* Printf.printf "%i, %i, %i %i\n%!" m n k batch_size; *)
            let xr = A.reshape x [| batch_size; m; k |] in
            let yr = A.reshape y [| batch_size; k; n |] in
            let zr = A.reshape z [| batch_size; m; n |] in
            for j = 0 to batch_size - 1 do
              let x1 = Bigarray.Genarray.sub_left xr j 1 in
              let x2 = Bigarray.Genarray.sub_left yr j 1 in
              let x3 = Bigarray.Genarray.sub_left zr j 1 in
              let alpha = 1. in
              let beta = 0. in
              let a = A.reshape x1 [| -1 |] |> Bigarray.array1_of_genarray in
              let b = A.reshape x2 [| -1 |] |> Bigarray.array1_of_genarray in
              let c = A.reshape x3 [| -1 |] |> Bigarray.array1_of_genarray in
              let layout = Owl_cblas_basic.CblasRowMajor in
              let transa = Owl_cblas_basic.CblasNoTrans in
              let transb = Owl_cblas_basic.CblasNoTrans in
              Owl_cblas_basic.gemm layout transa transb m n k alpha a k b n beta c n
            done;
            z)


      let bchol ?(upper = true) x =
        let shp = A.shape x in
        let ndims = Array.length shp in
        if ndims < 2
        then failwith "bchol: dimension must be greater than 2"
        else if ndims = 2
        then unpack_arr (Linalg.chol ~upper (pack_arr x))
        else (
          let x = A.copy x in
          let m = shp.(ndims - 1) in
          let n = shp.(ndims - 2) in
          if m <> n
          then failwith "bchol: last two dimensions do not match"
          else (
            let batch_size = Array.fold_left ( * ) 1 shp / m / n in
            let xr = A.reshape x [| batch_size; m; n |] in
            for j = 0 to batch_size - 1 do
              let a =
                Bigarray.Genarray.sub_left xr j 1 |> fun x -> A.reshape x [| m; m |]
              in
              if upper
              then (
                Owl_lapacke.potrf ~uplo:'U' ~a |> ignore;
                for t = 1 to m - 1 do
                  for s = 0 to t - 1 do
                    A.set a [| t; s |] 0.
                  done
                done)
              else (
                Owl_lapacke.potrf ~uplo:'L' ~a |> ignore;
                for t = 1 to m - 1 do
                  for s = 0 to t - 1 do
                    A.set a [| s; t |] 0.
                  done
                done)
            done;
            x))
    end

    let print_dim x =
      let shp = A.shape (unpack_arr x) in
      Array.iter (fun s -> Printf.printf "%i " s) shp;
      print_newline ();
      flush_all ()


    let rec _bmm =
      lazy
        (Builder.build_piso
           (module struct
             let label = "bmm"
             let ff_aa _ _ = raise Owl_exception.(NOT_IMPLEMENTED "bmm")
             let ff_ab _ _ = raise Owl_exception.(NOT_IMPLEMENTED "bmm")
             let ff_ba _ _ = raise Owl_exception.(NOT_IMPLEMENTED "bmm")
             let ff_bb a b = pack_arr (Primal.bmm a b)
             let df_da _cp _ap at bp = bmm at bp
             let df_db _cp ap _bp bt = bmm ap bt
             let df_dab _cp ap at bp bt = Maths.(bmm ap bt + bmm at bp)

             let dr_ab a b _cp ca =
               let dim_a = Array.length (shape a) in
               let dim_b = Array.length (shape b) in
               ( bmm !ca (Maths.swap (dim_b - 1) (dim_b - 2) (primal b))
               , bmm (Maths.swap (dim_a - 1) (dim_a - 2) (primal a)) !ca )


             let dr_a _a b _cp ca =
               let dim = Array.length (shape b) in
               bmm !ca (Maths.swap (dim - 1) (dim - 2) (primal b))


             let dr_b a _b _cp ca =
               let dim = Array.length (shape a) in
               bmm (Maths.swap (dim - 1) (dim - 2) (primal a)) !ca
           end : Builder.Piso))


    and bmm x = Stdlib.Lazy.force _bmm x
  end

  include AD.Primal
end
