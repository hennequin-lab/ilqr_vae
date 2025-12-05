type t = float * float

let increase (delta, mu) =
  let delta = max 2. (2. *. delta) in
  let mu = max 1E-6 (mu *. delta) in
  delta, mu


let decrease (delta, mu) =
  let delta = min (1. /. 2.) (delta /. 2.) in
  let mu = max 1E-6 (mu *. delta) in
  delta, mu
