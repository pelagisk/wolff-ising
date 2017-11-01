(* Implementation of Wolff cluster algorithm for 2D spin lattice, in OCaml *)
(* Written in 2014 by Axel Gagge *)

let ($) f g = function x -> f (g x);;

let rec in_list p = function
  | [] -> false
  | x::l -> if p=x then true else in_list p l;;

let pos_mod a b = let c = a mod b in if c>=0 then c else c+b;;

type spin = Up | Down;;
type lattice = spin array;;

let random_spin () = let c = Random.bool () in if c then Up else Down;;

let critical_ratio = (log (1.0 +. sqrt 2.0)) /. 2.0;;

let p_cluster_coeff temp =
  1.0 -. exp (- 2.0 *. critical_ratio /. temp);;

let probability temp =
  let c = Random.float 1.0 in
  if c < (p_cluster_coeff temp) then true else false;;

let init_lattice n f =
  let fmod i j = f () in
  let row_init i = Array.init n (fmod i)
  in Array.init n row_init;;

let min_img_d n (i,j) (k,l) =
  let image a b =
    let c = a - b in
    if c > (n/2) then abs (c - n)
    else if (-c) > (n/2) then abs (c + n)
    else abs c
  in let fimage a b = float_of_int (image a b)
  in let a,b = (fimage i k), (fimage j l)
  in sqrt (a**2.0 +. b**2.0);;

let lattice_iterij f lattice =
  let row_iter i elt = Array.iteri (f i) elt in
  Array.iteri row_iter lattice;;

(* todo check the type of this! *)
let lattice_fold_left f x lattice =
  let row_fold y row = Array.fold_left f y row in
  Array.fold_left row_fold x lattice;;

let sum_spins lattice =
  let add_up i spin = if spin=Up then i+1 else i-1 in
  lattice_fold_left add_up 0 lattice;;

let sample_spin_corr lattice =
  let n_lattice = Array.length lattice in
  let n_hist = 100 in
  let hist = Array.make n_hist 0 in
  (* todo: all indices seems to be one, why? *)
  let index (i,j) (m,n) = 2 * n_hist * int_of_float (
                                           (min_img_d n_lattice (i,j) (m,n)) /.
                                             (float_of_int n_lattice)) in
  let eachspin (i,j) m n elt =
    if (i,j) = (m,n) then ()
    else
      let cur_index = index (i,j) (m,n) in
      hist.(cur_index) <- hist.(cur_index) + 1
  in let sample_one_spin i j elt = lattice_iterij (eachspin (i,j)) lattice
     in let () = lattice_iterij sample_one_spin lattice
        in hist;;

let print_spin_corr hist =
  Array.iteri (Printf.printf "%d: %d\n") hist;;

let print_lattice lattice printfunc =
  let print_spin i j elt = printfunc (if elt=Up then "0 " else "1 ") in
  let row_iter i elt =
    begin
      Array.iteri (print_spin i) elt;
      printfunc "\n";
    end
  in Array.iteri row_iter lattice;;

let write_lattice lattice temp i =
  let file = Printf.sprintf "lattice_%f_%d.dat" temp i in
  begin
    let oc = open_out file in
    let fprint s = let _ = Printf.fprintf oc "%s" s in () in
    print_lattice lattice fprint;
    close_out oc;
  end;;

let wrapped_pos lattice (i,j) =
  let n = Array.length lattice in
  let corr i = pos_mod i n in
  ((corr i), (corr j));;

let spin_from_pos lattice (i,j) =
  let (ci,cj) = wrapped_pos lattice (i,j) in
  lattice.(ci).(cj);;

let flip_cluster lattice cluster =
  let flip_spin (i,j) =
    let (ci,cj) = wrapped_pos lattice (i,j) in
    let spin = if ((spin_from_pos lattice (i,j))=Down)
               then Up
               else Down
    in lattice.(ci).(cj) <- spin in
  let flipping (i,j) = flip_spin (i,j) in
  List.iter flipping cluster;;

(*
Add neighbor candidates of same spin to list by some probability
function, continue recursively. When done flip cluster (with
probability function we initially set to 1.0)
*)

let rec find_cluster lattice cluster (i,j) temp =
  let this_spin = spin_from_pos lattice (i,j) in
  let in_cluster pos = in_list pos cluster in
  let query_add pos = ((this_spin=(spin_from_pos lattice pos)) &&
                         ((in_cluster pos)=false) &&
                         (probability temp)) in
  begin
    for m=0 to 1
    do
      let new_pos = (i+(2*m-1),j) in
      if (query_add new_pos) = true
      then find_cluster lattice (new_pos::cluster) new_pos temp
      else ()
    done;
    for n=0 to 1
    do
      let new_pos = (i,j+(2*n-1)) in
      if (query_add new_pos) = true
      then find_cluster lattice (new_pos::cluster) new_pos temp
      else ()
    done;
    flip_cluster lattice cluster;
  end;;

let rec monte_carlo_step lattice i c c_samp temp =
  let select_random_spin lattice =
    let n = Array.length lattice in
    let i = Random.int n in
    let j = Random.int n in
    (i,j) in
  let pos = select_random_spin lattice in
  if i<c then
    begin
      Printf.printf "Iteration: %d of %d\n" i c;
      if (pos_mod i (c / c_samp))=0
      then
        begin
          (* print_lattice lattice print_string; *)
          write_lattice lattice temp i;
        end
      else ();
      find_cluster lattice [] pos temp;
      monte_carlo_step lattice (i+1) c c_samp temp;
    end
  else ();;

let rec temp_iter to_run t t_range t_samp =
  let (min_t,max_t) = t_range in
  if t >= max_t then ()
  else
    begin
      to_run t;
      let t_step = (max_t -. min_t) /. (float_of_int t_samp) in
      temp_iter to_run (t +. t_step) t_range t_samp;
    end;;

let run spin_init ~size ~cycles ~c_samp ~t_samp ~t_range =
  let (min_t,_) = t_range in
  let ising_lattice = init_lattice size spin_init in
  let mcstep = monte_carlo_step ising_lattice 0 cycles c_samp in
  begin
    temp_iter mcstep min_t t_range t_samp;
    Printf.printf "Now printing spin corr:\n";
    print_spin_corr (sample_spin_corr ising_lattice);
  end;;

(* todo: write functions to compute spin correlations and magnetisation! *)

run random_spin ~size:100 ~cycles:1000 ~c_samp:10 ~t_samp:10 ~t_range:(0.8,1.2);;
