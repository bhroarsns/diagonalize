extern crate nalgebra as na;

mod lanczos;
mod power;
pub mod qr;
mod rand_matrix;

use na::{Dynamic, OMatrix};
use std::env;

#[allow(unused_imports)]
use std::fs::File;

#[allow(unused_imports)]
use std::io::{BufWriter, Read, Write};

#[allow(dead_code)]
type DMatrixf64 = OMatrix<f64, Dynamic, Dynamic>;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args: Vec<String> = env::args().collect();

    // let nsize = &args[1].trim().parse::<usize>()?;
    // let repeat = &args[2].trim().parse::<usize>()?;

    // let nsize = &args[1].trim().parse::<usize>()?;
    // let times = &args[2].trim().parse::<usize>()?;
    // let para = &args[3].trim().parse::<usize>()?;

    let nsize = &args[1].trim().parse::<usize>()?;
    let times = &args[2].trim().parse::<usize>()?;
    let repeat = &args[3].trim().parse::<usize>()?;
    let para = &args[4].trim().parse::<usize>()?;

    // let mut f = File::open("./matrix.dat")?;
    // let mut matrix = String::new();
    // f.read_to_string(&mut matrix)?;
    // let a = DMatrixf64::from_vec(*nsize, *nsize, matrix.split('\n').flat_map(|str| {
    //     str.split_whitespace().map(|s| s.trim().parse().unwrap())
    // }).collect());
    // let a = rand_matrix::get_rand_with_eigens((0..10).map(|v| v as f64).collect());
    let a = rand_matrix::get_sparse_and_sym(*nsize, *para);
    // let adense = nalgebra_sparse::convert::serial::convert_csr_dense(&a);
    // let a = rand_matrix::that_matrix(*nsize);

    lanczos::eigenvalues_csr(&a, *times, &format!("./lanczos_{}_{}.dat", nsize, para))?;
    // lanczos::eigenvalues(&a, *times, "./lanczos2.dat")?;
    println!("lanczos method finished");

    // qr::eigenvalues_record(&a, *repeat, "./qr.dat")?;
    // println!("qr method finished");

    // power::max_eigenvalue_csr(&a, *repeat, &format!("./power_{}.dat", *nsize))?;
    power::max_eigenvalue_csr(&a, *repeat, &format!("./power_{}_{}.dat", nsize, para))?;
    // power::max_eigenvalue(&a, *repeat, "./power2.dat")?;
    println!("power method finished");

    // println!("right: {:?}", a.eigenvalues());

    Ok(())
}