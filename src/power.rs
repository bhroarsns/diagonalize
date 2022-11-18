use na::{U1, Dynamic, OMatrix};
use rand::Rng;
use nalgebra_sparse::CsrMatrix;
use std::fs::File;
use std::io::{BufWriter, Write};

type DMatrixf64 = OMatrix<f64, Dynamic, Dynamic>;
type DVectorf64 = OMatrix<f64, Dynamic, U1>;

#[allow(dead_code)]
pub fn max_eigenvalue(mat: &DMatrixf64, times: usize, dataname: &str) -> Result<(), Box<dyn std::error::Error>> {

    let nsize = mat.ncols();

    let mut rng = rand::thread_rng();
    let rand_vec: Vec<f64> = (0..nsize).map(|_| rng.gen_range(-1.0..1.0)).collect();
    let mut pre_vec = DVectorf64::from_vec(rand_vec);
    let mut cur_vec = mat * &pre_vec;

    let mut cand = &cur_vec.norm().powi(2) / &cur_vec.dot(&pre_vec);

    let mut out = String::new();
    out.push_str(&format!("0 {}", cand));

    for i in 1..(times + 1) {
        let tmp = mat * &cur_vec;

        cand = &cur_vec.norm().powi(2) / &cur_vec.dot(&pre_vec);
        out.push_str(&format!("\n{} {}", i, cand));
        pre_vec = cur_vec.clone();
        cur_vec = tmp;
    }

    let mut writer = BufWriter::new(File::create(dataname)?);
    writer.write_all(out.as_bytes())?;
    writer.flush()?;

    Ok(())
}

#[allow(dead_code)]
pub fn max_eigenvalue_csr(mat: &CsrMatrix<f64>, times: usize, dataname: &str) -> Result<(), Box<dyn std::error::Error>> {

    let nsize = mat.ncols();

    let mut rng = rand::thread_rng();
    let rand_vec: Vec<f64> = (0..nsize).map(|_| rng.gen_range(-1.0..1.0)).collect();
    let mut pre_vec = DVectorf64::from_vec(rand_vec);
    let mut cur_vec = mat * &pre_vec;
    let mut cand;

    let mut out = String::new();

    for i in 0..times {
        cand = &cur_vec.norm().powi(2) / &cur_vec.dot(&pre_vec);
        out.push_str(&format!("\n{} {}", i, cand));

        cur_vec = &cur_vec / cur_vec.norm();
        let tmp = mat * &cur_vec;
        pre_vec = cur_vec.clone();
        cur_vec = tmp;
    }

    let mut writer = BufWriter::new(File::create(dataname)?);
    writer.write_all(out.as_bytes())?;
    writer.flush()?;

    Ok(())
}