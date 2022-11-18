use na::{U1, Dynamic, OMatrix};
use rand::Rng;
use nalgebra_sparse::csr::CsrMatrix;
use std::fs::File;
use std::io::{BufWriter, Write};

#[allow(unused_imports)]
use super::qr;

type DMatrixf64 = OMatrix<f64, Dynamic, Dynamic>;
type DVectorf64 = OMatrix<f64, Dynamic, U1>;

#[allow(dead_code)]
pub fn eigenvalues(mat: &DMatrixf64, times: usize, filename: &str) -> Result<(), Box<dyn std::error::Error>> {

    let nsize = mat.ncols();
    let mut rng = rand::thread_rng();
    let rand_vec: Vec<f64> = (0..nsize).map(|_| rng.gen_range(-1.0..1.0)).collect();
    let mut init_vec = DVectorf64::from_vec(rand_vec);
    let norm = &init_vec.norm();
    init_vec = init_vec / *norm;

    let mut alpha = init_vec.dot(&(mat * &init_vec));
    let mut beta = (mat * &init_vec - alpha * &init_vec).norm();
    let mut beta_pre;
    let mut cur_vec = (mat * &init_vec - alpha * &init_vec) / beta;
    let mut pre_vec = init_vec;

    let mut res = DMatrixf64::from_vec(2, 2, vec![alpha, beta, beta, 0.0]);
    let mut out = String::new();

    for i in 1..times {
        alpha = cur_vec.dot(&(mat * &cur_vec));
        res[(i + 1) * (i + 1) - 1] = alpha;
        out.push_str("\n");
        out.push_str(&format!("{} ", i));
        out.push_str(&qr::eigenvalues(&res, 1.0e-9));

        beta_pre = beta;
        beta = (mat * &cur_vec - beta_pre * &pre_vec  - alpha * &cur_vec).norm();
        let tmp = (mat * &cur_vec - beta_pre * &pre_vec - alpha * &cur_vec) / beta;
        pre_vec = cur_vec.clone();
        cur_vec = tmp;

        res = res.resize(i + 2, i + 2, 0.0);
        res[(i + 2) * (i + 1) - 1] = beta;
        res[(i + 2) * (i + 2) - 2] = beta;
    }
    alpha = cur_vec.dot(&(mat * &cur_vec));
    res[(times + 1) * (times + 1) - 1] = alpha;
    out.push_str("\n");
    out.push_str(&format!("{} ", times));
    out.push_str(&qr::eigenvalues(&res, 1.0e-9));

    let mut writer = BufWriter::new(File::create(filename)?);
    writer.write_all(out.as_bytes())?;
    writer.flush()?;

    Ok(())
}

#[allow(dead_code)]
pub fn eigenvalues_csr(mat: &CsrMatrix<f64>, times: usize, filename: &str) -> Result<(), Box<dyn std::error::Error>> {

    let nsize = mat.ncols();
    let mut rng = rand::thread_rng();
    let rand_vec: Vec<f64> = (0..nsize).map(|_| rng.gen_range(-1.0..1.0)).collect();
    let mut init_vec = DVectorf64::from_vec(rand_vec);
    let norm = &init_vec.norm();
    init_vec = init_vec / *norm;

    let mut alpha = init_vec.dot(&(mat * &init_vec));
    let mut beta = (mat * &init_vec - alpha * &init_vec).norm();
    let mut beta_pre;
    let mut cur_vec = (mat * &init_vec - alpha * &init_vec) / beta;
    let mut pre_vec = init_vec;

    let mut res = DMatrixf64::from_vec(2, 2, vec![alpha, beta, beta, 0.0]);
    let mut out = String::from(&format!("0 {}", alpha));

    for i in 1..times {
        alpha = cur_vec.dot(&(mat * &cur_vec));
        res[(i + 1) * (i + 1) - 1] = alpha;
        out.push_str("\n");
        out.push_str(&format!("{} ", i));
        out.push_str(&qr::eigenvalues(&res, 1.0e-9));

        beta_pre = beta;
        beta = (mat * &cur_vec - beta_pre * &pre_vec  - alpha * &cur_vec).norm();
        let tmp = (mat * &cur_vec - beta_pre * &pre_vec - alpha * &cur_vec) / beta;
        pre_vec = cur_vec.clone();
        cur_vec = tmp;

        res = res.resize(i + 2, i + 2, 0.0);
        res[(i + 2) * (i + 1) - 1] = beta;
        res[(i + 2) * (i + 2) - 2] = beta;
    }
    alpha = cur_vec.dot(&(mat * &cur_vec));
    res[(times + 1) * (times + 1) - 1] = alpha;
    out.push_str("\n");
    out.push_str(&format!("{} ", times));
    out.push_str(&qr::eigenvalues(&res, 1.0e-9));

    let mut writer = BufWriter::new(File::create(filename)?);
    writer.write_all(out.as_bytes())?;
    writer.flush()?;

    Ok(())
}