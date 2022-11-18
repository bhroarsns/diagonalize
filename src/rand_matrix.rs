use na::{Dynamic, OMatrix, DVector};
use rand::Rng;
use nalgebra_sparse::{coo::CooMatrix, csr::CsrMatrix};

type DMatrixf64 = OMatrix<f64, Dynamic, Dynamic>;

#[allow(dead_code)]
pub fn get_rand_with_eigens(data: Vec<f64>) -> DMatrixf64 {
    let diag = DVector::from_vec(data);
    let diag = DMatrixf64::from_diagonal(&diag);

    let size = diag.ncols();

    let mut prop: Vec<f64> = Vec::new();
    let mut rng = rand::thread_rng();
    for _ in 0..size {
        for _ in 0..size {
            prop.push(rng.gen_range(-1.0..1.0));
        }
    }

    let p = DMatrixf64::from_vec(size, size, prop);
    let lu = p.clone().lu();
    let pinv = lu.try_inverse().unwrap();
    println!("matrix created");

    pinv * diag * p
}

#[allow(dead_code)]
pub fn get_sparse_and_sym(size: usize, param: usize) -> CsrMatrix<f64> {

    let mut res = CooMatrix::new(size, size);
    let mut rng = rand::thread_rng();

    for i in 0..size {
        for j in 0..size {
            if i < j {
                if rng.gen_range(0..size) < param {
                    let val = rng.gen_range(-1.0..1.0);
                    res.push(i, j, val);
                    res.push(j, i, val);
                }
            } else if i == j {
                if rng.gen_range(0..size) < param {
                    let val = rng.gen_range(-1.0..1.0);
                    res.push(i, j, val);
                } 
            }
        }
    }

    println!("matrix created");
    CsrMatrix::from(&res)
}

#[allow(dead_code)]
pub fn get_sparse_and_non_sym(size: usize, param: usize) -> CsrMatrix<f64> {

    let mut res = CooMatrix::new(size, size);
    let mut rng = rand::thread_rng();

    for i in 0..size {
        for j in 0..size {
            if rng.gen_range(0..size) < param {
                let val = rng.gen_range(-1.0..1.0);
                res.push(i, j, val);
            }
        }
    }

    println!("matrix created");
    CsrMatrix::from(&res)
}

#[allow(dead_code)]
pub fn that_matrix(size: usize) -> CsrMatrix<f64> {
    let mut res = CooMatrix::new(size, size);

    for i in 0..(size - 1) {
        res.push(i, i, 2.0);
        res.push(i, i + 1, -1.0);
        res.push(i + 1, i, -1.0);
    }

    res.push(size - 1, size - 1, 1.0);

    CsrMatrix::from(&res)
}