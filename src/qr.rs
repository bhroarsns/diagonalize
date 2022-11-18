use na::{Dynamic, OMatrix, DVector};
use std::fs::File;
use std::io::{BufWriter, Write};

type DMatrixf64 = OMatrix<f64, Dynamic, Dynamic>;

pub fn eigenvalues_record(mat: &DMatrixf64, times: usize, dataname: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut cur_mat = mat.clone();
    let mut writer = BufWriter::new(File::create(dataname)?);
    let mut out = String::new();
    out.push_str(&format!("0 {}", cur_mat.diagonal().into_iter().map(|v| format!("{}", v)).collect::<Vec<String>>().join(" ")));

    for i in 0..times {
        let qr = cur_mat.clone().qr();
        cur_mat = qr.r() * qr.q();
        out.push_str(&format!("\n{} {}", i + 1, cur_mat.diagonal().into_iter().map(|v| format!("{}", v)).collect::<Vec<String>>().join(" ")));
    }

    writer.write_all(out.as_bytes())?;
    writer.flush()?;

    Ok(())
}

pub fn eigenvalues(mat: &DMatrixf64, epsilon: f64) -> String {

    let nsize = mat.ncols();

    let mut pre_diag: DVector<f64> = DVector::repeat(nsize, 0.0);
    let mut cur_mat = mat.clone();

    let mut rem = epsilon + 1.0;

    let mut counter = 0;
    while rem > epsilon {
        counter += 1;
        let qr = cur_mat.clone().qr();
        cur_mat = qr.r() * qr.q();
        rem = (&pre_diag[0] - &cur_mat.diagonal()[0]).abs();
        pre_diag = cur_mat.diagonal()
    }

    println!("repeat time: {}", counter);

    cur_mat.diagonal().into_iter().map(|v| format!("{}", v)).collect::<Vec<String>>().join(" ")
}