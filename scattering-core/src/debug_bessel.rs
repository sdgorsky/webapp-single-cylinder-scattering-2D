use complex_bessel_rs::bessel_j::bessel_j;
use num_complex::Complex64;

fn main() {
    let z = Complex64::new(2.5, 0.0);
    
    // Test positive and negative orders
    let j2 = bessel_j(2.0, z).unwrap();
    let j_neg2 = bessel_j(-2.0, z).unwrap();
    let j3 = bessel_j(3.0, z).unwrap();
    let j_neg3 = bessel_j(-3.0, z).unwrap();
    
    println!("J_2(2.5) = {:?}", j2);
    println!("J_-2(2.5) = {:?}", j_neg2);
    println!("Expected J_-2 = J_2 (since (-1)^2 = 1): {:?}", j2);
    println!();
    println!("J_3(2.5) = {:?}", j3);
    println!("J_-3(2.5) = {:?}", j_neg3);
    println!("Expected J_-3 = -J_3 (since (-1)^3 = -1): {:?}", -j3);
    println!();
    
    // Check if the relationship holds
    println!("J_-2 == J_2 ? diff = {}", (j2 - j_neg2).norm());
    println!("J_-3 == -J_3 ? diff = {}", (j3 + j_neg3).norm());
}
