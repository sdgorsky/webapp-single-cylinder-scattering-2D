use criterion::{criterion_group, criterion_main, Criterion};

use scattering_core::field::{compute_field, FieldParams};
use scattering_core::performance_testing::{
    create_coord_array, create_coord_array_boxed, create_coord_vec, exterior_field_naive,
    exterior_field_optimized, exterior_field_spline, exterior_field_spline_symmetric, pw_vec_1,
    pw_vec_2, pw_vec_3,
};
use scattering_core::scattering::{calculate_scattering, Material, Polarization, ScatteringParams};
fn test_field_calculation() {
    let params = ScatteringParams {
        wavelength: 1.0,
        material: Material {
            permittivity_real: 4.0,
            permittivity_imag: 0.0,
            permeability_real: 1.0,
            permeability_imag: 0.0,
        },
        polarization: Polarization::TM,
        max_order: 10,
    };

    let scattering = calculate_scattering(&params);

    // Build field params
    let field_params = FieldParams {
        wavelength: params.wavelength,
        permittivity_real: params.material.permittivity_real,
        permittivity_imag: params.material.permittivity_imag,
        permeability_real: params.material.permeability_real,
        permeability_imag: params.material.permeability_imag,
        incident_coeffs_real: scattering
            .incident_coefficients
            .iter()
            .map(|c| c.re)
            .collect(),
        incident_coeffs_imag: scattering
            .incident_coefficients
            .iter()
            .map(|c| c.im)
            .collect(),
        scattering_coeffs_real: scattering
            .scattering_coefficients
            .iter()
            .map(|c| c.re)
            .collect(),
        scattering_coeffs_imag: scattering
            .scattering_coefficients
            .iter()
            .map(|c| c.im)
            .collect(),
        internal_coeffs_real: scattering
            .internal_coefficients
            .iter()
            .map(|c| c.re)
            .collect(),
        internal_coeffs_imag: scattering
            .internal_coefficients
            .iter()
            .map(|c| c.im)
            .collect(),
        orders: scattering.orders,
    };

    let _result = compute_field(&field_params);
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Field Benchmark Group");

    group.sample_size(100);

    //group.bench_function("field calc", |b| b.iter(|| test_field_calculation()));
    //group.bench_function("coord_vec", |b| b.iter(|| create_coord_vec()));
    //group.bench_function("coord_array", |b| b.iter(|| create_coord_array()));
    //group.bench_function("coord_array_boxed", |b| b.iter(|| create_coord_array_boxed()));

    group.bench_function("pw_vec_1", |b| b.iter(|| pw_vec_1()));
    group.bench_function("pw_vec_2", |b| b.iter(|| pw_vec_2()));
    group.bench_function("pw_vec_3", |b| b.iter(|| pw_vec_3()));

    // Exterior field benchmarks need fewer samples due to long runtime
    group.sample_size(10);
    group.bench_function("exterior_field_naive", |b| {
        b.iter(|| exterior_field_naive())
    });
    group.bench_function("exterior_field_optimized", |b| {
        b.iter(|| exterior_field_optimized())
    });
    group.bench_function("exterior_field_spline", |b| {
        b.iter(|| exterior_field_spline())
    });
    group.bench_function("exterior_field_spline_sym", |b| {
        b.iter(|| exterior_field_spline_symmetric())
    });

    // group.bench_function("pw_array", |b| b.iter(|| pw_array()));
    // group.bench_function("pw_array_boxed", |b| b.iter(|| pw_array_boxed()));

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
