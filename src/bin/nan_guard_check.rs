//! Standalone check for the non-finite accel guard in `kick_and_calc_accel`.
//!
//! Runs via `cargo run --release --bin nan_guard_check`. The crate's
//! `src/tests/` modules are stale (old API), so this binary is the clean way to
//! exercise the guard against the CURRENT API end-to-end:
//!   1. build a tiny 3-atom system (bonded pair + isolated atom),
//!   2. one healthy step -> energy must be finite,
//!   3. poison atom 0's coordinate to NaN -> step must not panic, must zero the
//!      atom's velocity (not NaN it), and must flag potential_energy = NaN so
//!      callers report `crashed` on this step.

use bio_files::{AtomGeneric, BondGeneric, BondType};
use dynamics::{
    ComputationDevice, FfMolType, HydrogenConstraint, MdConfig, MdOverrides, MdState, MolDynamics,
    SimBoxInit, Solvent,
    params::FfParamSet,
};
use lin_alg::f32::Vec3;
use na_seq::Element;

fn main() {
    let param_set = FfParamSet::new_amber().expect("amber params");
    let dev = ComputationDevice::Cpu;

    let cfg = MdConfig {
        sim_box: SimBoxInit::Fixed((Vec3::new(0., 0., 0.), Vec3::new(60., 60., 60.))),
        solvent: Solvent::None,
        barostat_cfg: None,
        hydrogen_constraint: HydrogenConstraint::Flexible,
        recenter_sim_box: false,
        max_init_relaxation_iters: None,
        overrides: MdOverrides {
            skip_water_relaxation: true,
            skip_counterion_insertion: true,
            ..Default::default()
        },
        ..Default::default()
    };

    let c = 30.0_f32;
    let atom_0 = AtomGeneric {
        serial_number: 1,
        posit: Vec3::new(c - 1.5, c, c).into(),
        force_field_type: Some("ca".to_string()),
        element: Element::Carbon,
        // partial_charge is in elementary charge here; 1.0 amber-unit charge
        // corresponds to 1/18.2223 ≈ 0.0549 e.
        partial_charge: Some(0.0549),
        ..Default::default()
    };
    let atom_1 = AtomGeneric {
        serial_number: 2,
        posit: Vec3::new(c + 1.5, c, c).into(),
        force_field_type: Some("ca".to_string()),
        element: Element::Carbon,
        partial_charge: Some(0.0549),
        ..Default::default()
    };
    let mol_a = MolDynamics {
        ff_mol_type: FfMolType::SmallOrganic,
        atoms: vec![atom_0, atom_1],
        bonds: vec![BondGeneric {
            atom_0_sn: 1,
            atom_1_sn: 2,
            bond_type: BondType::Aromatic,
        }],
        ..Default::default()
    };

    let mol_b = MolDynamics {
        ff_mol_type: FfMolType::SmallOrganic,
        atoms: vec![AtomGeneric {
            serial_number: 3,
            posit: Vec3::new(c, c - 8.0, c).into(),
            force_field_type: Some("ca".to_string()),
            element: Element::Carbon,
            partial_charge: Some(0.0549),
            ..Default::default()
        }],
        ..Default::default()
    };

    let (mut md, _explicit_solvent) =
        MdState::new(&dev, &cfg, &[mol_a, mol_b], &param_set).expect("MdState::new");

    // 1) Healthy step: populates nb_pairs, energy finite.
    md.step(&dev, 0.001, None);
    assert!(
        md.potential_energy.is_finite(),
        "healthy step produced non-finite energy: {:?}",
        md.potential_energy
    );
    println!("healthy step OK, u = {}", md.potential_energy);

    // 2) Poison atom 0's coordinate -> force recomputation yields NaN.
    md.atoms[0].posit = Vec3::new(f32::NAN, 0.0, 0.0);
    md.step(&dev, 0.001, None); // must NOT panic

    assert!(
        md.potential_energy.is_nan(),
        "guard should flag system as crashed (potential_energy = NaN), got {:?}",
        md.potential_energy
    );
    let v = &md.atoms[0].vel;
    assert!(
        v.x.is_finite() && v.y.is_finite() && v.z.is_finite(),
        "poisoned atom velocity must be finite (zeroed), got {v:?}"
    );
    println!("NAN_GUARD_OK: non-finite accel -> zeroed vel + potential_energy=NaN (crashed)");
}
