//! Verify the ff19SB 3C H-C-H valence-angle patch:
//!   HC-3C-HC, H1-3C-H1, HC-3C-H1
//! These are the CT-equivalent angles for the 3C atom type (ff19SB's "sp3
//! aliphatic C with three heavy atoms"), used by rebuilt sidechains whose
//! carbons get typed 3C. Without them, topology assignment fails with
//! "Missing valence angle params for HC-3C-HC".
//!
//! Run from the dynamics crate root:
//!   cargo run --example check_3c_angles --release

use std::path::Path;

use bio_files::md_params::ForceFieldParams;
use dynamics::merge_params;

fn main() {
    // Mirror the engine's real merged protein FF table (FfParamSet::new_amber:
    // merge_params(&parm19, &peptide_frcmod)).
    let parm19 = ForceFieldParams::load_dat(Path::new("param_data/parm19.dat"))
        .expect("load parm19.dat");
    let frcmod = ForceFieldParams::load_frcmod(Path::new("param_data/frcmod.ff19SB"))
        .expect("load frcmod.ff19SB");
    let merged = merge_params(&parm19, &frcmod);

    let checks: &[(&str, &str, &str)] = &[
        // patched entries (was missing before the fix):
        ("HC", "3C", "HC"),
        ("H1", "3C", "H1"),
        ("HC", "3C", "H1"),
        // pre-existing controls (from parm19.dat):
        ("HC", "CT", "HC"),
        ("H1", "CT", "H1"),
        ("HC", "2C", "HC"),
        ("HC", "3C", "2C"),
        ("XC", "3C", "HC"),
    ];

    let mut all_ok = true;
    for &(a, b, c) in checks {
        let hit = merged.get_valence_angle(
            &(a.to_string(), b.to_string(), c.to_string()),
            true,
        );
        let status = if hit.is_some() { "FOUND" } else { "MISSING" };
        println!("{a}-{b}-{c}: {status}");
        if hit.is_none() {
            all_ok = false;
        }
    }

    if all_ok {
        println!("\nALL 3C H-C-H ANGLES PRESENT ✔");
    } else {
        println!("\nMISSING ANGLES DETECTED ✘");
        std::process::exit(1);
    }
}
