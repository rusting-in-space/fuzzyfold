
use kinetics::ViennaRNA;
use kinetics::energy_tables::Base;
use kinetics::energy_tables::basify;
use kinetics::coaxial_stacking::eval_multibranch_loop;
use kinetics::coaxial_stacking::eval_exterior_loop;

fn main() {
    // turn strings into Vec<Base>
    let seg1 = basify("UUG");
    let seg2 = basify("CUG");
    let seg3 = basify("CG");
    let seg4 = basify("CUG");

    // create a Vec<&[Base]> by borrowing each Vec<Base>
    let binding: Vec<&[Base]> = vec![&seg1, &seg2, &seg3, &seg4];

    // load ViennaRNA parameter file
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/rna_turner2004.par");
    let model = ViennaRNA::from_parameter_file(path).unwrap();

    let mfe = eval_multibranch_loop(&binding, &model);
    println!("Minimum free energy: {}", mfe);

    let mfe = eval_exterior_loop(&binding, &model);
    println!("Minimum free energy: {}", mfe);


    //// construct DP with &[&[Base]]
    //let mut dp = CoaxialStackingDP::from_exterior_loop(&binding, &model);

    //println!("Minimum free energy: {}", dp.compute());
    //println!("{:?} {:?} {:?}", seg1, seg2, seg3);

    //// construct DP with &[&[Base]]
    //let mut dp = CoaxialStackingDP::from_multibranch_loop(&binding, &model);

    //println!("Minimum free energy: {}", dp.compute());
    //println!("{:?} {:?} {:?}", seg1, seg2, seg3);
}


