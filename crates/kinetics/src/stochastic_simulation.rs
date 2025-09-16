use std::fmt;
use rand::Rng;
use nohash_hasher::IntMap;
use energy::EnergyModel;
use crate::LoopStructure;

const KB: f64 = 0.001987204285; // kcal/(mol*K)
const K0: f64 = 273.15;

pub trait KineticModel {
    /// Given ΔE (in kcal/mol) and maybe other info, return the rate constant
    fn rate(&self, delta_e: i32) -> f64;
}

#[derive(Debug, Clone, Copy)]
pub struct Metropolis {
    kt: f64, // k_B * T in kcal/mol
    k0: f64,
}

impl Metropolis {
    pub fn new(celsius: f64, k0: f64) -> Self {
        if k0 <= 0. {
            panic!("k0 must be positive!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            kt: KB * t_kelvin,
            k0,
        }
    }
}

impl KineticModel for Metropolis {
    fn rate(&self, delta_e: i32) -> f64 {
        if delta_e <= 0 {
            self.k0
        } else {
            self.k0 * ((-delta_e as f64 / 100.) / self.kt).exp()
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Reaction {
    Add {
        i: usize,
        j: usize,
        delta_e: i32,
        rate: f64,
    },
    Del {
        i: usize,
        j: usize,
        delta_e: i32,
        rate: f64,
    },
}

impl Reaction {
    pub fn new_add<K: KineticModel>(model: &K, 
        i: usize, j: usize, delta_e: i32
) -> Self {
        let rate = model.rate(delta_e);
        Reaction::Add { i, j, delta_e, rate }
    }

    pub fn new_del<K: KineticModel>(model: &K, 
        i: usize, j: usize, delta_e: i32) -> Self {
        let rate = model.rate(delta_e);
        Reaction::Del { i, j, delta_e, rate }
    }

    pub fn ij(&self) -> (usize, usize) {
        match self {
            Reaction::Add { i, j, .. } => (*i, *j),
            Reaction::Del { i, j, .. } => (*i, *j),
        }
    }

    pub fn rate(&self) -> f64 {
        match self {
            Reaction::Add { rate, .. } => *rate,
            Reaction::Del { rate, .. } => *rate,
        }
    }

    pub fn delta_e(&self) -> i32 {
        match self {
            Reaction::Add { delta_e, .. } => *delta_e,
            Reaction::Del { delta_e, .. } => *delta_e,
        }
    }

}

pub struct LoopStructureSSA<'a, M: EnergyModel, K: KineticModel> {
    loopstructure: LoopStructure<'a, M>, // owns the RNA folding state
    ratemodel: &'a K,
    flux: f64,
    pair_flux: f64,

    loop_flux: IntMap<usize, f64>,
    loop_rxns: IntMap<usize, Vec<Reaction>>,
    pair_rxns: IntMap<usize, Reaction>
}

impl<'a, M, K> fmt::Debug for LoopStructureSSA<'a, M, K>
where
    M: EnergyModel,
    K: KineticModel + fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopStructureSSA")
            .field("ratemodel", &self.ratemodel)   // prints Debug for kinetic model
            .field("loopstructure", &format!("{}", self.loopstructure))
            .field("flux", &self.flux)
            //.field("num_reactions", &self.reactions.len())
            .finish()
    }
}

impl<'a, M: EnergyModel, K: KineticModel> From<(LoopStructure<'a, M>, &'a K)>
    for LoopStructureSSA<'a, M, K>
{
    fn from((loopstructure, ratemodel): (LoopStructure<'a, M>, &'a K)) -> Self {
        let mut flux = 0.0;
        let mut loop_flux: IntMap<usize, f64> = IntMap::default();
        let mut loop_rxns: IntMap<usize, Vec<Reaction>> = IntMap::default();

        //let loop_lookup = loopstructure.loop_lookup().clone();

        for (lli, add_neighbors) in loopstructure.loop_neighbors().iter() {
            let mut lflux = 0.;
            let mut lrxns = Vec::new();
            for &(i, j, delta) in add_neighbors {
                let rxn = Reaction::new_add(ratemodel, i, j, delta);
                lflux += rxn.rate();
                //println!("Loop: {} flux: {} ({})", lli, lflux, rxn.delta_e());
                lrxns.push(rxn);
            }
            loop_flux.insert(*lli, lflux);
            loop_rxns.insert(*lli, lrxns);
            flux += lflux;
            //println!("{} {}", flux, lflux);
        }

        let mut pair_rxns: IntMap<usize, Reaction> = IntMap::default();
        let mut pair_flux = 0.0;
        for (i, j, delta) in loopstructure.get_del_neighbors() {
            let rxn = Reaction::new_del(ratemodel, i, j, delta);
            flux += rxn.rate();
            pair_flux += rxn.rate();
            pair_rxns.insert(i, rxn);
            //println!("Flux: {} {}", flux, pair_flux);
        }

        Self {
            ratemodel,
            loopstructure,
            flux,
            pair_flux,
            loop_flux,
            loop_rxns,
            pair_rxns,
        }
    }
}

impl<'a, M: EnergyModel, K: KineticModel> LoopStructureSSA<'a, M, K> {
    pub fn current_structure(&self) -> String {
        format!("{}", self.loopstructure)
    }
   
    pub fn remove_loop_reaction(&mut self, lli: usize) {
        self.flux -= self.loop_flux.remove(&lli).expect("The lflux to be removed.");
        self.loop_rxns.remove(&lli);
    }

    pub fn remove_pair_reaction(&mut self, pli: usize) {
        let rxn = self.pair_rxns.remove(&pli).expect("The reaction to be removed.");
        self.flux -= rxn.rate();
        self.pair_flux -= rxn.rate();
        let (i, j) = rxn.ij();
        let &lli_outer = self.loopstructure.loop_lookup().get(&i).expect("i -> lli outer");
        let &lli_inner = self.loopstructure.loop_lookup().get(&j).expect("j -> lli inner");
        self.remove_loop_reaction(lli_inner);
        self.remove_loop_reaction(lli_outer);
    }

    pub fn insert_loop_reactions(&mut self, lli: usize, 
        add_neighbors: Vec<(usize, usize, i32)>
    ) {
        let mut lflux = 0.;
        let mut lrxns = Vec::new();
        for (i, j, delta) in add_neighbors {
            let rxn = Reaction::new_add(self.ratemodel, i, j, delta);
            lflux += rxn.rate();
            lrxns.push(rxn);
        }
        self.loop_flux.insert(lli, lflux);
        self.loop_rxns.insert(lli, lrxns);
        self.flux += lflux;
    }

    pub fn update_pair_reactions(&mut self, change: Vec<(usize, usize, i32)>) {
        for (i, j, delta) in change {
            if let Some(old) = self.pair_rxns.remove(&i) {
                self.flux -= old.rate();
                self.pair_flux -= old.rate();
            }
            let rxn = Reaction::new_del(self.ratemodel, i, j, delta);
            self.flux += rxn.rate();
            self.pair_flux += rxn.rate();
            self.pair_rxns.insert(i, rxn);
        }
    }

    pub fn simulate<R, F>(
        &mut self,
        rng: &mut R,
        t_max: f64,
        mut callback: F,
    )
    where
        R: Rng + ?Sized,
        F: FnMut(f64, &LoopStructure<'a, M>),
    {
        let mut t = 0.;

        while t < t_max {
            assert!(self.flux > 0.0, "Flux vanished, no reactions possible");

            // sample waiting time ~ Exp(flux)
            t += -rng.random::<f64>().ln() / self.flux;

            // Callback after waiting time, before applying move.
            callback(t, &self.loopstructure);

            // sample reaction, probably the bottleneck for now
            let mut r = rng.random::<f64>() * self.flux;

            let chosen = if r <= self.pair_flux {
                let mut found = None;
                for (pli, rxn) in self.pair_rxns.iter() {
                    r -= rxn.rate();
                    if r <= 0.0 {
                        found = Some((*pli, rxn.clone()));
                        break;
                    }
                }
                found
            } else {
                r -= self.pair_flux;
                let mut found = None;
                'outer: for (lli, flux) in self.loop_flux.iter() {
                    r -= flux;
                    if r <= 0.0 {
                        r += flux;
                        for rxn in self.loop_rxns[lli].iter() {
                            r -= rxn.rate();
                            if r <= 0.0 {
                                found = Some((*lli, rxn.clone()));
                                break 'outer;
                            }
                        }
                    }
                }
                found
            };

            if let Some((idx, rxn)) = chosen {
                match rxn {
                    Reaction::Add { i, j, .. } => {
                        self.remove_loop_reaction(idx);
                        let ((lli, ami), (llj, amj), pair_changes) = self
                            .loopstructure.apply_add_move(i, j);
                        self.insert_loop_reactions(lli, ami);
                        self.insert_loop_reactions(llj, amj);
                        self.update_pair_reactions(pair_changes);
                    },
                    Reaction::Del { i, j, .. } => {
                        self.remove_pair_reaction(idx);
                        let ((lli, neighbors), pair_changes) = self
                            .loopstructure.apply_del_move(i, j);
                        self.insert_loop_reactions(lli, neighbors);
                        self.update_pair_reactions(pair_changes);
                    },
                }
            } else {
                panic!("No reaction chosen despite positive flux");
            }

        }
    }
}

