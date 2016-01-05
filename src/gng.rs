use rand::random;
use std::cell::Cell;
use std::rc::Rc;

use vec_map::VecMap;

#[derive(Debug)]
struct Unit {
    v: Vec<f32>,
    e: f32,
    a: f32,
    n: VecMap<Rc<Cell<usize>>>,
    pos: (f32, f32, f32)
}

impl Unit {
    fn new(v: Vec<f32>) -> Unit {
        Unit {
            v: v,
            e: 0.0,
            a: 0.0,
            n: VecMap::new(),
            pos: (random::<f32>()-0.5, random::<f32>()-0.5, random::<f32>()-0.5),
        }
    }

    fn distance2(&self, vals: &[f32]) -> f32 {
        self.v.iter().zip(vals.iter()).map(|(x, y)| { let z = *x-*y; z*z }).fold(0., |x, y| x+y) / self.v.len() as f32
    }

    fn distance1(&self, vals: &[f32]) -> f32 {
        self.v.iter().zip(vals.iter()).map(|(x, y)| { let z = *x-*y; z.abs() }).fold(0., |x, y| x+y)
    }

    fn distancee(&self, vals: &[f32]) -> f32 {
        distancee(&self.v, vals)
    }

    fn drag(&mut self, vals: &[f32], eps: f32) {
        let dist = (-self.distance2(vals)).exp();
         for (w, d) in self.v.iter_mut().zip(vals.iter()) {
            *w += eps * (*d - *w) * dist;
        }
    }

    fn age_edges(&mut self) {
        for a in self.n.values() {
            a.set(a.get()+1)
        }
    }
}

pub struct GNGParams {
    pub w_size: usize,
    pub split_error_decrease_rate: f32,
    pub error_decrease_rate: f32,
    pub creation_period: usize,
    pub drag_rate_direct: f32,
    pub drag_rate_indirect: f32,
    pub max_age: usize,
    pub error_thresold: f32,
    pub activity_factor: f32,
}

pub struct GNG {
    params: GNGParams,
    next_id: usize,
    units: VecMap<Unit>
}

fn distancee(x: &[f32], y: &[f32]) -> f32 {
    let mut numer = 0.0;
    let mut denom = 0.0;
    for (u, v) in x.iter().zip(y.iter()) {
        let mut z = (*u - *v);
        z = z*z;
        let ez = z.exp();
        numer += z*ez;
        denom += ez;
    }
    numer/denom
}

impl GNG {
    pub fn new(params: GNGParams) -> GNG {
        let mut h = VecMap::new();
        h.insert(0, Unit::new((0..params.w_size).map(|_| random()).collect()));
        h.insert(1, Unit::new((0..params.w_size).map(|_| random()).collect()));
        let mut gng = GNG {
            next_id: 2,
            units: h,
            params: params,
        };
        gng.connect_reset(0,1);
        gng
    }

    fn connect_reset(&mut self, i: usize, j: usize) {
        let rc = Rc::new(Cell::new(0));
        self.units.get_mut(&i).map(|u| u.n.insert(j, rc.clone()));
        self.units.get_mut(&j).map(|u| u.n.insert(i, rc));
    }

    fn unconnect(&mut self, i: usize, j:usize) {
        self.units.get_mut(&i).map(|u| u.n.remove(&j));
        self.units.get_mut(&j).map(|u| u.n.remove(&i));
    }

    fn add_node(&mut self, values: Vec<f32>, pos: (f32, f32, f32)) -> usize {
        self.units.insert(self.next_id, { let mut u = Unit::new(values); u.pos = pos; u });
        let id = self.next_id;
        self.next_id += 1;
        id
    }

    fn split_edge(&mut self, i: usize, j: usize) {
        let (new_values, new_pos) = if let (Some(unit_i), Some(unit_j)) = (self.units.get(&i), self.units.get(&j)) {
            (
                unit_i.v.iter().zip(unit_j.v.iter()).map(|(x, y)| 0.5*(*x + *y)).collect(),
                ((unit_i.pos.0 + unit_j.pos.0)/2., (unit_i.pos.1 + unit_j.pos.1)/2., (unit_i.pos.2 + unit_j.pos.2)/2.)
            )
        } else { return };
        let k = self.add_node(new_values, new_pos);
        self.unconnect(i, j);
        self.connect_reset(i, k);
        self.connect_reset(j, k);
        let mut err = 0.;
        let alpha = self.params.split_error_decrease_rate;
        self.units.get_mut(&i).map(|u| {u.e *= alpha; err = u.e });
        self.units.get_mut(&j).map(|u| u.e *= alpha);
        self.units.get_mut(&k).map(|u| u.e = err);
    }

    pub fn handle_point(&mut self, datapoint: &[f32]) -> usize {
        // find nearest neighbor
        let (nearest, min_dist, sec_nearest, _sec_min_dist) = self.units.iter().map(|(i, u)| (i, u.distance2(datapoint))).fold(
            (0, ::std::f32::INFINITY, 1, ::std::f32::INFINITY),
            |(c1, c1m, c2, c2m), (n, nm)|
                if      nm < c1m { ( n,  nm, c1, c1m) }
                else if nm < c2m { (c1, c1m,  n,  nm) }
                else             { (c1, c1m, c2, c2m) }
        );
        // drag points
        let eps = self.params.drag_rate_direct;
        let mut neighbors = Vec::new();
        if let Some(w) = self.units.get_mut(&nearest).map(|u| { u.drag(datapoint, eps); neighbors.extend(u.n.keys()); u.e += min_dist; u.a += 1.; u.age_edges(); u.v.clone() } ) {
            let eps = self.params.drag_rate_indirect;
            for n in neighbors {
                self.units.get_mut(&n).map(|u| { u.drag(datapoint, eps); u.drag(&w, -eps / 2.0) });
            }
        }

        // update connections
        self.connect_reset(nearest, sec_nearest);
        let max_age = self.params.max_age;
        let to_prune = self.units.get_mut(&nearest).unwrap()
                                 .n.iter()
                                 .filter_map(|(k, v)| if v.get() > max_age { Some(k) } else { None })
                                 .collect::<Vec<usize>>();
        for j in to_prune {
            self.unconnect(nearest, j);
        }

        let to_delete = self.units.iter().filter_map(|(i, u)| { if u.n.len() == 0 { Some(i) } else { None } }).collect::<Vec<usize>>();
        for i in to_delete {
            self.units.remove(&i);
        }

        // create new unit if needed
        let (q, mut err) = self.units.iter().fold((0, 0.), |(m, em), (n, u)| if u.e > em { (n, u.e) } else { (m, em) });
        if err > self.params.error_thresold {
            let neigs = self.units.get(&q).map(|u| u.n.keys().collect()).unwrap_or(Vec::new());
            let (mut r, mut max_e) = (q, 0.);
            for n in neigs {
                self.units.get(&n).map(|u| if u.e < max_e { max_e = u.e; r = n} );
            }
            self.split_edge(q, r);
        }

        let d = self.params.error_decrease_rate;
        let r = self.params.activity_factor;
        for (_, v) in self.units.iter_mut() { v.e *= d; v.a *= r; }
        // result
        nearest
    }

    pub fn size(&self) -> usize {
        self.units.len()
    }

    pub fn weigths(&self, i: usize) -> Option<(f32, &[f32])> {
        self.units.get(&i).map(|u| (u.a, &u.v[..]))
    }

    pub fn move_all(&mut self) {
        let mut map = VecMap::new();
        for (i, u) in &self.units {
            let e = map.entry(i).or_insert((0.0, 0.0, 0.0));
            let mut force = {
                let dist = (u.pos.0*u.pos.0 + u.pos.1 * u.pos.1 + u.pos.2 * u.pos.2).sqrt() + 0.01;
                (-0.03*u.pos.0, -0.03*u.pos.1, -0.03*u.pos.2)
            };
            for (n,a) in u.n.iter() {
                let other = match self.units.get(&n) {
                    Some(u) => u.pos,
                    None => continue
                };
                let vec = (other.0 - u.pos.0, other.1 - u.pos.1, other.2 - u.pos.2);
                let f = (-(a.get() as f32)/50.0).exp();
                force.0 += f*0.005*vec.0;
                force.1 += f*0.005*vec.1;
                force.2 += f*0.005*vec.2;
            }
            for n in self.units.values() {
                let vec = (n.pos.0 - u.pos.0, n.pos.1 - u.pos.1, n.pos.2 - u.pos.2);
                let mut dist = (vec.0*vec.0 + vec.1 * vec.1 + vec.2 * vec.2).sqrt();
                if dist == 0.0 {
                    // explode !
                    force.0 += 0.001*(random::<f32>()-0.5);
                    force.1 += 0.001*(random::<f32>()-0.5);
                    force.2 += 0.001*(random::<f32>()-0.5);
                    dist = 0.01;
                };
                force.0 -= 0.001*vec.0/dist/dist;
                force.1 -= 0.001*vec.1/dist/dist;
                force.2 -= 0.001*vec.2/dist/dist;
            }
            *e = force;
        }
        for (i, u) in &mut self.units {
            let e = map.entry(i).or_insert((0.0, 0.0, 0.0));
            u.pos.0 += e.0;
            u.pos.1 += e.1;
            u.pos.2 += e.2;
        }
    }

    pub fn positions(&self) -> Vec<(f32, (f32, f32, f32))> {
        self.units.values().map(|u| (u.a, u.pos)).collect()
    }

    pub fn edges(&self) -> Vec<(f32, (f32, f32, f32), (f32, f32, f32))> {
        self.units.values().flat_map(|u| {
            let p = u.pos;
            u.n.iter().filter_map(move |(i, a)| self.units.get(&i).map(|v| (a.get() as f32, p, v.pos)))
        }).collect()
    }
}