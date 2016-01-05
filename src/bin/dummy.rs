extern crate byteorder;
extern crate rand;

use byteorder::{WriteBytesExt, NativeEndian};

use rand::thread_rng;
use rand::distributions::Sample;
use rand::distributions::normal::Normal;

use std::io::Write;

fn main() {
    let stdout = ::std::io::stdout();
    let mut outlock = stdout.lock();
    let mut stderr = ::std::io::stderr();

    outlock.write_u32::<NativeEndian>(2).unwrap();

    let mut norm = Normal::new(0.0, 3.0);
    let mut rng = thread_rng();

    let mut derive = 1.0;
    //let mut delta = 0.00001;

    loop {
        let mut f = norm.sample(&mut rng) + derive;
        if rand::random::<bool>() { f = -f };
        f = f.tanh();
        outlock.write_f32::<NativeEndian>(f as f32).unwrap();
        /*if derive >= 10.0 {
            delta = -0.00001;
        } else if derive <= 0.0 {
            delta = 0.00001;
        }
        derive += delta;
        writeln!(stderr, "{}", derive);*/
    }
}
