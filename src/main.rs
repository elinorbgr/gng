extern crate byteorder;
extern crate rand;
#[macro_use]
extern crate glium;
extern crate vec_map;

use byteorder::{ReadBytesExt, NativeEndian};
use std::io::Read;


mod gng;
mod gl;

const DISPLAY_PERIOD: usize = 1000;
const ERR_T: f32 = 0.045;

fn main() {
    let mut stdin = ::std::io::stdin();

    let feature_count = match stdin.read_u32::<NativeEndian>() {
        Ok(u) => u as usize,
        Err(e) => panic!("I/O Error: {:?}", e)
    };

    println!("feature count: {}", feature_count);

    let mut gng = gng::GNG::new(
        gng::GNGParams {
            w_size: feature_count,
            split_error_decrease_rate: 0.3,
            error_decrease_rate: 0.8,
            creation_period: 5000,
            drag_rate_direct: 0.15,
            drag_rate_indirect: 0.05,
            max_age: 50,
            error_thresold: ERR_T,
            activity_factor: 0.99999,
        }
    );

    let mut drawer = gl::Drawer::new();

    let mut show = DISPLAY_PERIOD;
    let mut print = 10;

    let mut v = vec![0f32; feature_count];

    let mut lastpos: (i32, i32) = (0, 0);
    let mut moving = false;

    let input_scale = 0.005;

    'outer : loop {

        let mut offset = (0f32, 0f32, 0f32);

        let mut moved = false;

        for ev in drawer.display.poll_events() {
            match ev {
                glium::glutin::Event::Closed => break 'outer,
                glium::glutin::Event::MouseInput(state, glium::glutin::MouseButton::Left) => {
                    moved = true;
                    moving = match state {
                        glium::glutin::ElementState::Pressed => true,
                        glium::glutin::ElementState::Released => false,
                    };
                },
                glium::glutin::Event::MouseMoved((x, y)) => {
                    if moving {
                        moved = true;
                        offset.0 += (x - lastpos.0) as f32 * input_scale;
                        offset.1 += (y - lastpos.1) as f32 * input_scale;
                    }
                    lastpos = (x, y);
                },
                glium::glutin::Event::MouseWheel(delta) => {
                    moved = true;
                    offset.2 -= match delta {
                        glium::glutin::MouseScrollDelta::LineDelta(_, z) => z * 10.0,
                        glium::glutin::MouseScrollDelta::PixelDelta(_, z) => z as f32
                    } * input_scale * 10.;
                }
                _ => ()
            }
        }
        if moved {
            drawer.rot_hori(offset.0);
            drawer.rot_vert(offset.1);
            drawer.zoom(offset.2);
        }

        for f in &mut v {
            match stdin.read_f32::<NativeEndian>() {
                Ok(n) => { *f = n; },
                Err(e) => panic!("I/O Error: {:?}", e)
            }
        }

        let class = gng.handle_point(&v);
        show -= 1;
        if show <= 0 {
            /*
            print!("{:>3} -- {:>3}", gng.size(), class);
            if let Some((a,w)) = gng.weigths(class) {
                print!(" {:>7.2} [", a);
                for v in w {
                    print!("{:.2}, ", v);
                }
                print!("]");
            }
            println!("");
            */
            // redraw

            gng.move_all();

            print -= 1;
            if print <= 0 {
                println!("{:?}", gng.size());
                print = 10;
            }

            drawer.redraw(&gng);

            show = DISPLAY_PERIOD;
        }
    }

}
