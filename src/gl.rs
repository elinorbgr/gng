use gng::GNG;

use glium::{Surface, VertexBuffer, Program};

const VERTEX_SHADER_SRC : &'static str = r#"
    #version 140

    in vec3 position;

    uniform float u;
    uniform float v;
    uniform float w;
    uniform float r;
    uniform mat4 perspective;
    uniform mat4 view;

    void main() {
        vec3 pos = r * position;
        pos.x += u;
        pos.y += v;
        pos.z += w;
        gl_Position = perspective * view * vec4(pos, 1.0);
    }
"#;

const FRAGMENT_SHADER_SRC : &'static str = r#"
    #version 140

    out vec4 color;

    uniform float t;
    uniform float s;

    void main() {
        color = vec4(s*1.0*(t*0.9 + 0.1) + t*(1.0-s), t*(1.0-s), t*(1.0-s), 0.9*s + 0.1);
    }
"#;

#[derive(Copy, Clone)]
struct Vertex {
    position: [f32; 3],
}

implement_vertex!(Vertex, position);

pub struct Drawer {
    pub display: ::glium::backend::glutin_backend::GlutinFacade,
    vbuffer: VertexBuffer<Vertex>,
    indices: ::glium::index::NoIndices,
    program: Program,
    position: [f32; 3],
    direction: [f32; 3],
    up: [f32; 3]
}

impl Drawer {
    pub fn new() -> Drawer {
        let display = {
            use glium::DisplayBuild;
            ::glium::glutin::WindowBuilder::new()
                                .with_depth_buffer(24)
                                .build_glium().unwrap()
        };

        let vbuffer = {
            let shape = vec![
                //
                Vertex { position: [1.0, 0.0, 0.0] },
                Vertex { position: [0.0, 1.0, 0.0] },
                Vertex { position: [0.0, 0.0, -1.0] },
                //
                Vertex { position: [1.0, 0.0, 0.0] },
                Vertex { position: [0.0, -1.0, 0.0] },
                Vertex { position: [0.0, 0.0, -1.0] },
                //
                Vertex { position: [1.0, 0.0, 0.0] },
                Vertex { position: [0.0, 1.0, 0.0] },
                Vertex { position: [0.0, 0.0, 1.0] },
                //
                Vertex { position: [1.0, 0.0, 0.0] },
                Vertex { position: [0.0, -1.0, 0.0] },
                Vertex { position: [0.0, 0.0, 1.0] },
                //
                Vertex { position: [-1.0, 0.0, 0.0] },
                Vertex { position: [0.0, 1.0, 0.0] },
                Vertex { position: [0.0, 0.0, -1.0] },
                //
                Vertex { position: [-1.0, 0.0, 0.0] },
                Vertex { position: [0.0, -1.0, 0.0] },
                Vertex { position: [0.0, 0.0, -1.0] },
                //
                Vertex { position: [-1.0, 0.0, 0.0] },
                Vertex { position: [0.0, 1.0, 0.0] },
                Vertex { position: [0.0, 0.0, 1.0] },
                //
                Vertex { position: [-1.0, 0.0, 0.0] },
                Vertex { position: [0.0, -1.0, 0.0] },
                Vertex { position: [0.0, 0.0, 1.0] },
            ];
            VertexBuffer::new(&display, &shape).unwrap()
        };

        let indices = ::glium::index::NoIndices(::glium::index::PrimitiveType::TrianglesList);

        let program = Program::from_source(&display, VERTEX_SHADER_SRC, FRAGMENT_SHADER_SRC, None).unwrap();

        Drawer {
            display: display,
            vbuffer: vbuffer,
            indices: indices,
            program: program,
            position: [0., 0., -2.],
            direction: [0., 0., 1.],
            up: [0., 1., 0.]
        }
    }

    pub fn zoom(&mut self, ampl: f32) {
        self.position[0] *= ampl.exp();
        self.position[1] *= ampl.exp();
        self.position[2] *= ampl.exp();
    }

    pub fn rot_hori(&mut self, ampl: f32) {
        self.position = rotate(self.position, self.up, ampl);
        // renormalize
        let ampl = self.position.iter().cloned().fold(0f32, |s, x| s + x*x).sqrt();
        for (d, p) in self.direction.iter_mut().zip(self.position.iter_mut()) {
            *d = - (*p) / ampl;
        }
    }

    pub fn rot_vert(&mut self, ampl: f32) {

    }

    pub fn redraw(&self, gng: &GNG) {
        let mut target = self.display.draw();

        // println!("{:?} - {:?} - {:?}", self.position, self.direction, self.up);

        let perspective = {
            let (width, height) = target.get_dimensions();
            let aspect_ratio = height as f32 / width as f32;

            let fov: f32 = 3.141592 / 3.0;
            let zfar = 1024.0;
            let znear = 0.1;

            let f = 1.0 / (fov / 2.0).tan();

            [
                [f *   aspect_ratio   ,    0.0,              0.0              ,   0.0],
                [         0.0         ,     f ,              0.0              ,   0.0],
                [         0.0         ,    0.0,  (zfar+znear)/(zfar-znear)    ,   1.0],
                [         0.0         ,    0.0, -(2.0*zfar*znear)/(zfar-znear),   0.0],
            ]
        };

        target.clear_color_and_depth((0.0, 0.0, 0.0, 1.0), 1.0);

        let params = ::glium::DrawParameters {
            depth: ::glium::Depth {
                test: ::glium::draw_parameters::DepthTest::IfLess,
                write: true,
                .. Default::default()
            },
            blend: ::glium::Blend::alpha_blending(),
            .. Default::default()
        };

        let view = view_matrix(&self.position, &self.direction, &self.up);

        let maxc = 1500.0;

        for (a, (u, v, w)) in gng.positions() {
            let t = if a > maxc { 1.0 } else { a / maxc };
            target.draw(&self.vbuffer, &self.indices, &self.program, &uniform!{ perspective: perspective, u: u, v: v, w: w, r: 0.02f32, t: t, s: 1.0f32, view: view }, &params).unwrap();
        }

        for (a, (x,y,z),(u,v,w)) in gng.edges() {
            let vbuffer = VertexBuffer::new(&self.display, &[Vertex { position: [x, y, z] }, Vertex { position: [u, v, w] }]).unwrap();
            let indices = ::glium::index::NoIndices(::glium::index::PrimitiveType::LinesList);
            target.draw(&vbuffer, &indices, &self.program, &uniform!{ perspective: perspective, u: 0.0f32, v: 0.0f32, w: 0.0f32, r: 1.0f32, t: (-a/50.0).exp(), s:0.0f32, view: view}, &params).unwrap();
        }

        target.finish().unwrap();
    }
}

fn view_matrix(position: &[f32; 3], direction: &[f32; 3], up: &[f32; 3]) -> [[f32; 4]; 4] {
    let f = {
        let f = direction;
        let len = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
        let len = len.sqrt();
        [f[0] / len, f[1] / len, f[2] / len]
    };

    let s = [up[1] * f[2] - up[2] * f[1],
             up[2] * f[0] - up[0] * f[2],
             up[0] * f[1] - up[1] * f[0]];

    let s_norm = {
        let len = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
        let len = len.sqrt();
        [s[0] / len, s[1] / len, s[2] / len]
    };

    let u = [f[1] * s_norm[2] - f[2] * s_norm[1],
             f[2] * s_norm[0] - f[0] * s_norm[2],
             f[0] * s_norm[1] - f[1] * s_norm[0]];

    let p = [-position[0] * s_norm[0] - position[1] * s_norm[1] - position[2] * s_norm[2],
             -position[0] * u[0] - position[1] * u[1] - position[2] * u[2],
             -position[0] * f[0] - position[1] * f[1] - position[2] * f[2]];

    [
        [s[0], u[0], f[0], 0.0],
        [s[1], u[1], f[1], 0.0],
        [s[2], u[2], f[2], 0.0],
        [p[0], p[1], p[2], 1.0],
    ]
}

fn rotate(vec: [f32; 3], axis: [f32; 3], angle: f32) -> [f32; 3] {
    let c = angle.cos();
    let s = angle.sin();
    let (ux, uy, uz) = (axis[0], axis[1], axis[2]);
    let (xx, yy, zz) = (vec[0], vec[1], vec[2]);
    [
        (c + ux*ux*(1. - c))  * xx + (ux*uy*(1.-c) - uz*s) * yy + (ux*uz*(1.-c) + uy*s) * zz,
        (uy*ux*(1.-c) + uz*s) * xx + (c + uy*uy*(1.-c))    * yy + (uy*uz*(1.-c) - ux*s) * zz,
        (uz*ux*(1.-c) - uy*s) * xx + (uz*uy*(1.-c) + ux*s) * yy + (c + uz*uz*(1.-c))    * zz
    ]
}