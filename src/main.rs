#[macro_use]
extern crate nalgebra as na;
use na::{MatrixNM, U51};
extern crate colored;
use colored::*;

type StorMat = MatrixNM<f64, U51, U51>;


struct BasicContext {
    xmin : f64,
    xmax : f64,
    ymin : f64,
    ymax : f64,

    nx : u16,  
    ny : u16,  

    cp     : f64,  // heat capacity     
    H      : f64,  // ummmm...? TODO
    dt     : f64,  // timestep 
    refrho : f64,  // reference density
    g      : f64,  // gravity
}

struct Context {
    u   : StorMat,  // vel in x
    un  : StorMat,  // vel in x
    v   : StorMat,  // vel in y
    vn  : StorMat,  // vel in y
    p   : StorMat,  // pressure
    pn  : StorMat,  // pressure
    rho : StorMat,  // density
    nu  : StorMat,  // viscosity
    t   : StorMat,  // temp
    tn  : StorMat,  // temp
    k   : StorMat,  // conductivity
}


fn apply_vel_bc(bc : &BasicContext, c : & mut Context) {
    let mut idx : usize;
    let mut off : usize;

    for j in 0..bc.ny {
        idx = (bc.nx * j) as usize;

        // left wall free-slip
        c.u[idx] = 0.;
        c.v[idx] = c.v[idx + 1];

        // right wall free-slip
        off = idx + (bc.nx - 1) as usize;
        c.u[off] = 0.;
        c.v[off] = c.v[off - 1];
    }

    for i in 0..bc.nx {
        // bottom wall free-slip
        c.u[i as usize] = c.u[(bc.nx + i) as usize];
        c.v[i as usize] = 0.;

        // top wall free-slip
        idx = (bc.nx * (bc.ny - 1) + i) as usize;
        c.u[idx] = c.u[(bc.nx * (bc.ny - 2) + i) as usize];
        c.v[idx] = 0.;
    }
}


fn apply_thermal_bc(bc : &BasicContext, c : & mut Context) {
    let mut idx : usize;
    let mut off : usize;

    for j in 0..bc.ny {
        idx = (bc.nx * j) as usize;

        // left wall, no temp escape
        c.t[idx] = c.t[idx + 1];

        // right wall, no temp escape
        off = idx + (bc.nx - 1) as usize;
        c.t[off] = c.t[off - 1];
    }

    for i in 0..bc.nx {
        // bottom wall 
        c.t[i as usize] = 1000.;

        // top wall
        idx = (bc.nx * (bc.ny - 1) + i) as usize;
        c.t[idx] = 0.;
    }
}

fn solve_advection_diffusion(bc : &BasicContext, c : & mut Context) {
    let mut kx : f64;
    let mut ky : f64;

    let dx = ( bc.xmax - bc.xmin ) / (bc.nx as f64);  // TODO put this in the context
    let dy = ( bc.ymax - bc.ymin ) / (bc.ny as f64);
    let dx2 = dx * dx;
    let dy2 = dy * dy;
    let twodx = 2. * dx;
    let twody = 2. * dy;

    for i in 0..(bc.nx * bc.ny) {
        c.tn[i as usize] = c.t[i as usize];
    }

    let mut m : usize;  // middle
    let mut n : usize;  // north
    let mut e : usize;
    let mut s : usize;
    let mut w : usize;  // west

    for j in 1..bc.ny-1 {
        for i in 1..bc.nx-1 {
            m = (bc.nx * j + i) as usize;
            n = (bc.nx * (j-1) + i) as usize;
            s = (bc.nx * (j+1) + i) as usize;
            e = (bc.nx * j + (i+1)) as usize;
            w = (bc.nx * j + (i-1)) as usize;

            kx = c.k[m] + (c.tn[e] - 2. * c.tn[m] + c.tn[w]) / dx2;
            ky = c.k[m] + (c.tn[s] - 2. * c.tn[m] + c.tn[n]) / dx2;

            c.t[m] = c.tn[m] + bc.dt * ((bc.H + kx + ky) / (c.rho[m] * bc.cp) 
                                        - (c.u[m] * ( (c.tn[e] - c.tn[w]) / twodx ))
                                        - (c.v[m] * ( (c.tn[s] - c.tn[n]) / twody )) 
                                       );
        }
    }

    apply_thermal_bc(bc, c);
}


fn print_field(bc : &BasicContext, field : & StorMat,
               min : f64, max : f64) {

    let mut idx : usize;
    let mut norm : f64;
    let bchar = "██";

    for jj in 0..bc.ny {
        let j = bc.ny - jj - 1;
        for i in 0..bc.nx {
            idx = (bc.nx * j + i) as usize;
            norm = (field[idx] - min) / (max - min);

            match norm {
                0f64       ... 0.1429f64  => print!("{}", bchar.blue()),
                0.1429f64  ... 0.28571f64 => print!("{}", bchar.cyan()),
                0.28571f64 ... 0.42857f64 => print!("{}", bchar.green()),
                0.42857f64 ... 0.57128f64 => print!("{}", bchar.white()),
                0.57128f64 ... 0.71429f64 => print!("{}", bchar.yellow()),
                0.71429f64 ... 0.85714f64 => print!("{}", bchar.red()),
                0.85714f64 ... 1f64       => print!("{}", bchar.magenta()),
                _                         => print!("{}", bchar.black()) // Error - this may an issue in the future
            }
        }
        print!("{}", "\n".normal());
    }

}


fn main() {
    let bc = BasicContext{ xmin   : 0.,
                           xmax   : 2.,
                           ymin   : 0.,
                           ymax   : 2.,
                           nx     : 51,
                           ny     : 51, 
 
                           cp     : 60.,   
                           H      : 0.,    
                           dt     : 9e-2,  
                           refrho : 100.,  
                           g      : 9.81,  
                         };

    let ele = bc.nx * bc.ny;
    let dx = ( bc.xmax - bc.xmin ) / (bc.nx as f64);
    let dy = ( bc.ymax - bc.ymin ) / (bc.ny as f64);

    

    println!("dx = {}\ndy = {}\nele = {}", dx, dy, ele);

    let mut c = Context{ u   : StorMat::from_element(0.),
                         un  : StorMat::from_element(0.),
                         v   : StorMat::from_element(0.),
                         vn  : StorMat::from_element(0.),
                         p   : StorMat::from_element(0.),
                         pn  : StorMat::from_element(0.),
                         rho : StorMat::from_element(bc.refrho),
                         nu  : StorMat::from_element(1.  ),
                         t   : StorMat::from_element(500.),
                         tn  : StorMat::from_element(500.),
                         k   : StorMat::from_element(100.),
                        };

    let mut currenttime = 0.;
    let mut timestep = 0;

    // Initial conditions
    for j in 0..bc.ny {
        for i in 0..bc.nx {
            let idx = (bc.nx * j + i) as usize;
            c.p[idx] = bc.refrho * dy as f64 * bc.g * (bc.ny - j) as f64;  // Approx lithostatic pressure
            // Apply thermal initial state
            if j < 4 && i < bc.nx/2 {
                c.t[idx] = 1000.;
            }
            else if j > bc.ny-5 && i > bc.nx/2 {
                c.t[idx] = 0.;
            }
        }
    }

    apply_vel_bc(&bc, & mut c);
    apply_thermal_bc(&bc, & mut c);

    // Main loop
    loop {

        if timestep > 100000 {
            break
        }

        solve_advection_diffusion(&bc, & mut c);

        timestep += 1;
        currenttime += bc.dt;
    }

    println!("Final time: {}", currenttime);
    print_field(&bc, &c.t, 0., 1000.);
}
