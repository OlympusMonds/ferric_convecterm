#[macro_use]
extern crate nalgebra as na;
use na::{MatrixNM, U101, U51};
extern crate colored;
use colored::*;

type StorMat = MatrixNM<f64, U101, U51>;


struct BasicContext {
    xmin : f64,
    xmax : f64,
    ymin : f64,
    ymax : f64,

    nx : u16,  
    ny : u16,  

    cp      : f64,  // heat capacity     
    H       : f64,  // ummmm...? TODO
    dt      : f64,  // timestep 
    refrho  : f64,  // reference density
    reftemp : f64,  // reference temp
    refnu   : f64,  // reference visc
    refk    : f64,  // reference conductivity
    a       : f64,  // thermal expansivity
    theta   : f64,  // Frank-Kam theta
    g       : f64,  // gravity
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
    b   : StorMat,  // b
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
            ky = c.k[m] + (c.tn[s] - 2. * c.tn[m] + c.tn[n]) / dy2;

            c.t[m] = c.tn[m] + bc.dt * ((bc.H + kx + ky) / (c.rho[m] * bc.cp) 
                                        - (c.u[m] * ( (c.tn[e] - c.tn[w]) / twodx ))
                                        - (c.v[m] * ( (c.tn[s] - c.tn[n]) / twody )) 
                                       );
        }
    }

    apply_thermal_bc(bc, c);
}


fn solve_pressure_poisson(bc : &BasicContext, c : & mut Context) {
    let mut idx : usize;
    let mut off : usize;

    let dx = ( bc.xmax - bc.xmin ) / (bc.nx as f64);  // TODO put this in the context
    let dy = ( bc.ymax - bc.ymin ) / (bc.ny as f64);
    let dx2 = dx * dx;
    let dy2 = dy * dy;
    let twodx = 2. * dx;
    let twody = 2. * dy;
    let inv_dt = 1. / bc.dt;

    let mut m : usize;  // middle
    let mut n : usize;  // north
    let mut e : usize;
    let mut s : usize;
    let mut w : usize;  // west

    // Pre-solve b term.
    for j in 1..bc.ny-1 {
        for i in 1..bc.nx-1 {
            m = (bc.nx * j + i) as usize;
            n = (bc.nx * (j-1) + i) as usize;
            s = (bc.nx * (j+1) + i) as usize;
            e = (bc.nx * j + (i+1)) as usize;
            w = (bc.nx * j + (i-1)) as usize;

            c.b[m] = (( c.rho[m] * dx2 * dy2 ) / ( 2. * (dx2 + dy2))) *
                (
                    inv_dt *
                    (
                        (( c.u[e] - c.u[w] ) / twodx ) +
                        (( c.v[s] - c.v[n] ) / twody )
                    ) -
                    (( c.u[e] - c.u[w] ) / twodx).powi(2) -
                    2. *
                    ((( c.u[s] - c.u[n] ) / twody ) *
                     (( c.v[e] - c.v[w] ) / twodx )) -
                    (( c.v[s] - c.v[n] ) / twody ).powi(2)
                );

        }
    }

    let mut stepcount = 0_usize;
    let mut diff = 1_f64;
    let mut pd  : f64;
    let mut pnt : f64;
    loop {
        if stepcount >= 1_000_000 {
            panic!("Unable to solve poisson: sc = {}, diff = {}", stepcount, diff);
        }
        if diff < 1e-5 && stepcount > 4 {
            break;
        }
        if stepcount % 10000 == 0 {
            println!("    Pressure diff: {}", diff);
        }

        for i in 0..(bc.nx * bc.ny) {
            c.pn[i as usize] = c.p[i as usize];
        }

        for j in 1..bc.ny-1 {
            for i in 1..bc.nx-1 {
                m = (bc.nx * j + i) as usize;
                n = (bc.nx * (j-1) + i) as usize;
                s = (bc.nx * (j+1) + i) as usize;
                e = (bc.nx * j + (i+1)) as usize;
                w = (bc.nx * j + (i-1)) as usize;

                c.p[m] = (( c.pn[e] + c.pn[w] ) * dy2 + (c.pn[s] + c.pn[n]) * dx2 ) /
                          ( 2. * (dx2 + dy2) ) - c.b[m];
            }
        }

        // Pressure BCs
        for j in 0..bc.ny {
            idx = (bc.nx * j) as usize;
            c.p[idx] = c.p[idx + 1];

            off = idx + (bc.nx - 1) as usize;
            c.p[off] = c.p[off - 1];
        }

        // Check for steady state
        pd = 0_f64;
        pnt = 0_f64;

        for i in 0..bc.nx {
            c.p[i as usize] = c.p[(bc.nx + i) as usize];

            idx = (bc.nx * (bc.ny - 1) + i) as usize;
            c.p[idx] = c.p[(bc.nx * (bc.ny - 2) + i) as usize];
        }

        for i in 0..(bc.nx * bc.ny) {
            c.tn[i as usize] = c.t[i as usize];
            pd += (c.p[i as usize].abs() - c.pn[i as usize].abs()).abs();
            pnt += c.pn[i as usize].abs();
        }

        if pnt > 0. {
            diff = (pd / pnt).abs();
        } else {
            diff = 1_f64;
        }

        stepcount += 1;
    }
}


fn solve_stokes_momentum(bc : &BasicContext, c : & mut Context) {
    let dx = ( bc.xmax - bc.xmin ) / (bc.nx as f64);  // TODO put this in the context
    let dy = ( bc.ymax - bc.ymin ) / (bc.ny as f64);
    let dtodx2 = bc.dt / (dx * dx);
    let dtody2 = bc.dt / (dy * dy);

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

            c.u[m] = c.un[m] - ( bc.dt / (c.rho[m] * 2. * dx) ) * (c.p[e] - c.p[w]) 
                     + c.nu[m] * (
                         (dtodx2 * (c.un[e] - 2.*c.un[m] + c.un[w])) +
                         (dtody2 * (c.un[s] - 2.*c.un[m] + c.un[n])) 
                         );

            c.v[m] = c.vn[m] - ( bc.dt / (c.rho[m] * 2. * dy) ) * (c.p[e] - c.p[w]) 
                     + c.nu[m] * (
                         (dtodx2 * (c.vn[e] - 2.*c.vn[m] + c.vn[w])) +
                         (dtody2 * (c.vn[s] - 2.*c.vn[m] + c.vn[n])) 
                         ) - (bc.g * c.rho[m]);

        }
    }
}


fn solve_flow(bc : &BasicContext, c : & mut Context) {
    let mut stepcount = 0_usize;
    let mut diff = 1_f64;

    let mut ud  : f64;
    let mut vd  : f64;
    let mut unt : f64;
    let mut vnt : f64;

    loop {
        if stepcount >= 50000 {
            panic!("Unable to solve: sc = {}, diff = {}", stepcount, diff);
        }
        if diff < 1e-5 && stepcount >= 2 {
            break;
        }
        if stepcount % 1000 == 0 {
            println!("Stokes diff: {}", diff);
        }

        for i in 0..(bc.nx * bc.ny) {
            c.un[i as usize] = c.u[i as usize];
            c.vn[i as usize] = c.v[i as usize];
        }

        solve_pressure_poisson(bc, c);
        solve_stokes_momentum(bc, c);

        apply_vel_bc(bc, c);

        ud  = 0_f64;
        vd  = 0_f64;
        unt = 0_f64;
        vnt = 0_f64;

        for i in 0..(bc.nx * bc.ny) {
            ud += (c.u[i as usize].abs() - c.un[i as usize].abs()).abs();
            unt += c.un[i as usize].abs();

            vd += (c.v[i as usize].abs() - c.vn[i as usize].abs()).abs();
            vnt += c.vn[i as usize].abs();
        }

        if unt > 0. && vnt > 0. {
            diff = ( (ud / unt).abs() - (vd / vnt).abs() ) / 2.;
        } else {
            diff = 1_f64;
        }

        stepcount += 1;
    }
}


fn update_k(bc : &BasicContext, c : & mut Context) {
    let mut idx : usize;

    for j in 0..bc.ny {
        for i in 0..bc.nx {
            idx = (bc.nx * j + i) as usize;
            c.k[idx] = bc.refk * (1. - (bc.a * (c.t[idx] - bc.reftemp)));
        }
    }
}


fn update_nu(bc : &BasicContext, c : & mut Context) {
    let mut idx : usize;

    for j in 0..bc.ny {
        for i in 0..bc.nx {
            idx = (bc.nx * j + i) as usize;
            c.nu[idx] = bc.refnu * (-bc.theta * ( (c.t[idx] - bc.reftemp) / bc.reftemp)).exp();
        }
    }
}


fn update_rho(bc : &BasicContext, c : & mut Context) {
    let mut idx : usize;

    for j in 0..bc.ny {
        for i in 0..bc.nx {
            idx = (bc.nx * j + i) as usize;
            c.rho[idx] = bc.refrho * (1. - (bc.a * (c.t[idx] - bc.reftemp)));
        }
    }
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
            if norm < 0.{
                norm = 0f64;
            }
            if norm > 1. {
                norm = 1f64;
            }

            match norm {
                0f64       ... 0.1429f64  => print!("{}", bchar.blue()),
                0.1428f64  ... 0.28571f64 => print!("{}", bchar.cyan()),
                0.28570f64 ... 0.42857f64 => print!("{}", bchar.green()),
                0.42856f64 ... 0.57128f64 => print!("{}", bchar.white()),
                0.57127f64 ... 0.71429f64 => print!("{}", bchar.yellow()),
                0.71428f64 ... 0.85714f64 => print!("{}", bchar.red()),
                0.85713f64 ... 1f64       => print!("{}", bchar.magenta()),
                _                         => print!("{}", bchar.black()) // Error - this may an issue in the future
            }
        }
        print!("{}", "\n".normal());
    }
    println!("\n");

}


fn main() {
    let bc = BasicContext{ xmin    : 0.,
                           xmax    : 5.,
                           ymin    : 0.,
                           ymax    : 2.5,
                           nx      : 101,
                           ny      : 51, 
 
                           cp      : 60.,   
                           H       : 0.,    
                           dt      : 9e-5,  
                           refrho  : 100.,  
                           reftemp : 100.,  
                           refnu   : 1.,  
                           refk    : 100.,  
                           a       : 0.001,  
                           theta   : 1.5,  
                           g       : 0.0003,  
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
                         b   : StorMat::from_element(0.),
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

        if timestep > 10000 {
            break
        }

        solve_advection_diffusion(&bc, & mut c);
        update_rho(&bc, & mut c);
        update_nu(&bc, & mut c);
        update_k(&bc, & mut c);
        solve_flow(&bc, & mut c);

        if timestep % 1000 == 0 {
            print_field(&bc, &c.t, 0., 1000.);
        }

        timestep += 1;
        currenttime += bc.dt;
    }

    println!("Final time: {}", currenttime);
    print_field(&bc, &c.t, 0., 1000.);
}
