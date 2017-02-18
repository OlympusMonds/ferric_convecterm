
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
    u   : Vec<f64>,  // vel in x
    un  : Vec<f64>,  // vel in x
    v   : Vec<f64>,  // vel in y
    vn  : Vec<f64>,  // vel in y
    p   : Vec<f64>,  // pressure
    pn  : Vec<f64>,  // pressure
    rho : Vec<f64>,  // density
    nu  : Vec<f64>,  // viscosity
    t   : Vec<f64>,  // temp
    tn  : Vec<f64>,  // temp
    k   : Vec<f64>,  // conductivity
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

    for j in 0..bc.ny {
        idx = (bc.nx * j) as usize;

        // left wall, no temp escape
        c.t[idx] = c.t[idx + 1];

        // right wall, no temp escape
        let off = idx + (bc.nx - 1) as usize;
        c.t[off] = c.t[off - 1];
    }

    for i in 0..bc.nx {
        // bottom wall 
        c.v[i as usize] = 1000.;

        // top wall
        idx = (bc.nx * (bc.ny - 1) + i) as usize;
        c.v[idx] = 0.;
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
    let mut n : u32;  // north
    let mut e : usize;
    let mut s : usize;
    let mut w : usize;  // west

    for j in 0..bc.ny {
        for i in 0..bc.nx {
            m = (bc.nx * j + i) as usize;
            n = (bc.nx * (j-1) + i);
            s = (bc.nx * (j+1) + i) as usize;
            e = (bc.nx * j + (i+1)) as usize;
            w = (bc.nx * j + (i-1)) as usize;

            kx = c.k[m] + (c.tn[e] - 2. * c.tn[m] + c.tn[w]) / dx2;
            ky = c.k[m] + (c.tn[s] - 2. * c.tn[m] + c.tn[n]) / dx2;

        }
    }

    apply_thermal_bc(bc, c);
}


fn print_field(bc : &BasicContext, field : & Vec<f64>) {
    for jj in 0..bc.ny {
        let j = bc.ny - jj - 1;
        for i in 0..bc.nx {
            let idx = (bc.nx * j + i) as usize;

            print!("{} ", field[idx]);
        }
        print!("\n");
    }

}


fn main() {
    let bc = BasicContext{ xmin   : 0.,
                           xmax   : 2.,
                           ymin   : 0.,
                           ymax   : 2.,
                           nx     : 21,
                           ny     : 21, 
 
                           cp     : 60.,   
                           H      : 0.,    
                           dt     : 9e-5,  
                           refrho : 100.,  
                           g      : 9.81,  
                         };

    let ele = bc.nx * bc.ny;
    let dx = ( bc.xmax - bc.xmin ) / (bc.nx as f64);
    let dy = ( bc.ymax - bc.ymin ) / (bc.ny as f64);

    println!("dx = {}\ndy = {}\nele = {}", dx, dy, ele);

    let mut c = Context{ u   : vec![0.      as f64; ele as usize],  // vel in x
                         un  : vec![0.      as f64; ele as usize],
                         v   : vec![0.      as f64; ele as usize],
                         vn  : vec![0.      as f64; ele as usize],
                         p   : vec![0.      as f64; ele as usize],
                         pn  : vec![0.      as f64; ele as usize],
                         rho : vec![bc.refrho;      ele as usize],
                         nu  : vec![1.      as f64; ele as usize],
                         t   : vec![500.    as f64; ele as usize],
                         tn  : vec![500.    as f64; ele as usize],
                         k   : vec![100.    as f64; ele as usize],
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
    print_field(&bc, &c.p);

    apply_vel_bc(&bc, & mut c);
    apply_thermal_bc(&bc, & mut c);

    // Main loop
    loop {

        if timestep > 100 {
            break
        }

        solve_advection_diffusion(&bc, & mut c);

        timestep += 1;
    }
}
