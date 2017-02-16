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
    v   : Vec<f64>,  // vel in y
    p   : Vec<f64>,  // pressure
    rho : Vec<f64>,  // density
    nu  : Vec<f64>,  // viscosity
    t   : Vec<f64>,  // temp
    k   : Vec<f64>,  // conductivity
}


fn print_field(bc : BasicContext, field : Vec<f64>) {
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
                         v   : vec![0.      as f64; ele as usize],
                         p   : vec![0.      as f64; ele as usize],
                         rho : vec![bc.refrho;      ele as usize],
                         nu  : vec![1.      as f64; ele as usize],
                         t   : vec![500.    as f64; ele as usize],
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
            /*
            if j < 4 && i < nx/2 {
                t[idx] = 1000.;
            }
            else if j > ny-5 && i > nx/2 {
                t[idx] = 0.;
            }
            */
        }
    }

    print_field(bc, c.p);
}
